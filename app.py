import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math

# ==========================================
# 0. 预设数据 (写入代码，永久保存)
# ==========================================
# 数据来源：kdp_user_calibration.json
CALIB_PRESETS = {
    "自定义输入": {
        "k_ab": 1.0, "b_ab": 0.0, "k_h": 1.0, "b_h": 0.0
    },
    "刘冰": {
        "k_ab": 0.96224, "b_ab": 0.22727, "k_h": 0.78007, "b_h": 0.35606
    },
    "刘乃慎": {
        "k_ab": 0.979371, "b_ab": 0.0454546, "k_h": 0.755245, "b_h": 0.454545
    },
    "不缩放 (标准)": {
        "k_ab": 1.0, "b_ab": 0.0, "k_h": 1.0, "b_h": 0.0
    }
}

# ==========================================
# 核心物理常数与计算公式
# ==========================================

# KDP 物理参数
RHO_KDP_SOLID = 2.33  # g/cm³
RHO_KDP_SOLID_TAB4 = 2.34
DEFAULT_PHI = 47.0    # 锥面角
C_OVER_A = 0.936      # Tab1 cap ratio

# 密度拟合系数
RHO_A = 0.96286258
RHO_B = 0.00888087

# 溶解度表
SAT_SOLUBILITY_REF = {
    20: 22.6, 25: 26.3, 30: 29.8, 35: 31.7,
    40: 33.5, 45: 37.3, 50: 41.0, 55: 45.6, 60: 50.1
}
WATER_DENSITY_REF = {
    20: 0.9982, 25: 0.9970, 30: 0.9957, 35: 0.9940,
    40: 0.9922, 45: 0.9903, 50: 0.9880, 55: 0.9857, 60: 0.9832
}

def get_wt_percent_from_T(T):
    return 10.68 + 0.3616 * T

def get_T_from_wt_percent(wt):
    return (wt - 10.68) / 0.3616

def get_rho_sat(wt_percent):
    return RHO_A + RHO_B * wt_percent

def get_wt_from_rho(rho):
    if RHO_B == 0: return 0
    return (rho - RHO_A) / RHO_B

def wt_percent_from_solubility(S):
    return 100.0 * S / (100.0 + S)

def pyramid_height_from_phi(a, phi_deg):
    return (a / 2.0) * math.tan(math.radians(phi_deg))

def crystal_masses_g(a_mm, density_g_cm3, h_over_a, phi_deg):
    a_cm = a_mm / 10.0
    h_prism_cm = a_cm * h_over_a
    h_pyr_cm = (a_cm / 2.0) * math.tan(math.radians(phi_deg))
    V_total = (a_cm ** 2) * h_prism_cm + (1.0/3.0) * (a_cm ** 2) * h_pyr_cm
    return V_total * density_g_cm3, V_total

# ==========================================
# 页面 UI 构建
# ==========================================

st.set_page_config(page_title="KDP 综合计算平台 Web", layout="wide", page_icon="💎")
st.title("💎 KDP 综合计算平台 v4.2 (内置标定版)")

# 创建四个选项卡
tab1, tab2, tab3, tab4 = st.tabs([
    "① 晶体质量/3D", 
    "② 饱和溶液换算", 
    "③ 配液方案 (Pro)", 
    "④ 生长控制"
])

# ==========================================
# Tab 1: 晶体质量 / 3D
# ==========================================
with tab1:
    st.header("1. 晶体质量计算 & 3D 可视化")
    
    col1_1, col1_2 = st.columns([1, 1.5])
    
    with col1_1:
        st.subheader("观测数据输入")
        with st.form("tab1_form"):
            obs_a = st.number_input("a 观测值 (mm)", value=26.5)
            obs_b = st.number_input("b 观测值 (mm)", value=26.5)
            obs_h = st.number_input("h 观测值 (mm)", value=35.0)
            
            st.markdown("---")
            st.markdown("**标定参数选择**")
            
            # --- 新增：标定参数选择器 ---
            # 1. 也是放在 form 里的 selectbox
            user_select = st.selectbox("选择使用者 / 预设参数", options=list(CALIB_PRESETS.keys()), index=1) 
            # 默认为index=1 (刘冰)，方便直接用
            
            # 2. 根据选择获取默认值
            current_vals = CALIB_PRESETS[user_select]
            
            st.caption("公式: obs = k · real + b  (下方数值可手动修改)")
            c1, c2 = st.columns(2)
            # 注意：这里 value 设置为从字典取出的值，这样切换下拉菜单时，输入框会自动变
            k_ab = c1.number_input("k (ab方向)", value=current_vals["k_ab"], format="%.5f")
            b_ab = c2.number_input("b (ab方向)", value=current_vals["b_ab"], format="%.5f")
            
            c3, c4 = st.columns(2)
            k_h = c3.number_input("k (h方向)", value=current_vals["k_h"], format="%.5f")
            b_h = c4.number_input("b (h方向)", value=current_vals["b_h"], format="%.5f")
            # ---------------------------
            
            h_mode = st.radio("h 的含义", ["总高度 Htot", "柱体高度 Hp"])
            
            submit_t1 = st.form_submit_button("计算并生成 3D 模型", type="primary")

    if submit_t1:
        try:
            # 1. 反算真实尺寸
            real_a = (obs_a - b_ab) / k_ab
            real_b = (obs_b - b_ab) / k_ab
            real_h = (obs_h - b_h) / k_h
            
            # 2. 计算几何
            def calc_geom(a, b, h, mode):
                a_eq = (a + b) / 2.0
                Hy = 0.5 * a_eq * C_OVER_A
                if mode == "总高度 Htot":
                    Htot = h
                    Hp = h - Hy
                else:
                    Hp = h
                    Htot = h + Hy
                
                a_cm = a_eq / 10.0
                Hp_cm = Hp / 10.0
                Hy_cm = Hy / 10.0
                V_total = (a_cm**2 * Hp_cm) + (1.0/3.0 * a_cm**2 * Hy_cm)
                mass = V_total * RHO_KDP_SOLID
                return a, b, Hp, Hy, Htot, V_total, mass

            res_raw = calc_geom(obs_a, obs_b, obs_h, h_mode)
            res_cal = calc_geom(real_a, real_b, real_h, h_mode)
            
            with col1_1:
                st.success(f"校正后质量: **{res_cal[6]:.2f} g**")
                df_res = pd.DataFrame({
                    "参数": ["a (mm)", "b (mm)", "柱高 Hp (mm)", "总高 H (mm)", "体积 (cm³)", "质量 (g)"],
                    "校正前 (Raw)": [f"{x:.2f}" for x in [res_raw[0], res_raw[1], res_raw[2], res_raw[4], res_raw[5], res_raw[6]]],
                    "校正后 (Real)": [f"{x:.2f}" for x in [res_cal[0], res_cal[1], res_cal[2], res_cal[4], res_cal[5], res_cal[6]]]
                })
                st.dataframe(df_res, hide_index=True)

            with col1_2:
                # 3D 绘图
                def get_mesh_data(a, b, hp, hy, color, opacity):
                    dx, dy = a/2, b/2
                    x_p = [-dx, dx, dx, -dx, -dx, dx, dx, -dx]
                    y_p = [-dy, -dy, dy, dy, -dy, -dy, dy, dy]
                    z_p = [0, 0, 0, 0, hp, hp, hp, hp]
                    i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
                    j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
                    k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]
                    prism = go.Mesh3d(x=x_p, y=y_p, z=z_p, i=i, j=j, k=k, color=color, opacity=opacity, name='Prism')
                    x_c = [-dx, dx, dx, -dx, 0]
                    y_c = [-dy, -dy, dy, dy, 0]
                    z_c = [hp, hp, hp, hp, hp+hy]
                    i_c = [0, 1, 2, 3]
                    j_c = [1, 2, 3, 0]
                    k_c = [4, 4, 4, 4]
                    cap = go.Mesh3d(x=x_c, y=y_c, z=z_c, i=i_c, j=j_c, k=k_c, color=color, opacity=opacity, name='Cap')
                    return [prism, cap]

                fig = go.Figure()
                meshes_raw = get_mesh_data(res_raw[0], res_raw[1], res_raw[2], res_raw[3], 'gray', 0.15)
                for m in meshes_raw: fig.add_trace(m)
                meshes_cal = get_mesh_data(res_cal[0], res_cal[1], res_cal[2], res_cal[3], '#0078D4', 0.6)
                for m in meshes_cal: fig.add_trace(m)
                
                fig.update_layout(scene=dict(aspectmode='data'), height=500, margin=dict(l=0, r=0, b=0, t=0))
                st.plotly_chart(fig, use_container_width=True)

        except Exception as e:
            st.error(f"计算错误: {e}")

# ==========================================
# Tab 2: 饱和溶液换算
# ==========================================
with tab2:
    st.header("2. 饱和溶液 质量↔体积")
    t2_col1, t2_col2 = st.columns([1, 1])
    with t2_col1:
        st.subheader("A. 快速查询")
        temp_query = st.slider("选择温度 (°C)", 20, 60, 40, step=5)
        S_ref = SAT_SOLUBILITY_REF.get(temp_query, 0)
        rho_w_ref = WATER_DENSITY_REF.get(temp_query, 0)
        wt_ref = wt_percent_from_solubility(S_ref)
        rho_sat_ref = get_rho_sat(wt_ref)
        st.info(f"**{temp_query}°C 参考数据:**\n* S: {S_ref} g/100g水\n* wt%: {wt_ref:.2f}%\n* ρ_sat: {rho_sat_ref:.4f} g/mL")
        
        st.subheader("B. 计算器")
        calc_mode = st.selectbox("计算模式", ["体积 → 质量", "质量 → 体积"])
        if calc_mode == "体积 → 质量":
            vol_in = st.number_input("输入溶液体积 (mL)", value=1000.0)
            if st.button("计算质量"):
                st.success(f"溶液总质量: {vol_in * rho_sat_ref:.2f} g")
        else:
            mass_in = st.number_input("输入溶液质量 (g)", value=1200.0)
            if st.button("计算体积"):
                st.success(f"溶液总体积: {mass_in / rho_sat_ref:.2f} mL")
                
    with t2_col2:
        st.subheader("参考数据表")
        data_rows = []
        for t in range(20, 61, 5):
            s = SAT_SOLUBILITY_REF[t]
            w = wt_percent_from_solubility(s)
            r = get_rho_sat(w)
            data_rows.append([t, s, f"{w:.2f}", f"{r:.4f}", WATER_DENSITY_REF[t]])
        df_ref = pd.DataFrame(data_rows, columns=["温度(°C)", "S (g/100g)", "wt%", "ρ_sat", "ρ_H2O"])
        st.dataframe(df_ref, hide_index=True, use_container_width=True)

# ==========================================
# Tab 3: 配液方案 (Pro)
# ==========================================
with tab3:
    st.header("3. 饱和溶液配制 & 生长预测 (Pro)")
    st.markdown("#### A. 配液计算 (任意选2算3)")
    calc_type = st.selectbox("已知条件组合", ["已知: 总重M & 温度T", "已知: 溶剂水W & 溶质S", "已知: 体积V & 温度T", "已知: 溶质S & 温度T"])
    c3_1, c3_2, c3_3 = st.columns(3)
    res_M, res_T, res_W, res_S, res_V = None, None, None, None, None
    
    with c3_1:
        if "总重M" in calc_type: in_M = st.number_input("总重 M (g)", value=1200.0)
        if "温度T" in calc_type: in_T = st.number_input("饱和温度 T (°C)", value=40.0, key="t3_T")
        if "溶剂水W" in calc_type: in_W = st.number_input("溶剂水 W (g)", value=800.0)
        if "溶质S" in calc_type: in_S = st.number_input("溶质 S (g)", value=400.0)
        if "体积V" in calc_type: in_V = st.number_input("体积 V (mL)", value=1000.0)
            
    if st.button("计算配方", type="primary"):
        try:
            if calc_type == "已知: 总重M & 温度T":
                res_M, res_T = in_M, in_T
                wt = get_wt_percent_from_T(in_T)
                res_S, res_W, res_V = in_M*(wt/100.0), in_M*(1-wt/100.0), in_M/get_rho_sat(wt)
            elif calc_type == "已知: 溶剂水W & 溶质S":
                res_W, res_S = in_W, in_S
                res_M = res_W + res_S
                wt = (res_S/res_M)*100.0
                res_T, res_V = get_T_from_wt_percent(wt), res_M/get_rho_sat(wt)
            elif calc_type == "已知: 体积V & 温度T":
                res_V, res_T = in_V, in_T
                wt = get_wt_percent_from_T(in_T)
                res_M = in_V * get_rho_sat(wt)
                res_S, res_W = res_M*(wt/100.0), res_M*(1-wt/100.0)
            elif calc_type == "已知: 溶质S & 温度T":
                res_S, res_T = in_S, in_T
                wt = get_wt_percent_from_T(in_T)
                res_M = res_S/(wt/100.0)
                res_W, res_V = res_M-res_S, res_M/get_rho_sat(wt)

            st.session_state['recipe_res'] = (res_M, res_T, res_W, res_S, res_V)
            st.success("计算完成")
            cc1, cc2, cc3, cc4, cc5 = st.columns(5)
            cc1.metric("总重 M", f"{res_M:.1f}")
            cc2.metric("温度 T", f"{res_T:.1f}")
            cc3.metric("溶质 S", f"{res_S:.1f}")
            cc4.metric("溶剂 W", f"{res_W:.1f}")
            cc5.metric("体积 V", f"{res_V:.1f}")
        except Exception as e: st.error(f"出错: {e}")
            
    st.markdown("#### B. 生长预测 (需先完成上方计算)")
    if 'recipe_res' in st.session_state:
        mM, mT, mW, mS, mV = st.session_state['recipe_res']
        with st.form("growth_pred"):
            st.write(f"当前状态: {mT:.1f}°C, 总重 {mM:.1f}g")
            gc1, gc2 = st.columns(2)
            target_T = gc1.number_input("降温至 (°C)", value=20.0)
            seed_mode = st.radio("籽晶模式", ["点籽晶 (Mode A)", "片状籽晶 (Mode B)"], horizontal=True)
            if "Mode A" in seed_mode: ratio_ha = gc2.number_input("高宽比 H/L", value=1.0)
            else: seed_L = gc2.number_input("籽晶边长 L (mm)", value=20.0)
                
            if st.form_submit_button("预测生长"):
                wt_end = get_wt_percent_from_T(target_T)
                S_end_sat = mW / (1.0 - wt_end/100.0) * (wt_end/100.0)
                dS = mS - S_end_sat
                if dS <= 0: st.warning(f"无析出 (dS={dS:.2f}g)")
                else:
                    st.success(f"析出质量: {dS:.2f} g")
                    V_crys = (dS / RHO_KDP_SOLID) * 1000.0
                    tan_g = math.tan(math.radians(DEFAULT_PHI))
                    if "Mode A" in seed_mode:
                        factor = ratio_ha - 1.0/(3.0*tan_g)
                        if factor <= 0: st.error("H/L太小")
                        else:
                            L_fin = (V_crys/factor)**(1/3.0)
                            st.info(f"L = {L_fin:.1f} mm, H = {L_fin*ratio_ha:.1f} mm")
                    else:
                        h_cap = seed_L / (2.0 * tan_g)
                        V_cap = (1.0/3.0)*(seed_L**2)*h_cap
                        if V_crys < V_cap: st.info("未长满锥帽")
                        else: st.info(f"生长总高: {h_cap + (V_crys-V_cap)/(seed_L**2):.1f} mm")

# ==========================================
# Tab 4: 生长控制
# ==========================================
with tab4:
    st.header("4. 生长过程表 & 5天温控")
    with st.expander("参数设置", expanded=True):
        c4_1, c4_2, c4_3, c4_4 = st.columns(4)
        M0 = c4_1.number_input("初始溶液 M0 (g)", value=2000.0)
        T0 = c4_2.number_input("初始饱和 T0 (°C)", value=55.0)
        a_min = c4_3.number_input("a 最小 (mm)", value=1.0)
        a_max = c4_4.number_input("a 最大 (mm)", value=80.0)
        c4_5, c4_6, c4_7, c4_8 = st.columns(4)
        step_val = c4_5.number_input("步长 (mm)", value=2.0)
        ha_ratio = c4_6.number_input("h/a 比例", value=1.0)
        phi_val = c4_7.number_input("φ 角度", value=DEFAULT_PHI)
        
    if st.button("生成生长过程表", type="primary"):
        C0 = get_wt_percent_from_T(T0)
        solute0 = M0 * (C0/100.0)
        growth_data = []
        curr_a = a_min
        while curr_a <= a_max + 1e-9:
            m_crys, V_crys = crystal_masses_g(curr_a, RHO_KDP_SOLID_TAB4, ha_ratio, phi_val)
            sol_m = max(M0 - m_crys, 1e-9)
            sol_s = max(solute0 - m_crys, 0)
            T_sat = get_T_from_wt_percent(100.0 * (sol_s / sol_m))
            if T_sat < 10.0: break
            h_prism = ha_ratio * curr_a
            growth_data.append({
                "a (mm)": round(curr_a, 2), "H (mm)": round(h_prism+pyramid_height_from_phi(curr_a, phi_val), 2),
                "晶体重(g)": round(m_crys, 2), "T_sat (°C)": round(T_sat, 2)
            })
            curr_a += step_val
        st.session_state['df_growth'] = pd.DataFrame(growth_data)
        st.success(f"已生成 {len(growth_data)} 条数据")

    if 'df_growth' in st.session_state:
        st.dataframe(st.session_state['df_growth'], use_container_width=True, height=250)
        st.subheader("速率计算 & 温控")
        rc1, rc2, rc3 = st.columns(3)
        df_s = st.session_state['df_growth']
        idx1 = rc1.number_input("起始Index", 0, len(df_s)-1, 0)
        idx2 = rc2.number_input("结束Index", 0, len(df_s)-1, min(5, len(df_s)-1))
        dt_h = rc3.number_input("时间间隔 (h)", value=24.0)
        
        if st.button("计算速率"):
            r1, r2 = df_s.iloc[idx1], df_s.iloc[idx2]
            v_a = (r2["a (mm)"] - r1["a (mm)"]) / (dt_h/24.0)
            st.info(f"a方向速度 = **{v_a:.2f} mm/天**")
            st.session_state['calc_va'] = v_a
            st.session_state['calc_state'] = r2
            
    if 'calc_va' in st.session_state:
        with st.form("plan_5day"):
            st.write(f"基准速度: {st.session_state['calc_va']:.2f} mm/天")
            T_curr = st.number_input("当前T_now (°C)", value=40.0)
            offset = st.number_input("过冷度偏移 ΔT", value=0.0)
            if st.form_submit_button("生成方案"):
                start_row = st.session_state['calc_state']
                C_start = get_wt_percent_from_T(T0) # 近似使用初始配方
                solute_tot = M0 * (C_start/100.0)
                plan_res = []
                for d in range(1, 6):
                    a_new = start_row["a (mm)"] + st.session_state['calc_va'] * d
                    m_new, _ = crystal_masses_g(a_new, RHO_KDP_SOLID_TAB4, ha_ratio, phi_val)
                    sol_m = M0 - m_new
                    sol_s = solute_tot - m_new
                    Tsat_new = get_T_from_wt_percent(100.0 * (sol_s/sol_m))
                    plan_res.append({"Day": d, "a": f"{a_new:.2f}", "T_sat": f"{Tsat_new:.2f}", "T_set": f"{Tsat_new-(start_row['T_sat (°C)']-T_curr)-offset:.2f}"})
                st.table(pd.DataFrame(plan_res))
