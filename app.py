import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math

# ==========================================
# 0. 预设数据 (写入代码，永久保存)
# ==========================================
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
# 核心物理常数与计算公式 (源自 kdp_common / tab2)
# ==========================================
RHO_KDP_SOLID = 2.33   # Tab1/3 使用
RHO_KDP_SOLID_TAB4 = 2.34 # Tab4 使用 (保持源码差异)
DEFAULT_PHI = 47.0
C_OVER_A = 0.936

# 密度拟合系数
RHO_A = 0.96286258
RHO_B = 0.00888087

# 溶解度参考表
SAT_SOLUBILITY_REF = {20:22.6, 25:26.3, 30:29.8, 35:31.7, 40:33.5, 45:37.3, 50:41.0, 55:45.6, 60:50.1}
WATER_DENSITY_REF = {20:0.9982, 25:0.9970, 30:0.9957, 35:0.9940, 40:0.9922, 45:0.9903, 50:0.9880, 55:0.9857, 60:0.9832}

# --- 公式函数 ---
def get_wt_percent_from_T(T):
    return 10.68 + 0.3616 * T

def get_T_from_wt_percent(wt):
    return (wt - 10.68) / 0.3616

def get_rho_sat(wt_percent):
    return RHO_A + RHO_B * wt_percent

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
st.set_page_config(page_title="KDP 综合计算平台", layout="wide", page_icon="💎")

st.title("💎 KDP 综合计算平台 v4.5 (全功能复原版)")

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
    st.header("1. 晶体质量计算 & 3D 可视化 (含个人标定)")
    
    col1_L, col1_R = st.columns([1, 1.2]) 
    
    with col1_L:
        st.subheader("A. 观测数据输入")
        c1, c2, c3 = st.columns(3)
        obs_a = c1.number_input("a 观测值 (mm)", value=26.5)
        obs_b = c2.number_input("b 观测值 (mm)", value=26.5)
        obs_h = c3.number_input("h 观测值 (mm)", value=35.0)
        
        st.markdown("---")
        st.subheader("B. 线性标定参数")
        # 标定选择
        user_select = st.selectbox("选择使用者 / 预设参数", list(CALIB_PRESETS.keys()), index=1)
        vals = CALIB_PRESETS[user_select]
        
        st.info(f"当前公式：观测值 = k × 真实值 + b (即 Real = (Obs - b) / k)")
        
        ck1, ck2 = st.columns(2)
        k_ab = ck1.number_input("k (ab方向)", value=vals["k_ab"], format="%.5f")
        b_ab = ck2.number_input("b (ab方向)", value=vals["b_ab"], format="%.5f")
        
        ck3, ck4 = st.columns(2)
        k_h = ck3.number_input("k (h方向)", value=vals["k_h"], format="%.5f")
        b_h = ck4.number_input("b (h方向)", value=vals["b_h"], format="%.5f")
        
        h_mode = st.radio("h 的含义", ["总高度 Htot (包含锥帽)", "柱体高度 Hp (仅柱面部分)"], horizontal=True)

        if st.button("计算并绘制 3D 模型", type="primary", use_container_width=True):
            try:
                # 1. 核心计算逻辑 (源自 kdp_tab1)
                real_a = (obs_a - b_ab) / k_ab
                real_b = (obs_b - b_ab) / k_ab
                real_h = (obs_h - b_h) / k_h
                
                # 防错
                if real_a <=0 or real_b <=0 or real_h <=0:
                    st.error("校正后尺寸出现负值或零，请检查标定参数。")
                    st.stop()

                def calc_details(a, b, h, mode):
                    a_eq = (a + b) / 2.0
                    Hy = 0.5 * a_eq * C_OVER_A
                    if "总高度" in mode:
                        Htot = h
                        Hp = h - Hy
                    else:
                        Hp = h
                        Htot = h + Hy
                    
                    if Hp < 0: Hp = 0 # 物理保护
                    
                    a_cm = a_eq / 10.0
                    Hp_cm = Hp / 10.0
                    Hy_cm = Hy / 10.0
                    
                    V_prism = a_cm**2 * Hp_cm
                    V_pyr = (1.0/3.0) * a_cm**2 * Hy_cm
                    V_total = V_prism + V_pyr
                    mass = V_total * RHO_KDP_SOLID
                    return a, b, a_eq, Hp, Hy, Htot, V_prism, V_pyr, V_total, mass

                raw = calc_details(obs_a, obs_b, obs_h, h_mode)
                cal = calc_details(real_a, real_b, real_h, h_mode)
                
                # 2. 结果展示 (复原表格列)
                st.success(f"校正后晶体质量: **{cal[9]:.2f} g**")
                
                # 构造详细对比表
                df_res = pd.DataFrame({
                    "参数指标": [
                        "a (mm)", "b (mm)", "等效边长 a_eq (mm)", 
                        "柱体高度 Hp (mm)", "锥帽高度 Hy (mm)", "总高度 Htot (mm)",
                        "柱体体积 (cm³)", "锥帽体积 (cm³)", "总体积 (cm³)", "质量 (g)"
                    ],
                    "校正前 (Raw)": [f"{x:.2f}" for x in raw],
                    "校正后 (Real)": [f"{x:.2f}" for x in cal]
                })
                st.table(df_res) # 使用 table 展示全部内容
                
                # 存入 Session 供右侧绘图
                st.session_state['t1_res'] = (obs_a, obs_b, raw[3], raw[4], real_a, real_b, cal[3], cal[4])
                
            except Exception as e:
                st.error(f"计算发生错误: {e}")

    with col1_R:
        # 3D 绘图区域
        st.subheader("C. 3D 模型预览")
        if 't1_res' in st.session_state:
            oa, ob, ohp, ohy, ra, rb, rhp, rhy = st.session_state['t1_res']
            
            # Plotly 绘图函数
            def get_mesh(a, b, hp, hy, color, opac, name):
                dx, dy = a/2, b/2
                # 柱体
                xp = [-dx, dx, dx, -dx, -dx, dx, dx, -dx]
                yp = [-dy, -dy, dy, dy, -dy, -dy, dy, dy]
                zp = [0, 0, 0, 0, hp, hp, hp, hp]
                i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
                j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
                k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]
                prism = go.Mesh3d(x=xp, y=yp, z=zp, i=i, j=j, k=k, color=color, opacity=opac, name=name+'-柱')
                # 锥帽
                xc = [-dx, dx, dx, -dx, 0]
                yc = [-dy, -dy, dy, dy, 0]
                zc = [hp, hp, hp, hp, hp+hy]
                ic, jc, kc = [0, 1, 2, 3], [1, 2, 3, 0], [4, 4, 4, 4]
                cap = go.Mesh3d(x=xc, y=yc, z=zc, i=ic, j=jc, k=kc, color=color, opacity=opac, name=name+'-帽')
                return [prism, cap]

            fig = go.Figure()
            # 绘制灰色虚影 (校正前)
            for m in get_mesh(oa, ob, ohp, ohy, 'gray', 0.15, 'Raw'): fig.add_trace(m)
            # 绘制蓝色实体 (校正后)
            for m in get_mesh(ra, rb, rhp, rhy, '#0078D4', 0.7, 'Real'): fig.add_trace(m)
            
            fig.update_layout(
                scene=dict(aspectmode='data', xaxis_title='X', yaxis_title='Y', zaxis_title='Z'), 
                margin=dict(l=0, r=0, b=0, t=20),
                height=600
            )
            st.plotly_chart(fig, use_container_width=True)
            st.caption("注：灰色半透明为校正前(原始观测)，蓝色实体为校正后(真实尺寸)。")
        else:
            st.info("👈 请在左侧输入数据并点击计算按钮")

# ==========================================
# Tab 2: 饱和溶液换算
# ==========================================
with tab2:
    st.header("2. 饱和溶液 质量↔体积 与 预混估算")
    
    col2_L, col2_R = st.columns([1, 1])
    
    with col2_L:
        st.subheader("A. 饱和点参考数据")
        temp_query = st.slider("选择温度 (°C)", 20, 60, 40, step=5)
        
        # 实时查询数据
        S_ref = SAT_SOLUBILITY_REF.get(temp_query, 0)
        wt_ref = wt_percent_from_solubility(S_ref)
        rho_sat_ref = get_rho_sat(wt_ref)
        rho_w_ref = WATER_DENSITY_REF.get(temp_query, 0)
        
        # 使用卡片展示详情
        st.info(f"""
        **{temp_query}°C 下的饱和 KDP 溶液性质：**
        * 溶解度 S: **{S_ref}** g/100g水
        * 质量浓度 wt%: **{wt_ref:.2f}%**
        * 饱和密度 ρ_sat: **{rho_sat_ref:.4f}** g/mL
        * 纯水密度 ρ_H2O: **{rho_w_ref:.4f}** g/mL
        """)
        
        st.markdown("---")
        st.subheader("B. 换算计算器")
        mode = st.radio("计算模式", ["已知体积 V → 求质量 M", "已知质量 M → 求体积 V"], horizontal=True)
        
        val_in = st.number_input("请输入数值 (mL 或 g)", value=1000.0, format="%.2f")
        
        if st.button("开始换算", use_container_width=True):
            if "已知体积" in mode:
                vol_in = val_in
                mass_out = vol_in * rho_sat_ref
                txt_res = f"溶液总质量: {mass_out:.2f} g"
                total_mass = mass_out
                total_vol = vol_in
            else:
                mass_in = val_in
                vol_out = mass_in / rho_sat_ref
                txt_res = f"溶液总体积: {vol_out:.2f} mL"
                total_mass = mass_in
                total_vol = vol_out
            
            st.success(f"计算结果：**{txt_res}**")
            
            # 找回“预混体积”计算 (源自 kdp_tab2.py)
            wt_frac = wt_ref / 100.0
            m_solute = total_mass * wt_frac
            m_water = total_mass * (1.0 - wt_frac)
            
            V_powder = m_solute / RHO_KDP_SOLID
            V_water_vol = m_water / rho_w_ref
            V_premix = V_powder + V_water_vol
            
            st.markdown("#### 详细组成与预混估算")
            st.write(f"1. **溶质 (KDP粉末)**: 质量 {m_solute:.2f} g | 对应粉末体积约 {V_powder:.2f} mL")
            st.write(f"2. **溶剂 (纯水)**: 质量 {m_water:.2f} g | 对应水体积 {V_water_vol:.2f} mL")
            st.warning(f"👉 **预混体积估算**: {V_powder:.2f} + {V_water_vol:.2f} = **{V_premix:.2f} mL** (通常小于最终饱和体积)")

    with col2_R:
        st.subheader("C. 20-60°C 参考数据表")
        # 复原完整表格
        data_rows = []
        for t in range(20, 65, 5):
            s = SAT_SOLUBILITY_REF.get(t, "-")
            if s != "-":
                w = wt_percent_from_solubility(s)
                r = get_rho_sat(w)
                rw = WATER_DENSITY_REF.get(t, "-")
                data_rows.append([t, s, f"{w:.2f}", f"{r:.4f}", rw])
        
        df_ref = pd.DataFrame(data_rows, columns=["温度(°C)", "溶解度 S", "浓度 wt%", "饱和密度 ρ_sat", "水密度 ρ_H2O"])
        st.dataframe(df_ref, hide_index=True, use_container_width=True, height=600)

# ==========================================
# Tab 3: 配液方案 (Pro)
# ==========================================
with tab3:
    st.header("3. 饱和溶液配制 & 生长预测 (Pro)")
    
    st.subheader("A. 配液计算 (任意勾选2个已知量)")
    
    # 使用直观的单选按钮
    calc_type = st.radio(
        "请选择已知条件组合：", 
        ["已知: 溶液总重 M & 饱和温度 T", 
         "已知: 溶剂水重 W & 溶质干重 S", 
         "已知: 溶液体积 V & 饱和温度 T", 
         "已知: 溶质干重 S & 饱和温度 T"],
        horizontal=True
    )
    
    # 动态输入区
    c3_1, c3_2, c3_3 = st.columns(3)
    in_val1, in_val2 = 0.0, 0.0
    
    with c3_1:
        if "总重 M" in calc_type: in_val1 = st.number_input("溶液总重 M (g)", value=1200.0)
        elif "水重 W" in calc_type: in_val1 = st.number_input("溶剂水重 W (g)", value=800.0)
        elif "体积 V" in calc_type: in_val1 = st.number_input("溶液体积 V (mL)", value=1000.0)
        elif "干重 S" in calc_type: in_val1 = st.number_input("溶质干重 S (g)", value=400.0)
        
    with c3_2:
        if "温度 T" in calc_type: in_val2 = st.number_input("饱和温度 T (°C)", value=40.0)
        elif "干重 S" in calc_type and "水重 W" in calc_type: in_val2 = st.number_input("溶质干重 S (g)", value=400.0)
        
    with c3_3:
        st.write("##") # 占位符，对齐按钮
        btn_calc_recipe = st.button("计算配方详情", type="primary", use_container_width=True)

    if btn_calc_recipe:
        try:
            rM, rT, rW, rS, rV = 0,0,0,0,0
            # 逻辑分支 (源自 kdp_tab3.py)
            if "总重 M" in calc_type:
                wt = get_wt_percent_from_T(in_val2)
                rM, rT = in_val1, in_val2
                rS = rM * wt/100.0; rW = rM - rS; rV = rM/get_rho_sat(wt)
            elif "水重 W" in calc_type: 
                rW, rS = in_val1, in_val2; rM = rW+rS
                wt = rS/rM*100.0; rT = get_T_from_wt_percent(wt); rV = rM/get_rho_sat(wt)
            elif "体积 V" in calc_type:
                rV, rT = in_val1, in_val2; wt = get_wt_percent_from_T(rT)
                rM = rV * get_rho_sat(wt); rS = rM*wt/100.0; rW = rM-rS
            elif "干重 S" in calc_type:
                rS, rT = in_val1, in_val2; wt = get_wt_percent_from_T(rT)
                rM = rS/(wt/100.0); rW = rM-rS; rV = rM/get_rho_sat(wt)
                
            st.session_state['t3_res'] = (rM, rT, rW, rS, rV)
            
            st.success("计算完成！详细配方如下：")
            cc1, cc2, cc3, cc4, cc5 = st.columns(5)
            cc1.metric("溶液总重 M", f"{rM:.1f} g")
            cc2.metric("饱和温度 T", f"{rT:.1f} °C")
            cc3.metric("溶质干重 S", f"{rS:.1f} g")
            cc4.metric("溶剂水重 W", f"{rW:.1f} g")
            cc5.metric("溶液体积 V", f"{rV:.1f} mL")
            
        except Exception as e:
            st.error(f"计算出错: {e}")

    st.markdown("---")
    st.subheader("B. 生长预测 (基于上方计算结果)")
    
    if 't3_res' in st.session_state:
        mM, mT, mW, mS, mV = st.session_state['t3_res']
        st.caption(f"当前溶液状态: 饱和点 {mT:.1f}°C, 总溶质 {mS:.1f}g")
        
        g_c1, g_c2, g_c3 = st.columns(3)
        t_final = g_c1.number_input("1. 降温至目标温度 (°C)", value=20.0)
        seed_mode = g_c2.selectbox("2. 籽晶模式", ["模式A：点籽晶 (3D生长)", "模式B：片状籽晶 (仅Z向)"])
        
        param_val = 0.0
        if "模式A" in seed_mode:
            param_val = g_c3.number_input("3. 设定高宽比 (H/L)", value=1.0)
        else:
            param_val = g_c3.number_input("3. 籽晶边长 L (mm)", value=20.0)
            
        if st.button("预测晶体尺寸", use_container_width=True):
            # 1. 析出量计算
            wt_end = get_wt_percent_from_T(t_final)
            S_end_sat = mW / (1.0 - wt_end/100.0) * (wt_end/100.0)
            dS = mS - S_end_sat
            
            if dS <= 0:
                st.warning(f"无法生长：目标温度下的溶解度更高，或析出量为负 (dS={dS:.2f}g)")
            else:
                st.success(f"理论析出晶体质量: **{dS:.2f} g**")
                
                # 2. 几何尺寸推算
                V_crys_mm3 = (dS / RHO_KDP_SOLID) * 1000.0
                tan_g = math.tan(math.radians(DEFAULT_PHI))
                
                if "模式A" in seed_mode:
                    ratio = param_val
                    factor = ratio - 1.0/(3.0*tan_g)
                    if factor <= 0:
                        st.error("错误：高宽比太小，无法形成完整四棱锥结构。")
                    else:
                        L_final = (V_crys_mm3 / factor)**(1/3.0)
                        H_final = L_final * ratio
                        st.info(f"预测最终尺寸: 边长 L = **{L_final:.1f} mm**, 总高 H = **{H_final:.1f} mm**")
                else:
                    L_plate = param_val
                    h_cap_full = L_plate / (2.0 * tan_g)
                    V_cap_full = (1.0/3.0)*(L_plate**2)*h_cap_full
                    
                    if V_crys_mm3 < V_cap_full:
                        st.info("晶体尚未长满锥帽部分。")
                    else:
                        V_prism = V_crys_mm3 - V_cap_full
                        h_prism_added = V_prism / (L_plate**2)
                        H_total = h_cap_full + h_prism_added
                        st.info(f"预测生长总高: **{H_total:.1f} mm** (其中柱面增高 {h_prism_added:.1f} mm)")
    else:
        st.info("请先完成上方 A 部分的配液计算。")

# ==========================================
# Tab 4: 生长控制
# ==========================================
with tab4:
    st.header("4. 生长过程表 & 5天温控方案")
    
    st.subheader("A. 初始参数设置")
    # 不折叠，直接展示所有参数
    col4_1, col4_2, col4_3, col4_4 = st.columns(4)
    M0 = col4_1.number_input("初始溶液总重 M0 (g)", value=2000.0)
    T0 = col4_2.number_input("初始饱和温度 T0 (°C)", value=55.0)
    a_min = col4_3.number_input("起始尺寸 a_min (mm)", value=1.0)
    a_max = col4_4.number_input("结束尺寸 a_max (mm)", value=80.0)
    
    col4_5, col4_6, col4_7, col4_8 = st.columns(4)
    step_val = col4_5.number_input("步长 step (mm)", value=2.0)
    ha_ratio = col4_6.number_input("高宽比 h/a", value=1.0)
    phi_val = col4_7.number_input("锥面角 φ (°)", value=DEFAULT_PHI)
    rho_val = col4_8.number_input("密度 ρ (g/cm³)", value=RHO_KDP_SOLID_TAB4)
    
    if st.button("生成详细生长过程表", type="primary"):
        # 计算初始状态
        C0 = get_wt_percent_from_T(T0)
        solute0 = M0 * (C0 / 100.0)
        
        growth_data = []
        curr_a = a_min
        
        while curr_a <= a_max + 1e-9:
            # 1. 晶体几何与质量
            m_crys, V_crys = crystal_masses_g(curr_a, rho_val, ha_ratio, phi_val)
            
            # 2. 溶液状态更新
            sol_m = max(M0 - m_crys, 1e-9)
            sol_s = max(solute0 - m_crys, 0.0)
            C_now = 100.0 * (sol_s / sol_m)
            T_sat = get_T_from_wt_percent(C_now)
            
            if T_sat < 10.0: break # 温度过低停止
            
            # 3. 几何高度
            h_prism = ha_ratio * curr_a
            h_pyr = pyramid_height_from_phi(curr_a, phi_val)
            H_total = h_prism + h_pyr
            
            # 找回所有列 (源自 kdp_tab4)
            growth_data.append({
                "a (mm)": round(curr_a, 2),
                "柱高 h (mm)": round(h_prism, 2),
                "总高 H (mm)": round(H_total, 2),
                "晶体体积 (cm³)": round(V_crys, 2),
                "晶体质量 (g)": round(m_crys, 2),
                "溶液质量 (g)": round(sol_m, 2),
                "剩余溶质 (g)": round(sol_s, 2),
                "浓度 C (g/100g)": round(C_now, 2),
                "饱和点 Tsat (°C)": round(T_sat, 2)
            })
            
            curr_a += step_val
            
        st.session_state['df_growth'] = pd.DataFrame(growth_data)
        st.success(f"生成成功！共 {len(growth_data)} 条数据")

    # 展示表格
    if 'df_growth' in st.session_state:
        df_show = st.session_state['df_growth']
        st.markdown("### 生长过程数据表")
        st.dataframe(df_show, use_container_width=True, height=300)
        
        st.markdown("---")
        st.subheader("B. 速率计算 (选择两行)")
        
        r_c1, r_c2, r_c3 = st.columns(3)
        # 使用索引选择行，模拟原来的“点击两行”操作
        idx_start = r_c1.number_input("起始行索引 (Index 1)", min_value=0, max_value=len(df_show)-1, value=0)
        idx_end = r_c2.number_input("结束行索引 (Index 2)", min_value=0, max_value=len(df_show)-1, value=min(5, len(df_show)-1))
        delta_hour = r_c3.number_input("生长间隔时间 (小时)", value=24.0)
        
        if st.button("计算区间平均速率"):
            row1 = df_show.iloc[idx_start]
            row2 = df_show.iloc[idx_end]
            
            da = row2["a (mm)"] - row1["a (mm)"]
            dm = row2["晶体质量 (g)"] - row1["晶体质量 (g)"]
            dt_day = delta_hour / 24.0
            
            v_a = da / dt_day
            v_m = dm / dt_day
            
            st.info(f"计算结果: a方向生长速度 = **{v_a:.2f} mm/天**, 质量生长速度 = **{v_m:.2f} g/天**")
            
            # 存入 Session 供 5天方案使用
            st.session_state['calc_va'] = v_a
            st.session_state['calc_row_end'] = row2 # 以结束行作为温控起点

    st.markdown("---")
    st.subheader("C. 生成 5 天温控方案")
    
    if 'calc_va' in st.session_state:
        # 恢复手动修改速度的功能
        col_p1, col_p2, col_p3 = st.columns(3)
        v_a_input = col_p1.number_input("a方向速度 (mm/天)", value=float(st.session_state['calc_va']), format="%.2f")
        T_now_input = col_p2.number_input("当前溶液温度 T_now (°C)", value=40.0)
        offset_input = col_p3.number_input("过冷度偏移量 ΔT (°C)", value=0.0)
        
        if st.button("生成 5 天降温表"):
            start_row = st.session_state['calc_row_end']
            a_start = start_row["a (mm)"]
            
            # 为了计算准确，这里需要重新反推当时的溶液状态
            # 简化逻辑：利用 M0 和 T0 计算初始总溶质，减去生长掉的质量
            C_init = get_wt_percent_from_T(T0)
            Solute_init = M0 * (C_init / 100.0)
            
            plan_rows = []
            T_sat_current_real = start_row["饱和点 Tsat (°C)"]
            delta_hold = T_sat_current_real - T_now_input
            
            st.write(f"基准状态: a={a_start}mm, Tsat={T_sat_current_real}°C, 当前过冷度={delta_hold:.2f}°C")
            
            for day in range(1, 6):
                a_new = a_start + v_a_input * day
                m_new, _ = crystal_masses_g(a_new, rho_val, ha_ratio, phi_val)
                
                sol_m_new = M0 - m_new
                sol_s_new = Solute_init - m_new
                
                # 计算新饱和点
                C_new = 100.0 * (sol_s_new / sol_m_new)
                Tsat_new = get_T_from_wt_percent(C_new)
                
                # 计算目标温度
                target_delta = delta_hold + offset_input
                T_end = Tsat_new - target_delta
                
                plan_rows.append({
                    "天次 (Day)": day,
                    "尺寸 a (mm)": f"{a_new:.2f}",
                    "晶体质量 (g)": f"{m_new:.2f}",
                    "溶液质量 (g)": f"{sol_m_new:.2f}",
                    "浓度 C (g/100g)": f"{C_new:.2f}",
                    "饱和点 Tsat (°C)": f"{Tsat_new:.2f}",
                    "目标温度 T_end (°C)": f"{T_end:.2f}",
                    "目标过冷度 ΔT": f"{target_delta:.2f}"
                })
                
            st.table(pd.DataFrame(plan_rows))
    else:
        st.info("请先在上方 B 部分计算速率，才能生成 5 天方案。")
