import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import math

# ==========================================
# 0. 预设数据 (写入代码)
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
# 物理公式
# ==========================================
RHO_KDP_SOLID = 2.33
RHO_KDP_SOLID_TAB4 = 2.34
DEFAULT_PHI = 47.0
C_OVER_A = 0.936
RHO_A = 0.96286258
RHO_B = 0.00888087
SAT_SOLUBILITY_REF = {20:22.6, 25:26.3, 30:29.8, 35:31.7, 40:33.5, 45:37.3, 50:41.0, 55:45.6, 60:50.1}
WATER_DENSITY_REF = {20:0.9982, 25:0.9970, 30:0.9957, 35:0.9940, 40:0.9922, 45:0.9903, 50:0.9880, 55:0.9857, 60:0.9832}

def get_wt_percent_from_T(T): return 10.68 + 0.3616 * T
def get_T_from_wt_percent(wt): return (wt - 10.68) / 0.3616
def get_rho_sat(wt_percent): return RHO_A + RHO_B * wt_percent
def wt_percent_from_solubility(S): return 100.0 * S / (100.0 + S)
def pyramid_height_from_phi(a, phi_deg): return (a / 2.0) * math.tan(math.radians(phi_deg))
def crystal_masses_g(a_mm, density_g_cm3, h_over_a, phi_deg):
    a_cm = a_mm / 10.0
    h_prism_cm = a_cm * h_over_a
    h_pyr_cm = (a_cm / 2.0) * math.tan(math.radians(phi_deg))
    V_total = (a_cm ** 2) * h_prism_cm + (1.0/3.0) * (a_cm ** 2) * h_pyr_cm
    return V_total * density_g_cm3, V_total

# ==========================================
# 页面 UI 构建
# ==========================================
st.set_page_config(page_title="KDP 计算平台", layout="wide", page_icon="💎")

# --- [关键修改] 注入 CSS 样式：强制缩小间距和字体 ---
st.markdown("""
    <style>
        /* 缩小顶部空白 */
        .block-container {
            padding-top: 1.5rem !important;
            padding-bottom: 1rem !important;
            padding-left: 2rem !important;
            padding-right: 2rem !important;
        }
        /* 缩小标题字体 */
        h1 { font-size: 1.6rem !important; margin-bottom: 0.5rem !important; }
        h2 { font-size: 1.3rem !important; margin-top: 0.5rem !important; margin-bottom: 0.5rem !important; }
        h3 { font-size: 1.1rem !important; margin-bottom: 0.2rem !important; }
        /* 缩小输入框的 Label 字体 */
        .stNumberInput label, .stSelectbox label, .stTextInput label {
            font-size: 0.85rem !important;
        }
        /* 调整 Tab 的内边距 */
        .stTabs [data-baseweb="tab"] {
            padding-top: 4px;
            padding-bottom: 4px;
        }
        /* 让 Plotly 图表贴紧一点 */
        .js-plotly-plot {
            margin-top: -20px;
        }
    </style>
""", unsafe_allow_html=True)

st.title("💎 KDP 综合计算平台 v4.3 (紧凑版)")

tab1, tab2, tab3, tab4 = st.tabs(["① 晶体3D", "② 饱和换算", "③ 配液Pro", "④ 生长控制"])

# ==========================================
# Tab 1: 晶体质量 / 3D
# ==========================================
with tab1:
    col1_L, col1_R = st.columns([1, 1.2]) # 左侧输入，右侧画图
    
    with col1_L:
        st.subheader("观测与标定")
        # [修改] 使用 3 列布局，让输入框变小且并排
        c1, c2, c3 = st.columns(3)
        with c1: obs_a = st.number_input("a 观测 (mm)", value=26.5)
        with c2: obs_b = st.number_input("b 观测 (mm)", value=26.5)
        with c3: obs_h = st.number_input("h 观测 (mm)", value=35.0)
        
        # [修改] 使用 expander 折叠不常用的标定参数，节省空间
        with st.expander("🛠️ 标定参数设置 (点击展开)", expanded=False):
            user_select = st.selectbox("选择预设", list(CALIB_PRESETS.keys()), index=1)
            vals = CALIB_PRESETS[user_select]
            
            st.caption("公式: obs = k · real + b")
            ck1, ck2, ck3, ck4 = st.columns(4)
            k_ab = ck1.number_input("k (ab)", value=vals["k_ab"], format="%.5f")
            b_ab = ck2.number_input("b (ab)", value=vals["b_ab"], format="%.5f")
            k_h = ck3.number_input("k (h)", value=vals["k_h"], format="%.5f")
            b_h = ck4.number_input("b (h)", value=vals["b_h"], format="%.5f")
            
            h_mode = st.radio("h 含义", ["总高度 Htot", "柱体高度 Hp"], horizontal=True)

        if st.button("计算并建模", type="primary", use_container_width=True):
            try:
                real_a = (obs_a - b_ab) / k_ab
                real_b = (obs_b - b_ab) / k_ab
                real_h = (obs_h - b_h) / k_h
                
                a_eq = (real_a + real_b) / 2.0
                Hy = 0.5 * a_eq * C_OVER_A
                if h_mode == "总高度 Htot":
                    Hp = real_h - Hy
                else:
                    Hp = real_h
                
                # 简单防错
                if Hp < 0: Hp = 0.1 
                
                a_cm = a_eq / 10.0
                V_total = (a_cm**2 * (Hp/10.0)) + (1.0/3.0 * a_cm**2 * (Hy/10.0))
                mass = V_total * RHO_KDP_SOLID
                
                st.success(f"校正后质量: **{mass:.2f} g**")
                
                # 存入 session_state 供画图使用
                st.session_state['t1_res'] = (obs_a, obs_b, obs_h, h_mode, real_a, real_b, Hp, Hy, mass)
                
            except Exception as e:
                st.error(f"计算错误: {e}")

    with col1_R:
        # 只有计算过才画图
        if 't1_res' in st.session_state:
            oa, ob, oh, mode, ra, rb, Hp, Hy, m = st.session_state['t1_res']
            
            # 这里为了省代码空间，稍微简化画图逻辑
            def get_mesh(a, b, hp, hy, color, opac):
                dx, dy = a/2, b/2
                # Prism
                xp = [-dx, dx, dx, -dx, -dx, dx, dx, -dx]
                yp = [-dy, -dy, dy, dy, -dy, -dy, dy, dy]
                zp = [0, 0, 0, 0, hp, hp, hp, hp]
                i = [7, 0, 0, 0, 4, 4, 6, 6, 4, 0, 3, 2]
                j = [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3]
                k = [0, 7, 2, 3, 6, 7, 1, 1, 5, 5, 7, 6]
                prism = go.Mesh3d(x=xp, y=yp, z=zp, i=i, j=j, k=k, color=color, opacity=opac, name='Prism')
                # Cap
                xc = [-dx, dx, dx, -dx, 0]
                yc = [-dy, -dy, dy, dy, 0]
                zc = [hp, hp, hp, hp, hp+hy]
                ic, jc, kc = [0, 1, 2, 3], [1, 2, 3, 0], [4, 4, 4, 4]
                cap = go.Mesh3d(x=xc, y=yc, z=zc, i=ic, j=jc, k=kc, color=color, opacity=opac, name='Cap')
                return [prism, cap]

            fig = go.Figure()
            # Raw (Gray) - 简单估算 raw 的几何用于显示
            raw_hp = oh - (0.5*(oa+ob)/2.0*C_OVER_A) if "总" in mode else oh
            if raw_hp < 0: raw_hp = 0
            raw_hy = 0.5*(oa+ob)/2.0*C_OVER_A
            
            for m in get_mesh(oa, ob, raw_hp, raw_hy, 'gray', 0.15): fig.add_trace(m)
            for m in get_mesh(ra, rb, Hp, Hy, '#0078D4', 0.7): fig.add_trace(m)
            
            fig.update_layout(
                scene=dict(aspectmode='data'), 
                margin=dict(l=0, r=0, b=0, t=10), # 紧凑边距
                height=400 # 限制高度
            )
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"蓝色实体: 真实尺寸 (L={ra:.1f}x{rb:.1f}, H={Hp+Hy:.1f})")
        else:
            st.info("👈 请在左侧输入数据并计算")

# ==========================================
# Tab 2: 饱和溶液
# ==========================================
with tab2:
    c1, c2 = st.columns([1, 1])
    with c1:
        st.subheader("快速计算")
        # [修改] 横向滑块
        temp_query = st.slider("温度 (°C)", 20, 60, 40, step=5)
        S = SAT_SOLUBILITY_REF.get(temp_query, 0)
        rho = get_rho_sat(wt_percent_from_solubility(S))
        
        st.markdown(f"**参考**: S={S}, ρ={rho:.4f}")
        
        # [修改] 紧凑的两列
        cc1, cc2 = st.columns(2)
        mode = cc1.selectbox("模式", ["体积→质量", "质量→体积"], label_visibility="collapsed")
        val = cc2.number_input("数值 (mL/g)", value=1000.0, label_visibility="collapsed")
        
        if st.button("换算", use_container_width=True):
            if "体积" in mode: st.success(f"质量: {val * rho:.2f} g")
            else: st.success(f"体积: {val / rho:.2f} mL")
            
    with c2:
        # 表格本身就很紧凑，不需要大改
        data = [[t, SAT_SOLUBILITY_REF[t], f"{get_rho_sat(wt_percent_from_solubility(SAT_SOLUBILITY_REF[t])):.4f}"] for t in range(20, 65, 5)]
        st.dataframe(pd.DataFrame(data, columns=["温度", "S", "密度"]), height=250, use_container_width=True)

# ==========================================
# Tab 3: 配液 Pro
# ==========================================
with tab3:
    # [修改] 布局优化
    calc_type = st.selectbox("已知条件", ["已知: 总重M & 温度T", "已知: 水W & 溶质S", "已知: 体积V & 温度T", "已知: S & T"])
    
    # 动态生成输入框，一行显示
    cols = st.columns(3)
    in_1, in_2 = 0, 0
    with cols[0]:
        if "M" in calc_type: in_1 = st.number_input("总重 M (g)", value=1200.0)
        elif "W" in calc_type: in_1 = st.number_input("水 W (g)", value=800.0)
        elif "V" in calc_type: in_1 = st.number_input("体积 V (mL)", value=1000.0)
        elif "S" in calc_type: in_1 = st.number_input("溶质 S (g)", value=400.0)
    with cols[1]:
        if "T" in calc_type: in_2 = st.number_input("温度 T (°C)", value=40.0)
        elif "S" in calc_type and "W" in calc_type: in_2 = st.number_input("溶质 S (g)", value=400.0)
    with cols[2]:
        btn_calc = st.button("计算配方", type="primary", use_container_width=True)

    if btn_calc:
        try:
            rM, rT, rW, rS, rV = 0,0,0,0,0
            # 简化逻辑展示
            if "M" in calc_type:
                wt = get_wt_percent_from_T(in_2)
                rM, rT = in_1, in_2
                rS = rM * wt/100; rW = rM - rS; rV = rM/get_rho_sat(wt)
            elif "W" in calc_type: # W & S
                rW, rS = in_1, in_2; rM = rW+rS
                wt = rS/rM*100; rT = get_T_from_wt_percent(wt); rV = rM/get_rho_sat(wt)
            elif "V" in calc_type:
                rV, rT = in_1, in_2; wt = get_wt_percent_from_T(rT)
                rM = rV * get_rho_sat(wt); rS = rM*wt/100; rW = rM-rS
            elif "S" in calc_type: # S & T
                rS, rT = in_1, in_2; wt = get_wt_percent_from_T(rT)
                rM = rS/(wt/100); rW = rM-rS; rV = rM/get_rho_sat(wt)
                
            st.session_state['t3_res'] = (rM, rT, rW, rS, rV)
            # 结果显示：一行 Metric
            cc1, cc2, cc3, cc4, cc5 = st.columns(5)
            cc1.metric("M (g)", f"{rM:.0f}")
            cc2.metric("T (°C)", f"{rT:.1f}")
            cc3.metric("S (g)", f"{rS:.0f}")
            cc4.metric("W (g)", f"{rW:.0f}")
            cc5.metric("V (mL)", f"{rV:.0f}")
            
        except: st.error("计算错")

    st.divider()
    
    if 't3_res' in st.session_state:
        mM, mT, mW, mS, mV = st.session_state['t3_res']
        st.caption(f"基于当前配方 (M={mM:.0f}g, T={mT:.1f}°C) 预测生长:")
        
        c_g1, c_g2, c_g3 = st.columns([1, 1, 1])
        t_end = c_g1.number_input("目标温度", value=20.0)
        mode = c_g2.selectbox("模式", ["点籽晶 (Mode A)", "片状 (Mode B)"])
        param = c_g3.number_input("H/L 或 边长L", value=1.0 if "A" in mode else 20.0)
        
        if st.button("预测析出量", use_container_width=True):
            wt_e = get_wt_percent_from_T(t_end)
            S_sat = mW / (1-wt_e/100) * (wt_e/100)
            dS = mS - S_sat
            if dS>0: st.success(f"析出: {dS:.2f} g")
            else: st.warning("无析出")

# ==========================================
# Tab 4: 生长控制
# ==========================================
with tab4:
    # [修改] 使用 expander 隐藏密密麻麻的参数
    with st.expander("⚙️ 初始参数设置 (M0, T0, 几何)", expanded=False):
        ec1, ec2, ec3, ec4 = st.columns(4)
        M0 = ec1.number_input("M0", 2000.0)
        T0 = ec2.number_input("T0", 55.0)
        a_min = ec3.number_input("a_min", 1.0)
        a_max = ec4.number_input("a_max", 80.0)
        ec5, ec6, ec7 = st.columns(3)
        step = ec5.number_input("step", 2.0)
        ha = ec6.number_input("h/a", 1.0)
        phi = ec7.number_input("phi", 47.0)
        
    if st.button("生成生长表", use_container_width=True):
        C0 = get_wt_percent_from_T(T0); sol0 = M0 * C0/100
        dat = []
        cur = a_min
        while cur <= a_max:
            mc, vc = crystal_masses_g(cur, RHO_KDP_SOLID_TAB4, ha, phi)
            sm = max(M0-mc, 1e-9); ss = max(sol0-mc, 0)
            ts = get_T_from_wt_percent(100*ss/sm)
            if ts<10: break
            dat.append({"a":round(cur,2), "H":round(ha*cur+pyramid_height_from_phi(cur,phi),2), "m":round(mc,2), "Ts":round(ts,2)})
            cur += step
        st.session_state['df4'] = pd.DataFrame(dat)

    if 'df4' in st.session_state:
        df = st.session_state['df4']
        st.dataframe(df, height=200, use_container_width=True)
        
        st.caption("温控方案生成:")
        xc1, xc2, xc3, xc4 = st.columns(4)
        i1 = xc1.number_input("Start Idx", 0, len(df)-1, 0)
        i2 = xc2.number_input("End Idx", 0, len(df)-1, min(5, len(df)-1))
        dh = xc3.number_input("Interval(h)", 24.0)
        off = xc4.number_input("Offset", 0.0)
        
        if st.button("生成方案"):
            r1, r2 = df.iloc[i1], df.iloc[i2]
            va = (r2['a']-r1['a'])/(dh/24)
            st.info(f"速度: {va:.2f} mm/d")
            # 简单生成5天
            plans = []
            for d in range(1,6):
                anew = r2['a'] + va*d
                mc, _ = crystal_masses_g(anew, RHO_KDP_SOLID_TAB4, ha, phi)
                C0 = get_wt_percent_from_T(T0); sol0 = M0 * C0/100
                sm = M0-mc; ss = sol0-mc
                ts = get_T_from_wt_percent(100*ss/sm)
                plans.append({"Day":d, "a":f"{anew:.2f}", "Ts":f"{ts:.2f}", "Tset":f"{ts-off:.2f}"})
            st.table(pd.DataFrame(plans))
