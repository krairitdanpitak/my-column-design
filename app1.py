import streamlit as st
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import numpy as np
import io
import base64
import streamlit.components.v1 as components

# ==========================================
# 1. SETUP & CSS
# ==========================================
st.set_page_config(page_title="RC Column Design (Auto)", layout="wide")

st.markdown("""
<style>
    /* ปุ่มพิมพ์ */
    .print-btn-internal {
        background-color: #008CBA; border: none; color: white !important;
        padding: 12px 28px; text-align: center; text-decoration: none;
        display: inline-block; font-size: 16px; margin: 10px 0px;
        cursor: pointer; border-radius: 5px; font-family: 'Sarabun', sans-serif;
        font-weight: bold; box-shadow: 0 2px 5px rgba(0,0,0,0.2);
    }
    .print-btn-internal:hover { background-color: #005f7f; }

    /* ตาราง */
    .report-table {width: 100%; border-collapse: collapse; font-family: sans-serif; font-size: 14px;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px;}
    .report-table th {background-color: #f2f2f2; text-align: center; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 15px;}
    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .load-value {color: #D32F2F !important; font-weight: bold;}

    /* รูปภาพ */
    .drawing-container { display: flex; justify-content: center; gap: 20px; flex-wrap: wrap; }
    .drawing-box { border: 1px solid #ddd; padding: 10px; background-color: #fff; text-align: center; min-width: 300px; }
</style>
""", unsafe_allow_html=True)

# ==========================================
# 2. DATABASE & HELPER
# ==========================================
BAR_INFO = {
    'RB6': {'A_cm2': 0.283, 'd_mm': 6},
    'RB9': {'A_cm2': 0.636, 'd_mm': 9},
    'DB10': {'A_cm2': 0.785, 'd_mm': 10},
    'DB12': {'A_cm2': 1.131, 'd_mm': 12},
    'DB16': {'A_cm2': 2.011, 'd_mm': 16},
    'DB20': {'A_cm2': 3.142, 'd_mm': 20},
    'DB25': {'A_cm2': 4.909, 'd_mm': 25},
    'DB28': {'A_cm2': 6.158, 'd_mm': 28},
    'DB32': {'A_cm2': 8.042, 'd_mm': 32}
}


def fmt(n, digits=3):
    try:
        val = float(n)
        if math.isnan(val): return "-"
        return f"{val:,.{digits}f}"
    except:
        return "-"


def beta1FromFc(fc_MPa):
    if fc_MPa <= 28: return 0.85
    b1 = 0.85 - 0.05 * ((fc_MPa - 28) / 7)
    return max(0.65, b1)


# ==========================================
# 3. CALCULATION LOGIC (ACI 318-19)
# ==========================================
def calculate_interaction_curve(b, h, cover, main_db, nx, ny, fc, fy):
    """Generate P-M Interaction Diagram Points"""
    d_prime = cover + 10 + main_db / 2
    d = h - d_prime

    total_bars = 2 * nx + 2 * max(0, ny - 2)
    Ast = total_bars * (math.pi * (main_db / 2) ** 2)
    As_face = Ast / 2.0  # Simplified 2-layer model for curve generation

    points = []

    # Po (Pure Compression)
    Ag = b * h
    Po = 0.85 * fc * (Ag - Ast) + fy * Ast
    Pn_max = 0.80 * Po  # Tied Column

    # Iterate Neutral Axis (c)
    c_values = np.linspace(1.5 * h, 0.1 * h, 30)

    for c in c_values:
        eps_cu = 0.003
        beta1 = beta1FromFc(fc)
        a = beta1 * c

        Cc = 0.85 * fc * b * min(a, h)

        # Steel Layer 1 (Compression side)
        eps_s1 = eps_cu * (c - d_prime) / c
        fs1 = min(fy, 200000 * eps_s1);
        fs1 = max(-fy, fs1)
        Fs1 = As_face * fs1

        # Steel Layer 2 (Tension side)
        eps_s2 = eps_cu * (c - d) / c
        fs2 = min(fy, 200000 * eps_s2);
        fs2 = max(-fy, fs2)
        Fs2 = As_face * fs2

        Pn = Cc + Fs1 + Fs2

        # Moment about Plastic Centroid
        Mc = Cc * (h / 2 - a / 2)
        Ms1 = Fs1 * (h / 2 - d_prime)
        Ms2 = -Fs2 * (d - h / 2)
        Mn = Mc + Ms1 + Ms2

        # Phi Factor
        eps_t = abs(eps_cu * (d - c) / c)
        if eps_t <= 0.002:
            phi = 0.65
        elif eps_t >= 0.005:
            phi = 0.90
        else:
            phi = 0.65 + (eps_t - 0.002) * (250 / 3)

        # Cap at Pn_max
        phiPn = phi * Pn
        limit_top = 0.65 * Pn_max
        if phiPn > limit_top: phiPn = limit_top

        points.append({'P': phiPn, 'M': phi * Mn})

    points.append({'P': 0, 'M': points[-1]['M']})  # Pure Bending

    return points, Ag, Ast, limit_top


def check_capacity(curve_points, Pu_target, Mu_target):
    """Check if (Mu, Pu) is inside the curve"""
    # 1. Check Axial Limit
    max_P = curve_points[0]['P']
    if Pu_target > max_P: return False

    # 2. Check Moment at Pu level
    # Interpolate M_cap at Pu_target
    m_cap = 0
    found = False
    for i in range(len(curve_points) - 1):
        p1 = curve_points[i]['P']
        p2 = curve_points[i + 1]['P']
        # Curve goes from High P to Low P
        if p2 <= Pu_target <= p1:
            ratio = (Pu_target - p2) / (p1 - p2 + 1e-9)
            m1 = curve_points[i]['M']
            m2 = curve_points[i + 1]['M']
            m_cap = m2 + ratio * (m1 - m2)
            found = True
            break

    if not found:
        # If Pu is negative (tension) or very low, assume pure bending cap
        if Pu_target <= 0:
            m_cap = curve_points[-1]['M']
        else:
            return False

    return Mu_target <= m_cap


def auto_design_reinforcement(inputs):
    """Loop to find min bars that satisfy loads"""
    b = inputs['b'] * 10;
    h = inputs['h'] * 10
    cover = inputs['cover'] * 10
    fc = inputs['fc'] * 0.0980665;
    fy = inputs['fy'] * 0.0980665
    main_key = inputs['mainBar']
    db_main = BAR_INFO[main_key]['d_mm']

    Pu_N = inputs['Pu'] * 9806.65
    Mu_Nmm = inputs['Mu'] * 9806650.0

    best_nx, best_ny = 0, 0
    min_steel_area = 999999
    found_solution = False

    # Loop ranges (e.g., 2 to 10 bars per face)
    for nx in range(2, 10):
        for ny in range(2, 10):
            total_bars = 2 * nx + 2 * max(0, ny - 2)
            Ast = total_bars * (math.pi * (db_main / 2) ** 2)
            rho = Ast / (b * h)

            # Constraints: 1% <= rho <= 8%
            if rho > 0.08: continue

            # Calculate Curve
            curve, _, _, _ = calculate_interaction_curve(b, h, cover, db_main, nx, ny, fc, fy)

            # Check Capacity
            if check_capacity(curve, Pu_N, Mu_Nmm):
                if Ast < min_steel_area:
                    min_steel_area = Ast
                    best_nx = nx
                    best_ny = ny
                    found_solution = True

    return found_solution, best_nx, best_ny


def process_column_calculation(inputs):
    rows = []

    def sec(title):
        rows.append(["SECTION", title, "", "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        rows.append([item, formula, subs, result, unit, status])

    # Inputs
    b = inputs['b'] * 10;
    h = inputs['h'] * 10
    cover = inputs['cover'] * 10
    fc = inputs['fc'] * 0.0980665
    fy = inputs['fy'] * 0.0980665
    main_key = inputs['mainBar']
    nx = int(inputs['nx'])
    ny = int(inputs['ny'])
    db_main = BAR_INFO[main_key]['d_mm']

    Pu_tf = inputs['Pu']
    Mu_tfm = inputs['Mu']

    # 1. Properties
    sec("1. PROPERTIES")
    row("Materials", "fc', fy", f"{fmt(fc, 2)}, {fmt(fy, 0)}", "-", "MPa")
    Ag = b * h
    row("Gross Area", "Ag = b·h", f"{b:.0f}·{h:.0f}", f"{Ag:.0f}", "mm²")

    # 2. Reinforcement
    sec("2. REINFORCEMENT")
    total_bars = 2 * nx + 2 * max(0, ny - 2)
    bar_area = BAR_INFO[main_key]['A_cm2'] * 100
    Ast = total_bars * bar_area
    rho_g = Ast / Ag

    st_rho = "OK" if 0.01 <= rho_g <= 0.08 else "FAIL"
    row("Provided", f"{total_bars}-{main_key}", f"{total_bars} x {bar_area:.0f}", f"{Ast:.0f}", "mm²")
    row("Ratio ρg", "Ast / Ag", f"{Ast:.0f} / {Ag:.0f}", f"{rho_g * 100:.2f}", "%", st_rho)

    # 3. Capacity
    sec("3. CAPACITY CHECK")
    curve_points, _, _, limit_top = calculate_interaction_curve(b, h, cover, db_main, nx, ny, fc, fy)

    Pu_N = Pu_tf * 9806.65
    Mu_Nmm = Mu_tfm * 9806650.0

    # Axial Check
    phiPn_max_tf = (limit_top) / 9806.65
    st_ax = "PASS" if Pu_N <= limit_top else "FAIL"
    row("Load Pu", "-", "-", f"{fmt(Pu_tf, 2)}", "tf", "")
    row("Max Axial", "φPn,max", f"0.65·0.80·Po", f"{fmt(phiPn_max_tf, 2)}", "tf", st_ax)

    # Interaction Check
    # Interpolate M capacity at Pu
    m_cap_Nmm = 0
    found = False
    for i in range(len(curve_points) - 1):
        p1 = curve_points[i]['P'];
        p2 = curve_points[i + 1]['P']
        if p2 <= Pu_N <= p1:
            r = (Pu_N - p2) / (p1 - p2 + 1e-9)
            m_cap_Nmm = curve_points[i + 1]['M'] + r * (curve_points[i]['M'] - curve_points[i + 1]['M'])
            found = True;
            break
    if not found and Pu_N <= curve_points[0]['P']: m_cap_Nmm = curve_points[-1]['M']  # Low load

    m_cap_tfm = m_cap_Nmm / 9806650.0
    st_pm = "PASS" if Mu_tfm <= m_cap_tfm and st_ax == "PASS" else "FAIL"

    row("Load Mu", "-", "-", f"{fmt(Mu_tfm, 2)}", "tf-m", "")
    row("Capacity φMn", "at Pu", f"Interpolated", f"{fmt(m_cap_tfm, 2)}", "tf-m", st_pm)

    sec("4. TIES")
    tie_key = inputs['tieBar']
    db_tie = BAR_INFO[tie_key]['d_mm']
    s1 = 16 * db_main;
    s2 = 48 * db_tie;
    s3 = min(b, h)
    s_req = min(s1, s2, s3)
    s_prov = math.floor(s_req / 25) * 25
    row("Spacing", "min(16db, 48dt, dim)", f"min({s1:.0f},{s2:.0f},{s3:.0f})", f"@{s_prov / 10:.0f} cm", "-", "OK")

    sec("5. CONCLUSION")
    overall = "OK" if st_rho == "OK" and st_pm == "PASS" else "NOT OK"
    row("Status", "-", "-", overall, "-", overall)

    return rows, curve_points, total_bars, s_prov


# ==========================================
# 4. PLOTTING
# ==========================================
def fig_to_base64(fig):
    buf = io.BytesIO();
    fig.savefig(buf, format='png', bbox_inches='tight');
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def plot_sect(b, h, cover, main_db, tie_db, nx, ny, tie_s):
    fig, ax = plt.subplots(figsize=(4, 4))
    rect = patches.Rectangle((0, 0), b, h, lw=2, ec='#333', fc='#eee')
    ax.add_patch(rect)
    margin = cover + tie_db / 2
    rect_tie = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, lw=2, ec='#1976D2', fc='none')
    ax.add_patch(rect_tie)

    start_x = margin + main_db / 2;
    end_x = b - margin - main_db / 2
    xs = np.linspace(start_x, end_x, nx) if nx > 1 else [b / 2]
    start_y = margin + main_db / 2;
    end_y = h - margin - main_db / 2
    ys = np.linspace(start_y, end_y, ny) if ny > 1 else [h / 2]

    for x in xs:
        ax.add_patch(patches.Circle((x, end_y), main_db / 2, ec='k', fc='#D32F2F'))
        ax.add_patch(patches.Circle((x, start_y), main_db / 2, ec='k', fc='#D32F2F'))
    if ny > 2:
        for y in ys[1:-1]:
            ax.add_patch(patches.Circle((start_x, y), main_db / 2, ec='k', fc='#D32F2F'))
            ax.add_patch(patches.Circle((end_x, y), main_db / 2, ec='k', fc='#D32F2F'))

    ax.set_xlim(-50, b + 50);
    ax.set_ylim(-50, h + 50);
    ax.axis('off');
    ax.set_aspect('equal')
    ax.set_title(f"{2 * nx + 2 * max(0, ny - 2)}-DB{main_db:.0f}", fontweight='bold')
    return fig


def plot_pm(curve, Pu, Mu):
    fig, ax = plt.subplots(figsize=(5, 5))
    Ms = [p['M'] / 9806650.0 for p in curve]
    Ps = [p['P'] / 9806.65 for p in curve]
    ax.plot(Ms, Ps, 'b-', lw=2, label='Capacity')
    ax.plot([0, Ms[-1]], [0, Ps[-1]], 'b--')
    ax.plot([0, 0], [0, Ps[0]], 'b-')
    ax.plot(Mu, Pu, 'ro', ms=8, label='Load')
    ax.set_xlabel('Moment (tf-m)');
    ax.set_ylabel('Axial (tf)')
    ax.grid(True, ls='--', alpha=0.6);
    ax.legend()
    return fig


# ==========================================
# 5. UI & REPORT
# ==========================================
st.title("RC Column Design (Auto)")

with st.sidebar.form("inputs"):
    st.header("Project Info")
    project = st.text_input("Project", "Office Building")
    col_id = st.text_input("Column ID", "C-01")
    engineer = st.text_input("Engineer", "Mr. Engineer")

    st.header("1. Properties")
    c1, c2 = st.columns(2)
    fc = c1.number_input("fc' (ksc)", 240);
    fy = c2.number_input("fy (ksc)", 4000)
    b = c1.number_input("b (cm)", 25);
    h = c2.number_input("h (cm)", 25)
    cover = st.number_input("Cover (cm)", 3.0)

    st.header("2. Reinforcement")
    # Auto Selection Mode
    design_mode = st.radio("Design Mode", ["Manual", "Auto-Design"])

    mainBar = st.selectbox("Main Bar", list(BAR_INFO.keys()), index=4)
    tieBar = st.selectbox("Tie Bar", ['RB6', 'RB9', 'DB10'], index=0)

    if design_mode == "Manual":
        nx = st.number_input("Nx (bars along X)", 2)
        ny = st.number_input("Ny (bars along Y)", 2)
    else:
        st.info(f"Program will find optimal Nx, Ny for {mainBar}")
        nx = 2  # Placeholder
        ny = 2  # Placeholder

    st.header("3. Loads")
    Pu = st.number_input("Pu (tf)", 40.0)
    Mu = st.number_input("Mu (tf-m)", 2.0)

    run_btn = st.form_submit_button("Run Design")

if run_btn:
    # Prepare inputs dictionary
    data = {
        'project': project, 'col_id': col_id, 'engineer': engineer,
        'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cover': cover,
        'mainBar': mainBar, 'tieBar': tieBar,
        'Pu': Pu, 'Mu': Mu
    }

    if design_mode == "Auto-Design":
        found, best_nx, best_ny = auto_design_reinforcement(data)
        if found:
            data['nx'] = best_nx
            data['ny'] = best_ny
            st.success(
                f"✅ Auto-Design Found: Use {2 * best_nx + 2 * max(0, best_ny - 2)}-{mainBar} (Nx={best_nx}, Ny={best_ny})")
        else:
            st.error("❌ Auto-Design Failed: Section too small or loads too high for this bar size.")
            st.stop()
    else:
        data['nx'] = nx
        data['ny'] = ny

    # Run Calculation
    rows, curve, total_bars, s_prov = process_column_calculation(data)

    # Plotting
    db_m = BAR_INFO[mainBar]['d_mm'];
    db_t = BAR_INFO[tieBar]['d_mm']
    img1 = fig_to_base64(plot_sect(b * 10, h * 10, cover * 10, db_m, db_t, int(data['nx']), int(data['ny']), s_prov))
    img2 = fig_to_base64(plot_pm(curve, Pu, Mu))

    # HTML Report
    t_rows = "".join([
                         f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td class='load-value'>{r[3]}</td><td>{r[4]}</td><td class='{('pass-ok' if 'PASS' in r[5] or 'OK' in r[5] else 'pass-no')}'>{r[5]}</td></tr>" if
                         r[0] != "SECTION" else f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>" for r in rows])

    html = f"""
    <div style="font-family: Sarabun, sans-serif; padding: 20px;">
        <div style="text-align:center; border-bottom: 2px solid #333; margin-bottom: 20px;">
            <div style="float:right; border:2px solid #333; padding:5px 10px; font-weight:bold;">{col_id}</div>
            <h2>ENGINEERING DESIGN REPORT</h2>
            <h4>RC Column Design (ACI 318-19)</h4>
        </div>
        <div style="display:flex; justify-content:space-between; margin-bottom:20px;">
            <div style="border:1px solid #ddd; padding:10px; width:48%;">
                <strong>Project:</strong> {project}<br><strong>Engineer:</strong> {engineer}
            </div>
            <div style="border:1px solid #ddd; padding:10px; width:48%;">
                <strong>Section:</strong> {b}x{h} cm<br><strong>Loads:</strong> Pu={Pu}tf, Mu={Mu}tf-m
            </div>
        </div>
        <div class="drawing-container">
            <div class="drawing-box"><img src="{img1}" style="max-width:100%;"></div>
            <div class="drawing-box"><img src="{img2}" style="max-width:100%;"></div>
        </div>
        <br>
        <table class="report-table">
            <thead><tr><th width="20%">Item</th><th width="25%">Formula</th><th width="30%">Substitution</th><th>Result</th><th>Unit</th><th>Status</th></tr></thead>
            <tbody>{t_rows}</tbody>
        </table>
        <div style="margin-top:40px; text-align:center;">
            <div style="display:inline-block; width:250px; text-align:left;">
                <strong>Designed by:</strong><br><br><div style="border-bottom:1px solid #000;"></div>
                <div style="text-align:center; margin-top:5px;">({engineer})<br>วิศวกรโครงสร้าง</div>
            </div>
        </div>
    </div>
    """
    st.components.v1.html(html, height=1200, scrolling=True)
else:
    st.info("👈 กรุณากรอกข้อมูลและกด Run Design")