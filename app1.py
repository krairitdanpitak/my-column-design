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
st.set_page_config(page_title="RC Column Design SDM (ACI 318-19)", layout="wide")

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Sarabun:wght@400;700&display=swap');

    /* CSS ปุ่มพิมพ์ */
    .print-btn-internal {
        background-color: #008CBA;
        border: none;
        color: white !important;
        padding: 12px 28px;
        text-align: center;
        text-decoration: none;
        display: inline-block;
        font-size: 16px;
        margin: 10px 0px;
        cursor: pointer;
        border-radius: 5px;
        font-family: 'Sarabun', sans-serif;
        font-weight: bold;
    }
    .print-btn-internal:hover { background-color: #005f7f; }

    /* CSS ตาราง */
    .report-table {width: 100%; border-collapse: collapse; font-family: 'Sarabun', sans-serif; font-size: 14px;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px;}
    .report-table th {background-color: #f2f2f2; text-align: center; font-weight: bold;}

    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 15px;}
    .load-value {color: #D32F2F !important; font-weight: bold;}
</style>
""", unsafe_allow_html=True)

# ==========================================
# 2. DATABASE & HELPER FUNCTIONS
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


def fmt(n, digits=2):
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
# 3. SHEAR CHECK FUNCTION (NEW)
# ==========================================
def check_column_shear_aci318_19(pu_tf, vu_tf, fc_ksc, fy_tie_ksc, b_cm, h_cm, d_mm, av_mm2, s_mm):
    """
    ตรวจสอบกำลังรับแรงเฉือนของเสาตาม ACI 318-19
    คืนค่าเป็น list ของ rows เพื่อนำไปต่อกับตารางหลัก
    """
    rows = []

    # Helper to add row
    def row(i, f, sub, r, u, st=""):
        rows.append([i, f, sub, r, u, st])

    # 1. Conversions
    Nu = pu_tf * 9806.65  # Ton -> Newton
    Vu = vu_tf * 9806.65  # Ton -> Newton
    fc_mpa = fc_ksc * 0.0980665
    fyt_mpa = fy_tie_ksc * 0.0980665
    bw = b_cm * 10
    Ag = bw * (h_cm * 10)

    lambda_val = 1.0
    phi_v = 0.75

    rows.append(["SECTION", "4. SHEAR CAPACITY CHECK (ACI 318-19)", "", "", "", ""])

    # 2. Concrete Capacity (Vc)
    # Vc = 0.17 * (1 + Nu / (14 * Ag)) * lambda * sqrt(fc) * bw * d
    nu_term = Nu / (14 * Ag)
    vc_val = 0.17 * (1 + nu_term) * lambda_val * math.sqrt(fc_mpa) * bw * d_mm
    phi_vc = phi_v * vc_val

    sub_vc = f"0.17(1 + {Nu / 1000:.0f}/(14·{Ag / 100:.0f}))√{fc_mpa:.0f}·{bw}·{d_mm:.0f}"
    row("Concrete Vc", "0.17(1 + Nu/14Ag)λ√fc bw d", sub_vc, f"{vc_val / 9806.65:.2f}", "tf", "")

    # 3. Steel Capacity (Vs)
    # Vs = (Av * fyt * d) / s
    vs_val = (av_mm2 * fyt_mpa * d_mm) / s_mm
    phi_vs = phi_v * vs_val

    sub_vs = f"({av_mm2:.0f} · {fyt_mpa:.0f} · {d_mm:.0f}) / {s_mm:.0f}"
    row("Steel Vs", "Av fyt d / s", sub_vs, f"{vs_val / 9806.65:.2f}", "tf", "")

    # Check Max Vs Limit (ACI 318-19)
    vs_max = 0.66 * math.sqrt(fc_mpa) * bw * d_mm
    if vs_val > vs_max:
        row("Vs Max Limit", "0.66√fc bw d", f"Vs {vs_val / 9806.65:.2f} > Max {vs_max / 9806.65:.2f}", "FAIL", "-",
            "FAIL")

    # 4. Total Capacity
    phi_vn = phi_vc + phi_vs
    status = "PASS" if phi_vn >= Vu else "FAIL"

    row("Total Capacity φVn", "φ(Vc + Vs)", f"{phi_vc / 9806.65:.2f} + {phi_vs / 9806.65:.2f}",
        f"{phi_vn / 9806.65:.2f}", "tf", status)
    row("Shear Demand Vu", "Input Load", "-", f"{vu_tf:.2f}", "tf", status)

    # 5. Spacing Check (Shear)
    limit_threshold = 0.33 * math.sqrt(fc_mpa) * bw * d_mm
    if vs_val <= limit_threshold:
        max_s_shear = min(d_mm / 2, 600)
        cond_text = "Vs ≤ 0.33√fc bw d"
    else:
        max_s_shear = min(d_mm / 4, 300)
        cond_text = "Vs > 0.33√fc bw d"

    status_s = "OK" if s_mm <= max_s_shear else "FAIL"
    row("Shear Spacing Check", f"Limit ({cond_text})", f"Limit: {max_s_shear:.0f} mm", f"Use: {s_mm:.0f}", "mm",
        status_s)

    return rows, status


# ==========================================
# 4. CALCULATION LOGIC (MAIN)
# ==========================================
def calculate_interaction_curve(b, h, cover, main_db, nx, ny, fc, fy):
    """สร้างจุดบนกราฟ P-M Interaction Diagram"""
    d_prime = cover + 10 + main_db / 2
    d = h - d_prime

    Ast = (2 * nx + 2 * max(0, ny - 2)) * (math.pi * (main_db / 2) ** 2)
    As_face = Ast / 2.0

    points = []
    Ag = b * h
    Po = 0.85 * fc * (Ag - Ast) + fy * Ast

    c_values = np.linspace(1.5 * h, 0.1 * h, 40)

    for c in c_values:
        eps_cu = 0.003
        beta1 = beta1FromFc(fc)
        a = beta1 * c
        Cc = 0.85 * fc * b * min(a, h)

        eps_s1 = eps_cu * (c - d_prime) / c
        fs1 = max(-fy, min(fy, 200000 * eps_s1))
        Fs1 = As_face * fs1

        eps_s2 = eps_cu * (c - d) / c
        fs2 = max(-fy, min(fy, 200000 * eps_s2))
        Fs2 = As_face * fs2

        Pn = Cc + Fs1 + Fs2
        Mc = Cc * (h / 2 - a / 2)
        Ms1 = Fs1 * (h / 2 - d_prime)
        Ms2 = -Fs2 * (d - h / 2)
        Mn = Mc + Ms1 + Ms2

        eps_t = abs(eps_cu * (d - c) / c)
        if eps_t <= 0.002:
            phi = 0.65
        elif eps_t >= 0.005:
            phi = 0.90
        else:
            phi = 0.65 + (eps_t - 0.002) * (250 / 3)

        phiPn_val = phi * Pn
        limit_top = 0.65 * 0.80 * Po
        if phiPn_val > limit_top: phiPn_val = limit_top

        points.append({'P': phiPn_val, 'M': phi * Mn, 'phi': phi})

    points.append({'P': 0, 'M': points[-1]['M']})
    return points, Ag, Ast, 0.65 * 0.80 * Po


def check_capacity(curve_points, max_load, Pu_target, Mu_target):
    if Pu_target > max_load: return False
    m_cap = 0
    found = False
    if Pu_target < curve_points[-1]['P']:
        m_cap = curve_points[-1]['M']
        found = True
    else:
        for i in range(len(curve_points) - 1):
            p1 = curve_points[i]['P']
            p2 = curve_points[i + 1]['P']
            if p2 <= Pu_target <= p1:
                ratio = (Pu_target - p2) / (p1 - p2 + 1e-9)
                m1 = curve_points[i]['M']
                m2 = curve_points[i + 1]['M']
                m_cap = m2 + ratio * (m1 - m2)
                found = True
                break
    if not found: return False
    return Mu_target <= m_cap


def auto_design_reinforcement(inputs):
    b = inputs['b'] * 10;
    h = inputs['h'] * 10
    cover = inputs['cover'] * 10
    fc = inputs['fc'] * 0.0980665;
    fy = inputs['fy'] * 0.0980665
    main_key = inputs['mainBar']
    db_main = BAR_INFO[main_key]['d_mm']
    Pu_N = inputs['Pu'] * 9806.65
    Mu_Nmm = inputs['Mu'] * 9806650.0

    valid_designs = []
    for nx in range(2, 9):
        for ny in range(2, 9):
            total_bars = 2 * nx + 2 * max(0, ny - 2)
            Ast = total_bars * (math.pi * (db_main / 2) ** 2)
            Ag = b * h
            rho = Ast / Ag
            if not (0.01 <= rho <= 0.08): continue

            curve, _, _, p_max = calculate_interaction_curve(b, h, cover, db_main, nx, ny, fc, fy)
            if check_capacity(curve, p_max, Pu_N, Mu_Nmm):
                valid_designs.append({'nx': nx, 'ny': ny, 'ast': Ast})

    if not valid_designs: return False, 2, 2
    valid_designs.sort(key=lambda x: x['ast'])
    return True, valid_designs[0]['nx'], valid_designs[0]['ny']


def process_column_calculation(inputs):
    rows = []

    def sec(title):
        rows.append(["SECTION", title, "", "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        rows.append([item, formula, subs, result, unit, status])

    # Inputs Conversion
    b = inputs['b'] * 10;
    h = inputs['h'] * 10  # mm
    cover = inputs['cover'] * 10
    fc_mpa = inputs['fc'] * 0.0980665
    fy_mpa = inputs['fy'] * 0.0980665

    main_key = inputs['mainBar'];
    tie_key = inputs['tieBar']
    nx = int(inputs['nx']);
    ny = int(inputs['ny'])
    Pu_tf = inputs['Pu'];
    Mu_tfm = inputs['Mu'];
    Vu_tf = inputs['Vu']

    # --- 1. MATERIAL & GEOMETRY ---
    sec("1. MATERIAL & SECTION PROPERTIES")
    row("Concrete Strength", "fc' (Input)", f"{inputs['fc']:.0f} ksc", f"{fmt(fc_mpa, 2)}", "MPa")
    row("Section Size", "b x h", f"{inputs['b']:.0f} x {inputs['h']:.0f}", f"{fmt(b, 0)}x{fmt(h, 0)}", "mm")

    Ag = b * h
    total_bars = 2 * nx + 2 * max(0, ny - 2)
    bar_area_one = BAR_INFO[main_key]['A_cm2'] * 100
    Ast = total_bars * bar_area_one
    rho_g = Ast / Ag

    row("Main Reinforcement", f"Total {total_bars}-{main_key}", f"{total_bars} x {fmt(bar_area_one, 2)}",
        f"{fmt(Ast, 0)}", "mm²")
    status_rho = "OK" if 0.01 <= rho_g <= 0.08 else "FAIL"
    row("Reinforcement Ratio", "ρg = Ast / Ag", f"{fmt(Ast, 0)} / {fmt(Ag, 0)}", f"{fmt(rho_g * 100, 2)}", "%",
        status_rho)

    # --- 2. AXIAL CAPACITY ---
    sec("2. AXIAL LOAD CAPACITY")
    term1 = 0.85 * fc_mpa * (Ag - Ast);
    term2 = fy_mpa * Ast
    Po_N = term1 + term2
    phiPn_max_N = 0.65 * 0.80 * Po_N
    phiPn_max_tf = phiPn_max_N / 9806.65

    row("Max Design Axial", "φPn,max = 0.52·Po", "-", f"{fmt(phiPn_max_tf, 2)}", "tf")
    status_axial = "PASS" if Pu_tf <= phiPn_max_tf else "FAIL"
    row("Axial Check", "Pu ≤ φPn,max", f"{fmt(Pu_tf, 2)} ≤ {fmt(phiPn_max_tf, 2)}", status_axial, "-", status_axial)

    # --- 3. TIE DESIGN ---
    sec("3. TIE (STIRRUP) DETAILING")
    db_main = BAR_INFO[main_key]['d_mm']
    db_tie = BAR_INFO[tie_key]['d_mm']

    s1 = 16 * db_main;
    s2 = 48 * db_tie;
    s3 = min(b, h)
    s_req = min(s1, s2, s3)

    row("Spacing Limit 1", "16 · db(main)", f"16 · {db_main}", f"{s1:.0f}", "mm")
    row("Spacing Limit 3", "Least Dimension", f"min({b},{h})", f"{s3:.0f}", "mm")

    s_prov = math.floor(s_req / 25.0) * 25.0
    if s_prov < 50: s_prov = 50
    row("Provide Ties", f"Use {tie_key}", f"min limits", f"@{s_prov / 10:.0f} cm", "-", "OK")

    # --- 4. SHEAR CAPACITY CHECK (NEW) ---
    # Prepare data for shear check
    d = h - cover - db_main / 2 - db_tie  # Effective depth
    tie_area_mm2 = 2 * BAR_INFO[tie_key]['A_cm2'] * 100  # 2 Legs

    shear_rows, status_shear = check_column_shear_aci318_19(
        Pu_tf, Vu_tf, inputs['fc'], inputs['fyt'], inputs['b'], inputs['h'], d, tie_area_mm2, s_prov
    )
    rows.extend(shear_rows)

    # --- 5. MOMENT CAPACITY CHECK ---
    sec("5. MOMENT CAPACITY CHECK")
    curve_points, _, _, _ = calculate_interaction_curve(b, h, cover, db_main, nx, ny, fc_mpa, fy_mpa)
    Pu_N = Pu_tf * 9806.65
    m_cap_Nmm = 0;
    found = False

    # Simple search for M capacity at Pu
    if Pu_N > curve_points[0]['P']:
        m_cap_Nmm = 0
    elif Pu_N < curve_points[-1]['P']:
        m_cap_Nmm = curve_points[-1]['M']
    else:
        for i in range(len(curve_points) - 1):
            p1 = curve_points[i]['P'];
            p2 = curve_points[i + 1]['P']
            if p2 <= Pu_N <= p1:
                ratio = (Pu_N - p2) / (p1 - p2 + 1e-9)
                m1 = curve_points[i]['M'];
                m2 = curve_points[i + 1]['M']
                m_cap_Nmm = m2 + ratio * (m1 - m2)
                found = True;
                break

    m_cap_tfm = m_cap_Nmm / 9806650.0
    status_pm = "PASS" if Mu_tfm <= m_cap_tfm else "FAIL"
    row("Interaction Check", "Mu ≤ φMn", f"{fmt(Mu_tfm, 2)} ≤ {fmt(m_cap_tfm, 2)}", status_pm, "-", status_pm)

    # --- FINAL STATUS ---
    sec("6. FINAL STATUS")
    overall = "OK" if (
                status_rho == "OK" and status_axial == "PASS" and status_shear == "PASS" and status_pm == "PASS") else "NOT OK"
    row("Overall", "-", "-", "DESIGN COMPLETE", "-", overall)

    return rows, curve_points, total_bars, s_prov


# ==========================================
# 5. PLOTTING & REPORT
# ==========================================
def fig_to_base64(fig):
    buf = io.BytesIO();
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0);
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def plot_column_section(b, h, cover, main_db, tie_db, nx, ny, tie_s, title="Column Section"):
    fig, ax = plt.subplots(figsize=(4, 4))
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#eee')
    ax.add_patch(rect)
    margin = cover + tie_db / 2
    rect_tie = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, linewidth=2, edgecolor='#1976D2',
                                 facecolor='none')
    ax.add_patch(rect_tie)

    start_x = margin + main_db / 2;
    end_x = b - margin - main_db / 2
    xs = np.linspace(start_x, end_x, nx) if nx > 1 else [b / 2]
    start_y = margin + main_db / 2;
    end_y = h - margin - main_db / 2
    ys = np.linspace(start_y, end_y, ny) if ny > 1 else [h / 2]

    for x in xs:
        ax.add_patch(patches.Circle((x, end_y), radius=main_db / 2, edgecolor='black', facecolor='#D32F2F'))
        ax.add_patch(patches.Circle((x, start_y), radius=main_db / 2, edgecolor='black', facecolor='#D32F2F'))
    if ny > 2:
        for y in ys[1:-1]:
            ax.add_patch(patches.Circle((start_x, y), radius=main_db / 2, edgecolor='black', facecolor='#D32F2F'))
            ax.add_patch(patches.Circle((end_x, y), radius=main_db / 2, edgecolor='black', facecolor='#D32F2F'))

    ax.set_xlim(-50, b + 50);
    ax.set_ylim(-50, h + 50);
    ax.axis('off');
    ax.set_aspect('equal')
    info = f"Size: {b / 10:.0f}x{h / 10:.0f} cm\nMain: {2 * nx + 2 * max(0, ny - 2)}-DB{main_db:.0f}\nTies: RB{tie_db:.0f}@{tie_s / 10:.0f}cm"
    ax.text(b / 2, -h * 0.2, info, ha='center', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
    return fig


def plot_interaction_diagram(curve_points, Pu_tf, Mu_tfm):
    fig, ax = plt.subplots(figsize=(5, 5))
    Ms = [p['M'] / 9806650.0 for p in curve_points]
    Ps = [p['P'] / 9806.65 for p in curve_points]
    ax.plot(Ms, Ps, 'b-', linewidth=2, label='Capacity φPn-φMn')
    ax.plot([0, Ms[-1]], [0, Ps[-1]], 'b--');
    ax.plot([0, 0], [0, Ps[0]], 'b-')
    ax.plot(Mu_tfm, Pu_tf, 'ro', markersize=8, label='Design Load')
    ax.set_xlabel('Moment φMn (tf-m)');
    ax.set_ylabel('Axial Load φPn (tf)')
    ax.set_title('P-M Interaction Diagram', fontweight='bold');
    ax.grid(True, linestyle='--', alpha=0.6);
    ax.legend()
    return fig


def generate_column_report(inputs, rows, img_sect, img_pm):
    table_rows = ""
    for r in rows:
        if r[0] == "SECTION":
            table_rows += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            status_cls = "pass-ok" if "OK" in r[5] or "PASS" in r[5] else "pass-no"
            val_cls = "load-value" if "Load Input" in str(r[0]) else ""
            table_rows += f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td class='{val_cls}'>{r[3]}</td><td>{r[4]}</td><td class='{status_cls}'>{r[5]}</td></tr>"

    html = f"""
    <!DOCTYPE html>
    <html lang="th">
    <head>
        <meta charset="UTF-8">
        <title>Column Design Report</title>
        <link href="https://fonts.googleapis.com/css2?family=Sarabun:wght@400;700&display=swap" rel="stylesheet">
        <style>
            body {{ font-family: 'Sarabun', sans-serif; padding: 20px; color: black; }}
            h1, h3 {{ text-align: center; margin: 5px; }}
            .header {{ margin-bottom: 20px; border-bottom: 2px solid #333; padding-bottom: 10px; position: relative; }}
            .col-id-box {{ position: absolute; top: 0; right: 0; border: 2px solid #333; padding: 5px 15px; font-weight: bold; font-size: 18px; }}
            .info-container {{ display: flex; justify-content: space-between; margin-bottom: 20px; }}
            .info-box {{ width: 48%; border: 1px solid #ddd; padding: 10px; }}
            .images {{ display: flex; justify-content: space-around; margin: 20px 0; align-items: center; }}
            .images img {{ width: 40%; border: 1px solid #ddd; padding: 5px; }}
            table {{ width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 12px; }}
            th, td {{ border: 1px solid #444; padding: 6px; }}
            th {{ background-color: #eee; }}
            .sec-row {{ background-color: #ddd; font-weight: bold; }}
            .pass-ok {{ color: green; font-weight: bold; text-align: center; }}
            .pass-no {{ color: red; font-weight: bold; text-align: center; }}
            .load-value {{ color: #D32F2F !important; font-weight: bold; }}
            .footer-section {{ margin-top: 40px; page-break-inside: avoid; }}
            .signature-block {{ width: 300px; text-align: center; }}
            .sign-line {{ border-bottom: 1px solid #000; margin: 40px 0 10px 0; }}
            @media print {{ .no-print {{ display: none !important; }} body {{ padding: 0; }} }}
            .print-btn-internal {{ background-color: #4CAF50; color: white; padding: 12px 24px; border: none; border-radius: 5px; cursor: pointer; }}
        </style>
    </head>
    <body>
        <div class="no-print" style="text-align: center;">
            <button onclick="window.print()" class="print-btn-internal">🖨️ Print This Page / พิมพ์หน้านี้</button>
        </div>
        <div class="header">
            <div class="col-id-box">{inputs['col_id']}</div>
            <h1>ENGINEERING DESIGN REPORT</h1>
            <h3>RC Column Design SDM (ACI 318-19)</h3>
        </div>
        <div class="info-container">
            <div class="info-box">
                <strong>Project:</strong> {inputs['project']}<br>
                <strong>Engineer:</strong> {inputs['engineer']}<br>
                <strong>Date:</strong> 15/12/2568
            </div>
            <div class="info-box">
                <strong>Materials:</strong> fc'={inputs['fc']} ksc, fy={inputs['fy']} ksc<br>
                <strong>Section:</strong> {inputs['b']} x {inputs['h']} cm<br>
                <strong>Rebar:</strong> Main {inputs['mainBar']}, Tie {inputs['tieBar']}
            </div>
        </div>
        <h3>Design Summary</h3>
        <div class="images">
            <img src="{img_sect}" />
            <img src="{img_pm}" />
        </div>
        <br><br><br>
        <h3>Calculation Details</h3>
        <table>
            <thead><tr><th width="20%">Item</th><th width="30%">Formula</th><th width="25%">Substitution</th><th>Result</th><th>Unit</th><th>Status</th></tr></thead>
            <tbody>{table_rows}</tbody>
        </table>
        <div class="footer-section">
            <div class="signature-block">
                <div style="text-align: left; font-weight: bold;">Designed by:</div>
                <div class="sign-line"></div>
                <div>({inputs['engineer']})</div>
                <div>วิศวกรโครงสร้าง</div>
            </div>
        </div>
    </body>
    </html>
    """
    return html


# ==========================================
# 6. MAIN UI
# ==========================================
st.title("RC Column Design SDM")

if 'calc_done' not in st.session_state: st.session_state['calc_done'] = False

with st.sidebar.form("inputs"):
    st.header("Project Info")
    project = st.text_input("Project Name", "อาคารสำนักงาน 2 ชั้น")
    col_id = st.text_input("Column Number", "C-01")
    engineer = st.text_input("Engineer Name", "นายไกรฤทธิ์ ด่านพิทักษ์")

    st.header("1. Material & Geometry")
    c1, c2 = st.columns(2)
    fc = c1.number_input("fc' (ksc)", 240)
    fy = c2.number_input("fy (ksc)", 4000)
    fyt = st.number_input("fyt (Tie) (ksc)", 2400)

    c1, c2, c3 = st.columns(3)
    b = c1.number_input("b (cm)", 25)
    h = c2.number_input("h (cm)", 25)
    cover = c3.number_input("Cover (cm)", 3.0)

    st.header("2. Reinforcement")
    design_mode = st.radio("Mode", ["Manual", "Auto-Design"])
    c1, c2 = st.columns(2)
    mainBar = c1.selectbox("Main Bar", list(BAR_INFO.keys()), index=4)
    tieBar = c2.selectbox("Tie Bar", ['RB6', 'RB9', 'DB10'], index=0)

    nx, ny = 2, 2
    if design_mode == "Manual":
        st.write("Number of bars per face:")
        c1, c2 = st.columns(2)
        nx = c1.number_input("Nx (bars along X)", 2)
        ny = c2.number_input("Ny (bars along Y)", 2)

    st.header("3. Loads (Factored)")
    Pu = st.number_input("Axial Load Pu (tf)", 40.0)
    Mu = st.number_input("Moment Mu (tf-m)", 2.0)
    # Added Shear Input
    Vu = st.number_input("Shear Load Vu (tf)", 1.5, min_value=0.0)

    run_btn = st.form_submit_button("Run Design")

if run_btn:
    inputs = {
        'project': project, 'col_id': col_id, 'engineer': engineer,
        'fc': fc, 'fy': fy, 'fyt': fyt,
        'b': b, 'h': h, 'cover': cover,
        'mainBar': mainBar, 'tieBar': tieBar,
        'nx': nx, 'ny': ny,
        'Pu': Pu, 'Mu': Mu, 'Vu': Vu
    }

    if design_mode == "Auto-Design":
        found, best_nx, best_ny = auto_design_reinforcement(inputs)
        if found:
            inputs['nx'] = best_nx;
            inputs['ny'] = best_ny
            st.success(f"✅ Auto-Design Found: Use {2 * best_nx + 2 * max(0, best_ny - 2)}-{mainBar}")
        else:
            st.error("❌ Auto-Design Failed: Section too small or load too high.");
            st.stop()

    rows, curve, total_bars, s_prov = process_column_calculation(inputs)

    main_db = BAR_INFO[mainBar]['d_mm'];
    tie_db = BAR_INFO[tieBar]['d_mm']
    img_sect = fig_to_base64(
        plot_column_section(b * 10, h * 10, cover * 10, main_db, tie_db, int(inputs['nx']), int(inputs['ny']), s_prov))
    img_pm = fig_to_base64(plot_interaction_diagram(curve, Pu, Mu))

    html_report = generate_column_report(inputs, rows, img_sect, img_pm)
    st.components.v1.html(html_report, height=1200, scrolling=True)

else:
    st.info("👈 กรุณากรอกข้อมูลและกด Run Design")
