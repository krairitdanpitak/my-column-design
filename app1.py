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
st.set_page_config(page_title="RC Column Design SDM", layout="wide")

st.markdown("""
<style>
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
        box-shadow: 0 2px 5px rgba(0,0,0,0.2);
    }
    .print-btn-internal:hover { background-color: #005f7f; }

    /* CSS ตาราง */
    .report-table {width: 100%; border-collapse: collapse; font-family: sans-serif; font-size: 14px;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px;}
    .report-table th {background-color: #f2f2f2; text-align: center; font-weight: bold;}

    .pass-ok {color: green; font-weight: bold;}
    .pass-no {color: red; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 15px;}
    .load-value {color: #D32F2F !important; font-weight: bold;}
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
# 3. CALCULATION LOGIC (ACI 318-19 COLUMN)
# ==========================================
def calculate_interaction_curve(b, h, cover, main_db, nx, ny, fc, fy):
    """สร้างจุดบนกราฟ P-M Interaction Diagram"""
    d_prime = cover + 10 + main_db / 2  # approx d'
    d = h - d_prime

    Ast = (2 * nx + 2 * (ny - 2)) * (math.pi * (main_db / 2) ** 2)
    As_face = Ast / 2.0  # Simplify to 2 layers for curve generation

    points = []

    # Pure Compression Point (Po)
    Ag = b * h
    Po = 0.85 * fc * (Ag - Ast) + fy * Ast
    Pn_max = 0.80 * Po  # Tied Column

    # Generate points by varying Neutral Axis (c)
    c_values = np.linspace(1.5 * h, 0.1 * h, 40)

    for c in c_values:
        eps_cu = 0.003
        beta1 = beta1FromFc(fc)
        a = beta1 * c

        # Concrete Force
        Cc = 0.85 * fc * b * min(a, h)

        # Steel Forces
        # Layer 1 (Compression side)
        eps_s1 = eps_cu * (c - d_prime) / c
        fs1 = min(fy, 200000 * eps_s1)
        fs1 = max(-fy, fs1)
        Fs1 = As_face * fs1

        # Layer 2 (Tension side)
        eps_s2 = eps_cu * (c - d) / c
        fs2 = min(fy, 200000 * eps_s2)
        fs2 = max(-fy, fs2)
        Fs2 = As_face * fs2

        Pn = Cc + Fs1 + Fs2

        # Moment about Plastic Centroid (h/2)
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

        # Cap Pn at Pn_max/0.65 (Theoretical limit before phi)
        # But we plot phi*Pn directly
        phiPn = min(phi * Pn, 0.65 * Pn_max / 0.65)  # Apply 0.80 factor logic roughly
        # Better: apply phi to Pn, but cap at phi*Pn_max(0.80Po)
        # Note: phi for compression controlled is 0.65

        phiPn_val = phi * Pn
        # Cap at phi(0.65) * 0.80 * Po
        limit_top = 0.65 * 0.80 * Po
        if phiPn_val > limit_top: phiPn_val = limit_top

        points.append({'P': phiPn_val, 'M': phi * Mn, 'phi': phi})

    # Pure Bending Point (approx)
    points.append({'P': 0, 'M': points[-1]['M']})

    return points, Ag, Ast


def process_column_calculation(inputs):
    rows = []

    def sec(title):
        rows.append(["SECTION", title, "", "", "", "", ""])

    def row(item, formula, subs, result, unit, status=""):
        rows.append([item, formula, subs, result, unit, status])

    # 1. Inputs & Conversions
    b = inputs['b'] * 10;
    h = inputs['h'] * 10  # mm
    cover = inputs['cover'] * 10
    fc = inputs['fc'] * 0.0980665  # MPa
    fy = inputs['fy'] * 0.0980665
    fyt = inputs['fyt'] * 0.0980665

    main_key = inputs['mainBar']
    tie_key = inputs['tieBar']
    nx = int(inputs['nx'])
    ny = int(inputs['ny'])

    Pu_tf = inputs['Pu']
    Mu_tfm = inputs['Mu']

    # --- 1. MATERIAL & GEOMETRY ---
    sec("1. MATERIAL & SECTION PROPERTIES")
    row("Concrete & Steel", "fc', fy", f"{fmt(fc, 2)}, {fmt(fy, 0)}", "-", "MPa")
    row("Section Size", "b x h", f"{fmt(b, 0)} x {fmt(h, 0)}", "-", "mm")

    beta1 = beta1FromFc(fc)
    row("β1 Factor", "0.85 - 0.05(fc'-28)/7", f"fc'={fmt(fc, 2)}", f"{fmt(beta1, 2)}", "-")

    Ag = b * h
    row("Gross Area (Ag)", "b · h", f"{fmt(b, 0)}·{fmt(h, 0)}", f"{fmt(Ag, 0)}", "mm²")

    # Rebar Calculation
    total_bars = 2 * nx + 2 * max(0, ny - 2)
    bar_area_one = BAR_INFO[main_key]['A_cm2'] * 100
    Ast = total_bars * bar_area_one

    row("Main Reinforcement", f"Total {total_bars}-{main_key}",
        f"{total_bars} x {fmt(bar_area_one, 0)}",
        f"{fmt(Ast, 0)}", "mm²")

    rho_g = Ast / Ag
    status_rho = "OK" if 0.01 <= rho_g <= 0.08 else "FAIL"
    row("Reinforcement Ratio", "ρg = Ast / Ag", f"{fmt(Ast, 0)} / {fmt(Ag, 0)}", f"{fmt(rho_g * 100, 2)}", "%",
        status_rho)

    # --- 2. AXIAL CAPACITY (Compression) ---
    sec("2. AXIAL LOAD CAPACITY")

    # Po
    # Po = 0.85 fc (Ag - Ast) + fy Ast
    Po_N = 0.85 * fc * (Ag - Ast) + fy * Ast
    Po_tf = Po_N / 9806.65

    row("Nominal Axial (Po)", "0.85fc'(Ag-Ast) + fy·Ast",
        f"0.85·{fmt(fc, 1)}·({fmt(Ag, 0)}-{fmt(Ast, 0)}) + ...",
        f"{fmt(Po_tf, 2)}", "tf")

    # Phi Pn Max
    # phi = 0.65, factor = 0.80
    phi_c = 0.65
    phiPn_max_N = phi_c * 0.80 * Po_N
    phiPn_max_tf = phiPn_max_N / 9806.65

    row("Max Design Axial", "φPn,max = 0.65·0.80·Po",
        f"0.52 · {fmt(Po_tf, 2)}",
        f"{fmt(phiPn_max_tf, 2)}", "tf")

    # Check Axial
    row("Load Input (Pu)", "-", "-", f"{fmt(Pu_tf, 3)}", "tf", "", )  # Red color handled by CSS

    status_axial = "PASS" if Pu_tf <= phiPn_max_tf else "FAIL"
    row("Axial Check", "Pu ≤ φPn,max", f"{fmt(Pu_tf, 2)} ≤ {fmt(phiPn_max_tf, 2)}", status_axial, "-", status_axial)

    # --- 3. TIE DESIGN ---
    sec("3. TIE (STIRRUP) DESIGN")
    db_main = BAR_INFO[main_key]['d_mm']
    db_tie = BAR_INFO[tie_key]['d_mm']

    s1 = 16 * db_main
    s2 = 48 * db_tie
    s3 = min(b, h)
    s_req = min(s1, s2, s3)

    row("Spacing Limit 1", "16 · db(main)", f"16 · {db_main}", f"{s1:.0f}", "mm")
    row("Spacing Limit 2", "48 · db(tie)", f"48 · {db_tie}", f"{s2:.0f}", "mm")
    row("Spacing Limit 3", "Least Dimension", f"min({b},{h})", f"{s3:.0f}", "mm")

    s_prov = math.floor(s_req / 25.0) * 25.0
    if s_prov < 50: s_prov = 50

    row("Provide Ties", f"Use {tie_key}", f"min({s1:.0f},{s2:.0f},{s3:.0f})", f"@{s_prov / 10:.0f} cm", "-", "OK")

    # --- 4. INTERACTION CHECK ---
    sec("4. MOMENT CAPACITY CHECK")

    curve_points, _, _ = calculate_interaction_curve(b, h, cover, db_main, nx, ny, fc, fy)

    # Find M capacity at Pu (Interpolation)
    Pu_N = Pu_tf * 9806.65
    m_cap_Nmm = 0

    # Simple search
    found = False
    for i in range(len(curve_points) - 1):
        p1 = curve_points[i]['P']
        p2 = curve_points[i + 1]['P']
        # Check if Pu is within range [p2, p1] (Descending P)
        if p2 <= Pu_N <= p1:
            ratio = (Pu_N - p2) / (p1 - p2 + 1e-9)
            m1 = curve_points[i]['M']
            m2 = curve_points[i + 1]['M']
            m_cap_Nmm = m2 + ratio * (m1 - m2)
            found = True
            break

    if not found:
        if Pu_N > curve_points[0]['P']:
            m_cap_Nmm = 0  # Above max
        else:
            m_cap_Nmm = curve_points[-1]['M']  # Very low load

    m_cap_tfm = m_cap_Nmm / 9806650.0

    row("Load Input (Mu)", "-", "-", f"{fmt(Mu_tfm, 3)}", "tf-m", "")
    row("Moment Capacity", "φMn @ Pu", f"Interpolated from Curve", f"{fmt(m_cap_tfm, 2)}", "tf-m")

    status_pm = "PASS" if Mu_tfm <= m_cap_tfm else "FAIL"
    row("Interaction Check", "Mu ≤ φMn", f"{fmt(Mu_tfm, 2)} ≤ {fmt(m_cap_tfm, 2)}", status_pm, "-", status_pm)

    sec("5. FINAL STATUS")
    overall = "OK" if (status_rho == "OK" and status_axial == "PASS" and status_pm == "PASS") else "NOT OK"
    row("Overall", "-", "-", "DESIGN COMPLETE", "-", overall)

    return rows, curve_points, total_bars, s_prov


# ==========================================
# 4. PLOTTING
# ==========================================
def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def plot_column_section(b, h, cover, main_db, tie_db, nx, ny, tie_s, title="Column Section"):
    fig, ax = plt.subplots(figsize=(4, 4))
    rect = patches.Rectangle((0, 0), b, h, linewidth=2, edgecolor='#333', facecolor='#eee')
    ax.add_patch(rect)
    margin = cover + tie_db / 2
    rect_tie = patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, linewidth=2, edgecolor='#1976D2',
                                 facecolor='none', linestyle='-')
    ax.add_patch(rect_tie)

    start_x = margin + main_db / 2;
    end_x = b - margin - main_db / 2
    xs = np.linspace(start_x, end_x, nx) if nx > 1 else [b / 2]
    start_y = margin + main_db / 2;
    end_y = h - margin - main_db / 2
    ys = np.linspace(start_y, end_y, ny) if ny > 1 else [h / 2]

    # Plot Bars
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
    ax.set_title(title, fontweight='bold')

    info = f"Size: {b / 10:.0f}x{h / 10:.0f} cm\nMain: {2 * nx + 2 * max(0, ny - 2)}-DB{main_db:.0f}\nTies: RB{tie_db:.0f}@{tie_s / 10:.0f}cm"
    ax.text(b / 2, -h * 0.2, info, ha='center', va='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
    return fig


def plot_interaction_diagram(curve_points, Pu_tf, Mu_tfm):
    fig, ax = plt.subplots(figsize=(5, 5))
    Ms = [p['M'] / 9806650.0 for p in curve_points]
    Ps = [p['P'] / 9806.65 for p in curve_points]

    ax.plot(Ms, Ps, 'b-', linewidth=2, label='Capacity φPn-φMn')
    ax.plot([0, Ms[-1]], [0, Ps[-1]], 'b--')
    ax.plot([0, 0], [0, Ps[0]], 'b-')
    ax.plot(Mu_tfm, Pu_tf, 'ro', markersize=8, label='Design Load')

    ax.set_xlabel('Moment φMn (tf-m)');
    ax.set_ylabel('Axial Load φPn (tf)')
    ax.set_title('P-M Interaction Diagram', fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6);
    ax.legend()
    return fig


# ==========================================
# 5. REPORT GENERATOR
# ==========================================
def generate_column_report(inputs, rows, img_sect, img_pm):
    table_rows = ""
    for r in rows:
        if r[0] == "SECTION":
            table_rows += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            status_cls = "pass-ok" if "OK" in r[5] or "PASS" in r[5] else "pass-no"
            val_cls = "load-value" if "Load Input" in str(r[0]) else ""
            table_rows += f"""
            <tr>
                <td>{r[0]}</td>
                <td>{r[1]}</td>
                <td>{r[2]}</td>
                <td class='{val_cls}'>{r[3]}</td>
                <td>{r[4]}</td>
                <td class='{status_cls}'>{r[5]}</td>
            </tr>
            """

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
            .header {{ position: relative; margin-bottom: 20px; border-bottom: 2px solid #333; padding-bottom: 10px; }}
            .beam-box {{
                position: absolute; top: 0; right: 0;
                border: 2px solid #333; padding: 5px 15px;
                font-size: 18px; font-weight: bold;
            }}
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

            @media print {{
                .no-print {{ display: none !important; }}
                body {{ padding: 0; }}
            }}
            .print-btn-internal {{
                background-color: #4CAF50; color: white; padding: 12px 24px;
                border: none; border-radius: 5px; cursor: pointer; font-size: 16px; margin-bottom: 20px;
            }}
        </style>
    </head>
    <body>
        <div class="no-print" style="text-align: center;">
            <button onclick="window.print()" class="print-btn-internal">🖨️ Print This Page / พิมพ์หน้านี้</button>
        </div>

        <div class="header">
            <div class="beam-box">{inputs['col_id']}</div>
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

        <br><br><br><br>

        <h3>Calculation Details</h3>
        <table>
            <thead>
                <tr>
                    <th width="20%">Item</th>
                    <th width="30%">Formula</th>
                    <th width="25%">Substitution</th>
                    <th>Result</th>
                    <th>Unit</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
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

if 'calc_done' not in st.session_state:
    st.session_state['calc_done'] = False

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
    c1, c2 = st.columns(2)
    mainBar = c1.selectbox("Main Bar", list(BAR_INFO.keys()), index=4)  # DB16
    tieBar = c2.selectbox("Tie Bar", ['RB6', 'RB9', 'DB10'], index=0)

    st.write("Number of bars per face:")
    c1, c2 = st.columns(2)
    nx = c1.number_input("Nx (bars along X)", 2, help="จำนวนเส้นในแนวแกน X (รวมมุม)")
    ny = c2.number_input("Ny (bars along Y)", 2, help="จำนวนเส้นในแนวแกน Y (รวมมุม)")

    st.header("3. Loads (Factored)")
    Pu = st.number_input("Axial Load Pu (tf)", 40.0)
    Mu = st.number_input("Moment Mu (tf-m)", 2.0)

    run_btn = st.form_submit_button("Run Design")

if run_btn:
    inputs = {
        'project': project, 'col_id': col_id, 'engineer': engineer,
        'fc': fc, 'fy': fy, 'fyt': fyt,
        'b': b, 'h': h, 'cover': cover,
        'mainBar': mainBar, 'tieBar': tieBar,
        'nx': nx, 'ny': ny,
        'Pu': Pu, 'Mu': Mu
    }

    # 1. Calculate
    rows, curve, total_bars, s_prov = process_column_calculation(inputs)

    # 2. Draw
    # Section
    main_db = BAR_INFO[mainBar]['d_mm']
    tie_db = BAR_INFO[tieBar]['d_mm']
    fig_sect = plot_column_section(b * 10, h * 10, cover * 10, main_db, tie_db, int(nx), int(ny), s_prov)
    img_sect = fig_to_base64(fig_sect)

    # P-M Curve
    fig_pm = plot_interaction_diagram(curve, Pu, Mu)
    img_pm = fig_to_base64(fig_pm)

    # 3. Report
    html_report = generate_column_report(inputs, rows, img_sect, img_pm)

    st.success(f"✅ Design Complete: Use {total_bars}-{mainBar}, Ties {tieBar}@{s_prov / 10:.0f}cm")
    st.components.v1.html(html_report, height=800, scrolling=True)

else:
    st.info("👈 กรุณากรอกข้อมูลเสาด้านซ้ายแล้วกด 'Run Design'")