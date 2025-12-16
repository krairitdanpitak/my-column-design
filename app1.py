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
# 1. SETUP & UTILS
# ==========================================
st.set_page_config(page_title="RC Design Suite Pro", layout="wide", page_icon="🏗️")

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


def fig_to_base64(fig):
    buf = io.BytesIO();
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=100);
    buf.seek(0)
    plt.close(fig)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


def draw_dim(ax, p1, p2, text, offset=30, color='black'):
    x1, y1 = p1;
    x2, y2 = p2
    angle = math.atan2(y2 - y1, x2 - x1);
    perp = angle + math.pi / 2
    ox = offset * math.cos(perp);
    oy = offset * math.sin(perp)
    p1o = (x1 + ox, y1 + oy);
    p2o = (x2 + ox, y2 + oy)
    ax.plot([x1, p1o[0]], [y1, p1o[1]], color=color, lw=0.5)
    ax.plot([x2, p2o[0]], [y2, p2o[1]], color=color, lw=0.5)
    ax.annotate('', xy=p1o, xytext=p2o, arrowprops=dict(arrowstyle='<->', color=color, lw=0.8))
    mx = (p1o[0] + p2o[0]) / 2;
    my = (p1o[1] + p2o[1]) / 2
    deg = math.degrees(angle)
    if 90 < deg <= 270:
        deg -= 180
    elif -270 <= deg < -90:
        deg += 180
    tx = mx + 10 * math.cos(perp);
    ty = my + 10 * math.sin(perp)
    ax.text(tx, ty, text, ha='center', va='center', rotation=deg, fontsize=8, color=color,
            bbox=dict(fc='white', ec='none', alpha=0.8))


st.markdown("""
<style>
    .report-table {width: 100%; border-collapse: collapse; font-family: sans-serif; font-size: 13px;}
    .report-table th, .report-table td {border: 1px solid #ddd; padding: 8px;}
    .report-table th {background-color: #f2f2f2; text-align: center; font-weight: bold;}
    .sec-row {background-color: #e0e0e0; font-weight: bold; font-size: 14px;}
    .pass-ok {color: green; font-weight: bold; text-align: center;}
    .pass-no {color: red; font-weight: bold; text-align: center;}
    .load-val {color: #D32F2F !important; font-weight: bold;}
    .drawing-container {display: flex; justify-content: center; gap: 20px; flex-wrap: wrap; margin-top: 20px;}
    .drawing-box {border: 1px solid #ddd; padding: 10px; background: white; text-align: center; min-width: 300px;}
</style>
""", unsafe_allow_html=True)


# ==========================================
# 2. MODULE: BEAM DESIGN
# ==========================================
def process_beam_detailed(inputs):
    rows = []

    def sec(t):
        rows.append(["SECTION", t, "", "", "", ""])

    def row(i, f, s, r, u, st=""):
        rows.append([i, f, s, r, u, st])

    b = inputs['b'] * 10;
    h = inputs['h'] * 10;
    cov = inputs['cov'] * 10
    fc = inputs['fc'] * 0.0981;
    fy = inputs['fy'] * 0.0981
    db_s = BAR_INFO[inputs['s_bar']]['d_mm']
    db_m = BAR_INFO[inputs['m_bar']]['d_mm']
    d = h - cov - db_s - db_m / 2

    sec("1. MATERIAL & SECTION PARAMETERS")
    row("Materials", "fc', fy", f"{inputs['fc']:.0f}, {inputs['fy']:.0f}", "-", "ksc", "")
    row("Section Size", "b x h", f"{b:.0f} x {h:.0f}", "-", "mm", "")
    row("Effective Depth d", "h-cov-ds-dm/2", f"{h:.0f}-{cov}-{db_s}-{db_m / 2:.1f}", f"{d:.1f}", "mm", "")
    as_min = max(0.25 * math.sqrt(fc) / fy, 1.4 / fy) * b * d
    row("As,min", "max(...)bd", f"max(...)*{b}*{d:.0f}", f"{as_min:.0f}", "mm²", "")

    locs = [
        ('Left (Top) Mu-', inputs['mu_Ln']), ('Left (Bot) Mu+', inputs['mu_Lp']),
        ('Mid (Top) Mu-', inputs['mu_Mn']), ('Mid (Bot) Mu+', inputs['mu_Mp']),
        ('Right (Top) Mu-', inputs['mu_Rn']), ('Right (Bot) Mu+', inputs['mu_Rp'])
    ]
    bar_res = {}
    sec("2. FLEXURAL DESIGN")

    for name, mu in locs:
        mu_nmm = mu * 9806650
        req_as = mu_nmm / (0.9 * fy * 0.9 * d) if mu > 0 else 0
        des_as = max(req_as, as_min) if mu > 0 else 0
        n = max(math.ceil(des_as / (BAR_INFO[inputs['m_bar']]['A_cm2'] * 100)), 2)
        if mu <= 0.01: n = 2

        prov_as = n * BAR_INFO[inputs['m_bar']]['A_cm2'] * 100
        a = (prov_as * fy) / (0.85 * fc * b);
        phi_mn = 0.9 * prov_as * fy * (d - a / 2)
        stt = "PASS" if phi_mn >= mu_nmm else "FAIL"

        key = name.split()[0] + ("_Top" if "Top" in name else "_Bot")
        bar_res[key] = f"{n}-{inputs['m_bar']}"

        if mu > 0.01:
            row(f"{name}: Mu", "-", "-", f"{mu:.2f}", "tf-m", "")
            row(f"{name}: As,req", "Mn≥Mu", f"As_min={as_min:.0f}", f"{req_as:.0f}", "mm²", "")
            row(f"{name}: Provide", f"Use {inputs['m_bar']}", f"Req:{des_as:.0f}",
                f"{n}-{inputs['m_bar']} ({prov_as:.0f})", "mm²", "OK")
            row(f"{name}: Strength", "φMn ≥ Mu", f"{phi_mn / 9.8e6:.2f} ≥ {mu:.2f}", "PASS", "tf-m", stt)

    sec("3. SHEAR DESIGN")
    vc = 0.17 * math.sqrt(fc) * b * d;
    phi_vc = 0.75 * vc
    row("Capacity φVc", "0.75·0.17√fc·bd", f"0.75·0.17·√{fc:.1f}·{b}·{d:.0f}", f"{phi_vc / 9806:.2f}", "tf", "")

    shear_locs = [("Left", inputs['vu_L']), ("Mid", inputs['vu_M']), ("Right", inputs['vu_R'])]
    stir_res = {}
    av = 2 * BAR_INFO[inputs['s_bar']]['A_cm2'] * 100

    for loc, vu in shear_locs:
        vu_n = vu * 9806
        if vu_n > phi_vc:
            vs_req = (vu_n / 0.75) - vc
            s_req = (av * fy * d) / vs_req
            s_prov = math.floor(min(s_req, d / 2, 600) / 10) * 10
            stir_res[loc] = f"@{s_prov / 10:.0f}cm"
            row(f"{loc}: Vu", "-", "-", f"{vu:.2f}", "tf", "")
            row(f"{loc}: Vs,req", "Vu/φ - Vc", f"{vu:.2f}/0.75 - {vc / 9806:.2f}", f"{vs_req / 9806:.2f}", "tf", "")
            row(f"{loc}: Provide", f"Use {inputs['s_bar']}", f"min({s_req:.0f}, d/2)", f"@{s_prov / 10:.0f} cm", "-",
                "OK")
        else:
            s_prov = math.floor(min(d / 2, 600) / 10) * 10
            stir_res[loc] = f"@{s_prov / 10:.0f}cm"
            row(f"{loc}: Vu", "-", "-", f"{vu:.2f}", "tf", "")
            row(f"{loc}: Check", "Vu ≤ φVc", f"{vu:.2f} ≤ {phi_vc / 9806:.2f}", "Min Stirrup", "-", "PASS")

    return rows, bar_res, stir_res


def plot_beam_elevation(bar_res, stir_res, b, h):
    fig, ax = plt.subplots(figsize=(10, 3.5))
    L = 1000
    ax.add_patch(patches.Rectangle((0, 0), L, h, ec='black', fc='white', lw=2))
    ax.plot([L / 3, L / 3], [0, h], 'k--', lw=0.5);
    ax.plot([2 * L / 3, 2 * L / 3], [0, h], 'k--', lw=0.5)

    ax.text(L / 6, h + 20, "Left Support", ha='center', fontweight='bold')
    ax.text(L / 2, h + 20, "Mid Span", ha='center', fontweight='bold')
    ax.text(5 * L / 6, h + 20, "Right Support", ha='center', fontweight='bold')

    cols = ['blue', 'red']
    for i, zone in enumerate(['Left', 'Mid', 'Right']):
        x_pos = L / 6 + i * (L / 3)
        ax.text(x_pos, h - 30, bar_res.get(f'{zone}_Top', '-'), ha='center', color=cols[0], fontweight='bold')
        ax.text(x_pos, 30, bar_res.get(f'{zone}_Bot', '-'), ha='center', color=cols[1], fontweight='bold')
        ax.text(x_pos, h / 2, f"Stir: {stir_res.get(zone, '-')}", ha='center', fontsize=9,
                bbox=dict(fc='white', ec='none', alpha=0.8))

    ax.set_xlim(-50, L + 50);
    ax.set_ylim(-50, h + 50);
    ax.axis('off')
    return fig


# ==========================================
# 3. MODULE: COLUMN DESIGN (Auto + Detail)
# ==========================================
def calculate_pm_curve(b, h, cover, db, nx, ny, fc, fy):
    points = []
    d_prime = cover + 10 + db / 2
    ast = (2 * nx + 2 * max(0, ny - 2)) * (math.pi * (db / 2) ** 2)
    po = 0.85 * fc * (b * h - ast) + fy * ast
    pn_max = 0.8 * po

    c_vals = np.linspace(1.5 * h, 0.1 * h, 30)
    for c in c_vals:
        a = 0.85 * c;
        cc = 0.85 * fc * b * min(a, h)
        fs1 = min(fy, 200000 * 0.003 * (c - d_prime) / c);
        fs1 = max(-fy, fs1)
        fs2 = min(fy, 200000 * 0.003 * (c - (h - d_prime)) / c);
        fs2 = max(-fy, fs2)
        pn = cc + (ast / 2) * fs1 + (ast / 2) * fs2
        mn = cc * (h / 2 - a / 2) + (ast / 2) * fs1 * (h / 2 - d_prime) - (ast / 2) * fs2 * (h / 2 - d_prime)

        phi = 0.65
        phi_pn = phi * pn;
        phi_mn = phi * mn
        if phi_pn > 0.65 * pn_max: phi_pn = 0.65 * pn_max
        points.append({'P': phi_pn, 'M': phi_mn})
    points.append({'P': 0, 'M': points[-1]['M']})
    return points, b * h, ast, po, 0.65 * pn_max


def process_column_detailed(inputs):
    rows = []

    def sec(t):
        rows.append(["SECTION", t, "", "", "", ""])

    def row(i, f, s, r, u, st=""):
        rows.append([i, f, s, r, u, st])

    b = inputs['b'] * 10;
    h = inputs['h'] * 10;
    cov = inputs['cov'] * 10
    fc = inputs['fc'] * 0.0981;
    fy = inputs['fy'] * 0.0981

    nx, ny = (2, 2);
    db = BAR_INFO[inputs['m_bar']]['d_mm']

    # --- AUTO DESIGN LOGIC START ---
    if inputs['mode'] == 'Auto':
        found = False
        candidates = []
        # Try configurations from 2x2 up to 8x8
        for i_x in range(2, 9):
            for i_y in range(2, 9):
                # Calculate Properties
                curr_ast = (2 * i_x + 2 * max(0, i_y - 2)) * (math.pi * (db / 2) ** 2)
                rho = curr_ast / (b * h)

                # Check Rho Limits (1% - 8%)
                if not (0.01 <= rho <= 0.08): continue

                # Check Capacity
                curve, _, _, _, pmax = calculate_pm_curve(b, h, cov, db, i_x, i_y, fc, fy)
                pu_n = inputs['pu'] * 9806

                # Check if Pu <= Pmax (Basic check first)
                if pu_n > pmax: continue

                # Check if point (Mu, Pu) is inside Curve
                # Simple interpolation check
                is_safe = False
                mu_nmm = inputs['mu'] * 9806650
                m_cap_at_pu = 0
                for i in range(len(curve) - 1):
                    if curve[i + 1]['P'] <= pu_n <= curve[i]['P']:
                        r = (pu_n - curve[i + 1]['P']) / (curve[i]['P'] - curve[i + 1]['P'] + 1e-9)
                        m_cap_at_pu = curve[i + 1]['M'] + r * (curve[i]['M'] - curve[i + 1]['M'])
                        break

                if mu_nmm <= m_cap_at_pu:
                    candidates.append({'nx': i_x, 'ny': i_y, 'ast': curr_ast})

        # Select best candidate (Min Steel Area)
        if candidates:
            candidates.sort(key=lambda x: x['ast'])
            best = candidates[0]
            nx, ny = best['nx'], best['ny']
        else:
            # Fallback if no design found (e.g. section too small)
            nx, ny = 2, 2
            # We will still run calc to show it fails
    else:
        nx, ny = int(inputs['nx']), int(inputs['ny'])
    # --- AUTO DESIGN LOGIC END ---

    curve, ag, ast, po, pmax = calculate_pm_curve(b, h, cov, db, nx, ny, fc, fy)

    # 1. PROPERTIES
    sec("1. MATERIAL & SECTION PROPERTIES")
    row("Concrete & Steel", "fc', fy", f"{inputs['fc']:.0f}, {inputs['fy']:.0f}", "-", "ksc", "")
    row("Section Size", "b x h", f"{b} x {h}", "-", "mm", "")
    beta1 = 0.85 if fc <= 28 else max(0.65, 0.85 - 0.05 * (fc - 28) / 7)
    row("β1 Factor", "0.85-0.05(fc'-28)/7", f"fc'={fc:.1f} MPa", f"{beta1:.2f}", "-", "")
    row("Gross Area (Ag)", "b·h", f"{b}·{h}", f"{ag:.0f}", "mm²", "")
    row("Main Reinforcement", f"Total {2 * nx + 2 * max(0, ny - 2)}-{inputs['m_bar']}", f"Total Area", f"{ast:.0f}",
        "mm²", "")
    rho = ast / ag
    row("Reinforcement Ratio", "ρg = Ast/Ag", f"{ast:.0f}/{ag:.0f}", f"{rho * 100:.2f}", "%",
        "OK" if 0.01 <= rho <= 0.08 else "FAIL")

    # 2. AXIAL
    sec("2. AXIAL LOAD CAPACITY")
    row("Nominal Axial (Po)", "0.85fc'(Ag-Ast)+fyAst", f"0.85·{fc:.1f}·({ag:.0f}-{ast:.0f})+...", f"{po / 9806:.2f}",
        "tf", "")
    row("Max Design Axial", "φPn,max = 0.65·0.80·Po", f"0.52·{po / 9806:.2f}", f"{pmax / 9806:.2f}", "tf", "")
    pu_n = inputs['pu'] * 9806
    row("Load Input (Pu)", "-", "-", f"{inputs['pu']:.2f}", "tf", "")
    row("Axial Check", "Pu ≤ φPn,max", f"{inputs['pu']:.2f} ≤ {pmax / 9806:.2f}", "PASS" if pu_n <= pmax else "FAIL",
        "-", "PASS" if pu_n <= pmax else "FAIL")

    # 3. TIE DESIGN
    sec("3. TIE (STIRRUP) DESIGN")
    db_main = BAR_INFO[inputs['m_bar']]['d_mm']
    db_tie = BAR_INFO[inputs['t_bar']]['d_mm']
    s1 = 16 * db_main;
    s2 = 48 * db_tie;
    s3 = min(b, h)
    row("Spacing Limit 1", "16·db(main)", f"16·{db_main}", f"{s1:.0f}", "mm", "")
    row("Spacing Limit 2", "48·db(tie)", f"48·{db_tie}", f"{s2:.0f}", "mm", "")
    row("Spacing Limit 3", "Least Dimension", f"min({b},{h})", f"{s3:.0f}", "mm", "")
    s_prov = math.floor(min(s1, s2, s3) / 10) * 10
    row("Provide Ties", f"Use {inputs['t_bar']}", f"min({s1:.0f},{s2:.0f},{s3:.0f})", f"@{s_prov / 10:.0f} cm", "-",
        "OK")

    # 4. MOMENT
    sec("4. MOMENT CAPACITY CHECK")
    mu_nmm = inputs['mu'] * 9806650;
    m_cap = 0
    for i in range(len(curve) - 1):
        if curve[i + 1]['P'] <= pu_n <= curve[i]['P']:
            r = (pu_n - curve[i + 1]['P']) / (curve[i]['P'] - curve[i + 1]['P'] + 1e-9)
            m_cap = curve[i + 1]['M'] + r * (curve[i]['M'] - curve[i + 1]['M'])
            break

    row("Load Input (Mu)", "-", "-", f"{inputs['mu']:.2f}", "tf-m", "")
    row("Moment Capacity", "φMn @ Pu", "Interpolated from Curve", f"{m_cap / 9.8e6:.2f}", "tf-m", "")
    row("Interaction Check", "Mu ≤ φMn", f"{inputs['mu']:.2f} ≤ {m_cap / 9.8e6:.2f}",
        "PASS" if mu_nmm <= m_cap else "FAIL", "-", "PASS" if mu_nmm <= m_cap else "FAIL")

    return rows, curve, nx, ny, s_prov


def plot_col_sect_detailed(b, h, cov, nx, ny, m_bar, t_bar, s_tie):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.add_patch(patches.Rectangle((0, 0), b, h, ec='k', fc='#eee', lw=2))
    margin = cov + 6
    ax.add_patch(patches.Rectangle((margin, margin), b - 2 * margin, h - 2 * margin, ec='blue', fc='none'))

    xs = np.linspace(margin + 6, b - margin - 6, nx);
    ys = np.linspace(margin + 6, h - margin - 6, ny)
    for x in xs: ax.add_patch(patches.Circle((x, margin + 6), 6, fc='red')); ax.add_patch(
        patches.Circle((x, h - margin - 6), 6, fc='red'))
    for y in ys[1:-1]: ax.add_patch(patches.Circle((margin + 6, y), 6, fc='red')); ax.add_patch(
        patches.Circle((b - margin - 6, y), 6, fc='red'))

    text = f"Size: {b / 10}x{h / 10} cm\nMain: {2 * nx + 2 * max(0, ny - 2)}-{m_bar}\nTies: {t_bar}@{s_tie / 10:.0f}cm"
    ax.text(b / 2, -h * 0.2, text, ha='center', va='top', fontsize=10, bbox=dict(fc='white', ec='black'))
    ax.set_xlim(-b * 0.2, b * 1.2);
    ax.set_ylim(-h * 0.5, h * 1.2);
    ax.axis('off');
    ax.set_aspect('equal')
    return fig


def plot_pm_curve(curve, pu, mu):
    fig, ax = plt.subplots(figsize=(5, 5))
    ms = [p['M'] / 9.8e6 for p in curve];
    ps = [p['P'] / 9806 for p in curve]
    ax.plot(ms, ps, 'b-', lw=2, label='Capacity φPn-φMn')
    ax.plot([0, ms[-1]], [0, ps[-1]], 'b--');
    ax.plot([0, 0], [0, ps[0]], 'b-')
    ax.plot(mu, pu, 'ro', ms=8, label='Design Load')
    ax.set_xlabel('Moment φMn (tf-m)');
    ax.set_ylabel('Axial Load φPn (tf)')
    ax.grid(True, ls='--', alpha=0.5);
    ax.legend();
    ax.set_title("P-M Interaction Diagram")
    return fig


# ==========================================
# 4. MODULE: FOOTING DESIGN
# ==========================================
def process_footing_detailed(inputs):
    rows = []

    def sec(t):
        rows.append(["SECTION", t, "", "", "", ""])

    def row(i, f, s, r, u, st=""):
        rows.append([i, f, s, r, u, st])

    fc = inputs['fc'] * 0.0981;
    fy = inputs['fy'] * 0.0981
    pu = inputs['pu'];
    n_pile = int(inputs['n_pile'])
    s = inputs['s'] * 1000;
    edge = inputs['edge'] * 1000
    dp = inputs['dp'] * 1000;
    col = 300
    h_final = inputs['h'] * 1000;
    cover = 75
    db = BAR_INFO[inputs['m_bar']]['d_mm']
    d = h_final - cover - db

    coords = []
    if n_pile == 1:
        coords = [(0, 0)]
    elif n_pile == 2:
        coords = [(-s / 2, 0), (s / 2, 0)]
    elif n_pile == 3:
        coords = [(-s / 2, -s * 0.288), (s / 2, -s * 0.288), (0, s * 0.577)]
    elif n_pile == 4:
        coords = [(-s / 2, -s / 2), (s / 2, -s / 2), (-s / 2, s / 2), (s / 2, s / 2)]
    elif n_pile == 5:
        coords = [(-s / 2, -s / 2), (s / 2, -s / 2), (-s / 2, s / 2), (s / 2, s / 2), (0, 0)]

    bx = (max([abs(x) for x, _ in coords]) * 2) + dp + 2 * edge if n_pile > 1 else dp + 2 * edge
    by = (max([abs(y) for _, y in coords]) * 2) + dp + 2 * edge if n_pile > 1 else dp + 2 * edge
    bx = max(bx, col + 2 * edge);
    by = max(by, col + 2 * edge)

    sec("1. GEOMETRY")
    row("Size", "B x L", f"{bx:.0f}x{by:.0f}", f"h={h_final:.0f}", "mm", "")
    lambda_s = math.sqrt(2 / (1 + 0.004 * d))
    row("Size Effect λs", "√(2/(1+0.004d))", f"√(2/(1+0.004*{d:.0f}))", f"{lambda_s:.3f}", "-", "≤1.0")

    sec("2. PILE REACTION")
    p_avg = pu / n_pile
    row("Load per Pile", "Ru = Pu / N", f"{pu}/{n_pile}", f"{p_avg:.2f}", "tf",
        "PASS" if p_avg <= inputs['cap'] else "FAIL")

    sec("3. FLEXURAL DESIGN (X & Y)")
    p_n = p_avg * 9806;
    mx = 0;
    my = 0
    for x, y in coords:
        lx = abs(x) - col / 2;
        ly = abs(y) - col / 2
        if lx > 0: mx += p_n * lx
        if ly > 0: my += p_n * ly

    dirs = [('X-Dir (Long)', mx, by), ('Y-Dir (Short)', my, bx)]
    res_bars = {}

    for label, mom, width in dirs:
        req_as = mom / (0.9 * fy * 0.9 * d) if mom > 0 else 0
        min_as = 0.0018 * width * h_final
        des_as = max(req_as, min_as)
        n = math.ceil(des_as / (BAR_INFO[inputs['m_bar']]['A_cm2'] * 100))
        if n_pile == 1: n = max(n, 4)
        prov_as = n * BAR_INFO[inputs['m_bar']]['A_cm2'] * 100

        row(f"Moment Mu ({label})", "Σ P(arm)", "-", f"{mom / 9.8e6:.2f}", "tf-m", "")
        row(f"As Req ({label})", "Max(Calc, Min)", f"Max({req_as:.0f}, {min_as:.0f})", f"{des_as:.0f}", "mm²", "")
        row(f"Provide ({label})", f"{n}-{inputs['m_bar']}", f"As={prov_as:.0f}", "OK", "-", "")
        key = 'X-Dir' if 'X-Dir' in label else 'Y-Dir'
        res_bars[key] = n

    if n_pile > 1:
        sec("4. SHEAR (ACI 318-19)")
        bo = 4 * (col + d)
        vc_p = 0.33 * lambda_s * math.sqrt(fc) * bo * d
        vu_p = sum([p_n for x, y in coords if max(abs(x), abs(y)) > (col + d) / 2])
        phi_vc_p = 0.75 * vc_p

        row("Punching: Vu", "Sum Outside", "-", f"{vu_p / 9806:.2f}", "tf", "")
        row("Punching: φVc", "0.75·0.33λs√fc·bo·d", f"0.75·0.33·{lambda_s:.2f}...", f"{phi_vc_p / 9806:.2f}", "tf",
            "PASS" if vu_p <= phi_vc_p else "FAIL")

        prov_as_x = res_bars.get('X-Dir', 4) * BAR_INFO[inputs['m_bar']]['A_cm2'] * 100
        rho_w = prov_as_x / (by * d);
        rho_term = math.pow(rho_w, 1 / 3)
        vc_b = 0.66 * lambda_s * rho_term * math.sqrt(fc) * by * d
        vu_b = sum([p_n for x, y in coords if abs(x) > col / 2 + d])
        phi_vc_b = 0.75 * vc_b

        row("Beam: Vu", "Sum Outside d", "-", f"{vu_b / 9806:.2f}", "tf", "")
        row("Beam: φVc", "0.75·0.66λs(ρ)^1/3√fc·bd", f"ρ^1/3={rho_term:.2f}...", f"{phi_vc_b / 9806:.2f}", "tf",
            "PASS" if vu_b <= phi_vc_b else "FAIL")

    return rows, coords, bx, by, res_bars


def plot_foot_combined(coords, bx, by, n_x, n_y, bar, h_mm, col_mm):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # PLAN
    ax1.set_title("PLAN VIEW", fontweight='bold')
    ax1.add_patch(patches.Rectangle((-bx / 2, -by / 2), bx, by, ec='k', fc='#f9f9f9', lw=2))
    ys = np.linspace(-by / 2 + 100, by / 2 - 100, min(n_y, 8))
    for y in ys: ax1.plot([-bx / 2 + 50, bx / 2 - 50], [y, y], 'b-', lw=1, alpha=0.5)
    xs = np.linspace(-bx / 2 + 100, bx / 2 - 100, min(n_x, 8))
    for x in xs: ax1.plot([x, x], [-by / 2 + 50, by / 2 - 50], 'r-', lw=1, alpha=0.5)
    ax1.add_patch(patches.Rectangle((-col_mm / 2, -col_mm / 2), col_mm, col_mm, ec='k', fc='#ddd', hatch='//'))
    for x, y in coords: ax1.add_patch(patches.Circle((x, y), 120, ec='k', ls='--'))

    draw_dim(ax1, (-bx / 2, -by / 2 - 200), (bx / 2, -by / 2 - 200), f"L={bx / 1000:.2f}m", 0)
    draw_dim(ax1, (-bx / 2 - 200, -by / 2), (-bx / 2 - 200, by / 2), f"B={by / 1000:.2f}m", 0)
    ax1.text(0, by / 2 + 150, f"{n_x}-{bar} (Y-Dir)", ha='center', color='red', fontweight='bold',
             bbox=dict(fc='white', ec='red'))
    ax1.text(bx / 2 + 150, 0, f"{n_y}-{bar} (X-Dir)", va='center', rotation=90, color='blue', fontweight='bold',
             bbox=dict(fc='white', ec='blue'))
    ax1.set_xlim(-bx / 1.1, bx / 1.1);
    ax1.set_ylim(-by / 1.1, by / 1.1);
    ax1.set_aspect('equal');
    ax1.axis('off')

    # SECTION
    ax2.set_title("SECTION DETAIL", fontweight='bold')
    ax2.plot([-bx, bx], [0, 0], 'k-', lw=0.5)
    ax2.add_patch(patches.Rectangle((-bx / 2, -h_mm), bx, h_mm, ec='k', fc='#f0f0f0', lw=2))
    ax2.add_patch(patches.Rectangle((-col_mm / 2, 0), col_mm, h_mm / 2, ec='k', fc='#fff', hatch='///'))

    cov = 75
    ax2.plot([-bx / 2 + cov, bx / 2 - cov], [-h_mm + cov, -h_mm + cov], 'r-', lw=3)
    ax2.plot([-bx / 2 + cov, -bx / 2 + cov], [-h_mm + cov, -h_mm + cov + h_mm * 0.6], 'r-', lw=3)
    ax2.plot([bx / 2 - cov, bx / 2 - cov], [-h_mm + cov, -h_mm + cov + h_mm * 0.6], 'r-', lw=3)

    draw_dim(ax2, (bx / 2 + 200, 0), (bx / 2 + 200, -h_mm), f"h={h_mm / 1000:.2f}m", 50)
    ax2.text(0, -h_mm + cov - 150, f"Main Reinforcement", ha='center', color='red', fontweight='bold')
    ax2.set_xlim(-bx / 1.1, bx / 1.1);
    ax2.set_ylim(-h_mm * 2, h_mm);
    ax2.set_aspect('equal');
    ax2.axis('off')

    return fig


# ==========================================
# 5. GENERATE REPORT
# ==========================================
def generate_report(title, rows, imgs, proj, eng):
    t_rows = ""
    for r in rows:
        if r[0] == "SECTION":
            t_rows += f"<tr class='sec-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            cls = "pass-ok" if "PASS" in r[5] or "OK" in r[5] else ("pass-no" if "FAIL" in r[5] else "")
            val_cls = "load-val" if "Mu" in r[0] or "Vu" in r[0] or "Pu" in r[0] else ""
            t_rows += f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td class='{val_cls}'>{r[3]}</td><td>{r[4]}</td><td class='{cls}'>{r[5]}</td></tr>"

    img_html = "".join([f"<div class='drawing-box'><img src='{i}' style='max-width:100%'></div>" for i in imgs])
    print_btn = """<div style="text-align:center; margin-bottom:20px;">
        <button onclick="window.print()" class="print-btn">🖨️ Print / Save PDF</button>
    </div>"""

    return f"""
    <div style="font-family: Sarabun, sans-serif; padding: 20px;">
        {print_btn}
        <h2 style="text-align:center; border-bottom: 2px solid #333;">{title}</h2>
        <div style="display:flex; justify-content:space-between; margin-bottom:15px;">
            <div><strong>Project:</strong> {proj}</div><div><strong>Engineer:</strong> {eng}</div>
        </div>
        <div class="drawing-container">{img_html}</div><br>
        <table class="report-table">
            <thead><tr><th width="20%">Item</th><th width="25%">Formula</th><th width="30%">Substitution</th><th>Result</th><th>Unit</th><th>Status</th></tr></thead>
            <tbody>{t_rows}</tbody>
        </table>
        <div style="margin-top:40px; text-align:center;">
            <div style="display:inline-block; width:250px; text-align:left;">
                <strong>Designed by:</strong><br><br><div style="border-bottom:1px solid #000;"></div>
                <div style="text-align:center;">({eng})<br>วิศวกรโครงสร้าง</div>
            </div>
        </div>
    </div>
    """


# ==========================================
# 6. UI ROUTER
# ==========================================
st.sidebar.title("🏗️ Design Suite")
mode = st.sidebar.radio("Select Module", ["Beam", "Column", "Footing"])
st.sidebar.markdown("---")
proj = st.sidebar.text_input("Project", "Project A", key='p')
eng = st.sidebar.text_input("Engineer", "Eng. A", key='e')

if mode == "Beam":
    st.header("Beam Design (Detailed)")
    with st.sidebar.form("b"):
        c1, c2 = st.columns(2)
        fc = c1.number_input("fc'", value=240);
        fy = c2.number_input("fy", value=4000)
        b = c1.number_input("b (cm)", value=25);
        h = c2.number_input("h (cm)", value=50)
        cov = st.number_input("Cover", value=3.0);
        m_bar = st.selectbox("Main", list(BAR_INFO.keys()), index=4)
        s_bar = st.selectbox("Stirrup", list(BAR_INFO.keys()));
        st.write("Moments (tf-m):")
        c1, c2, c3 = st.columns(3)
        mu_Ln = c1.number_input("L-Neg", value=8.0);
        mu_Lp = c1.number_input("L-Pos", value=4.0)
        mu_Mn = c2.number_input("M-Neg", value=0.0);
        mu_Mp = c2.number_input("M-Pos", value=8.0)
        mu_Rn = c3.number_input("R-Neg", value=8.0);
        mu_Rp = c3.number_input("R-Pos", value=4.0)
        st.write("Shear (tf):")
        c1, c2, c3 = st.columns(3)
        vu_L = c1.number_input("Vu-L", value=12.0);
        vu_M = c2.number_input("Vu-M", value=8.0);
        vu_R = c3.number_input("Vu-R", value=12.0)
        run = st.form_submit_button("Calculate")
    if run:
        d = {'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cov': cov, 'm_bar': m_bar, 's_bar': s_bar,
             'mu_Ln': mu_Ln, 'mu_Lp': mu_Lp, 'mu_Mn': mu_Mn, 'mu_Mp': mu_Mp, 'mu_Rn': mu_Rn, 'mu_Rp': mu_Rp,
             'vu_L': vu_L, 'vu_M': vu_M, 'vu_R': vu_R}
        rows, bar_res, stir_res = process_beam_detailed(d)
        img = fig_to_base64(plot_beam_elevation(bar_res, stir_res, b * 10, h * 10))
        st.components.v1.html(generate_report("Beam Calculation Report", rows, [img], proj, eng), height=1200,
                              scrolling=True)

elif mode == "Column":
    st.header("Column Design (Auto Mode)")
    with st.sidebar.form("c"):
        c1, c2 = st.columns(2)
        fc = c1.number_input("fc'", value=240);
        fy = c2.number_input("fy", value=4000)
        b = c1.number_input("b", value=25);
        h = c2.number_input("h", value=25)
        cov = st.number_input("Cover", value=3.0);
        opt = st.radio("Mode", ["Auto", "Manual"])
        nx = st.number_input("Nx", value=2);
        ny = st.number_input("Ny", value=2)
        m_bar = st.selectbox("Main", list(BAR_INFO.keys()), index=4);
        t_bar = st.selectbox("Tie", ['RB6', 'RB9'])
        pu = st.number_input("Pu", value=40.0);
        mu = st.number_input("Mu", value=2.0)
        run = st.form_submit_button("Calculate")
    if run:
        d = {'fc': fc, 'fy': fy, 'b': b, 'h': h, 'cov': cov, 'mode': opt, 'nx': nx, 'ny': ny, 'm_bar': m_bar,
             't_bar': t_bar, 'pu': pu, 'mu': mu}
        rows, curve, bnx, bny, s_tie = process_column_detailed(d)
        img1 = fig_to_base64(plot_col_sect_detailed(b * 10, h * 10, cov * 10, bnx, bny, m_bar, t_bar, s_tie))
        img2 = fig_to_base64(plot_pm_curve(curve, pu, mu))
        st.components.v1.html(generate_report("Column Calculation Report", rows, [img1, img2], proj, eng), height=1200,
                              scrolling=True)

elif mode == "Footing":
    st.header("Footing Design (Detailed)")
    with st.sidebar.form("f"):
        c1, c2 = st.columns(2)
        fc = c1.number_input("fc'", value=240);
        fy = c2.number_input("fy", value=4000)
        n = st.selectbox("Piles", [1, 2, 3, 4, 5], index=3);
        dp = st.number_input("Dia", value=0.22)
        s = st.number_input("Space", value=0.8);
        h = st.number_input("Thk", value=0.5)
        edge = st.number_input("Edge", value=0.25);
        m_bar = st.selectbox("Main", list(BAR_INFO.keys()), index=4)
        pu = st.number_input("Pu", value=60.0);
        cap = st.number_input("Cap", value=30.0)
        auto = st.checkbox("Auto-H", value=True);
        run = st.form_submit_button("Calculate")
    if run:
        d = {'fc': fc, 'fy': fy, 'n_pile': n, 'dp': dp, 's': s, 'h': h, 'edge': edge, 'm_bar': m_bar, 'pu': pu,
             'cap': cap, 'auto_h': auto}
        rows, coords, bx, by, res_bars = process_footing_detailed(d)
        nx = res_bars.get('X-Dir', 4);
        ny = res_bars.get('Y-Dir', 4)
        img = fig_to_base64(plot_foot_combined(coords, bx, by, nx, ny, m_bar, h * 1000, 300))
        st.components.v1.html(generate_report("Pile Cap Calculation Report", rows, [img], proj, eng), height=1200,
                              scrolling=True)
