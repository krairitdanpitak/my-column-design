import streamlit as st
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import math
import io
import base64
import streamlit.components.v1 as components
from datetime import date

# ==========================================
# 1. SETUP & CSS (REPORT STYLE MATCHING PDF)
# ==========================================
st.set_page_config(page_title="RC Two-Way Slab Design SDM", layout="wide")

st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Sarabun:wght@400;500;700&display=swap');

    /* Main settings */
    .report-container {
        font-family: 'Sarabun', sans-serif;
        max-width: 210mm; /* A4 Width */
        margin: 0 auto;
        padding: 20px;
        background-color: white;
        color: #333;
    }

    /* Header */
    .report-header { text-align: center; margin-bottom: 15px; border-bottom: 2px solid #ddd; padding-bottom: 10px; }
    .report-title { font-size: 22px; font-weight: bold; margin: 0; }
    .report-subtitle { font-size: 16px; color: #555; margin-top: 5px; }
    .id-badge { float: right; border: 2px solid #333; padding: 2px 8px; font-weight: bold; font-size: 14px; }

    /* Info Box (Grey Background) */
    .info-box {
        background-color: #f8f9fa;
        border: 1px solid #dee2e6;
        padding: 15px;
        margin-bottom: 20px;
        border-radius: 4px;
        display: flex;
        justify-content: space-between;
        font-size: 14px;
    }
    .info-col { width: 48%; }

    /* Table Styles */
    .calc-table { width: 100%; border-collapse: collapse; font-size: 13px; margin-bottom: 25px; }
    .calc-table th { background-color: #f1f3f5; border: 1px solid #ced4da; padding: 10px; text-align: center; font-weight: bold; }
    .calc-table td { border: 1px solid #ced4da; padding: 8px 10px; vertical-align: middle; }

    /* Section Headers in Table */
    .section-row { background-color: #e9ecef; font-weight: bold; text-align: left; padding-left: 15px !important; color: #495057; }

    /* Status Colors */
    .status-ok { color: #28a745; font-weight: bold; text-align: center; }
    .status-fail { color: #dc3545; font-weight: bold; text-align: center; }
    .status-warn { color: #ffc107; font-weight: bold; text-align: center; }
    .result-val { font-weight: bold; color: #000; text-align: center; }

    /* Print Button (Green) */
    .print-btn-internal {
        background-color: #28a745; color: white; padding: 10px 20px;
        border: none; border-radius: 5px; font-family: 'Sarabun', sans-serif;
        font-weight: bold; cursor: pointer; text-decoration: none;
        display: inline-block; margin-bottom: 20px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .print-btn-internal:hover { background-color: #218838; }

    /* Plot Image */
    .plot-container { text-align: center; margin-bottom: 20px; border: 1px solid #eee; padding: 10px; }
    .plot-img { max-width: 100%; height: auto; }

    /* Conclusion Box */
    .conclusion-box {
        border: 2px solid #28a745; padding: 15px; text-align: center; margin-top: 30px; border-radius: 5px;
    }

    @media print {
        .no-print { display: none !important; }
        .report-container { width: 100%; max-width: none; padding: 0; box-shadow: none; }
        body { margin: 0; padding: 0; background-color: white; }
    }
</style>
""", unsafe_allow_html=True)

# ==========================================
# 2. DATA & LOGIC (UPDATED FULL ACI TABLE)
# ==========================================
BAR_INFO = {
    'RB6': {'A_cm2': 0.283, 'd_mm': 6},
    'RB9': {'A_cm2': 0.636, 'd_mm': 9},
    'DB10': {'A_cm2': 0.785, 'd_mm': 10},
    'DB12': {'A_cm2': 1.131, 'd_mm': 12},
    'DB16': {'A_cm2': 2.011, 'd_mm': 16},
    'DB20': {'A_cm2': 3.142, 'd_mm': 20}
}

# ACI 318 Method 2 Coefficients Dictionary (Complete Table)
# Structure: Case ID -> m (La/Lb) -> [Ca_neg, Ca_dl, Ca_ll, Cb_neg, Cb_dl, Cb_ll] *Modified order for function return
# Standard Table Data is usually: [Ca_neg, Cb_neg, Ca_dl, Cb_dl, Ca_ll, Cb_ll]
# Let's use the standard format: [Ca_neg, Cb_neg, Ca_dl, Cb_dl, Ca_ll, Cb_ll]
ACI_METHOD2_DATA = {
    # Case 1: All 4 edges Discontinuous
    1: {
        1.00: [0.00, 0.00, 0.036, 0.036, 0.036, 0.036],
        0.95: [0.00, 0.00, 0.040, 0.033, 0.040, 0.033],
        0.90: [0.00, 0.00, 0.045, 0.029, 0.045, 0.029],
        0.85: [0.00, 0.00, 0.050, 0.026, 0.050, 0.026],
        0.80: [0.00, 0.00, 0.056, 0.023, 0.056, 0.023],
        0.75: [0.00, 0.00, 0.061, 0.019, 0.061, 0.019],
        0.70: [0.00, 0.00, 0.068, 0.016, 0.068, 0.016],
        0.65: [0.00, 0.00, 0.074, 0.013, 0.074, 0.013],
        0.60: [0.00, 0.00, 0.081, 0.010, 0.081, 0.010],
        0.55: [0.00, 0.00, 0.088, 0.008, 0.088, 0.008],
        0.50: [0.00, 0.00, 0.095, 0.006, 0.095, 0.006],
    },
    # Case 2: All 4 edges Continuous
    2: {
        1.00: [0.033, 0.033, 0.018, 0.018, 0.027, 0.027],
        0.95: [0.038, 0.030, 0.020, 0.016, 0.030, 0.025],
        0.90: [0.043, 0.026, 0.022, 0.014, 0.034, 0.022],
        0.85: [0.049, 0.023, 0.025, 0.012, 0.037, 0.019],
        0.80: [0.055, 0.019, 0.028, 0.010, 0.041, 0.017],
        0.75: [0.061, 0.016, 0.031, 0.009, 0.045, 0.015],
        0.70: [0.068, 0.013, 0.034, 0.007, 0.049, 0.013],
        0.65: [0.074, 0.010, 0.037, 0.006, 0.053, 0.010],
        0.60: [0.080, 0.007, 0.041, 0.005, 0.058, 0.008],
        0.55: [0.085, 0.005, 0.044, 0.004, 0.062, 0.006],
        0.50: [0.090, 0.004, 0.048, 0.003, 0.066, 0.005],
    },
    # Case 3: 1 Short edge Discontinuous (3 Continuous)
    3: {
        1.00: [0.033, 0.033, 0.018, 0.023, 0.027, 0.027],
        0.95: [0.038, 0.029, 0.020, 0.020, 0.030, 0.024],
        0.90: [0.043, 0.025, 0.022, 0.018, 0.034, 0.022],
        0.85: [0.049, 0.021, 0.025, 0.016, 0.037, 0.019],
        0.80: [0.055, 0.018, 0.028, 0.013, 0.041, 0.017],
        0.75: [0.061, 0.015, 0.031, 0.011, 0.045, 0.015],
        0.70: [0.068, 0.012, 0.034, 0.009, 0.049, 0.012],
        0.65: [0.074, 0.009, 0.037, 0.008, 0.053, 0.010],
        0.60: [0.080, 0.007, 0.041, 0.006, 0.058, 0.008],
        0.55: [0.085, 0.005, 0.044, 0.005, 0.062, 0.006],
        0.50: [0.090, 0.003, 0.048, 0.004, 0.066, 0.005],
    },
    # Case 4: 1 Long edge Discontinuous (3 Continuous)
    4: {
        1.00: [0.033, 0.033, 0.022, 0.018, 0.027, 0.027],
        0.95: [0.041, 0.027, 0.025, 0.016, 0.030, 0.025],
        0.90: [0.048, 0.023, 0.029, 0.014, 0.035, 0.022],
        0.85: [0.055, 0.019, 0.032, 0.012, 0.040, 0.019],
        0.80: [0.062, 0.016, 0.036, 0.010, 0.044, 0.017],
        0.75: [0.069, 0.013, 0.040, 0.009, 0.049, 0.015],
        0.70: [0.076, 0.010, 0.044, 0.007, 0.054, 0.013],
        0.65: [0.082, 0.008, 0.048, 0.006, 0.059, 0.011],
        0.60: [0.088, 0.006, 0.053, 0.005, 0.063, 0.009],
        0.55: [0.093, 0.004, 0.057, 0.004, 0.068, 0.007],
        0.50: [0.098, 0.003, 0.062, 0.003, 0.072, 0.005],
    },
    # Case 5: 2 Short edges Discontinuous (2 Long Continuous)
    5: {
        1.00: [0.033, 0.033, 0.018, 0.026, 0.027, 0.027],
        0.95: [0.038, 0.000, 0.020, 0.024, 0.030, 0.022],
        0.90: [0.043, 0.000, 0.022, 0.021, 0.034, 0.021],
        0.85: [0.049, 0.000, 0.025, 0.018, 0.037, 0.019],
        0.80: [0.055, 0.000, 0.028, 0.015, 0.041, 0.017],
        0.75: [0.061, 0.000, 0.031, 0.013, 0.045, 0.015],
        0.70: [0.068, 0.000, 0.034, 0.011, 0.049, 0.013],
        0.65: [0.074, 0.000, 0.037, 0.009, 0.053, 0.010],
        0.60: [0.080, 0.000, 0.041, 0.008, 0.058, 0.008],
        0.55: [0.085, 0.000, 0.044, 0.006, 0.062, 0.006],
        0.50: [0.090, 0.000, 0.048, 0.005, 0.066, 0.005],
    },
    # Case 6: 2 Long edges Discontinuous (2 Short Continuous)
    6: {
        1.00: [0.033, 0.033, 0.027, 0.018, 0.027, 0.027],
        0.95: [0.000, 0.029, 0.031, 0.016, 0.030, 0.025],
        0.90: [0.000, 0.025, 0.036, 0.014, 0.035, 0.022],
        0.85: [0.000, 0.021, 0.041, 0.012, 0.040, 0.019],
        0.80: [0.000, 0.017, 0.047, 0.010, 0.045, 0.017],
        0.75: [0.000, 0.014, 0.052, 0.009, 0.050, 0.015],
        0.70: [0.000, 0.011, 0.058, 0.007, 0.056, 0.013],
        0.65: [0.000, 0.008, 0.064, 0.006, 0.062, 0.011],
        0.60: [0.000, 0.006, 0.070, 0.005, 0.068, 0.009],
        0.55: [0.000, 0.004, 0.076, 0.004, 0.074, 0.007],
        0.50: [0.000, 0.003, 0.083, 0.003, 0.080, 0.005],
    },
    # Case 7: 2 Adjacent edges Discontinuous (Corner)
    7: {
        1.00: [0.033, 0.033, 0.022, 0.022, 0.028, 0.028],
        0.95: [0.041, 0.029, 0.025, 0.020, 0.031, 0.025],
        0.90: [0.048, 0.025, 0.029, 0.018, 0.036, 0.022],
        0.85: [0.055, 0.021, 0.032, 0.016, 0.041, 0.019],
        0.80: [0.062, 0.018, 0.036, 0.014, 0.045, 0.017],
        0.75: [0.069, 0.015, 0.040, 0.012, 0.050, 0.015],
        0.70: [0.076, 0.012, 0.044, 0.010, 0.055, 0.013],
        0.65: [0.082, 0.009, 0.048, 0.008, 0.060, 0.011],
        0.60: [0.088, 0.007, 0.053, 0.007, 0.064, 0.009],
        0.55: [0.093, 0.005, 0.057, 0.005, 0.069, 0.007],
        0.50: [0.098, 0.004, 0.062, 0.004, 0.074, 0.005],
    },
    # Case 8: 3 edges Discontinuous (1 Long Continuous)
    8: {
        1.00: [0.033, 0.000, 0.026, 0.022, 0.028, 0.028],
        0.95: [0.041, 0.000, 0.030, 0.020, 0.031, 0.025],
        0.90: [0.048, 0.000, 0.035, 0.018, 0.036, 0.022],
        0.85: [0.055, 0.000, 0.039, 0.016, 0.041, 0.019],
        0.80: [0.062, 0.000, 0.044, 0.014, 0.045, 0.017],
        0.75: [0.069, 0.000, 0.048, 0.012, 0.050, 0.015],
        0.70: [0.076, 0.000, 0.052, 0.010, 0.055, 0.013],
        0.65: [0.082, 0.000, 0.057, 0.008, 0.060, 0.011],
        0.60: [0.088, 0.000, 0.062, 0.007, 0.064, 0.009],
        0.55: [0.093, 0.000, 0.067, 0.005, 0.069, 0.007],
        0.50: [0.098, 0.000, 0.073, 0.004, 0.074, 0.005],
    },
    # Case 9: 3 edges Discontinuous (1 Short Continuous)
    9: {
        1.00: [0.000, 0.033, 0.022, 0.026, 0.028, 0.028],
        0.95: [0.000, 0.029, 0.025, 0.024, 0.031, 0.025],
        0.90: [0.000, 0.025, 0.029, 0.021, 0.036, 0.022],
        0.85: [0.000, 0.021, 0.032, 0.018, 0.041, 0.019],
        0.80: [0.000, 0.018, 0.036, 0.015, 0.045, 0.017],
        0.75: [0.000, 0.015, 0.040, 0.013, 0.050, 0.015],
        0.70: [0.000, 0.012, 0.044, 0.011, 0.055, 0.013],
        0.65: [0.000, 0.009, 0.048, 0.009, 0.060, 0.011],
        0.60: [0.000, 0.007, 0.053, 0.008, 0.064, 0.009],
        0.55: [0.000, 0.005, 0.057, 0.006, 0.069, 0.007],
        0.50: [0.000, 0.004, 0.062, 0.005, 0.074, 0.005],
    }
}


def get_moment_coefficients(case_id, m):
    """
    Returns [Ca_neg, Ca_dl, Ca_ll, Cb_neg, Cb_dl, Cb_ll]
    Uses linear interpolation for exact m values.
    """
    if m < 0.5: m = 0.5
    if m > 1.0: m = 1.0

    data = ACI_METHOD2_DATA.get(case_id)
    if not data: return 0.0, 0.05, 0.05, 0.0, 0.05, 0.05

    # Check for exact match
    if m in data:
        vals = data[m]
        # Data structure: [Ca_neg, Cb_neg, Ca_dl, Cb_dl, Ca_ll, Cb_ll]
        # Return format: [Ca_neg, Ca_dl, Ca_ll, Cb_neg, Cb_dl, Cb_ll] (to match unpacking in logic)
        return vals[0], vals[2], vals[4], vals[1], vals[3], vals[5]

    # Linear Interpolation
    ratios = sorted(data.keys())
    lower = max([r for r in ratios if r <= m])
    upper = min([r for r in ratios if r >= m])

    if lower == upper:
        vals = data[lower]
        return vals[0], vals[2], vals[4], vals[1], vals[3], vals[5]

    frac = (m - lower) / (upper - lower)
    res = []
    v_low = data[lower]
    v_up = data[upper]

    for i in range(6):
        val = v_low[i] + frac * (v_up[i] - v_low[i])
        res.append(val)

    # Return reordered for logic: Ca_neg, Ca_dl, Ca_ll, Cb_neg, Cb_dl, Cb_ll
    return res[0], res[2], res[4], res[1], res[3], res[5]


def fmt(n, digits=2):
    try:
        return f"{float(n):,.{digits}f}"
    except:
        return "-"


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    buf.seek(0)
    return f"data:image/png;base64,{base64.b64encode(buf.read()).decode()}"


# ==========================================
# 3. PLOTTING (SECTION VIEW)
# ==========================================
def plot_slab_section_detailed(Lx, h, cover, main_bar, s_main, top_s, case_id):
    fig, ax = plt.subplots(figsize=(8, 3.5))

    h_m = h / 100
    cov_m = cover / 100
    supp_w = 0.25

    # Concrete
    ax.add_patch(patches.Rectangle((0, 0), Lx, h_m, facecolor='white', edgecolor='black', linewidth=1.2))
    ax.add_patch(
        patches.Rectangle((-supp_w, -0.4), supp_w, 0.4 + h_m, facecolor='white', edgecolor='black', linewidth=1.2))
    ax.add_patch(patches.Rectangle((Lx, -0.4), supp_w, 0.4 + h_m, facecolor='white', edgecolor='black', linewidth=1.2))

    # Rebar - Bottom (Main)
    y_bot = cov_m
    ax.plot([-0.15, Lx + 0.15], [y_bot, y_bot], 'k-', linewidth=1.5)
    # Rebar - Top
    L_top = Lx / 3.5
    y_top = h_m - cov_m
    if top_s > 0:
        ax.plot([-0.15, L_top], [y_top, y_top], 'k-', linewidth=1.5)
        ax.plot([Lx - L_top, Lx + 0.15], [y_top, y_top], 'k-', linewidth=1.5)
        # Hooks
        ax.plot([L_top, L_top], [y_top, y_top - 0.08], 'k-', linewidth=1.5)
        ax.plot([Lx - L_top, Lx - L_top], [y_top, y_top - 0.08], 'k-', linewidth=1.5)

    # Dots (Cross section bars)
    n_dots = int(Lx / 0.2)
    for i in range(1, n_dots):
        ax.add_patch(patches.Circle((i * 0.2, y_bot + 0.015), 0.008, color='black'))

    # Dimensions
    # Span
    ax.annotate('', xy=(0, -0.15), xytext=(Lx, -0.15), arrowprops=dict(arrowstyle='|-|', linewidth=0.8))
    ax.text(Lx / 2, -0.2, f"{Lx:.2f} m.", ha='center', va='top', fontsize=10)

    # Labels
    top_txt = f"{main_bar}@{top_s:.2f}" if top_s > 0 else "Min Steel"
    bot_txt = f"{main_bar}@{s_main:.2f}"

    if top_s > 0:
        ax.text(L_top / 2, y_top + 0.1, top_txt, ha='center', fontsize=9, color='blue')
    ax.text(Lx / 2, y_bot - 0.1, bot_txt, ha='center', fontsize=9, color='blue')

    ax.axis('off')
    ax.set_ylim(-0.5, h_m + 0.4)
    ax.set_xlim(-0.5, Lx + 0.5)
    return fig


# ==========================================
# 4. REPORT GENERATOR
# ==========================================
def generate_report(inputs, rows, img_sec):
    table_html = ""
    for r in rows:
        if r[0].startswith("SECTION"):
            table_html += f"<tr class='section-row'><td colspan='6'>{r[1]}</td></tr>"
        else:
            item, formula, sub, res, unit, stat = r
            stat_cls = "status-ok" if stat in ["OK", "PASS", "COMPLETE"] else (
                "status-warn" if stat == "WARN" else "status-fail")
            if stat == "": stat_cls = ""

            table_html += f"""
            <tr>
                <td>{item}</td>
                <td style='color:#555;'>{formula}</td>
                <td style='color:#555;'>{sub}</td>
                <td class='result-val'>{res}</td>
                <td style='text-align:center;'>{unit}</td>
                <td class='{stat_cls}'>{stat}</td>
            </tr>"""

    today = date.today().strftime("%d/%m/%Y")

    html = f"""
    <div class="no-print" style="text-align: center;">
        <button onclick="window.print()" class="print-btn-internal">🖨️ Print This Page / พิมพ์หน้านี้</button>
    </div>

    <div class="report-container">
        <div class="report-header">
            <div class="id-badge">{inputs['slab_id']}</div>
            <div class="report-title">ENGINEERING DESIGN REPORT</div>
            <div class="report-subtitle">RC Two-Way Slab Design SDM (ACI 318 Method 2)</div>
        </div>

        <div class="info-box">
            <div class="info-col">
                <strong>Project:</strong> {inputs['project']}<br>
                <strong>Engineer:</strong> {inputs['engineer']}<br>
                <strong>Date:</strong> {today}
            </div>
            <div class="info-col">
                <strong>Panel Size:</strong> {inputs['Lx']} x {inputs['Ly']} m.<br>
                <strong>Thickness:</strong> {inputs['h']} cm (Cover {inputs['cover']} cm)<br>
                <strong>Materials:</strong> fc'={inputs['fc']} ksc, fy={inputs['fy']} ksc
            </div>
        </div>

        <div style="font-weight: bold; margin-bottom: 5px;">Design Visualization (Section View)</div>
        <div class="plot-container">
            <img src="{img_sec}" class="plot-img">
            <div style="font-size:12px; margin-top:5px;">Section showing Short Span Reinforcement</div>
        </div>

        <div style="font-weight: bold; margin-bottom: 5px;">Calculation Details</div>
        <table class="calc-table">
            <thead>
                <tr>
                    <th width="20%">Item</th>
                    <th width="20%">Formula</th>
                    <th width="25%">Substitution</th>
                    <th width="15%">Result</th>
                    <th width="10%">Unit</th>
                    <th width="10%">Status</th>
                </tr>
            </thead>
            <tbody>
                {table_html}
            </tbody>
        </table>

        <div class="conclusion-box">
            <div style="font-size: 14px; color: #555;">CONCLUSION</div>
            <div style="font-size: 24px; font-weight: bold; color: #28a745; margin: 5px 0;">DESIGN COMPLETE</div>
            <div style="font-size: 12px; color: #777;">Signed by: {inputs['engineer']}</div>
        </div>
    </div>
    """
    return html


# ==========================================
# 5. MAIN LOGIC
# ==========================================
st.title("RC Two-Way Slab Design SDM")
st.caption("ออกแบบพื้นคอนกรีตเสริมเหล็กสองทาง (ACI Method 2) | Generate Calculation Report")

with st.sidebar.form("input_form"):
    st.header("Project Info")
    project = st.text_input("Project Name", "อาคารสำนักงาน 2 ชั้น")
    slab_id = st.text_input("Slab Mark", "S-01")
    engineer = st.text_input("Engineer", "นายไกรฤทธิ์ ด่านพิทักษ์")

    st.header("1. Geometry")
    c1, c2 = st.columns(2)
    # Allowed 0.0 input, validation happens later
    Lx = c1.number_input("Short Span (m)", min_value=0.0, value=4.0, step=0.1)
    Ly = c2.number_input("Long Span (m)", min_value=0.0, value=5.0, step=0.1)

    c3, c4 = st.columns(2)
    h = c3.number_input("Thickness (cm)", min_value=0.0, value=15.0, step=1.0)
    cover = c4.number_input("Cover (cm)", min_value=0.0, value=2.5)

    st.header("2. Loads & Material")
    sdl = st.number_input("SDL (kg/m²)", 150.0)
    ll = st.number_input("LL (kg/m²)", 300.0)
    fc = st.number_input("fc' (ksc)", 240)
    fy = st.number_input("fy (ksc)", 4000)

    st.header("3. Design Parameters")
    case_opts = {1: "1. Simple", 2: "2. All Cont", 3: "3. One Short Discont", 4: "4. One Long Discont",
                 5: "5. Two Short Discont", 6: "6. Two Long Discont", 7: "7. Corner", 8: "8. One Long Cont",
                 9: "9. One Short Cont"}
    case_sel = st.selectbox("Support Case (ACI)", list(case_opts.values()))
    mainBar = st.selectbox("Main Rebar", list(BAR_INFO.keys()), index=1)

    run_btn = st.form_submit_button("Run Calculation")

if run_btn:
    # Validation for Zero inputs
    if Lx <= 0 or Ly <= 0 or h <= 0:
        st.error("❌ Error: Dimensions (Lx, Ly, Thickness) must be greater than 0.")
    else:
        # Prepare inputs
        inputs = {'project': project, 'slab_id': slab_id, 'engineer': engineer,
                  'Lx': Lx, 'Ly': Ly, 'h': h, 'cover': cover,
                  'fc': fc, 'fy': fy, 'sdl': sdl, 'll': ll,
                  'case': int(case_sel.split(".")[0]), 'mainBar': mainBar}

        # Calculation
        rows = []
        short, long_s = min(Lx, Ly), max(Lx, Ly)
        m = short / long_s

        # 1. Geometry
        rows.append(["SECTION", "1. GEOMETRY & SLAB TYPE", "", "", "", ""])
        rows.append(["Short Span (Lx)", "-", "-", fmt(short), "m", ""])
        rows.append(["Long Span (Ly)", "-", "-", fmt(long_s), "m", ""])
        rows.append(["Ratio Ly/Lx", "Long / Short", f"{fmt(long_s)}/{fmt(short)}", fmt(long_s / short), "-",
                     "OK" if (long_s / short) <= 2 else "WARN"])

        # 2. Loads
        w_d = 2400 * (h / 100) + sdl
        wu = 1.2 * w_d + 1.6 * ll
        rows.append(["SECTION", "2. LOAD ANALYSIS", "", "", "", ""])
        rows.append(["Factored Load (wu)", "1.2D + 1.6L", f"1.2({w_d:.0f})+1.6({ll})", fmt(wu), "kg/m²", ""])

        # 3. Moments & Steel
        Ca_neg, Ca_dl, Ca_ll, Cb_neg, Cb_dl, Cb_ll = get_moment_coefficients(inputs['case'], m)
        La_sq = short ** 2
        Ma_pos = (Ca_dl * 1.2 * w_d + Ca_ll * 1.6 * ll) * La_sq
        Ma_neg = Ca_neg * wu * La_sq

        rows.append(["SECTION", "3. SHORT SPAN DESIGN (MAIN)", "", "", "", ""])
        rows.append(["Design Moment (+)", "C_pos * wu * Lx²", "Method 2 Coeffs", fmt(Ma_pos), "kg-m", ""])

        # Rebar Logic
        db = BAR_INFO[mainBar]['d_mm']
        d = h - cover - db / 20


        def get_spacing(Mu):
            if Mu <= 1: return 0
            Rn = (Mu * 100) / (0.9 * 100 * d ** 2) * 1000
            rho = 0.85 * fc / fy * (1 - math.sqrt(max(0, 1 - 2 * Rn / (0.85 * fc))))
            As_req = max(rho * 100 * d, 0.0018 * 100 * h)
            s = math.floor(min((BAR_INFO[mainBar]['A_cm2'] * 100) / As_req, 3 * h, 45) * 2) / 2
            return min(s, 30.0)


        s_pos = get_spacing(Ma_pos)
        rows.append(["Effective Depth (d)", "h - cov - db/2", f"{h}-{cover}-{db / 10}/2", fmt(d), "cm", ""])
        rows.append(["Provide Bottom", f"Use {mainBar}", "-", f"@{s_pos} cm", "-", "OK"])

        s_neg = 0
        if Ca_neg > 0:
            s_neg = get_spacing(Ma_neg)
            rows.append(["Design Moment (-)", "C_neg * wu * Lx²", f"{Ca_neg:.3f}*{wu:.0f}...", fmt(Ma_neg), "kg-m", ""])
            rows.append(["Provide Top (Supp)", f"Use {mainBar}", "-", f"@{s_neg} cm", "-", "OK"])

        # 4. Long Span
        rows.append(["SECTION", "4. LONG SPAN DESIGN", "", "", "", ""])
        Mb_pos = (Cb_dl * 1.2 * w_d + Cb_ll * 1.6 * ll) * La_sq
        s_long = get_spacing(Mb_pos)
        rows.append(["Provide Bottom", f"Use {mainBar}", f"Mu={fmt(Mb_pos)}", f"@{s_long} cm", "-", "OK"])

        # 5. Checks
        rows.append(["SECTION", "5. CHECKS", "", "", "", ""])
        Vu = wu * short / 2
        phiVc = 0.85 * 0.53 * math.sqrt(fc) * 100 * d
        chk = "PASS" if phiVc >= Vu else "FAIL"
        rows.append(["Shear Capacity", "phiVc >= Vu", f"{fmt(phiVc)} >= {fmt(Vu)}", chk, "kg", chk])

        # Generate & Show
        img = fig_to_base64(plot_slab_section_detailed(Lx, h, cover, mainBar, s_pos, s_neg, inputs['case']))
        html_report = generate_report(inputs, rows, img)
        st.success("✅ Calculation Finished")
        components.html(html_report, height=1200, scrolling=True)
else:
    st.info("👈 Please enter slab dimensions and click Run.")