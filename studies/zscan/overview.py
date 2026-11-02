"""Overview of all z-scan datasets (from zscan_<dataset>.json): light speed per fibre type and dataset, the same channel in datasets that
share the DRS channel map, the attenuation results, and the test fibres (testfibres.py) fibre by fibre.
Usage: python3 studies/zscan/overview.py     Page: results/html/ZScan/ZScan_overview.html"""
import json, os
import numpy as np
from common import DATA, TYPES, TLAB
from datasets import PLOTS, HTML, DATASETS
import ROOT
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html

R = {d: json.load(open(f"{DATA}/zscan_{d}.json")) for d in DATASETS if os.path.exists(f"{DATA}/zscan_{d}.json") and not DATASETS[d].get("channel_map")}
TF = json.load(open(f"{DATA}/zscan_testfibres.json")) if os.path.exists(f"{DATA}/zscan_testfibres.json") else None
if TF: TF = dict(TF["selections"][TF["main_selection"]], selection=TF["main_selection"])
DRS = {d: J["selections"][J["main_selection"]]["drs"] for d, J in R.items()}
OUTD = f"{PLOTS}/overview"; os.makedirs(OUTD, exist_ok=True); OUT = f"{HTML}/ZScan_overview.html"
COL = {"Sci": ROOT.kGreen + 2, "CerQuartz": ROOT.kBlue + 1, "CerPlastic": ROOT.kRed + 1}; MK = {"Sci": 20, "CerQuartz": 21, "CerPlastic": 22}
SHORT = {"Sci": "scintillating", "CerQuartz": "quartz", "CerPlastic": "plastic"}
timed = lambda p: np.isfinite(p["time_slope"]) and p["time_resid"] < 0.1
fmt = lambda m, e: f"{m:+.3f} &plusmn; {e:.3f}" if np.isfinite(m) else "&ndash;"
TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
def lam(m, e):
    if not np.isfinite(m): return "&ndash;"
    return f"{1 / m:.1f}" if m - e > 0 else (f"&gt; {1 / (m + e):.1f}" if m + e > 0 else "none seen")

# light speed: median slope per type and dataset, with the 16-84% range of the channels
rows, G, L, C, M = [], [], [], [], []
slopes = {t: {} for t in TYPES}
for t in TYPES:
    row = [TLAB[t]]
    for d in R:
        v = np.array([p["time_slope"] for p in DRS[d]["channels"].values() if p["type"] == t and timed(p)])
        if not len(v): row += ["&ndash;"] * 3; continue
        m = np.median(v); slopes[t][d] = (m, np.percentile(v, 16), np.percentile(v, 84), len(v))
        row += [f"{m:.2f} ({len(v)} ch.)", f"{1 / abs(m) / 0.2998:.3f} c", f"{abs(m) * 0.2998:.2f}"]
    rows.append(row)
speed_tab = table(["Fibres"] + [f"{d}: {h}" for d in R for h in ("dt/dX [ns/m]", "speed", "n<sub>eff</sub>")], rows)
for t in TYPES:
    ds_ = [d for d in R if d in slopes[t]]
    g = ROOT.TGraphAsymmErrors(len(ds_))
    for i, d in enumerate(ds_):
        m, lo, hi, _ = slopes[t][d]; xi = list(R).index(d) + 0.1 * (TYPES.index(t) - 1)
        g.SetPoint(i, xi, m); g.SetPointError(i, 0, 0, m - lo, hi - m)
    G.append(g); L.append(TLAB[t]); C.append(COL[t]); M.append(MK[t])
extra = []
for i, d in enumerate(R):
    tx = ROOT.TLatex(i, -7.6, d); tx.SetTextAlign(21); tx.SetTextSize(0.04); extra.append(tx)
DrawHistos(G, L, -0.5, len(R) - 0.5, "dataset", -7.8, -4.5, "dt/dX [ns/m] (median, 16#minus84% of channels)", "speed_datasets", dology=False,
           drawoptions="P", legendoptions="P", markerstyles=M, mycolors=C, legendPos=[0.18, 0.74, 0.60, 0.89], extraToDraw=extra, outdir=OUTD,
           lumitext="z scans", usePDF=False)
plots = ["speed_datasets"]
# the same channel in two datasets with the same channel map: arrival-time slope
pair_rows = []
names = list(R)
for a in range(len(names)):
    for b in range(a + 1, len(names)):
        A_, B_ = DRS[names[a]]["channels"], DRS[names[b]]["channels"]
        common = [c for c in A_ if c in B_ and timed(A_[c]) and timed(B_[c])]
        if len(common) < 10: continue
        G, L, C = [], [], []
        for t in TYPES:
            cs = [c for c in common if A_[c]["type"] == t]
            if len(cs) < 3: continue
            x = np.array([A_[c]["time_slope"] for c in cs]); y = np.array([B_[c]["time_slope"] for c in cs])
            g = ROOT.TGraph(len(cs))
            for i in range(len(cs)): g.SetPoint(i, x[i], y[i])
            G.append(g); C.append(COL[t]); L.append(f"{SHORT[t]} ({len(cs)}): #Delta {np.median(y - x):+.2f}, rms {np.std(y - x):.2f} ns/m")
            pair_rows.append([f"{names[a]} vs {names[b]}", TLAB[t], len(cs), f"{np.median(y - x):+.3f}", f"{np.std(y - x):.3f}"])
        diag = ROOT.TLine(-7.5, -7.5, -4.5, -4.5); diag.SetLineStyle(2); diag.SetLineColor(ROOT.kGray + 2)
        nm = f"time_slope_{names[a]}_vs_{names[b]}"
        DrawHistos(G, L, -7.5, -4.5, f"dt/dX, {names[a]} [ns/m]", -7.5, -4.5, f"dt/dX, {names[b]} [ns/m]", nm, dology=False, drawoptions="P",
                   legendoptions="P", markerstyles=[20] * len(G), mycolors=C, legendPos=[0.15, 0.74, 0.80, 0.89], extraToDraw=[diag], outdir=OUTD,
                   lumitext="same channel, two datasets", usePDF=False)
        plots.append(nm)
pair_tab = table(["Datasets", "Fibres", "channels", "median difference [ns/m]", "rms [ns/m]"], pair_rows) if pair_rows else ""
# attenuation: FERS whole-module sums, DRS sums per category
rows = []
for d, J in R.items():
    if not J["fers"]: continue
    F = J["selections"][J["main_selection"]]["fers"]["sums"]["main"]
    rows.append([d, "FERS, whole module"] + [f"{fmt(F[t]['fit']['inv_lambda'], F[t]['fit']['err'])} ({lam(F[t]['fit']['inv_lambda'], F[t]['fit']['err'])} m)" for t in TYPES])
fers_tab = table(["Dataset", "Sum", *[TLAB[t] for t in TYPES]], rows) if rows else ""
cats = []
for d in R:
    for n in DRS[d]["sums"]["main"]["categories"]:
        if n not in cats: cats.append(n)
rows = []
for n in cats:
    row = [n]
    for d in R:
        c = DRS[d]["sums"]["main"]["categories"].get(n); f = c["fits"]["base"].get("adc") if c else None
        row.append(f"{fmt(f['common']['inv_lambda'], f['common']['err'])} ({len(c['channels'])} ch.)" if f else "&ndash;")
    rows.append(row)
drs_tab = table(["DRS sum (common run offsets), 1/&lambda; [/m]"] + list(R), rows)
nav = " &middot; ".join(f'<a href="ZScan_{d}.html">{d}: {DATASETS[d]["title"]}</a>' for d in list(R) + (["testfibres"] if TF else []))
nav += ' &middot; <a href="../Timing/Timing_resolution.html">timing resolution (studies/timing)</a>'
tf_tab = ""
if TF:
    from plotting.calox_plot_helper import fiber_base_type
    rows = []
    for l in sorted(TF["fibres"], key=lambda l: (fiber_base_type(l), l)):
        f = TF["fibres"][l]
        if "time_slope" not in f: continue
        c = f["fits"].get("lpc"); s_ = 1 / abs(f["time_slope"])
        rows.append([l, fmt(f["time_slope"], f["time_slope_err"]), f"{s_ / 0.2998:.3f} c", f"{0.2998 / s_:.2f}",
                     f"{fmt(c['inv_lambda'], c['err'])} ({lam(c['inv_lambda'], c['err'])} m)" if c else "&ndash;"])
    tf_tab = table(["Test fibre", "dt/dX [ns/m]", "speed", "n<sub>eff</sub>", "1/&lambda; [/m], light per crossing"], rows)
intro = f"""
<p>Z scans: 90&deg; positron runs in which the table moves the shower along the fibres (larger X = closer to the readout). Datasets: {nav}.</p>
<h3 style="margin:14px 0 4px 0">Light speed</h3>
{speed_tab}
<p>The arrival time depends only on how far the light travels to the readout: the speeds agree between the datasets although they sample
different parts of the module with different showers. The light depends on how much of the shower a channel sees.</p>
{pair_tab}
<h3 style="margin:14px 0 4px 0">Attenuation</h3>
<p>Only a sum over the whole module is independent of where the shower sits; with the fibres tilted relative to the table travel the
shower moves across the channels during a scan. Whole-module sums exist only where the calorimeter FERS was recorded:</p>
{fers_tab}
<p>DRS sums cover part of the module (the unsaturated channels of a category; different channels in different datasets) and are shown
for comparison:</p>
{drs_tab}
<h3 style="margin:14px 0 4px 0">Test fibres (same in-time hits, each fibre on its own; positrons with an MCP pulse, no hodoscope cuts)</h3>
{tf_tab}
"""
generate_html([f"{p}.png" for p in plots], OUTD, plots_per_row=4, output_html=OUT, title="Z scans: overview", intro_text=intro)
print(OUT)
