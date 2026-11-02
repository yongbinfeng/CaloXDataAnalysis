"""Per-channel galleries of a z-scan dataset (from zscan_<dataset>.json and the dumps), one row per channel:
  DRS:  in-time integral distributions at up to five positions; light per hit, hit fraction, light per event, all-event mean
        (bright channels) and arrival time vs X, with the fits (light: after the fitted run offsets, as in the fit)
  FERS: calibrated value distributions at up to five positions; mean per event vs X with the fit
Runs outside the main run set (2025: 1501) are drawn open and grey, not fitted.
Usage: python3 studies/zscan/gallery.py DATASET     Pages: results/html/ZScan/Gallery_<dataset>_<DRS|FERS>_<type>.html"""
import json, os, sys
import numpy as np
from common import Dataset, DATA, TYPES, TLAB, A_HIT, SAT_FRAC, MAX_RESID, load_drs, load_fers, matching
from datasets import PLOTS, HTML
import ROOT
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html

ds = Dataset(sys.argv[1])
J = json.load(open(f"{DATA}/zscan_{ds.name}.json")); MODE = J["main_selection"]; SEL = J["selections"][MODE]
if J["fers"]: matching(ds)
OUTD = f"{PLOTS}/{ds.name}/gallery"; os.makedirs(OUTD, exist_ok=True); os.makedirs(HTML, exist_ok=True)
D = SEL["drs"]; RUNS = D["runs"]; XM = np.array(D["x_m"]) * 1000; GRP = D["groups"]; MAIN = set(J["runsets"]["main"])
GCOL = [ROOT.kBlue + 1, ROOT.kRed + 1]; G_NAMES = sorted(set(GRP))
COLS5 = [ROOT.kBlue + 1, ROOT.kAzure + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1]
order = np.argsort(XM); first = [i for i in order if GRP[i] == G_NAMES[0]]
PICK = [first[j] for j in np.unique(np.linspace(0, len(first) - 1, min(5, len(first))).astype(int))]
short = lambda ch: ch.replace("DRS_Brg", "B").replace("DRS_Board", "B").replace("_Board", ".").replace("_Group", ".").replace("_Channel", ".")


def graph(xs, ys, es=None):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, float(x), float(y)); g.SetPointError(i, 0, float(es[i]) if es is not None else 0)
    return g


def vs_x(name, y, ey, ytit, lab, slope=None, line=False, note=""):
    """values vs X: one colour per group, runs outside the main run set open grey; exponential (or straight line) per group"""
    y = np.array(y, float); ey = np.array(ey, float) if ey is not None else None
    if not np.isfinite(y).sum() >= 2: return None
    G, C, M, extra = [], [], [], []
    for k, g in enumerate(G_NAMES):
        m = np.array([gg == g and r in MAIN and np.isfinite(yy) for gg, r, yy in zip(GRP, RUNS, y)])
        if not m.any(): continue
        G.append(graph(XM[m], y[m], ey[m] if ey is not None else None)); C.append(GCOL[k % 2]); M.append(20)
        if slope is not None and np.isfinite(slope) and m.sum() >= 2:
            b = slope / 1000.0
            if line:
                a0 = np.mean(y[m] - b * XM[m]); f = ROOT.TF1(f"f_{name}_{k}", "[0] + [1]*x", XM[m].min() - 10, XM[m].max() + 10)
            else:
                w = (y[m] / ey[m]) ** 2 if ey is not None else np.ones(m.sum())
                a0 = np.average(np.log(y[m]) - b * XM[m], weights=w); f = ROOT.TF1(f"f_{name}_{k}", "exp([0] + [1]*x)", XM[m].min() - 10, XM[m].max() + 10)
            f.SetParameters(a0, b); f.SetLineColor(GCOL[k % 2]); f.SetLineWidth(2); extra.append((f, "SAME"))
    out = np.array([r not in MAIN and np.isfinite(yy) for r, yy in zip(RUNS, y)])
    if out.any():
        g1 = graph(XM[out], y[out], ey[out] if ey is not None else None); g1.SetMarkerStyle(24); g1.SetMarkerColor(ROOT.kGray + 2)
        g1.SetLineColor(ROOT.kGray + 2); extra.append((g1, "P SAME"))
    if slope is not None and np.isfinite(slope):
        tx = ROOT.TLatex(0.18, 0.84, (f"{slope:+.2f} ns/m" if line else f"1/#lambda = {slope:+.2f} /m") + note); tx.SetNDC(); tx.SetTextSize(0.035)
        extra.append(tx)
    fin = y[np.isfinite(y)]; lo, hi = fin.min(), fin.max(); pad = 0.3 * (hi - lo) + 1e-3 * abs(hi) + 1e-6
    DrawHistos(G, [], XM.min() - 20, XM.max() + 20, "Table X [mm]", lo - pad, hi + 2.2 * pad, ytit, name, dology=False, drawoptions="P",
               markerstyles=M, mycolors=C, extraToDraw=extra, outdir=OUTD, lumitext=lab, usePDF=False)
    return name


def dists(name, vals, lab, xt):
    """normalised distributions of one channel at the picked positions"""
    allv = np.concatenate([v for v in vals if len(v)]) if any(len(v) for v in vals) else np.array([0.0])
    hi = np.percentile(allv, 99.5) if len(allv) > 10 else 1.0; lo = min(0.0, np.percentile(allv, 0.5))
    hs, lg = [], []
    for i, v in zip(PICK, vals):
        h = ROOT.TH1F(f"h_{name}_{i}", "", 30, lo, hi)
        for x in v: h.Fill(float(x))
        if h.Integral() > 0: h.Scale(1 / h.Integral())
        h.SetDirectory(0); hs.append(h); lg.append(f"X = {XM[i]:+.0f}: {len(v)}")
    DrawHistos(hs, lg, lo, hi, xt, 0, max(h.GetMaximum() for h in hs) * 1.6 + 1e-6, "Fraction", name, dology=False, drawoptions="HIST",
               legendoptions="L", mycolors=COLS5, legendPos=[0.50, 0.62, 0.93, 0.89], outdir=OUTD, lumitext=lab, usePDF=False)
    return name


def page(det, t, rows, heads, intro):
    names, cats = [], []
    for cat, ps in rows:
        names += [f"{p}.png" for p in ps if p] + ["NEWLINE"]; cats.append(cat)
    out = f"{HTML}/Gallery_{ds.name}_{det}_{t}.html"
    generate_html(names, OUTD, plots_per_row=6, output_html=out, title=f"Z scan {ds.name}, {det}, {TLAB[t]}: per-channel distributions and fits",
                  intro_text=intro)
    html = open(out).read(); parts = html.split('    <div class="grid">\n'); body = parts[0]
    for k, part in enumerate(parts[1:]):
        if k < len(cats) and (k == 0 or cats[k] != cats[k - 1]):
            body += f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{heads[cats[k]]}</h3>\n'
        body += '    <div class="grid">\n' + part
    open(out, "w").write(body)
    print(out, len(rows))


links = " &middot; ".join(f'<a href="Gallery_{ds.name}_{d}_{t}.html">{d} {TLAB[t]}</a>' for d in (["DRS", "FERS"] if J["fers"] else ["DRS"]) for t in TYPES)
back = f'<a href="ZScan_{ds.name}.html">Z scan {ds.name}</a>'
main_fit = J["main_fit"]
fit = lambda p, est: p["fit"][est]["main"][main_fit] if est in p["fit"] else (np.nan, np.nan, np.nan)
# ---------------- DRS ----------------
DC = D["channels"]; ev = {}
for i in PICK:
    r = RUNS[i]; L = load_drs(ds, r, MODE)
    ev[i] = {}
    for j, ch in enumerate(L["chans"]):
        if ch not in DC or not np.isfinite(DC[ch]["t0"][i]): ev[i][ch] = np.array([]); continue
        m = (np.abs(L["T"][:, j] - DC[ch]["t0"][i]) < ds.w_in(ch)) & (L["A"][:, j] > A_HIT[DC[ch]["is6mm"]])
        ev[i][ch] = L["E"][:, j][m]
cat_of = lambda p: "sat" if p["high_amp"]["1500"] >= SAT_FRAC else ("used" if fit(p, "lph")[2] < MAX_RESID else "resid")
CATS = ("used", "sat", "resid")
idx_main = [RUNS.index(r) for r in J["runsets"]["main"] if r in RUNS]
for t in TYPES:
    chans = sorted([c for c, p in DC.items() if p["type"] == t], key=lambda c: (CATS.index(cat_of(DC[c])), -np.nanmean(DC[c]["lpe"])))
    rows = []
    for c in chans:
        p = DC[c]; lab = short(c)
        cor = {}
        for est in ("lph", "lpe"):                      # shown after the fitted run offsets (main run set), as in the fit
            o = np.zeros(len(RUNS)); o[idx_main] = p["offsets"].get(est, [0.0] * len(idx_main)); cor[est] = np.exp(-o)
        ps = [dists(f"DRS_{c}_dist", [ev[i][c] for i in PICK], lab, "In-time CFD integral [ADC]"),
              vs_x(f"DRS_{c}_lph", np.array(p["lph"]) * cor["lph"], np.array(p["lph_err"]) * cor["lph"], "Light per hit [ADC]", lab, fit(p, "lph")[0],
                   note=f", res. {100 * fit(p, 'lph')[2]:.1f}%"),
              vs_x(f"DRS_{c}_hf", p["hf"], None, "Hit fraction", lab),
              vs_x(f"DRS_{c}_lpe", np.array(p["lpe"]) * cor["lpe"], np.array(p["lpe_err"]) * cor["lpe"], "Light per event [ADC]", lab, fit(p, "lpe")[0]),
              vs_x(f"DRS_{c}_all", p["all"], p["all_err"], "All-event mean [ADC] (no MCP)", lab, fit(p, "all")[0]) if "all" in p["fit"] else None,
              vs_x(f"DRS_{c}_time", p["time"], p["time_err"], "Arrival time [ns]", lab, p["time_slope"], line=True)]
        rows.append((cat_of(p), ps))
    n = {k: sum(1 for c in chans if cat_of(DC[c]) == k) for k in CATS}
    intro = (f"<p>{back} &middot; {links}</p><p>{len(chans)} channels with an in-time signal in every run: {n['used']} used, {n['sat']} saturated "
             f"(more than 1% of events above 1500 ADC peak), {n['resid']} with run-to-run residuals above 6%. Events with an MCP pulse. Row: in-time "
             f"integral distributions at up to five positions; light per hit, hit fraction, light per event (after the fitted run offsets, as in "
             f"the fit), all-event mean of all positrons (channels hit in more than half of the events), arrival time, vs table X. One colour per "
             f"beam-energy group; open grey: runs not in the main run set. Sorted by light per event.</p>")
    page("DRS", t, rows, {"used": "Used", "sat": "Saturated, not used", "resid": "Residuals above 6%, not averaged"}, intro)
# ---------------- FERS ----------------
if J["fers"]:
    FC = SEL["fers"]["channels"]; ev = {}
    for i in PICK:
        L = load_fers(ds, RUNS[i], MODE, False); ev[i] = {n: L["vals"][:, k] for k, n in enumerate(L["names"]) if n in FC}
    for t in TYPES:
        chans = sorted([n for n, p in FC.items() if p["type"] == t], key=lambda n: (fit(FC[n], "lpe")[2] >= MAX_RESID, -np.mean(FC[n]["signal"])))
        rows = []
        for n in chans:
            p = FC[n]; b, eb, res = fit(p, "lpe")
            rows.append(("used" if res < MAX_RESID else "resid",
                         [dists(f"FERS_{n}_dist", [ev[i][n] for i in PICK], f"{n} ({p['tx']:+.1f}, {p['ty']:+.2f})", "Calibrated FERS value"),
                          vs_x(f"FERS_{n}_mean", p["signal"], p["err"], "Mean per event", n, b, note=f", res. {100 * res:.1f}%")]))
        intro = (f"<p>{back} &middot; {links}</p><p>{len(chans)} FERS channels with a mean above 0.05 in every run (all positrons, "
                 f"{MODE} selection). Row: calibrated value per event at up to five positions; mean vs X with the fit. Sorted by mean signal.</p>")
        page("FERS", t, rows, {"used": "Residuals below 6%", "resid": "Residuals above 6% (not averaged)"}, intro)
