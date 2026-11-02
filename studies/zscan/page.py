"""Page of one z-scan dataset (from zscan_<dataset>.json): light speed, attenuation (FERS whole-module sums where FERS was recorded,
DRS sums over the unsaturated channels of a category), and the per-channel cross-checks in collapsed sections.
Test fibres (dataset testfibres): testfibres.page().
Usage: python3 studies/zscan/page.py DATASET     Page: results/html/ZScan/ZScan_<dataset>.html"""
import json, os, sys
import numpy as np
from common import DATA, TYPES, TLAB, SAT_FRAC, MAX_RESID, wmean
from datasets import PLOTS, HTML, DATASETS
if DATASETS[sys.argv[1]].get("channel_map"):
    import testfibres; testfibres.page(); sys.exit()
import ROOT
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html

NAME = sys.argv[1]
J = json.load(open(f"{DATA}/zscan_{NAME}.json")); MODE = J["main_selection"]; SEL = J["selections"][MODE]; D = SEL["drs"]
OUTD = f"{PLOTS}/{NAME}"; os.makedirs(OUTD, exist_ok=True); OUT = f"{HTML}/ZScan_{NAME}.html"
RUNS = D["runs"]; XM = np.array(D["x_m"]) * 1000; GRP = D["groups"]; G_NAMES = sorted(set(GRP)); MAIN = J["runsets"]["main"]
MF = J["main_fit"]; RS = list(J["runsets"])
COL = {"Sci": ROOT.kGreen + 2, "CerQuartz": ROOT.kBlue + 1, "CerPlastic": ROOT.kRed + 1}
MK = {"Sci": 20, "CerQuartz": 21, "CerPlastic": 22}; GMK = [20, 25, 22, 26]
SHORT = {"Sci": "scintillating", "CerQuartz": "quartz", "CerPlastic": "plastic"}
RSLAB = {"main": "main run set", "all": "all runs"}
FITLAB = {"plain": "no run offsets", "offsets": "with run offsets"}


def graph(xs, ys, es=None):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, float(x), float(y)); g.SetPointError(i, 0, float(es[i]) if es is not None else 0)
    return g


TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
fmt = lambda m, e: f"{m:+.3f} &plusmn; {e:.3f}" if np.isfinite(m) else "&ndash;"
def lam(m, e):
    if not np.isfinite(m): return "&ndash;"
    return f"{1 / m:.1f} ({1 / (m + e):.1f}&ndash;{1 / (m - e):.1f})" if m - e > 0 else (f"&gt; {1 / (m + e):.1f}" if m + e > 0 else "no attenuation seen")
details = lambda title, body: f'<details style="margin:10px 0"><summary><b>{title}</b></summary>{body}</details>'
timed = lambda p: np.isfinite(p["time_slope"]) and p["time_resid"] < 0.1
unsat = lambda p: p["high_amp"]["1500"] < SAT_FRAC
used = lambda p, est, rs="main", kind=MF: est in p["fit"] and unsat(p) and p["fit"][est][rs][kind][2] < MAX_RESID
main_plots, other_plots = [], []

# ---------------- light speed ----------------
G, L, Cc, M, extra = [], [], [], [], []
for t in TYPES:
    ps = [p for p in D["channels"].values() if p["type"] == t and timed(p)]
    if len(ps) < 3: continue
    sl = np.median([p["time_slope"] for p in ps])
    for k, g in enumerate(G_NAMES):
        m = np.array([gg == g for gg in GRP])
        d = np.nanmedian([np.array(p["time"], float)[m] - np.nanmean(np.array(p["time"], float)[m]) for p in ps], axis=0)
        G.append(graph(XM[m], d)); Cc.append(COL[t]); M.append(GMK[k % len(GMK)])
        L.append(f"{SHORT[t]} ({len(ps)} ch.): {sl:+.2f} ns/m" if k == 0 else f"{SHORT[t]}, {g}")
        fl = ROOT.TF1(f"ftx_{t}_{k}", f"{sl / 1000.0}*(x - ({XM[m].mean()}))", XM[m].min() - 10, XM[m].max() + 10)
        fl.SetLineColor(COL[t]); fl.SetLineWidth(2); extra.append((fl, "SAME"))
DrawHistos(G, L, XM.min() - 20, XM.max() + 20, "Table X [mm] (readout #rightarrow)", -1.5, 1.5, "Arrival time - channel mean [ns]", "time_vs_x",
           dology=False, drawoptions="P", legendoptions="P", markerstyles=M, mycolors=Cc, legendPos=[0.40, 0.66, 0.92, 0.89], extraToDraw=extra,
           outdir=OUTD, lumitext=f"{NAME}, median over channels", usePDF=False)
hs, lg, cl = [], [], []
for t in TYPES:
    v = [p["time_slope"] for p in D["channels"].values() if p["type"] == t and timed(p)]
    h = ROOT.TH1F(f"ht_{t}", "", 50, -8, -3)
    for x in v: h.Fill(x)
    h.SetDirectory(0); hs.append(h); cl.append(COL[t]); lg.append(f"{TLAB[t]} ({len(v)}): {np.median(v):.2f} ns/m" if v else TLAB[t])
DrawHistos(hs, lg, -8, -3, "d t_{arrival} / d X [ns/m]", 0, max(h.GetMaximum() for h in hs) * 1.6 + 1, "Channels", "time_slope", dology=False,
           drawoptions="HIST", legendoptions="L", mycolors=cl, legendPos=[0.17, 0.72, 0.85, 0.89], outdir=OUTD, lumitext=NAME, usePDF=True)
main_plots += [("Light speed", ["time_vs_x", "time_slope"])]
sp_rows = []
for t in TYPES:
    v = np.array([p["time_slope"] for p in D["channels"].values() if p["type"] == t and timed(p)])
    if not len(v): continue
    m = np.median(v); s = 1 / abs(m) / 0.2998
    sp_rows.append([TLAB[t], len(v), f"{m:.2f} ({np.percentile(v, 16):.2f}..{np.percentile(v, 84):.2f})", f"{s:.3f} c", f"{1000 / abs(m):.0f}", f"{1 / s:.2f}"])
speed_tab = table(["Fibres", "ch.", "dt/dX [ns/m] (16&ndash;84%)", "speed", "[mm/ns]", "n<sub>eff</sub> = c/v"], sp_rows)

# ---------------- DRS sums ----------------
SU = D["sums"]["main"]["categories"]
gs, labs, cols, extra = [], [], [], []
PAL = [ROOT.kGreen + 2, ROOT.kBlue + 1, ROOT.kAzure + 1, ROOT.kRed + 1, ROOT.kOrange + 7, ROOT.kBlack]
idx_main = [RUNS.index(r) for r in MAIN if r in RUNS]; gi = np.array([G_NAMES.index(GRP[i]) for i in idx_main])
for ic, (name, c) in enumerate(SU.items()):
    f = c["fits"]["base"].get("adc")
    if not f: continue
    y = np.array(c["values"]["base"]["adc"]) / np.exp(np.array(f["coef"])[gi]); ey = np.array(c["values"]["base"]["adc_err"]) / np.exp(np.array(f["coef"])[gi])
    gs.append(graph(XM[idx_main], y, ey)); labs.append(f"{name}: {f['inv_lambda']:+.3f}#pm{f['err']:.3f}/m"); cols.append(PAL[ic % len(PAL)])
    fl = ROOT.TF1(f"fs_{ic}", f"exp({f['inv_lambda']}*x/1000.)", XM.min(), XM.max()); fl.SetLineColor(PAL[ic % len(PAL)]); fl.SetLineWidth(2)
    extra.append((fl, "SAME"))
DrawHistos(gs, labs, XM.min() - 20, XM.max() + 20, "Table X [mm]", 0.8, 1.2, "Light per event / fit at X = 0", "sums_drs", dology=False, drawoptions="P",
           legendoptions="P", markerstyles=[20] * len(gs), mycolors=cols, legendPos=[0.18, 0.66, 0.75, 0.92], extraToDraw=extra, outdir=OUTD,
           lumitext=f"{NAME}, DRS sums (main run set)", usePDF=False)
for name, c in SU.items():
    f = c["fits"]["base"].get("pe")
    if not f: continue
    y, ey = np.array(c["values"]["base"]["pe"]), np.array(c["values"]["base"]["pe_err"])
    gs, labs, cl_, ex = [], [], [], []
    for k, g in enumerate(G_NAMES):
        m = gi == k
        if not m.any(): continue
        gs.append(graph(XM[idx_main][m], y[m], ey[m])); cl_.append([ROOT.kBlue + 1, ROOT.kRed + 1][k % 2])
        labs.append(g + (f": 1/#lambda = {f['inv_lambda']:+.3f}#pm{f['err']:.3f}/m" if k == 0 else ""))
        fl = ROOT.TF1(f"fpe_{name}_{k}".replace(" ", "_"), f"exp({f['coef'][k]} + {f['inv_lambda']}*x/1000.)", XM[idx_main][m].min(), XM[idx_main][m].max())
        fl.SetLineColor(cl_[-1]); fl.SetLineWidth(2); ex.append((fl, "SAME"))
    lo, hi = (y - ey).min(), (y + ey).max(); pad = 0.25 * (hi - lo)
    nm = f"sums_pe_{name.replace(' ', '_')}"
    DrawHistos(gs, labs, XM.min() - 20, XM.max() + 20, "Table X [mm]", lo - pad, hi + 1.5 * pad, "Mean p.e. per event", nm, dology=False,
               drawoptions="P", legendoptions="P", markerstyles=[20] * len(gs), mycolors=cl_, legendPos=[0.18, 0.76, 0.85, 0.89], extraToDraw=ex,
               outdir=OUTD, lumitext=f"{NAME}, {name}", usePDF=False)
    other_plots.append(("Sums in photoelectrons (3 mm channels)", nm))
CHK = ("sbpre", "thr6_15", "thr6_20", "thr6_30")
rows = []
for name, c in SU.items():
    f = c["fits"]["base"].get("adc")
    if not f: continue
    m, e = f["common"]["inv_lambda"], f["common"]["err"]
    alt = [c["fits"][k]["adc"]["common"]["inv_lambda"] for k in CHK if "adc" in c["fits"].get(k, {})]
    rows.append([name, len(c["channels"]), f"{m:+.3f} &plusmn; {e:.3f} ({f['common']['chi2ndf']:.1f})", lam(m, e), f"{max(abs(a - m) for a in alt):.3f}" if alt else "&ndash;"])
drs_sum_tab = table(["Sum (DRS)", "ch.", "1/&lambda; [/m] (&chi;<sup>2</sup>/ndf), common run offsets", "&lambda; [m]", "largest change in the checks [/m]"], rows)
rows = []
for rs in RS:
    for name, c in D["sums"][rs]["categories"].items():
        f = c["fits"]["base"].get("adc")
        if not f: continue
        rows.append([RSLAB.get(rs, rs) if name == next(iter(D["sums"][rs]["categories"])) else "", name, len(c["channels"]),
                     f"{100 * f['stat']:.1f} / {100 * f['resid_rms']:.1f}", f"{fmt(f['inv_lambda'], f['err'])} ({f['chi2ndf']:.1f})",
                     fmt(f["drift_corr"]["inv_lambda"], f["drift_corr"]["err"]), f"{fmt(f['common']['inv_lambda'], f['common']['err'])} ({f['common']['chi2ndf']:.1f})",
                     f"{f['drift_corr']['common']:+.3f}"])
all_sums_tab = table(["Run set", "Sum", "ch.", "stat. / scatter per run [%]", "1/&lambda; [/m] (&chi;<sup>2</sup>/ndf)", "drift-corrected",
                      "common run offsets", "common offsets + drift"], rows)
pe_rows = [[name, len(c["channels_pe"]), f"{c['spe_amp']:.0f} / {c['spe_int']:.0f}", f"{np.median(c['values']['base']['pe']):.1f}",
            f"{fmt(c['fits']['base']['pe']['inv_lambda'], c['fits']['base']['pe']['err'])} ({c['fits']['base']['pe']['chi2ndf']:.1f})"]
           for name, c in SU.items() if c["fits"]["base"].get("pe")]
pe_tab = table(["Sum (3 mm)", "ch.", "single p.e. height / integral [ADC]", "p.e. per event", "1/&lambda; [/m] (&chi;<sup>2</sup>/ndf)"], pe_rows)

# ---------------- per channel (cross-check) ----------------
for est, lab_ in (("lph", "light per hit"), ("lpe", "light per event")):
    hs, lg, cl = [], [], []
    for t in TYPES:
        v = [p["fit"][est]["main"][MF][0] for p in D["channels"].values() if p["type"] == t and used(p, est)]
        h = ROOT.TH1F(f"h_{est}_{t}", "", 32, -0.8, 0.8)
        for x in v: h.Fill(x)
        h.SetDirectory(0); hs.append(h); cl.append(COL[t]); lg.append(f"{TLAB[t]} ({len(v)}): median {np.median(v) if v else np.nan:+.2f} /m")
    DrawHistos(hs, lg, -0.8, 0.8, f"1 / #lambda [m^{{-1}}] ({lab_})", 0, max(h.GetMaximum() for h in hs) * 1.6 + 1, "Channels", f"inv_lambda_{est}",
               dology=False, drawoptions="HIST", legendoptions="L", mycolors=cl, legendPos=[0.17, 0.72, 0.85, 0.89], outdir=OUTD, lumitext=NAME, usePDF=False)
    other_plots.append(("Per channel (cross-check)", f"inv_lambda_{est}"))
for key, xr, xt, nm, logx in (("ty", (-3, 9), "Tower row t_{y} (vertical)", "lph_vs_ty", False),
                              ("amp", (1e-4, 1.0), "Fraction of events above 1500 ADC", "lph_vs_amp", True)):
    G, L, Cc = [], [], []
    for t in TYPES:
        ps = [p for p in D["channels"].values() if p["type"] == t and p["fit"]["lph"]["main"][MF][2] < MAX_RESID and (key == "amp" or unsat(p))]
        xv = (lambda p: p["ty"]) if key == "ty" else (lambda p: max(p["high_amp"]["1500"], 1.5e-4))
        if ps: G.append(graph([xv(p) for p in ps], [p["fit"]["lph"]["main"][MF][0] for p in ps], [p["fit"]["lph"]["main"][MF][1] for p in ps])); L.append(TLAB[t]); Cc.append(COL[t])
    DrawHistos(G, L, xr[0], xr[1], xt, -1.0, 1.0, "1 / #lambda [m^{-1}] (light per hit)", nm, dology=False, dologx=logx, drawoptions="P", legendoptions="P",
               markerstyles=[20] * len(G), mycolors=Cc, legendPos=[0.17, 0.75, 0.60, 0.89], outdir=OUTD, lumitext=NAME, usePDF=False)
    other_plots.append(("Per channel (cross-check)", nm))
rows = []
for est, elab in (("lph", "light per hit"), ("lpe", "light per event"), ("all", "all-event mean (bright ch.)")):
    for t in TYPES:
        row = [elab if t == "Sci" else "", TLAB[t]]
        for rs in RS:
            for kind in ("plain", "offsets"):
                s = D["summary"].get(est, {}).get(rs, {}).get(kind, {}).get(t)
                row.append(f"{fmt(s['inv_lambda'], s['err'])} ({s['n']})" if s and s["n"] else "&ndash;")
        rows.append(row)
per_tab = table(["Estimator", "Fibres"] + [f"{RSLAB.get(rs, rs)}, {FITLAB[k]}" for rs in RS for k in ("plain", "offsets")], rows)
SCUTS = [("none", lambda h: True), ("&lt;1% above 1500 (used)", lambda h: h["1500"] < 0.01), ("&lt;0.1% above 1500", lambda h: h["1500"] < 0.001),
         ("&lt;1% above 1000", lambda h: h["1000"] < 0.01), ("saturated only", lambda h: h["1500"] >= 0.01)]
rows = []
for t in TYPES:
    allp = [p for p in D["channels"].values() if p["type"] == t]
    sat = [p for p in allp if not unsat(p)]
    row = [TLAB[t], f"{len(sat)} / {len(allp)}"]
    for _, f in SCUTS:
        ps = [p for p in allp if f(p["high_amp"]) and p["fit"]["lph"]["main"][MF][2] < MAX_RESID]
        m, e, _ = wmean([p["fit"]["lph"]["main"][MF][0] for p in ps], [p["fit"]["lph"]["main"][MF][1] for p in ps])
        row.append(f"{fmt(m, e)} ({len(ps)})" if ps else "&ndash;")
    rows.append(row)
sat_tab = table(["Fibres", "saturated / all ch."] + [f"light per hit 1/&lambda;, cut: {c}" for c, _ in SCUTS], rows)
rows = []
for est, elab in (("lph", "light per hit"), ("lpe", "light per event")):
    for rs in RS:
        dr = D["drift"].get(est, {}).get(rs)
        if not dr: continue
        rows.append([elab, RSLAB.get(rs, rs), fmt(dr["ky"], dr["eky"])] + [fmt(dr["types"][t]["inv_lambda"], dr["types"][t]["err"]) if t in dr["types"] else "&ndash;" for t in TYPES])
drift_tab = table(["Estimator", "Run set", "k<sub>y</sub> [hodoscope cm/m]"] + [f"1/&lambda; {TLAB[t]}" for t in TYPES], rows)

# ---------------- FERS ----------------
fers_tab = fers_details = ""
if J["fers"]:
    F, FM = SEL["fers"], SEL["fers_mcp"]
    rows = []
    for t in TYPES:
        a = F["sums"]["main"][t]["fit"]
        rows.append([TLAB[t], f"<b>{fmt(a['inv_lambda'], a['err'])}</b> ({a['chi2ndf']:.1f})", f"<b>{lam(a['inv_lambda'], a['err'])}</b>"]
                    + [fmt(X_["sums"][rs][t]["fit"]["inv_lambda"], X_["sums"][rs][t]["fit"]["err"]) for X_, rs in ((F, "all"), (FM, "main"), (J["selections"]["spot"]["fers"], "main"))])
    fers_tab = table(["Fibres", "FERS whole-module sum 1/&lambda; [/m] (&chi;<sup>2</sup>/ndf)", "&lambda; [m]", "all runs", "MCP events only", "without position matching"], rows)
    # FERS sums vs X and the shower centroid
    G, L_, Cc, M, extra = [], [], [], [], []
    for t in TYPES:
        s = F["sums"]["main"][t]; f = s["fit"]; v = np.array(F["sums"]["all"][t]["values"]); e = np.array(F["sums"]["all"][t]["errs"])
        norm = np.exp(f["coef"][0] + f["inv_lambda"] * XM[idx_main].mean() / 1000.0)        # the fit at the mean position of the fitted runs
        inm = np.array([r in MAIN for r in RUNS])
        G.append(graph(XM[inm], (v / norm)[inm], (e / norm)[inm])); L_.append(f"{TLAB[t]}: 1/#lambda = {f['inv_lambda']:+.2f} #pm {f['err']:.2f} /m")
        Cc.append(COL[t]); M.append(MK[t])
        if (~inm).any():
            g1 = graph(XM[~inm], (v / norm)[~inm]); g1.SetMarkerStyle(MK[t]); g1.SetMarkerColor(ROOT.kGray + 1); extra.append((g1, "P SAME"))
        fl = ROOT.TF1(f"ff_{t}", f"exp({f['inv_lambda']}*(x-({XM[idx_main].mean()}))/1000.)", XM.min() - 10, XM.max() + 10)
        fl.SetLineColor(COL[t]); fl.SetLineWidth(2); extra.append((fl, "SAME"))
    DrawHistos(G, L_, XM.min() - 30, XM.max() + 30, "Table X [mm] (readout #rightarrow)", 0.8, 1.2, "Whole-module sum / fit", "sums_fers", dology=False,
               drawoptions="EP", legendoptions="EP", markerstyles=M, mycolors=Cc, legendPos=[0.22, 0.18, 0.72, 0.36], extraToDraw=extra, outdir=OUTD,
               lumitext=f"{NAME}, FERS (grey: not in the fit)", usePDF=False)
    Gc = [graph(XM, J["selections"]["spot"]["fers"]["centroid"]), graph(XM, F["centroid"])]
    DrawHistos(Gc, ["beam spot around the run median", "position-matched"], XM.min() - 30, XM.max() + 30, "Table X [mm] (readout #rightarrow)", -0.8, -0.2,
               "Shower centroid [tower rows]", "centroid_vs_x", dology=False, drawoptions="P", legendoptions="P", markerstyles=[20, 21],
               mycolors=[ROOT.kBlack, ROOT.kRed + 1], legendPos=[0.20, 0.75, 0.70, 0.89], outdir=OUTD, lumitext="FERS, scintillation-weighted, median", usePDF=False)
    main_plots += [("Attenuation: FERS whole-module sums, shower position in the module", ["sums_fers", "centroid_vs_x"])]
    hs, lg, cl = [], [], []
    for t in TYPES:
        v = [p["fit"]["lpe"]["main"]["plain"][0] for p in F["channels"].values() if p["type"] == t and p["fit"]["lpe"]["main"]["plain"][2] < MAX_RESID]
        h = ROOT.TH1F(f"hf_{t}", "", 40, -1, 1)
        for x in v: h.Fill(x)
        h.SetDirectory(0); hs.append(h); cl.append(COL[t]); lg.append(f"{TLAB[t]} ({len(v)}): median {np.median(v):+.2f} /m")
    DrawHistos(hs, lg, -1, 1, "1 / #lambda [m^{-1}] (FERS channel)", 0, max(h.GetMaximum() for h in hs) * 1.6 + 1, "Channels", "inv_lambda_fers", dology=False,
               drawoptions="HIST", legendoptions="L", mycolors=cl, legendPos=[0.17, 0.72, 0.85, 0.89], outdir=OUTD, lumitext=NAME, usePDF=False)
    other_plots.append(("Per channel (cross-check)", "inv_lambda_fers"))
    brows = []
    for t in TYPES:
        P_ = [p for p in F["channels"].values() if p["type"] == t and p["fit"]["lpe"]["main"]["plain"][2] < MAX_RESID]
        sg = np.array([np.mean(p["signal"]) for p in P_]); row = [TLAB[t]]
        for lo, hi in ((0.05, 0.15), (0.15, 0.4), (0.4, 1.0), (1.0, 99)):
            ps = [p for p, s_ in zip(P_, sg) if lo <= s_ < hi]
            m, e, _ = wmean([p["fit"]["lpe"]["main"]["plain"][0] for p in ps], [p["fit"]["lpe"]["main"]["plain"][1] for p in ps])
            row.append(f"{fmt(m, e)} ({len(ps)})" if ps else "&ndash;")
        brows.append(row)
    S_ = F["summary"]["lpe"]
    frows = [[TLAB[t]] + [f"{fmt(S_[rs][k][t]['inv_lambda'], S_[rs][k][t]['err'])} ({S_[rs][k][t]['n']})" for rs in RS for k in ("plain", "offsets")] for t in TYPES]
    M_ = J["matching"]
    fers_details = details("FERS: per channel, brightness, drift",
        f"<p>Shower drift: {M_['centroid_slope_rows_per_m']:+.2f} tower rows per m of table travel = {M_['hodo_cm_per_m']:+.2f} hodoscope cm per m "
        f"({M_['rows_per_hodo_cm']:.2f} rows per hodoscope cm). Raw LG below {F['lg_max']:.0f} of 8191 (99.9% of events): no FERS saturation.</p>"
        + table(["Fibres"] + [f"{RSLAB.get(rs, rs)}, {FITLAB[k]}" for rs in RS for k in ("plain", "offsets")], frows)
        + "<p>Per-channel slope vs channel brightness (mean signal, calibrated units):</p>"
        + table(["Fibres", "0.05&ndash;0.15", "0.15&ndash;0.4", "0.4&ndash;1", "&gt; 1 (shower core)"], brows)
        + table(["Run set", "k<sub>y</sub> [cm/m]"] + [TLAB[t] for t in TYPES],
                [[RSLAB.get(rs, rs), fmt(F["drift"][rs]["ky"], F["drift"][rs]["eky"])] + [fmt(F["drift"][rs]["types"][t]["inv_lambda"], F["drift"][rs]["types"][t]["err"])
                  if t in F["drift"][rs]["types"] else "&ndash;" for t in TYPES] for rs in RS if F["drift"].get(rs)]))
main_plots += [("Attenuation: DRS sums (part of the module)", ["sums_drs"])]

# ---------------- runs, checks ----------------
pc = D["position_check"]
run_rows = []
for i, r in enumerate(RUNS):
    row = [r, f"{XM[i]:+.1f}", GRP[i], "" if r in MAIN else "not in main set", D["n_sel"][i], f"{100 * D['n_sel'][i] / max(D['n_positrons'][i], 1):.0f}%",
           f"{pc['dx_mm'][i]:+.1f}" if pc else "&ndash;"]
    if J["fers"]:
        row += [f"{J['selections']['spot']['fers']['centroid'][i]:+.3f}", f"{J['matching']['shift_cm'][i]:+.2f}", f"{SEL['fers']['centroid'][i]:+.3f}"]
    run_rows.append(row)
run_tab = table(["Run", "X [mm]", "group", "", "DRS events (MCP)", "MCP fraction", "X from arrival time &minus; logged [mm]"]
                + (["centroid, beam spot [rows]", "hodoscope shift [cm]", "centroid, matched [rows]"] if J["fers"] else []), run_rows)
checks = ""
VAR = f"{DATA}/checks/variants_{NAME}.json"
if os.path.exists(VAR):
    VJ = json.load(open(VAR)); V = VJ["variants"]; rows = []
    for vname, chs in V.items():
        for t in TYPES:
            ps = [q for q in chs.values() if q["type"] == t and q["resid"] < 0.06]
            m, e, _ = wmean([q["inv_lambda"] for q in ps], [q["err"] for q in ps])
            rows.append([vname if t == "Sci" else "", TLAB[t], f"{100 * np.median([q['stat'] for q in ps]):.1f}", fmt(m, e),
                         f"{np.median([q['hit_fraction'] for q in ps]):.2f}"])
    checks += details("Check: integral variants (checks/variants.py" + (", before the MCP selection)" if VJ.get("selection") == "pre-MCP" else ")"),
        "<p>Short integral (3 ns) holds only ~40% of a scintillation pulse and is not proportional to the light: long integral for "
        "scintillating. Linear baseline and in-time search change nothing.</p>" + table(["Integral", "Fibres", "stat./point [%]", "1/&lambda; [/m]", "hit frac."], rows))
SH = f"{DATA}/checks/integral_shape_plots_{NAME}.json"
if os.path.exists(SH):
    for n in json.load(open(SH)): other_plots.append(("All-event integral distributions at five positions (checks/integral_shapes.py)", f"checks/{n}"))

# ---------------- page ----------------
other = [d for d in DATASETS if d != NAME]
gal = " &middot; ".join(f'<a href="Gallery_{NAME}_{d}_{t}.html">{d} {SHORT[t]}</a>' for d in (["DRS", "FERS"] if J["fers"] else ["DRS"]) for t in TYPES)
nav = " &middot; ".join([f'<a href="ZScan_overview.html">overview</a>'] + [f'<a href="ZScan_{d}.html">{d}</a>' for d in other])
fers_text = (f"<p><b>FERS whole-module sums</b> (all positrons, position-matched, main run set): the attenuation measurement. They do not "
             f"depend on the MCP requirement.</p>{fers_tab}") if J["fers"] else \
            "<p>No calorimeter FERS in this dataset: no whole-module sum.</p>"
intro = f"""
<p>{J['desc']}</p>
<p>Galleries: {gal}. Other pages: {nav}.</p>
<h3 style="margin:14px 0 4px 0">Light speed (DRS arrival time vs X, median over well-timed channels)</h3>
{speed_tab}
<p>The arrival time depends only on how far the light travels to the readout, not on how much of the shower a channel sees: it is the robust
z-scan measurement. n<sub>eff</sub> = c/v is an effective index (light reaching the readout travels at an angle to the fibre axis; a pulse
moves at the group velocity).</p>
<h3 style="margin:14px 0 4px 0">Attenuation</h3>
{fers_text}
<p><b>DRS sums</b> of light per event over the unsaturated channels of a category: a part of the module, so they depend on where the shower
sits; with the fibres tilted relative to the table travel the shower moves across the channels during the scan. Checks: randoms from the
sideband before the pulse only; 6 mm hit threshold 15 / 20 / 30 ADC instead of 25.</p>
{drs_sum_tab}
{details("Runs", run_tab)}
{details("Method", f'''<p><b>Events</b>: positrons (beam-PID flags of the period), 1&ndash;3 hodoscope fibres per axis, beam within 1 cm of the run
median ({"vertically: of a common centre after shifting each run by the measured drift of the shower in the module (position matching)" if MODE == "matched" else "both axes"}).
DRS: a pulse in {J['mcp_select']} above its configured cut (in 2026 the test-fibre <code>mcp_clean</code>) and a valid time reference
(MCP_DS_1).</p>
<p><b>DRS signal</b>: CFD integral, long (10 ns) for scintillating, short (3 ns) for Cherenkov. Hits: in time (&plusmn;4 ns scintillating,
&plusmn;1.2 ns Cherenkov) and peak above 40 ADC (3 mm) or 25 ADC (6 mm); randoms from sidebands subtracted. Light per event: events without a
hit count as zero; light per hit leaves them out and gives too flat a slope for few-p.e. channels. Channels with more than 1% of their
events above 1500 ADC peak (DRS clipping) are dropped.</p>
<p><b>Fits</b>: ln(signal) = a<sub>group</sub> + X/&lambda; (one normalisation per beam-energy group); run offsets
{"per DRS chip" if J["offsets"] == "chip" else "common to all channels"}, kept free of a constant and a linear trend in X. Sums: common
run offsets of all sums. Drift term: 1/&lambda;<sub>ch</sub> = 1/&lambda;<sub>type</sub> + k<sub>y</sub> g<sub>y,ch</sub>.</p>''')}
{details("All DRS sums (plain, drift-corrected, common run offsets) and photoelectrons", all_sums_tab + pe_tab)}
{details("DRS per channel (cross-check), saturation, lateral drift", per_tab + sat_tab + drift_tab)}
{fers_details}
{checks}
"""
pl = sum(([f"{n}.png" for n in g] + ["NEWLINE"] for _, g in main_plots), [])[:-1]
generate_html(pl, OUTD, plots_per_row=4, output_html=OUT, title=f"Z scan {NAME}: light speed and attenuation, {J['title']}", intro_text=intro)
html = open(OUT).read(); parts = html.split('    <div class="grid">\n')
html = parts[0] + "".join(f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{t}</h3>\n    <div class="grid">\n' + part
                          for (t, _), part in zip(main_plots, parts[1:]))
rel = os.path.relpath(OUTD, HTML)
more, last = "", None
for t, n in other_plots:
    if t != last: more += f"<h4>{t}</h4>"; last = t
    more += f'<a href="{rel}/{n}.png" target="_blank"><img src="{rel}/{n}.png" style="width:24%;min-width:220px"></a>'
html = html.replace("</body>", f'<details style="margin:20px"><summary><b>All other plots</b></summary>{more}</details>\n</body>')
open(OUT, "w").write(html)
print(OUT)
