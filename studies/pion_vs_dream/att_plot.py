"""Attenuation along the fibres from the tb2025 90-degree positron runs at Y ~ 1064 mm (1501, 1507, 1511, 1513, 15130):
type-summed FERS and DRS signals vs table X, fit without 1501 (shower in different towers). Sci, quartz, plastic separately."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json, os, sys
import numpy as np
REPO = REPO_DIR; SCR = DATA_DIR
os.chdir(REPO); sys.path.insert(0, REPO); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import ROOT
import plot_drs_mpv_scan as P
from plotting.my_function import DrawHistos
OUTD = "results/plots/FERS_PionResponse"
rl = json.load(open("data/Runlist.json"))
RUNS = [1513, 1501, 1507, 1511, 15130]; EXCL = 1501
X = {r: float(rl[str(r)]["Table position X"].replace("mm", "")) for r in RUNS}
EL = ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
TYPES = [("Sci", "scintillation", ROOT.kGreen + 2, 20), ("CerQuartz", "Cherenkov quartz", ROOT.kBlue + 1, 21), ("CerPlastic", "Cherenkov plastic", ROOT.kRed + 1, 22)]

def fers(run):
    d = np.load(f"{SCR}/scan90ch_{run}.npz"); sel = np.ones(len(d["vals"]), bool)
    for f in EL: sel &= d[f].astype(bool)
    v, m = d["vals"][sel], d["meta"]; mean = v.mean(0); ev = v.std(0) / np.sqrt(len(v))
    masks = {"Sci": m[:, 2] == 0, "CerQuartz": (m[:, 2] == 1) & (m[:, 3] == 1), "CerPlastic": (m[:, 2] == 1) & (m[:, 3] == 0)}
    out = {}
    for t, k in masks.items():
        s = v[:, k].sum(1); out[t] = (s.mean(), s.std() / np.sqrt(len(s)))
    return out
order = {o[0]: o[4] for o in P._build_channel_order(P.build_drs_boards(run_number=1499))}
def drs(run):
    path = f"{SCR}/mirror/results/root/Run{run}" if run >= 10000 else f"results/root/Run{run}"
    st = json.load(open(f"{path}/drs_energy_stats.json"))
    return {t: (sum(v[0] for ch, v in st.items() if order.get(ch) == t), 0.0) for t, *_ in TYPES}

def fit(x, y):
    """ln y = a + b x; slope error from the residual scatter (2 dof for 4 points)."""
    A = np.column_stack([np.ones_like(x), x]); c, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c; s2 = (r ** 2).sum() / max(len(x) - 2, 1); cov = s2 * np.linalg.inv(A.T @ A)
    return c, np.sqrt(cov[1, 1]), r

results = {}
for src, getter in (("FERS", fers), ("DRS", drs)):
    vals = {r: getter(r) for r in RUNS}
    graphs, gex, legs, cols, mks, extra = [], [], [], [], [], []
    for t, lab, col, mk in TYPES:
        use = [r for r in RUNS if r != EXCL]
        x = np.array([X[r] for r in use]) / 1000; y = np.log([vals[r][t][0] for r in use])
        c, eb, res = fit(x, y)
        lam, elam = 1 / c[1], eb / c[1] ** 2
        results.setdefault(src, {})[t] = dict(slope=c[1], eslope=eb, lam=lam, elam=elam, res_pct=list(100 * res), runs=use,
                                              r1501_pct=100 * (np.log(vals[EXCL][t][0]) - (c[0] + c[1] * X[EXCL] / 1000)))
        norm = np.exp(c[0] + c[1] * (-0.218))           # normalise to the fit at X = -218 mm
        g = ROOT.TGraphErrors(len(use)); g1 = ROOT.TGraphErrors(1)
        for i, r in enumerate(use):
            g.SetPoint(i, X[r], vals[r][t][0] / norm); g.SetPointError(i, 0, vals[r][t][1] / norm)
        g1.SetPoint(0, X[EXCL], vals[EXCL][t][0] / norm); g1.SetMarkerStyle(mk); g1.SetMarkerColor(ROOT.kGray + 1); g1.SetMarkerSize(1.3)
        f = ROOT.TF1(f"f_{src}_{t}", f"exp({c[1]}*(x+218)/1000.)", -420, -30); f.SetLineColor(col); f.SetLineWidth(2)
        graphs.append(g); legs.append(f"{lab}: #lambda = {lam:.1f} #pm {elam:.1f} m"); cols.append(col); mks.append(mk)
        extra += [(f, "SAME"), (g1, "P SAME")]
    t = ROOT.TLatex(0.18, 0.86, "grey: run 1501 (shower in other towers), not in the fit"); t.SetNDC(); t.SetTextSize(0.028)
    extra.append(t)
    name = f"att90_Y1064_{src}"
    DrawHistos(graphs, legs, -430, -20, "Table X [mm] (readout #rightarrow)", 0.80 if src == "FERS" else 0.75,
               1.15 if src == "FERS" else 1.20, "Type-summed signal / fit at X = -218 mm", name, dology=False, drawoptions="EP",
               legendoptions="EP", markerstyles=mks, mycolors=cols, legendPos=[0.30, 0.18, 0.80, 0.36], extraToDraw=extra,
               outdir=OUTD, lumitext=f"e^{{+}} 40 GeV, 90#circ, Y #approx 1064 mm, {src}", usePDF=True)
    print(src, {t: f"lambda {v['lam']:.2f} +- {v['elam']:.2f} m, res% {np.round(v['res_pct'], 2)}, 1501 off fit {v['r1501_pct']:+.1f}%" for t, v in results[src].items()})
json.dump(results, open(f"{SCR}/att90_Y1064.json", "w"), indent=1, default=float)
