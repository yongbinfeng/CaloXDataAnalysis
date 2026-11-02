"""Plots for the shower-position (depth) dependence of S and C and its correction, 80/100/120 GeV pions."""
import json
import os
import sys

import numpy as np

sys.path.insert(0, "/data1/yofeng/CaloX/CaloXDataAnalysis")
os.chdir("/data1/yofeng/CaloX/CaloXDataAnalysis")
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
OUT = "results/plots/FERS_PionResponse"
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"]
ALL_RUNS = [(1441, ROOT.kGray + 2), (1439, ROOT.kOrange + 7), (1438, ROOT.kMagenta + 1), (1437, ROOT.kCyan + 2),
            (1429, ROOT.kRed + 1), (1432, ROOT.kBlue + 1), (1433, ROOT.kGreen + 2), (1434, ROOT.kViolet + 1)]
MIN_EVENTS = 300        # the position fit needs a few hundred events


def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6:
            break
        m, s = m2, s2
    return m, s


def qs_energy(S, C, SO):
    Sc = S + SO
    b, a = 1 / 1.3, 0.2
    f = np.clip((a - C / Sc * b) / (C / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    return Sc / (f + (1 - f) * b)


def load(run):
    d = np.load(f"{SCR}/tail_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION:
        sel &= d[f].astype(bool)
    S = d["FERS_energyMix_Sci_calib_sum"][sel]; C = d["FERS_energyMix_Cer_calib_sum"][sel]
    SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    X = d["FERS_energyMix_Sci_calib_weighted_center_X"][sel]; Y = d["FERS_energyMix_Sci_calib_weighted_center_Y"][sel]
    k = S > 0.3 * E
    S, C, SO, X, Y = S[k], C[k], SO[k], X[k], Y[k]
    dx, dy = X - np.median(X), Y - np.median(Y)

    def fit(sig):
        y = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dx) < 3) & (np.abs(dy) < 3)
        A = np.column_stack([np.ones(len(y)), dx, dy])
        for _ in range(5):
            coef, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            res = y - A @ coef; s = 1.4826 * np.median(np.abs(res[m] - np.median(res[m])))
            m = m & (np.abs(res) < 2.5 * s)
        return coef
    cS, cC = fit(S), fit(C)
    wS = np.exp(-(cS[1] * dx + cS[2] * dy)); wC = np.exp(-(cC[1] * dx + cC[2] * dy))
    raw = {"S": S, "C": C, "DR": (S - 0.29 * C) / 0.71, "QS": qs_energy(S, C, SO)}
    cor = {"S": S * wS, "C": C * wC, "DR": (S * wS - 0.29 * C * wC) / 0.71, "QS": qs_energy(S * wS, C * wC, SO * wS)}
    return dict(E=E, dx=dx, dy=dy, raw=raw, cor=cor, slopes=dict(Sx=cS[1], Sy=cS[2], Cx=cC[1], Cy=cC[2]))


data, RUNS = {}, []
for run, col in ALL_RUNS:
    if not os.path.exists(f"{SCR}/tail_{run}.npz"):
        print(f"run {run}: no dump, skipped"); continue
    dd = load(run)
    if len(dd["raw"]["S"]) < MIN_EVENTS:
        print(f"run {run} ({dd['E']:.0f} GeV): {len(dd['raw']['S'])} events, below {MIN_EVENTS}; skipped"); continue
    data[run] = dd; RUNS.append((run, col))
names = []

# 1. position profiles: median S/E and C/E vs centre offset in X and in Y
edges = np.arange(-3.0, 3.01, 0.5)
for axis in ("dx", "dy"):
  for key, mk in (("S", 20), ("C", 21)):
        graphs, legs, cols, mks = [], [], [], []
        for run, col in RUNS:
            dd = data[run]; v = dd["raw"][key] / dd["E"]; off = dd[axis]
            pts = [(0.5 * (a + b), np.median(v[(off >= a) & (off < b)]), np.std(v[(off >= a) & (off < b)]) / np.sqrt(((off >= a) & (off < b)).sum()) * 1.25)
                   for a, b in zip(edges[:-1], edges[1:]) if ((off >= a) & (off < b)).sum() >= 25]
            g = ROOT.TGraphErrors(len(pts))
            for i, (x, y, e) in enumerate(pts):
                g.SetPoint(i, x, y); g.SetPointError(i, 0.0, e)
            graphs.append(g); legs.append(f"{dd['E']:.0f} GeV"); cols.append(col); mks.append(mk)
        nm = f"poscorr_profile_{key}_{'X' if axis == 'dx' else 'Y'}"
        ylo, yhi = (0.6, 1.35) if key == "S" else (0.4, 1.0)
        DrawHistos(graphs, legs, -3.2, 3.2, f"Shower centre {'X' if axis == 'dx' else 'Y'} - median [cm]",
                   ylo, yhi, f"Median {key} / beam energy", nm, dology=False, drawoptions="EP", legendoptions="EP",
                   markerstyles=mks, mycolors=cols, legendPos=[0.18, 0.70, 0.92, 0.88], legendNCols=4, outdir=OUT,
                   lumitext="#pi^{+}, pion-selected, hole veto", usePDF=True)
        names.append(nm)

# 2. distributions before / after the correction, S and Q/S, in GeV
for run, _ in RUNS:
    dd = data[run]; E = dd["E"]; nb = 80
    hs, legs = [], []
    for key, lab in (("S", "S"), ("DR", "DR"), ("QS", "Q/S")):
        for kind, tag in (("raw", "before"), ("cor", "after")):
            x = dd[kind][key]
            h = ROOT.TH1F(f"h_{key}_{kind}_{run}", "", nb, 0.0, 1.8 * E)
            for v in x:
                h.Fill(float(v))
            h.SetDirectory(0); hs.append(h)
            legs.append(f"{lab} {tag}: {x.std() / x.mean() * 100:.1f}%")
    ymax = max(h.GetMaximum() for h in hs) * 1.45
    nm = f"poscorr_dist_E{E:.0f}"
    DrawHistos(hs, legs, 0.0, 1.8 * E, "Energy [GeV]", 0.0, ymax, f"Events / {1.8 * E / nb:.2g} GeV", nm,
               dology=False, drawoptions="HIST", legendoptions="L",
               mycolors=[ROOT.kRed - 7, ROOT.kRed + 1, ROOT.kGreen - 6, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack],
               linestyles=[2, 1, 2, 1, 2, 1], legendPos=[0.17, 0.64, 0.92, 0.88], legendNCols=2, outdir=OUT, lumitext=f"#pi^{{+}} {E:.0f} GeV, run {run}",
               addOverflow=True, addUnderflow=True, usePDF=True)
    names.append(nm)

# 3. resolution (rms/mean) vs energy, before / after, with DREAM
dream = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
Es = [data[r]["E"] for r, _ in RUNS]


def gr(xs, ys):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, x, y); g.SetPointError(i, 0.0, 0.0)
    return g


for metric, ylab in (("rms", "Energy resolution #sigma_{rms}/mean [%]"), ("core", "Energy resolution, Gaussian core #sigma/mean [%]")):
    graphs, legs, cols, mks = [], [], [], []
    COL = {("S", "raw"): ROOT.kRed - 9, ("S", "cor"): ROOT.kRed + 1, ("DR", "raw"): ROOT.kGreen - 8, ("DR", "cor"): ROOT.kGreen + 2,
           ("QS", "raw"): ROOT.kGray + 1, ("QS", "cor"): ROOT.kBlack}
    MK = {"S": 20, "DR": 23, "QS": 22}
    for key, lab in (("S", "S"), ("DR", "DR"), ("QS", "Q/S")):
        for kind, tag in (("raw", ""), ("cor", ", pos. corr.")):
            ys = []
            for run, _ in RUNS:
                x = data[run][kind][key]
                if metric == "rms":
                    ys.append(x.std() / x.mean() * 100)
                else:
                    cm, cs = core(x); ys.append(cs / cm * 100)
            graphs.append(gr([e + (2.0 if kind == "cor" else -2.0) for e in Es], ys)); legs.append(f"{lab}{tag}")
            cols.append(COL[(key, kind)]); mks.append(MK[key])
    if metric == "rms":
        graphs += [gr(dream["energy_GeV"], dream["S_res_pct"]), gr(dream["energy_GeV"], dream["QS_res_pct"])]
        legs += ["DREAM S", "DREAM Q/S"]; cols += [ROOT.kOrange + 1, ROOT.kAzure + 1]; mks += [21, 21]
    nm = f"poscorr_resolution_{metric}"
    note = ROOT.TLatex(0.18, 0.20, "DR: #chi = 0.29 on raw S;  Q/S: DREAM method"); note.SetNDC(); note.SetTextSize(0.026)
    DrawHistos(graphs, legs, 0, 260, "Pion energy [GeV]", 0, 32, ylab, nm, dology=False, drawoptions="EP", legendoptions="P",
               markerstyles=mks, mycolors=cols, legendPos=[0.20, 0.66, 0.93, 0.89], legendNCols=2, extraToDraw=[note], outdir=OUT,
               lumitext="#pi^{+}, pion-selected, hole veto", usePDF=True)
    names.append(nm)
json.dump({"names": names, "slopes": {run: data[run]["slopes"] for run, _ in RUNS}}, open(f"{SCR}/poscorr_names.json", "w"),
          indent=1, default=float)
print("plots:", names)
