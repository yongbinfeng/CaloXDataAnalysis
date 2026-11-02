"""S vs C planes and chi before/after the event-by-event depth correction (hodoscope displacement), with DREAM's chi."""
import os
import sys

import numpy as np

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("edepth_test.py")).read()
exec(src[:src.index("cx, cy = fit_map")])
cx, cy = fit_map(load(1527, "electron"))
sys.path.insert(0, "/data1/yofeng/CaloX/CaloXDataAnalysis"); os.chdir("/data1/yofeng/CaloX/CaloXDataAnalysis")
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
OUT = "results/plots/FERS_PionResponse"


def fem_from(qs, ehS=1.3, ehC=5.0):
    a, b = 1 / ehC, 1 / ehS
    return (a - qs * b) / (qs * (1 - b) - (1 - a))


def versions(run):
    p = load(run, "pion"); E = p["E"]
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py; ix, iy = px - np.median(px), py - np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix * ix, iy * iy, ix * iy])

    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX) < 3) & (np.abs(dY) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            r = yy - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    Ai = A.copy(); Ai[:, 1:3] = 0
    Am = A.copy(); Am[:, 1] = dX - np.median(dX); Am[:, 2] = dY - np.median(dY)
    out = {}
    for lab, M in (("impact", Ai), ("median", Am), ("edepth", A)):
        out[lab] = (p["S"] * np.exp(-(M[:, 1:] @ cS[1:])) / E, p["C"] * np.exp(-(M[:, 1:] @ cC[1:])) / E)
    return E, out


def chis(s, c):
    chi_c = (1 - s.mean()) / (1 - c.mean())
    f = fem_from(c / s); ok = (f > 0) & (f < 1)
    p1, p0 = np.polyfit(f[ok], s[ok], 1); eh = 1 + p1 / p0
    grid = np.arange(-0.4, 0.81, 0.01)
    res = [((s - x * c) / (1 - x)).std() / ((s - x * c) / (1 - x)).mean() for x in grid]
    return chi_c, (1 - 1 / eh) / (1 - 0.2), grid[int(np.argmin(res))]


names = []
RUNS = [1438, 1437, 1429, 1432, 1433, 1434]
allv = {run: versions(run) for run in RUNS}

# 1. S/E vs C/E at 100 GeV, without and with the depth correction, with the chi = 0.29 line through the electron point
E, v = allv[1432]
for lab, title in (("impact", "impact corr. only"), ("edepth", "+ depth (e^{+} depth)")):
    s, c = v[lab]
    h = ROOT.TH2F(f"h_sc_{lab}", "", 90, 0.2, 1.1, 90, 0.5, 1.4)
    for a, b in zip(c, s):
        h.Fill(float(a), float(b))
    h.SetDirectory(0)
    chi_c, chi_f, chi_r = chis(s, c)
    # DR line: S = 1 - chi (1 - C) passes through the electron point (C, S) = (1, 1)
    l29 = ROOT.TLine(0.2, 1 - 0.29 * 0.8, 1.1, 1 - 0.29 * (-0.1)); l29.SetLineColor(ROOT.kRed + 1); l29.SetLineWidth(3)
    lfit = ROOT.TLine(0.2, 1 - chi_r * 0.8, 1.1, 1 - chi_r * (-0.1)); lfit.SetLineColor(ROOT.kBlack); lfit.SetLineWidth(3); lfit.SetLineStyle(2)
    ep = ROOT.TMarker(1.0, 1.0, 29); ep.SetMarkerSize(2.5); ep.SetMarkerColor(ROOT.kRed + 1)
    cm = ROOT.TMarker(float(c.mean()), float(s.mean()), 20); cm.SetMarkerSize(1.8); cm.SetMarkerColor(ROOT.kBlack)
    t = ROOT.TLatex(0.17, 0.86, f"{title}, best #chi = {chi_r:.2f}"); t.SetNDC(); t.SetTextSize(0.028)
    t2 = ROOT.TLatex(0.17, 0.82, "red: #chi = 0.29 through e^{+} (star); dashed: best #chi"); t2.SetNDC(); t2.SetTextSize(0.025)
    nm = f"chi_SvsC_E100_{lab}"
    DrawHistos([h], [], 0.2, 1.1, "C / beam energy", 0.5, 1.4, "S / beam energy", nm, dology=False, doth2=True,
               drawoptions="COLZ", zmin=0, zmax=h.GetMaximum(), extraToDraw=[l29, lfit, ep, cm, t, t2], outdir=OUT, lumitext="#pi^{+} 100 GeV, hodoscope sel.",
               usePDF=True)
    names.append(nm)

# 2. <S/E> vs f_em at 100 GeV (DREAM Fig. 27), before/after
graphs, legs, cols, mks = [], [], [], []
for lab, leg, col, mk in (("impact", "impact-position corrected only", ROOT.kGray + 1, 20), ("edepth", "+ depth (positron depth)", ROOT.kRed + 1, 21)):
    s, c = v[lab]; f = fem_from(c / s)
    edges = np.arange(0.1, 1.01, 0.075)
    pts = [(0.5 * (a + b), s[(f >= a) & (f < b)].mean(), s[(f >= a) & (f < b)].std() / np.sqrt(((f >= a) & (f < b)).sum()))
           for a, b in zip(edges[:-1], edges[1:]) if ((f >= a) & (f < b)).sum() >= 40]
    g = ROOT.TGraphErrors(len(pts))
    for i, (x, y, e) in enumerate(pts):
        g.SetPoint(i, x, y); g.SetPointError(i, 0, e)
    graphs.append(g); legs.append(leg); cols.append(col); mks.append(mk)
exp_line = ROOT.TF1("eh13", "[0]*(x + (1-x)/1.3)", 0, 1); exp_line.SetParameter(0, 1.0); exp_line.SetLineColor(ROOT.kBlue + 1); exp_line.SetLineStyle(2)
t = ROOT.TLatex(0.18, 0.20, "blue dashed: e/h_{S} = 1.3 (DREAM), electron scale"); t.SetNDC(); t.SetTextSize(0.028)
DrawHistos(graphs, legs, 0.0, 1.05, "f_{em} (from C/S)", 0.6, 1.25, "Mean S / beam energy", "chi_S_vs_fem_E100",
           dology=False, drawoptions="EP", legendoptions="EP", markerstyles=mks, mycolors=cols, legendPos=[0.18, 0.74, 0.80, 0.88],
           extraToDraw=[exp_line, t], outdir=OUT, lumitext="#pi^{+} 100 GeV, hodoscope sel.", usePDF=True)
names.append("chi_S_vs_fem_E100")

# 3. chi vs energy, three methods, without / with the depth correction
graphs, legs, cols, mks = [], [], [], []
Es = [allv[r][0] for r in RUNS]
for lab, tag, shade in (("impact", "impact", 0), ("edepth", "+depth", 1)):
    vals = np.array([chis(*allv[r][1][lab]) for r in RUNS])
    for k, (name, mk) in enumerate((("centroid", 20), ("Fig. 27", 21), ("best res.", 22))):
        g = ROOT.TGraphErrors(len(Es))
        for i, e in enumerate(Es):
            g.SetPoint(i, e + (2 * k - 2) * 1.5, vals[i, k]); g.SetPointError(i, 0, 0)
        graphs.append(g); legs.append(f"{name}, {tag}")
        cols.append([ROOT.kGray + 1, ROOT.kAzure - 9, ROOT.kRed - 9][k] if shade == 0 else [ROOT.kBlack, ROOT.kAzure + 2, ROOT.kRed + 1][k]); mks.append(mk)
band = ROOT.TBox(0, 0.25, 175, 0.29); band.SetFillColorAlpha(ROOT.kOrange, 0.35)
tb = ROOT.TLatex(0.18, 0.20, "orange: DREAM, #chi = 0.25-0.29 (e/h_{S} = 1.25-1.3, e/h_{C} = 5)"); tb.SetNDC(); tb.SetTextSize(0.026)
DrawHistos(graphs, legs, 0, 175, "Pion energy [GeV]", -0.1, 0.9, "#chi", "chi_vs_energy", dology=False, drawoptions="EP", legendoptions="P",
           markerstyles=mks, mycolors=cols, legendPos=[0.18, 0.66, 0.93, 0.89], legendNCols=2, extraToDraw=[band, tb], outdir=OUT,
           lumitext="#pi^{+}, hodoscope sel.", usePDF=True)
names.append("chi_vs_energy")
print("plots:", names)
