"""Uncorrected S/E vs the shower displacement (shower centre - hodoscope impact) in X and Y, relative to the positron depth:
deeper showers (larger displacement) give more scintillation light, the attenuation effect. Same selection as final_analysis.py."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final.py")).read(); exec(src[:src.index("def pion(")])          # cx, cy, x0, y0, SEL
RUNS = (1439, 1438, 1437, 1429, 1432, 1433, 1434)
data = {}
for run in RUNS:
    d = np.load(f"{SCR}/hodoq_{run}.npz"); E = float(d["energy"]); sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"]>=1)&(d["TTU_Hodo_nHitX"]<=3)&(d["TTU_Hodo_nHitY"]>=1)&(d["TTU_Hodo_nHitY"]<=3)&(d["FERS_energyMix_Sci_calib_sum"]>0.3*E)
    g = lambda k: d[k][sel].astype(float)
    px = cx[0]*g("hodoX_c")+cx[1]*g("hodoY_c")+cx[2]; py = cy[0]*g("hodoX_c")+cy[1]*g("hodoY_c")+cy[2]
    data[int(E)] = dict(dX=g("FERS_energyMix_Sci_calib_weighted_center_X")-px-x0, dY=g("FERS_energyMix_Sci_calib_weighted_center_Y")-py-y0,
                        s=g("FERS_energyMix_Sci_calib_sum") / E)
REPO = REPO_DIR; os.chdir(REPO); sys.path.insert(0, REPO)
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
COLS = [ROOT.kOrange + 7, ROOT.kRed + 1, ROOT.kMagenta + 1, ROOT.kViolet + 1, ROOT.kBlue + 1, ROOT.kAzure + 1, ROOT.kGreen + 2]
edges = np.arange(-1.5, 3.51, 0.25)
for axis in ("dX", "dY"):
    G, Lg = [], []
    for (E, v), col in zip(sorted(data.items()), COLS):
        pts = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (v[axis] >= lo) & (v[axis] < hi)
            if m.sum() >= 40: pts.append((0.5 * (lo + hi), v["s"][m].mean(), v["s"][m].std() / np.sqrt(m.sum())))
        g = ROOT.TGraphErrors(len(pts))
        for i, (x, y, e) in enumerate(pts): g.SetPoint(i, x, y); g.SetPointError(i, 0, e)
        G.append(g); Lg.append(f"{E} GeV")
    zero = ROOT.TLine(0, 0.55, 0, 1.25); zero.SetLineStyle(2); zero.SetLineColor(ROOT.kGray + 2)
    t = ROOT.TLatex(0.18, 0.86, "0 = positron depth; larger = deeper shower"); t.SetNDC(); t.SetTextSize(0.028)
    lab = "X" if axis == "dX" else "Y"
    DrawHistos(G, Lg, -1.5, 3.5, f"Shower centre {lab} #minus hodoscope impact {lab} [cm]", 0.55, 1.25, "Mean S / beam energy (uncorrected)",
               f"final_uncorr_S_vs_{axis}", dology=False, drawoptions="EP", legendoptions="EP", markerstyles=[20] * len(G),
               mycolors=COLS[:len(G)], legendPos=[0.62, 0.16, 0.92, 0.48], extraToDraw=[zero, t], outdir="results/plots/FERS_PionResponse",
               lumitext="#pi^{+}", usePDF=True)
    for E, v in sorted(data.items()):
        k = np.abs(v[axis]) < 3; p = np.polyfit(v[axis][k], v["s"][k], 1)
        print(f"{axis} {E:4d} GeV: N {len(v['s'])}, median {axis} {np.median(v[axis]):+.2f} cm, linear slope of S/E {100 * p[0] / np.polyval(p, 0):+.1f}% per cm")
