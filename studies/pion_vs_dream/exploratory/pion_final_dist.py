"""1D distributions for the pion_response_final correction: S, C, DR (chi = 0.29), Q/S per energy, on the electron scale,
dY treated as depth (solid) or not (dashed)."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final_res.py")).read()
exec(src[:src.index("d32 = json.load")])                     # pion(), ARR, L, qsfun (same correction as the response/resolution plots)
REPO = REPO_DIR; os.chdir(REPO); sys.path.insert(0, REPO)
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
SIG = (("S", "S", ROOT.kGreen + 2), ("C", "C", ROOT.kBlue + 1), ("DR", "DR", ROOT.kOrange + 7), ("QS", "Q/S", ROOT.kRed + 1))
names = []
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    hs, legs, cols, lss = [], [], [], []
    for yd, ls in ((1, 1), (0, 2)):
        ARR.clear(); E, _ = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], yd)
        S, C, So = ARR["final"]; vals = {"S": S / E, "C": C / E, "DR": (S - 0.29 * C) / 0.71 / E, "QS": qsfun(S, C, So) / E}
        for key, lab, col in SIG:
            x = vals[key]; nb = 90 if len(x) > 3000 else 45
            h = ROOT.TH1F(f"h_fin_{key}_{yd}_{run}", "", nb, 0.0, 1.8)
            for v in x: h.Fill(float(v))
            h.SetDirectory(0); hs.append(h); cols.append(col); lss.append(ls)
            legs.append(f"{lab}: {x.mean():.3f}, {100 * x.std() / x.mean():.1f}%")
    ymax = max(h.GetMaximum() for h in hs) * 2.0
    nm = f"pion_dist_final_E{E:.0f}"
    note = ROOT.TLatex(0.18, 0.58, "depth from dX+dY (solid) or dX only (dashed); mean, rms/mean"); note.SetNDC(); note.SetTextSize(0.026)
    DrawHistos(hs, legs, 0.0, 1.8, "Signal / (beam energy #times e^{+} response)", 0.0, ymax, "Events", nm, dology=False,
               drawoptions="HIST", legendoptions="L", mycolors=cols, linestyles=lss, legendPos=[0.17, 0.62, 0.93, 0.88], legendNCols=2, extraToDraw=[note],
               outdir="results/plots/FERS_PionResponse", lumitext=f"#pi^{{+}} {E:.0f} GeV, N = {len(S)}",
               addOverflow=True, addUnderflow=True, usePDF=True)
    names.append(nm)
print("plots:", names)
