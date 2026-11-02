"""Plots: resolutions for the pion_response_final correction (rms/mean and Gaussian core), both dY options, with DREAM
(S, C on the rms plot; Q/S, most likely Gaussian widths, on the core plot)."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json, os, sys
REPO = REPO_DIR; SCR = DATA_DIR
os.chdir(REPO); sys.path.insert(0, REPO)
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
R = json.load(open(f"{SCR}/pion_final_res.json")); d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
DREAM = {"S": list(zip(d32["energy_GeV"], d32["S_res_pct"])), "C": list(zip(d32["energy_GeV"], d32["C_res_pct"])),
         "QS": list(zip(d32["energy_GeV"], d32["QS_res_pct"]))}
SIG = (("S", "S", ROOT.kGreen + 2, ROOT.kGreen - 6), ("C", "C", ROOT.kBlue + 1, ROOT.kAzure - 4),
       ("DR", "DR", ROOT.kOrange + 7, None), ("QS", "Q/S", ROOT.kRed + 1, ROOT.kRed - 7))
Es = sorted(int(e) for e in R if int(e) >= 30)       # 20 GeV: 750 events, left out as in section 2
for idx, kind, dream_keys in ((0, "rms", ("S", "C")), (1, "core", ("QS",))):
    graphs, legs, cols, mks = [], [], [], []
    dgraphs = []
    for key, lab, col, light in SIG:
        for v, tag, mk, dx in (("final_dY0", "", 20, -1.5), ("final_dY1", ", dX+dY", 34, 1.5)):
            g = ROOT.TGraph(len(Es))
            for i, e in enumerate(Es): g.SetPoint(i, e + dx, R[str(e)][v][key][idx])
            graphs.append(g); legs.append(f"{lab}{tag}"); cols.append(col); mks.append(mk)
        if key in dream_keys:
            pts = [(e, y) for e, y in DREAM[key] if e >= 30]; g = ROOT.TGraph(len(pts))
            for i, (e, y) in enumerate(pts): g.SetPoint(i, e, y)
            dgraphs.append((g, f"DREAM {lab}", light))
    for g, lab, col in dgraphs:
        graphs.append(g); legs.append(lab); cols.append(col); mks.append(29)
    ylab = "#sigma_{rms} / mean [%]" if kind == "rms" else "Gaussian-core #sigma / mean [%]"
    DrawHistos(graphs, legs, 0, 175, "Pion energy [GeV]", 0, 62 if kind == "rms" else 45, ylab, f"pion_resolution_final_{kind}",
               dology=False, drawoptions="P", legendoptions="P", markerstyles=mks, mycolors=cols, legendPos=[0.20, 0.66, 0.93, 0.89],
               legendNCols=2, outdir="results/plots/FERS_PionResponse", lumitext="#pi^{+} at e^{+} depth", usePDF=True)
print("ok")
