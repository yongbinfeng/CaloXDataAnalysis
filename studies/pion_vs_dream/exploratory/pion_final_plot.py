"""Plot: pion S, C, Q/S on our electron scale at positron depth (measured attenuation), both dY options, vs DREAM."""
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
R = json.load(open(f"{SCR}/pion_final.json"))["rows"]; d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json"))
DREAM = {"S": list(zip(d11["energy_GeV"], d11["S_over_E"])), "C": [(100, 0.640)], "Q": list(zip(d11["energy_GeV"], d11["QS_over_E"]))}
graphs, legs, cols, mks = [], [], [], []
def add(pts, leg, col, mk):
    g = ROOT.TGraphErrors(len(pts))
    for i, (x, y, e) in enumerate(pts): g.SetPoint(i, x, y); g.SetPointError(i, 0, e)
    graphs.append(g); legs.append(leg); cols.append(col); mks.append(mk)
# one legend row per signal: CaloX (dY not depth), CaloX (dY = depth), DREAM
for key, lab, col, light in (("S", "S", ROOT.kGreen + 2, ROOT.kGreen - 6), ("C", "C", ROOT.kBlue + 1, ROOT.kAzure - 4),
                             ("Q", "Q/S", ROOT.kRed + 1, ROOT.kRed - 7)):
    for yd, tag, mk, dx in (("0", "", 20, -1.5), ("1", ", dX+dY", 34, 1.5)):
        add(sorted((float(E) + dx, v[key], v["e" + key]) for E, v in R[yd].items()), f"{lab}{tag}", col, mk)
    add([(e, y, 0) for e, y in DREAM[key]], f"DREAM {lab}", light, 29)
DrawHistos(graphs, legs, 0, 175, "Pion energy [GeV]", 0.45, 1.45, "Signal / (beam energy #times e^{+} response)", "pion_response_final",
           dology=False, drawoptions="EP", legendoptions="P", markerstyles=mks, mycolors=cols, legendPos=[0.17, 0.74, 0.93, 0.89],
           legendNCols=3, outdir="results/plots/FERS_PionResponse", lumitext="#pi^{+} at e^{+} depth", usePDF=True)
print("ok")
