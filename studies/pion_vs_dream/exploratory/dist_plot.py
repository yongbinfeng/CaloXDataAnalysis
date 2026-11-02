"""Per-energy distributions of S, C and raw-S DR (signal / beam energy), pion selection + hole veto."""
import json, os, sys
import numpy as np
sys.path.insert(0, "/data1/yofeng/CaloX/CaloXDataAnalysis"); os.chdir("/data1/yofeng/CaloX/CaloXDataAnalysis")
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
OUT = "results/plots/FERS_PionResponse"
ev = np.load(f"{SCR}/calox_events_holeveto.npz")
res = json.load(open(f"{SCR}/calox_res_holeveto.json"))
names = []
for run in ["1441", "1439", "1438", "1437", "1429", "1432", "1433", "1434"]:
    E = int(res[run]["E"])
    nb = 60 if res[run]["N"] < 500 else 120
    hs, legs = [], []
    for k, lab in (("S", "S"), ("C", "C"), ("DR", "DR (#chi = 0.29)")):
        x = ev[f"{run}_{k}"]
        h = ROOT.TH1F(f"h_{k}_{run}", "", nb, 0.0, 1.6 * E)
        for v in x: h.Fill(float(v) * E)              # x is signal / beam energy; plot the energy in GeV
        h.SetDirectory(0)
        hs.append(h)
        legs.append(f"{lab}:  {x.mean():.3f},  {x.std() / x.mean() * 100:.1f}%")
    ymax = max(h.GetMaximum() for h in hs) * 1.45
    name = f"pion_dist_E{E}"
    note = ROOT.TLatex(0.18, 0.66, "mean / beam energy,  rms / mean"); note.SetNDC(); note.SetTextSize(0.028)
    DrawHistos(hs, legs, 0.0, 1.6 * E, "Energy [GeV]", 0.0, ymax, f"Events / {1.6 * E / nb:.2g} GeV",
               name, dology=False, drawoptions="HIST", legendoptions="L", mycolors=[ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2],
               legendPos=[0.17, 0.72, 0.70, 0.88], extraToDraw=[note], outdir=OUT, lumitext=f"#pi^{{+}} {E} GeV, run {run}, N = {res[run]['N']}",
               addOverflow=True, addUnderflow=True, usePDF=True)
    names.append(name)
json.dump(names, open(f"{SCR}/dist_names.json", "w"))
print("plots in", OUT, names)
