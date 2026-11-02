"""Single-pion energy resolution vs energy: CaloX tb2025 vs DREAM 2005 (Fig. 32a)."""
import json, os, sys
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
ours = json.load(open(f"{SCR}/calox_res_holeveto.json"))
dream = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
def g(xs, ys):
    gr = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)): gr.SetPoint(i, x, y); gr.SetPointError(i, 0.0, 0.0)
    return gr
runs = [r for r in ours if ours[r]["E"] >= 30]          # 20 GeV: 104 events after the cuts
E = [ours[r]["E"] for r in runs]
graphs = [g(E, [ours[r]["S_rms"] for r in runs]), g(E, [ours[r]["C_rms"] for r in runs]), g(E, [ours[r]["DR_rms"] for r in runs]),
          g(E, [ours[r]["QS_rms"] for r in runs]),
          g(dream["energy_GeV"], dream["S_res_pct"]), g(dream["energy_GeV"], dream["C_res_pct"]), g(dream["energy_GeV"], dream["QS_res_pct"])]
legends = ["CaloX S", "CaloX C", "CaloX DR (#chi = 0.29, raw S)", "CaloX Q/S-corrected", "DREAM S", "DREAM C",
           "DREAM Q/S-corrected"]
markers = [20, 21, 23, 22, 20, 21, 22]
colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kBlack, ROOT.kRed - 9, ROOT.kBlue - 9, ROOT.kGray + 1]  # DREAM: light, filled (open markers render wrongly in this ROOT)
note = ROOT.TLatex(0.18, 0.20, "#sigma_{rms}/mean. DREAM: NIM A 537 (2005) 537, Fig. 32a"); note.SetNDC(); note.SetTextSize(0.026)
os.makedirs(OUT, exist_ok=True)
DrawHistos(graphs, legends, 0, 320, "Pion energy [GeV]", 0, 45, "Energy resolution [%]",
           "pion_resolution_vs_DREAM", dology=False, drawoptions="EP", legendoptions="P",
           markerstyles=markers, mycolors=colors, legendPos=[0.42, 0.55, 0.93, 0.88], legendNCols=1,
           extraToDraw=[note], outdir=OUT, lumitext="#pi^{+}, pion-selected, hole veto", usePDF=True)
print("plot in", OUT)
