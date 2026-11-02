"""Pion response vs energy (signal / beam energy, electron scale): CaloX tb2025 vs DREAM 2005."""
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
NO_TARGET = [1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434]
WITH_TARGET = [1469, 1471, 1465, 1470]
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]


def response(run):
    d = np.load(f"{SCR}/resp_{run}.npz")
    E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION:
        sel &= d[f].astype(bool)
    out = {"E": E, "N": int(sel.sum()), "N_all": len(sel)}
    for key, col in (("S", "FERS_energyMix_Sci_calib_sum"), ("C", "FERS_energyMix_Cer_calib_sum"),
                     ("S_lc", "FERS_energyMix_Sci_calib_sum_LeakCorr"), ("C_lc", "FERS_energyMix_Cer_calib_sum_LeakCorr")):
        x = d[col][sel] / E
        out[key] = (float(x.mean()), float(x.std() / np.sqrt(len(x))))
    # DREAM Q/S method (Section 5.3): S corrected for lateral leakage, S (1 + f_out); em fraction from the
    # leakage-corrected C/S with e/h = 1.3 (S) and 5 (C); S_final = S_corr / (f_em + (1 - f_em) / 1.3)
    S, C = d["FERS_energyMix_Sci_calib_sum"][sel], d["FERS_energyMix_Cer_calib_sum"][sel]
    Sc = S + d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    a, b = 1 / 5.0, 1 / 1.3
    qs = C / Sc
    f = np.clip((a - qs * b) / (qs * (1 - b) - (1 - a)), -0.5, 1.5)
    x = Sc / (f + (1 - f) * b) / E
    out["QS"] = (float(x.mean()), float(x.std() / np.sqrt(len(x))))
    # dual-readout energy as in variables/fers.py getFERSEnergyDR: raw S and C, chi = 0.29
    x = (S - 0.29 * C) / 0.71 / E
    out["DR"] = (float(x.mean()), float(x.std() / np.sqrt(len(x))))
    return out


def graph(points, key, xshift=0.0):
    g = ROOT.TGraphErrors(len(points))
    for i, p in enumerate(points):
        g.SetPoint(i, p["E"] * (1 + xshift), p[key][0])
        g.SetPointError(i, 0.0, p[key][1])
    return g


rows = {}
for run in NO_TARGET + WITH_TARGET:
    if os.path.exists(f"{SCR}/resp_{run}.npz"):
        rows[run] = response(run)
print(f"{'run':>5s} {'E':>5s} {'pions':>12s} {'S/E':>7s} {'C/E':>7s} {'S_lc/E':>7s} {'C_lc/E':>7s} {'C/S':>6s}")
for run, r in rows.items():
    print(f"{run:5d} {r['E']:5.0f} {r['N']:6d}/{r['N_all']:5d} {r['S'][0]:7.3f} {r['C'][0]:7.3f} {r['S_lc'][0]:7.3f} "
          f"{r['C_lc'][0]:7.3f} {r['C'][0] / r['S'][0]:6.3f}  {'(target)' if run in WITH_TARGET else ''}")

dream = json.load(open(f"{STUDY_DIR}/dream_fig11.json"))
gd = ROOT.TGraph(len(dream["energy_GeV"]))
for i, (e, v) in enumerate(zip(dream["energy_GeV"], dream["S_over_E"])):
    gd.SetPoint(i, e, v)
gdc = ROOT.TGraph(1); gdc.SetPoint(0, 100.0, dream["C_over_E_100GeV"])
gdq = ROOT.TGraphErrors(len(dream["energy_GeV"]))
for i, (e, v) in enumerate(zip(dream["energy_GeV"], dream["QS_over_E"])):
    gdq.SetPoint(i, e, v); gdq.SetPointError(i, 0.0, 0.0)

nt = [rows[r] for r in NO_TARGET if r in rows]
wt = [rows[r] for r in WITH_TARGET if r in rows]
graphs = [graph(nt, "S"), graph(nt, "C"), graph(nt, "DR"), graph(nt, "QS"), gd, gdc, gdq]
legends = ["CaloX S", "CaloX C", "CaloX DR (#chi = 0.29, raw S)", "CaloX Q/S-corrected", "DREAM S", "DREAM C (100 GeV)",
           "DREAM Q/S-corrected"]
markers = [20, 21, 23, 22, 20, 21, 22]
colors = [ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kBlack, ROOT.kRed - 9, ROOT.kBlue - 9, ROOT.kGray + 1]
line = ROOT.TLine(0, 1.0, 320, 1.0); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray + 2)
os.makedirs(OUT, exist_ok=True)
ref = ROOT.TLatex(0.18, 0.20, "DREAM: NIM A 537 (2005) 537, Figs. 11, 35 and text"); ref.SetNDC(); ref.SetTextSize(0.026)
DrawHistos(graphs, legends, 0, 320, "Pion energy [GeV]", 0.2, 1.35, "Mean signal / beam energy",
           "pion_response_vs_DREAM", dology=False, drawoptions="EP", legendoptions="EP",
           markerstyles=markers, mycolors=colors, legendPos=[0.40, 0.22, 0.93, 0.58], legendNCols=1,
           extraToDraw=[line, ref], outdir=OUT, noLumi=False, lumitext="#pi^{+}, pion-selected", usePDF=True)

# leakage-corrected version (outer ring added back with weight 1, as in DREAM Eq. 6 for S)
graphs_lc = [graph(nt, "S_lc"), graph(nt, "C_lc"), gd, gdc]
DrawHistos(graphs_lc, ["CaloX S, leakage-corr.", "CaloX C, leakage-corr.", "DREAM S", "DREAM C (100 GeV)"],
           0, 320, "Pion energy [GeV]", 0.2, 1.15, "Mean signal / beam energy",
           "pion_response_vs_DREAM_leakcorr", dology=False, drawoptions="EP", legendoptions="EP",
           markerstyles=[20, 21, 20, 21], mycolors=[colors[0], colors[1], colors[4], colors[5]], legendPos=[0.50, 0.27, 0.93, 0.46], legendNCols=1,
           extraToDraw=[line, ref], outdir=OUT, noLumi=False, lumitext="#pi^{+}, pion-selected", usePDF=True)
json.dump(rows, open(f"{OUT}/pion_response.json", "w"), indent=1)
print("plots in", OUT)
