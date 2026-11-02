"""CaloX tb2025 pions vs DREAM, final version (one correction, one tilt):
  displacement (shower centre - hodoscope impact) -> depth relative to the positrons, for a total tilt of 2.5 deg along the
  measured direction (dY/dX = 0.71) -> attenuation factor removed with lambda per fibre type (90-degree positron runs)
  -> divided by our positron response at the same energy.
Outputs: final_analysis.json, final_dists.npz, plots final_*.png in results/plots/FERS_PionResponse."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final_res.py")).read(); exec(src[:src.index("d32 = json.load")])     # pion() + ARR, L, e_at, qsfun, R_TILT
TILT = 2.5                                                   # deg, total; best Q/S resolution at 60-160 GeV (scan 2-4 deg)
T = np.tan(np.radians(TILT)) / np.sqrt(1 + R_TILT ** 2)      # X component used inside pion()
d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
DRES = {int(e): (s, c, q) for e, s, c, q in zip(d32["energy_GeV"], d32["S_res_pct"], d32["C_res_pct"], d32["QS_res_pct"])}
RUNS = (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434)

def coreres(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return 100 * s / m
def fem_from(q, ehS=1.3, ehC=5.0):
    a, b = 1/ehC, 1/ehS; return (a - q*b)/(q*(1-b) - (1-a))
def chis(s, c):
    """chi three ways: centroid of the (C, S) cloud relative to the electron point; DREAM Fig. 27 slope; best DR resolution."""
    chi_c = (1 - s.mean()) / (1 - c.mean())
    f = fem_from(c / s); ok = (f > 0) & (f < 1); p1, p0 = np.polyfit(f[ok], s[ok], 1); eh = 1 + p1 / p0
    grid = np.arange(-0.4, 0.81, 0.01)
    res = [((s - x*c)/(1 - x)).std() / ((s - x*c)/(1 - x)).mean() for x in grid]
    return chi_c, (1 - 1/eh) / (1 - 0.2), grid[int(np.argmin(res))], eh

lam0 = {t: L[t][0] for t in L}
out, dists = {"tilt_deg": TILT, "r_tilt": R_TILT, "lam": L, "esc": {str(k): v for k, v in esc.items()}, "ref0": [x0, y0], "runs": {}}, {}
for run in RUNS:
    ARR.clear(); E, nom = pion(run, lam0["Sci"], lam0["CerQuartz"], lam0["CerPlastic"], 1)
    (S, C, So), (S0, C0, So0) = ARR["final"], ARR["raw"]
    var = []                                                 # attenuation-length uncertainties
    for t in L:
        for sg in (1, -1):
            lam = dict(lam0); lam[t] = lam0[t] + sg * L[t][1]
            var.append(pion(run, lam["Sci"], lam["CerQuartz"], lam["CerPlastic"], 1)[1])
    err = {k: float(np.sqrt(sum((v[k] - nom[k]) ** 2 for v in var) / 2)) for k in ("S", "C", "Q")}
    sig = lambda a, b, c: {"S": a, "C": b, "DR": (a - 0.29*b)/0.71, "QS": qsfun(a, b, c)}
    fin, raw = sig(S, C, So), sig(S0, C0, So0)
    cc, cf, cb, eh = chis(S / E, C / E)
    out["runs"][str(int(E))] = dict(run=run, N=nom["N"], fq=nom["fq"],
        resp={k: float(v.mean() / E) for k, v in fin.items()}, eresp=err, CS=float(C.mean() / S.mean()),
        resp_raw={k: float(v.mean() / E) for k, v in raw.items()},
        res={k: (float(100 * v.std() / v.mean()), float(coreres(v))) for k, v in fin.items()},
        res_raw={k: (float(100 * v.std() / v.mean()), float(coreres(v))) for k, v in raw.items()},
        chi=dict(centroid=float(cc), fig27=float(cf), best=float(cb), eh_S=float(eh)))
    for k in fin: dists[f"{int(E)}_{k}_fin"] = fin[k] / E; dists[f"{int(E)}_{k}_raw"] = raw[k] / E
    r = out["runs"][str(int(E))]
    print(f"{E:4.0f} N {r['N']:5d} | S/e {r['resp']['S']:.3f} C/e {r['resp']['C']:.3f} C/S {r['CS']:.3f} Q/S/e {r['resp']['QS']:.3f} | "
          f"res rms S {r['res']['S'][0]:.1f} C {r['res']['C'][0]:.1f} DR {r['res']['DR'][0]:.1f} QS {r['res']['QS'][0]:.1f}/{r['res']['QS'][1]:.1f} "
          f"(uncorr. QS {r['res_raw']['QS'][0]:.1f}) | chi {cc:.2f} {cf:.2f} {cb:.2f}")
json.dump(out, open(f"{SCR}/final_analysis.json", "w"), indent=1)
np.savez(f"{SCR}/final_dists.npz", **dists)

# ---------------- plots
REPO = REPO_DIR; os.chdir(REPO); sys.path.insert(0, REPO)
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
OUTD = "results/plots/FERS_PionResponse"
R = out["runs"]; Es = sorted(int(e) for e in R)
COL = {"S": ROOT.kGreen + 2, "C": ROOT.kBlue + 1, "DR": ROOT.kOrange + 7, "QS": ROOT.kRed + 1}
LIGHT = {"S": ROOT.kGreen - 6, "C": ROOT.kAzure - 4, "QS": ROOT.kRed - 7}
LAB = {"S": "S", "C": "C", "DR": "DR", "QS": "Q/S"}
def graph(pts):
    g = ROOT.TGraphErrors(len(pts))
    for i, (x, y, e) in enumerate(pts): g.SetPoint(i, x, y); g.SetPointError(i, 0, e)
    return g
lumi = f"#pi^{{+}}, tilt {TILT}#circ"

# response
G, Lg, Cl, Mk = [], [], [], []
DRESP = {"S": list(zip(d11["energy_GeV"], d11["S_over_E"])), "C": [(100, 0.640)], "QS": list(zip(d11["energy_GeV"], d11["QS_over_E"]))}
for k in ("S", "C", "QS"):
    G.append(graph([(e, R[str(e)]["resp"][k], R[str(e)]["eresp"]["Q" if k == "QS" else k]) for e in Es])); Lg.append(f"CaloX {LAB[k]}"); Cl.append(COL[k]); Mk.append(20)
    G.append(graph([(e, y, 0) for e, y in DRESP[k]])); Lg.append(f"DREAM {LAB[k]}"); Cl.append(LIGHT[k]); Mk.append(29)
DrawHistos(G, Lg, 0, 175, "Pion energy [GeV]", 0.45, 1.3, "Signal / (beam energy #times e^{+} response)", "final_response",
           dology=False, drawoptions="EP", legendoptions="P", markerstyles=Mk, mycolors=Cl, legendPos=[0.20, 0.74, 0.90, 0.89],
           legendNCols=2, outdir=OUTD, lumitext=lumi, usePDF=True)
# resolution: rms (DREAM S, C) and Gaussian core (DREAM Q/S); 20 GeV (750 events) left out
for idx, kind, dk in ((0, "rms", ("S", "C")), (1, "core", ("QS",))):
    G, Lg, Cl, Mk = [], [], [], []
    for k in ("S", "C", "DR", "QS"):
        G.append(graph([(e, R[str(e)]["res"][k][idx], 0) for e in Es if e >= 30])); Lg.append(f"CaloX {LAB[k]}"); Cl.append(COL[k]); Mk.append(20)
    for k in dk:
        G.append(graph([(e, DRES[e][("S", "C", "QS").index(k)], 0) for e in sorted(DRES) if e >= 30])); Lg.append(f"DREAM {LAB[k]}"); Cl.append(LIGHT[k]); Mk.append(29)
    DrawHistos(G, Lg, 0, 175, "Pion energy [GeV]", 0, 50 if kind == "rms" else 36,
               "#sigma_{rms} / mean [%]" if kind == "rms" else "Gaussian-core #sigma / mean [%]", f"final_resolution_{kind}",
               dology=False, drawoptions="P", legendoptions="P", markerstyles=Mk, mycolors=Cl, legendPos=[0.40, 0.68, 0.93, 0.89],
               legendNCols=2, outdir=OUTD, lumitext=lumi, usePDF=True)
# distributions: corrected (solid) and uncorrected (dashed)
for e in Es:
    hs, lg, cl, ls = [], [], [], []
    for kind, style in (("fin", 1), ("raw", 2)):
        for k in ("S", "C", "DR", "QS"):
            x = dists[f"{e}_{k}_{kind}"]; nb = 90 if len(x) > 3000 else 45
            h = ROOT.TH1F(f"h_{e}_{k}_{kind}", "", nb, 0, 1.8)
            for v in x: h.Fill(float(v))
            h.SetDirectory(0); hs.append(h); cl.append(COL[k]); ls.append(style)
            lg.append(f"{LAB[k]}: {x.mean():.3f}, {100 * x.std() / x.mean():.1f}%")
    note = ROOT.TLatex(0.18, 0.58, "solid: corrected, dashed: uncorrected (mean, rms/mean)"); note.SetNDC(); note.SetTextSize(0.026)
    DrawHistos(hs, lg, 0, 1.8, "Signal / (beam energy #times e^{+} response)", 0, max(h.GetMaximum() for h in hs) * 2.0, "Events",
               f"final_dist_E{e}", dology=False, drawoptions="HIST", legendoptions="L", mycolors=cl, linestyles=ls,
               legendPos=[0.17, 0.62, 0.93, 0.88], legendNCols=2, extraToDraw=[note], outdir=OUTD,
               lumitext=f"#pi^{{+}} {e} GeV, N = {len(dists[f'{e}_S_fin'])}", addOverflow=True, addUnderflow=True, usePDF=True)
print("plots done")
