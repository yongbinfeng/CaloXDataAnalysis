"""Resolution vs 1/sqrt(E) with straight-line fits sigma/E = a/sqrt(E) + b (DREAM's linear sum), CaloX (30-160 GeV) and DREAM
(Fig. 32a); quadrature fits a/sqrt(E) (+) b for the table. Two sets: final correction (S, C, Q/S) and no correction (S, C, DR;
DREAM's DR is their Q/S). Writes res_fit.json, final_resolution_fit.png and final_uncorr_resolution_fit.png."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("res_fit.py")).read(); exec(src[:src.index("def show(")])          # fit(), F, d32
Es = [e for e in sorted(int(k) for k in F) if e >= 30]
SETS = {"corr": dict(res="res", N="N", sig=(("S", "S", "S_res_pct"), ("C", "C", "C_res_pct"), ("QS", "Q/S", "QS_res_pct")),
                     plot="final_resolution_fit", tag=""),
        "unc": dict(res="res_unc", N="N_unc", sig=(("S", "S", "S_res_pct"), ("C", "C", "C_res_pct"), ("DR", "DR", "QS_res_pct")),
                    plot="final_uncorr_resolution_fit", tag=", uncorr.")}
out = {"calox_E": Es, "dream_E": d32["energy_GeV"]}
for name, st in SETS.items():
    for key, lab, dkey in st["sig"]:
        y = [F[str(e)][st["res"]][key][0] for e in Es]; ey = [v / np.sqrt(2 * F[str(e)][st["N"]]) for v, e in zip(y, Es)]
        dy = d32[dkey]; dey = [0.03 * v for v in dy]
        d = {"calox": {"y": y, "lin": fit(Es, y, ey, False), "quad": fit(Es, y, ey, True)},
             "dream": {"y": dy, "lin": fit(d32["energy_GeV"], dy, dey, False), "quad": fit(d32["energy_GeV"], dy, dey, True)}}
        # stochastic term with the constant term fixed to DREAM's
        bD = d["dream"]["lin"][2]; xx = 1 / np.sqrt(np.array(Es, float)); yy, ee = np.array(y), np.array(ey); w = 1 / ee ** 2
        aF = (w * xx * (yy - bD)).sum() / (w * xx * xx).sum(); c2 = ((yy - bD - aF * xx) ** 2 * w).sum() / (len(Es) - 1)
        d["calox"]["fixb"] = [aF, np.sqrt(1 / (w * xx * xx).sum()) * max(1, np.sqrt(c2)), bD, c2]
        if name == "corr": out[key] = d
        else: out.setdefault("unc", {})[key] = d
json.dump(out, open(f"{SCR}/res_fit.json", "w"), indent=1, default=float)

REPO = REPO_DIR; os.chdir(REPO); sys.path.insert(0, REPO)
import ROOT
ROOT.gROOT.SetBatch(True)
from plotting.my_function import DrawHistos
COL = {"S": ROOT.kGreen + 2, "C": ROOT.kBlue + 1, "QS": ROOT.kRed + 1, "DR": ROOT.kOrange + 7}
LIGHT = {"S": ROOT.kGreen - 6, "C": ROOT.kAzure - 4, "QS": ROOT.kRed - 7, "DR": ROOT.kRed - 7}
for (name, st), mode in [(kv, m) for kv in SETS.items() for m in ("free", "fixb")]:
    src_ = out if name == "corr" else out["unc"]
    G, Lg, Cl, Mk, extra = [], [], [], [], []
    for key, lab, _ in st["sig"]:
        for who, E, y, col, mk, ls in (("calox", Es, src_[key]["calox"]["y"], COL[key], 20, 1),
                                       ("dream", d32["energy_GeV"], src_[key]["dream"]["y"], LIGHT[key], 29, 2)):
            g = ROOT.TGraph(len(E))
            for i, (e, v) in enumerate(zip(E, y)): g.SetPoint(i, 1 / np.sqrt(e), v)
            a, ea, b, eb, _c = src_[key][who]["lin"]
            if who == "calox" and mode == "fixb": a, ea, b, _c = src_[key][who]["fixb"]
            G.append(g); Cl.append(col); Mk.append(mk)
            dlab = "DR (Q/S)" if (who == "dream" and key in ("QS", "DR")) else lab
            Lg.append(f"{'CaloX ' + lab + st['tag'] if who == 'calox' else 'DREAM ' + dlab}: {a:.0f}%/#sqrt{{E}} + {b:.1f}%")
            f = ROOT.TF1(f"f_{name}_{mode}_{key}_{who}", "[0]*x+[1]", 0.0, 0.25); f.SetParameters(a, b); f.SetLineColor(col); f.SetLineStyle(ls); f.SetLineWidth(2)
            extra.append((f, "SAME"))
    note = ROOT.TLatex(0.42, 0.17, "#splitline{fits a/#sqrt{E} + b" + (", CaloX b = DREAM's" if mode == "fixb" else "") + "}{CaloX 30-160 GeV solid, DREAM dashed}"); note.SetNDC(); note.SetTextSize(0.026)
    extra.append(note)
    DrawHistos(G, Lg, 0, 0.25, "1 / #sqrt{E/GeV}", 0, 60, "#sigma_{rms} / mean [%]", st["plot"] + ("_fixb" if mode == "fixb" else ""),
               dology=False, drawoptions="P",
               legendoptions="P", markerstyles=Mk, mycolors=Cl, legendPos=[0.17, 0.60, 0.78, 0.90], legendNCols=1, extraToDraw=extra,
               outdir="results/plots/FERS_PionResponse", lumitext="#pi^{+}", usePDF=True)
for name, st in SETS.items():
  src_ = out if name == "corr" else out["unc"]
  for key, lab, _ in st["sig"]:
    for who in ("calox", "dream"):
        l, q = src_[key][who]["lin"], src_[key][who]["quad"]
        print(f"{name:4s} {lab:4s} {who:5s}: linear a {l[0]:5.1f}+-{l[1]:4.1f} b {l[2]:4.1f}+-{l[3]:3.1f} (chi2/ndf {l[4]:4.1f}) | quad a {q[0]:5.1f}+-{q[1]:4.1f} b {q[2]:4.1f}+-{q[3]:3.1f}")
