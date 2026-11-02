"""Pion response and resolution with the fibre-hodoscope position correction, vs DREAM."""
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
SEL = {"electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"],
       "pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]}
PION_RUNS = [(1441, ROOT.kGray + 2), (1439, ROOT.kOrange + 7), (1438, ROOT.kMagenta + 1), (1437, ROOT.kCyan + 2),
             (1429, ROOT.kRed + 1), (1432, ROOT.kBlue + 1), (1433, ROOT.kGreen + 2), (1434, ROOT.kViolet + 1)]
MAP_RUN = 1527
MIN_EVENTS = 300


def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6:
            break
        m, s = m2, s2
    return m, s


def load(run, part):
    d = np.load(f"{SCR}/hodo_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL[part]:
        sel &= d[f].astype(bool)
    nx, ny = d["TTU_Hodo_nHitX"], d["TTU_Hodo_nHitY"]
    sel &= (nx >= 1) & (nx <= 3) & (ny >= 1) & (ny <= 3)
    g = lambda k: d[k][sel]
    out = dict(E=E, hx=g("hodoX_c"), hy=g("hodoY_c"), ex=g("FERS_energyMix_Sci_calib_weighted_center_X"),
               ey=g("FERS_energyMix_Sci_calib_weighted_center_Y"), S=g("FERS_energyMix_Sci_calib_sum"),
               C=g("FERS_energyMix_Cer_calib_sum"), SO=g("FERS_energyMix_Sci_calib_sum_OuterRing"))
    if part == "pion":
        k = out["S"] > 0.3 * E
        out = {kk: (v[k] if isinstance(v, np.ndarray) else v) for kk, v in out.items()}
    return out


def fit_map(h):
    X = np.column_stack([h["hx"], h["hy"], np.ones(len(h["hx"]))]); m = np.ones(len(X), bool)
    for _ in range(6):
        cx, *_ = np.linalg.lstsq(X[m], h["ex"][m], rcond=None); cy, *_ = np.linalg.lstsq(X[m], h["ey"][m], rcond=None)
        rx, ry = h["ex"] - X @ cx, h["ey"] - X @ cy
        sx, sy = 1.4826 * np.median(np.abs(rx[m])), 1.4826 * np.median(np.abs(ry[m]))
        m = (np.abs(rx) < 3 * sx) & (np.abs(ry) < 3 * sy)
    return cx, cy, sx, sy


def energies(S, C, SO):
    Sc = S + SO; b, a = 1 / 1.3, 0.2
    f = np.clip((a - C / Sc * b) / (C / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    return {"S": S, "C": C, "DR": (S - 0.29 * C) / 0.71, "QS": Sc / (f + (1 - f) * b)}


def fit_corr(sig, vx, vy, ix, iy):
    """ln(signal) vs displacement (vx, vy; linear) and impact position (ix, iy; linear + quadratic), robust LSQ.
    Returns (coefficients, multiplicative correction factor per event)."""
    y = np.log(np.maximum(sig, 1e-3))
    A = np.column_stack([np.ones(len(y)), vx, vy, ix, iy, ix * ix, iy * iy, ix * iy])
    m = (np.abs(vx) < 3) & (np.abs(vy) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
    for _ in range(5):
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
        r = y - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
    return c, np.exp(-(A[:, 1:] @ c[1:]))


# ---- mapping from the positron run, cross-checked on another positron run
pos = load(MAP_RUN, "electron")
cx, cy, sx, sy = fit_map(pos)
print(f"map from run {MAP_RUN}: x = {cx[0]:.3f} hx + {cx[1]:+.3f} hy + {cx[2]:+.3f};  y = {cy[0]:+.3f} hx + {cy[1]:.3f} hy + {cy[2]:+.3f};"
      f"  residual {sx * 10:.2f} / {sy * 10:.2f} mm")
if os.path.exists(f"{SCR}/hodo_1526.npz"):
    chk = load(1526, "electron"); c2x, c2y, s2x, s2y = fit_map(chk)
    rx = chk["ex"] - (cx[0] * chk["hx"] + cx[1] * chk["hy"] + cx[2]); ry = chk["ey"] - (cy[0] * chk["hx"] + cy[1] * chk["hy"] + cy[2])
    print(f"check on run 1526 (60 GeV e+): own fit scale {c2x[0]:.3f}/{c2y[1]:.3f}, offset {c2x[2]:+.3f}/{c2y[2]:+.3f};"
          f" with the run-{MAP_RUN} map: median residual {np.median(rx) * 10:+.2f} / {np.median(ry) * 10:+.2f} mm")
names = []

# mapping check plot: positron shower centre vs mapped hodoscope position (x)
h2 = ROOT.TH2F("h_map_x", "", 80, -1.5, 2.5, 80, -1.5, 2.5)
for a, b in zip(cx[0] * pos["hx"] + cx[1] * pos["hy"] + cx[2], pos["ex"]):
    h2.Fill(float(a), float(b))
h2.SetDirectory(0)
DrawHistos([h2], [], -1.5, 2.5, "Hodoscope X (detector frame) [cm]", -1.5, 2.5, "Shower centre X (S) [cm]",
           "hodo_map_check_X", dology=False, doth2=True, dologz=True, drawoptions="COLZ", outdir=OUT,
           lumitext=f"e^{{+}} {pos['E']:.0f} GeV, run {MAP_RUN}", usePDF=True)
names.append("hodo_map_check_X")

# ---- pions
data, runs = {}, []
for run, col in PION_RUNS:
    if not os.path.exists(f"{SCR}/hodo_{run}.npz"):
        print(f"run {run}: no dump"); continue
    p = load(run, "pion")
    if len(p["S"]) < MIN_EVENTS:
        print(f"run {run}: {len(p['S'])} events, skipped"); continue
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py
    vx, vy = dX - np.median(dX), dY - np.median(dY)
    ix, iy = px - np.median(px), py - np.median(py)                 # impact point, from the hodoscope
    (cS, wS), (cC, wC) = fit_corr(p["S"], vx, vy, ix, iy), fit_corr(p["C"], vx, vy, ix, iy)
    p.update(dX=dX, dY=dY, raw=energies(p["S"], p["C"], p["SO"]), cor=energies(p["S"] * wS, p["C"] * wC, p["SO"] * wS),
             slopes=(float(cS[1]), float(cS[2]), float(cC[1]), float(cC[2])))
    data[run] = p; runs.append((run, col))
    print(f"run {run} {p['E']:.0f} GeV: {len(p['S'])} pions, median displacement ({np.median(dX) * 10:+.1f}, {np.median(dY) * 10:+.1f}) mm,"
          f" S slopes {cS[1] * 100:+.1f} / {cS[2] * 100:+.1f} %/cm, C {cC[1] * 100:+.1f} / {cC[2] * 100:+.1f} %/cm")

# S/E vs displacement, per energy
edges = np.arange(-2.0, 3.01, 0.5)
for axis, lab in (("dX", "X"), ("dY", "Y")):
    graphs, legs, cols = [], [], []
    for run, col in runs:
        p = data[run]; v = p["S"] / p["E"]; off = p[axis]
        pts = [(0.5 * (a + b), np.median(v[(off >= a) & (off < b)]), 1.25 * np.std(v[(off >= a) & (off < b)]) / np.sqrt(((off >= a) & (off < b)).sum()))
               for a, b in zip(edges[:-1], edges[1:]) if ((off >= a) & (off < b)).sum() >= 25]
        g = ROOT.TGraphErrors(len(pts))
        for i, (x, y, e) in enumerate(pts):
            g.SetPoint(i, x, y); g.SetPointError(i, 0.0, e)
        graphs.append(g); legs.append(f"{p['E']:.0f} GeV"); cols.append(col)
    nm = f"hodo_profile_S_d{lab}"
    DrawHistos(graphs, legs, -2.2, 3.2, f"Shower centre {lab} - hodoscope impact {lab} [cm]", 0.6, 1.35, "Median S / beam energy", nm,
               dology=False, drawoptions="EP", legendoptions="EP", markerstyles=[20] * len(graphs), mycolors=cols,
               legendPos=[0.18, 0.72, 0.92, 0.88], legendNCols=4, outdir=OUT, lumitext="#pi^{+}, pion-selected", usePDF=True)
    names.append(nm)

Es = [data[r]["E"] for r, _ in runs]
dream11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json")); dream32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))


def gr(xs, ys, es=None):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, x, y); g.SetPointError(i, 0.0, es[i] if es else 0.0)
    return g


COL = {("S", "raw"): ROOT.kRed - 9, ("S", "cor"): ROOT.kRed + 1, ("C", "raw"): ROOT.kBlue - 9, ("C", "cor"): ROOT.kBlue + 1,
       ("DR", "raw"): ROOT.kGreen - 8, ("DR", "cor"): ROOT.kGreen + 2, ("QS", "raw"): ROOT.kGray + 1, ("QS", "cor"): ROOT.kBlack}
MK = {"S": 20, "C": 21, "DR": 23, "QS": 22}
LAB = {"S": "S", "C": "C", "DR": "DR", "QS": "Q/S"}

# response vs energy
graphs, legs, cols, mks = [], [], [], []
for key in ("S", "C", "DR", "QS"):
    for kind, tag in (("raw", ""), ("cor", ", hodo corr.")):
        ys = [data[r][kind][key].mean() / data[r]["E"] for r, _ in runs]
        es = [data[r][kind][key].std() / np.sqrt(len(data[r][kind][key])) / data[r]["E"] for r, _ in runs]
        graphs.append(gr([e + (1.5 if kind == "cor" else -1.5) for e in Es], ys, es)); legs.append(f"{LAB[key]}{tag}")
        cols.append(COL[(key, kind)]); mks.append(MK[key])
graphs += [gr(dream11["energy_GeV"], dream11["S_over_E"]), gr(dream11["energy_GeV"], dream11["QS_over_E"]), gr([100.0], [dream11["C_over_E_100GeV"]])]
legs += ["DREAM S", "DREAM Q/S", "DREAM C"]; cols += [ROOT.kOrange + 1, ROOT.kAzure + 1, ROOT.kAzure - 9]; mks += [21, 21, 21]
line = ROOT.TLine(0, 1.0, 320, 1.0); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray + 2)
DrawHistos(graphs, legs, 0, 320, "Pion energy [GeV]", 0.1, 1.45, "Mean signal / beam energy", "hodo_response", dology=False,
           drawoptions="EP", legendoptions="P", markerstyles=mks, mycolors=cols, legendPos=[0.50, 0.16, 0.94, 0.50], legendNCols=2,
           extraToDraw=[line], outdir=OUT, lumitext="#pi^{+}, hodoscope sel.", usePDF=True)
names.append("hodo_response")

# resolution vs energy
for metric, ylab in (("rms", "Energy resolution #sigma_{rms}/mean [%]"), ("core", "Energy resolution, Gaussian core #sigma/mean [%]")):
    graphs, legs, cols, mks = [], [], [], []
    for key in ("S", "C", "DR", "QS"):
        for kind, tag in (("raw", ""), ("cor", ", hodo corr.")):
            ys = []
            for r, _ in runs:
                x = data[r][kind][key]
                if metric == "rms":
                    ys.append(x.std() / x.mean() * 100)
                else:
                    m_, s_ = core(x); ys.append(s_ / m_ * 100)
            graphs.append(gr([e + (1.5 if kind == "cor" else -1.5) for e in Es], ys)); legs.append(f"{LAB[key]}{tag}")
            cols.append(COL[(key, kind)]); mks.append(MK[key])
    if metric == "rms":       # DREAM raw S and C are sigma_rms/mean (Section 5.1)
        graphs += [gr(dream32["energy_GeV"], dream32["S_res_pct"]), gr(dream32["energy_GeV"], dream32["C_res_pct"])]
        legs += ["DREAM S (rms)", "DREAM C (rms)"]; cols += [ROOT.kOrange + 1, ROOT.kAzure - 9]; mks += [21, 21]
    else:                     # DREAM corrected distributions are quoted from Gaussian fits (Figs. 16, 29, 33)
        graphs += [gr(dream32["energy_GeV"], dream32["QS_res_pct"])]
        legs += ["DREAM Q/S (Gauss.)"]; cols += [ROOT.kAzure + 1]; mks += [21]
    note = ROOT.TLatex(0.18, 0.20, "DR: #chi = 0.29 on raw S;  Q/S: DREAM method (leakage-corrected S)"); note.SetNDC(); note.SetTextSize(0.024)
    nm = f"hodo_resolution_{metric}"
    DrawHistos(graphs, legs, 0, 260, "Pion energy [GeV]", 0, 62, ylab, nm, dology=False, drawoptions="EP", legendoptions="P",
               markerstyles=mks, mycolors=cols, legendPos=[0.20, 0.60, 0.93, 0.89], legendNCols=2, extraToDraw=[note], outdir=OUT,
               lumitext="#pi^{+}, hodoscope sel.", usePDF=True)
    names.append(nm)

# distributions before / after, S, DR, Q/S
for run, _ in runs:
    p = data[run]; E = p["E"]; nb = 80
    hs, legs = [], []
    for key in ("S", "DR", "QS"):
        for kind, tag in (("raw", "before"), ("cor", "after")):
            x = p[kind][key]
            h = ROOT.TH1F(f"h_hodo_{key}_{kind}_{run}", "", nb, 0.0, 1.8 * E)
            for v in x:
                h.Fill(float(v))
            h.SetDirectory(0); hs.append(h); legs.append(f"{LAB[key]} {tag}: {x.std() / x.mean() * 100:.1f}%")
    ymax = max(h.GetMaximum() for h in hs) * 1.9
    nm = f"hodo_dist_E{E:.0f}"
    DrawHistos(hs, legs, 0.0, 1.8 * E, "Energy [GeV]", 0.0, ymax, f"Events / {1.8 * E / nb:.2g} GeV", nm, dology=False,
               drawoptions="HIST", legendoptions="L",
               mycolors=[ROOT.kRed - 7, ROOT.kRed + 1, ROOT.kGreen - 6, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack],
               linestyles=[2, 1, 2, 1, 2, 1], legendPos=[0.17, 0.64, 0.92, 0.88], legendNCols=2, outdir=OUT,
               lumitext=f"#pi^{{+}} {E:.0f} GeV, N = {len(p['S'])}", addOverflow=True, addUnderflow=True, usePDF=True)
    names.append(nm)

summary = {run: {"E": data[run]["E"], "N": int(len(data[run]["S"])), "slopes": data[run]["slopes"],
                 "median_displacement_mm": [float(np.median(data[run]["dX"]) * 10), float(np.median(data[run]["dY"]) * 10)],
                 **{f"{kind}_{key}_{q}": float(v) for kind in ("raw", "cor") for key in ("S", "C", "DR", "QS")
                    for q, v in (("mean", data[run][kind][key].mean() / data[run]["E"]),
                                 ("rms", data[run][kind][key].std() / data[run][kind][key].mean() * 100),
                                 ("core", (lambda m_, s_: s_ / m_ * 100)(*core(data[run][kind][key]))))}}
           for run, _ in runs}
json.dump({"names": names, "map": {"x": list(map(float, cx)), "y": list(map(float, cy)), "resid_mm": [float(sx * 10), float(sy * 10)]},
           "runs": summary}, open(f"{SCR}/hodo_summary.json", "w"), indent=1)
print("plots:", names)
