"""Stochastic and constant terms: sigma/E = a/sqrt(E) + b (linear sum, DREAM's choice) and a/sqrt(E) (+) b (quadrature),
weighted least squares, statistical errors on rms/mean ~ (rms/mean)/sqrt(2N), errors scaled by sqrt(chi2/ndf) when > 1."""
import json
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
F = json.load(open(f"{SCR}/final_analysis.json"))["runs"]; d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
def fit(E, y, ey, quad):
    E, y, ey = map(np.asarray, (E, y, ey))
    if quad:      # y^2 = A/E + B, linear in (A, B)
        X = np.column_stack([1 / E, np.ones_like(E)]); Y = y ** 2; eY = 2 * y * ey
    else:
        X = np.column_stack([1 / np.sqrt(E), np.ones_like(E)]); Y = y; eY = ey
    W = 1 / eY ** 2; cov = np.linalg.inv(X.T @ (X * W[:, None])); p = cov @ (X.T @ (W * Y))
    chi2 = (((Y - X @ p) / eY) ** 2).sum(); ndf = len(E) - 2; scale = max(1.0, np.sqrt(chi2 / ndf)); ep = np.sqrt(np.diag(cov)) * scale
    if quad:
        a, b = np.sqrt(max(p[0], 0)), np.sqrt(max(p[1], 0))
        ea = ep[0] / (2 * a) if a > 0 else np.nan; eb = ep[1] / (2 * b) if b > 0 else np.nan
        return a, ea, b, eb, chi2 / ndf
    return p[0], ep[0], p[1], ep[1], chi2 / ndf
def show(label, E, y, N=None, eyfix=None):
    ey = [v / np.sqrt(2 * n) for v, n in zip(y, N)] if N is not None else [eyfix * v for v in y]
    out = []
    for quad, tag in ((False, "a/sqrtE + b"), (True, "a/sqrtE (+) b")):
        a, ea, b, eb, c = fit(E, y, ey, quad)
        out.append(f"{tag}: a = {a:5.1f} +- {ea:4.1f} %, b = {b:4.1f} +- {eb:3.1f} % (chi2/ndf {c:5.1f})")
    print(f"  {label:26s} " + " | ".join(out))
Es = sorted(int(e) for e in F)
for rng, lo in (("20-160 GeV", 20), ("30-160 GeV", 30)):
    E = [e for e in Es if e >= lo]
    print(f"CaloX, {rng} (rms/mean):")
    for key, lab in (("S", "S"), ("C", "C"), ("DR", "DR (chi=0.29)"), ("QS", "Q/S")):
        show(f"{lab}, corrected", E, [F[str(e)]["res"][key][0] for e in E], [F[str(e)]["N"] for e in E])
    show("Q/S core, corrected", E, [F[str(e)]["res"]["QS"][1] for e in E], [F[str(e)]["N"] for e in E])
    for key, lab in (("S", "S"), ("C", "C"), ("DR", "DR")):
        show(f"{lab}, uncorrected", E, [F[str(e)]["res_unc"][key][0] for e in E], [F[str(e)]["N_unc"] for e in E])
print("DREAM (Fig. 32a as read off; 3% relative errors assumed):")
E = d32["energy_GeV"]
for key, lab in (("S_res_pct", "S"), ("C_res_pct", "C"), ("QS_res_pct", "Q/S")):
    show(lab, E, d32[key], eyfix=0.03)
print("energies used for DREAM:", E)
