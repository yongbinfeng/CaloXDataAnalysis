"""Normalise pion events to positron depth (displacement d = 0) and test against pions that shower near d = 0."""
import json
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
SEL = {"electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"],
       "pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]}


def load(run, part):
    d = np.load(f"{SCR}/hodo_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL[part]:
        sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"] >= 1) & (d["TTU_Hodo_nHitX"] <= 3) & (d["TTU_Hodo_nHitY"] >= 1) & (d["TTU_Hodo_nHitY"] <= 3)
    if part == "pion":
        sel &= d["FERS_energyMix_Sci_calib_sum"] > 0.3 * E
    g = lambda k: d[k][sel]
    return dict(E=E, hx=g("hodoX_c"), hy=g("hodoY_c"), ex=g("FERS_energyMix_Sci_calib_weighted_center_X"),
                ey=g("FERS_energyMix_Sci_calib_weighted_center_Y"), S=g("FERS_energyMix_Sci_calib_sum"), C=g("FERS_energyMix_Cer_calib_sum"))


def fit_map(h):
    X = np.column_stack([h["hx"], h["hy"], np.ones(len(h["hx"]))]); m = np.ones(len(X), bool)
    for _ in range(6):
        cx, *_ = np.linalg.lstsq(X[m], h["ex"][m], rcond=None); cy, *_ = np.linalg.lstsq(X[m], h["ey"][m], rcond=None)
        rx, ry = h["ex"] - X @ cx, h["ey"] - X @ cy
        sx, sy = 1.4826 * np.median(np.abs(rx[m])), 1.4826 * np.median(np.abs(ry[m]))
        m = (np.abs(rx) < 3 * sx) & (np.abs(ry) < 3 * sy)
    return cx, cy


cx, cy = fit_map(load(1527, "electron"))
d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json")); dS = dict(zip(d11["energy_GeV"], d11["S_over_E"]))
print("S/E and C/E: as measured; normalised to positron depth (d = 0) event by event; and the pions that already shower near d = 0")
print(f"{'E':>4s} {'N':>6s} | {'S/E meas':>8s} {'S/E @ d=0':>9s} {'N |d|<2mm':>9s} {'their S/E':>9s} | {'C/E meas':>8s} {'C/E @ d=0':>9s} {'their C/E':>9s} | "
      f"{'C/S @ d=0':>9s} | {'DREAM S/E':>9s}")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py                 # positrons: 0 on average, by construction of the map
    ix, iy = px - np.median(px), py - np.median(py)     # impact point relative to the beam-spot centre
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix * ix, iy * iy, ix * iy])

    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX) < 3) & (np.abs(dY) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            r = yy - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    S_e = p["S"] * np.exp(-(A[:, 1:] @ cS[1:])); C_e = p["C"] * np.exp(-(A[:, 1:] @ cC[1:]))
    near = np.hypot(dX, dY) < 0.2
    nS = p["S"][near].mean() / E if near.sum() >= 20 else np.nan; nC = p["C"][near].mean() / E if near.sum() >= 20 else np.nan
    print(f"{E:4.0f} {len(p['S']):6d} | {p['S'].mean() / E:8.3f} {S_e.mean() / E:9.3f} {near.sum():9d} {nS:9.3f} | "
          f"{p['C'].mean() / E:8.3f} {C_e.mean() / E:9.3f} {nC:9.3f} | {C_e.mean() / S_e.mean():9.3f} | {dS.get(int(E), float('nan')):9.3f}")
