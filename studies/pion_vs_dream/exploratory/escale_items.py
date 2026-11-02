"""Item 1: pion response relative to our own measured positron response, vs DREAM.
Item 2: positron response away from the module centre, and the electron scale as seen by pion showers."""
import json
import os

import numpy as np

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("edepth_res.py")).read()
exec(src[:src.index('print("resolution sigma_rms')])           # load, fit_map, core, qs, cx, cy
rl = json.load(open("/data1/yofeng/CaloX/CaloXDataAnalysis/data/Runlist.json"))
num = lambda s: float(str(s).replace("mm", "").strip())

# ---------------- item 1 ----------------
POS_CENTRE = {20: 1423, 30: 1416, 40: 1514, 60: 1526, 80: 1409, 100: 1527, 120: 1411}
eresp = {}
for E, run in POS_CENTRE.items():
    if not os.path.exists(f"{SCR}/hodo_{run}.npz"):
        continue
    e = load(run, "electron")
    s, c = e["S"] / e["E"], e["C"] / e["E"]
    k = (s > 0.7) & (s < 1.4)
    eresp[E] = (core(s[k])[0], core(c[k])[0])
print("item 1: positron response (Gaussian-core peak) at the central table position")
for E in sorted(eresp):
    print(f"   e+ {E:4d} GeV (run {POS_CENTRE[E]}): S/E {eresp[E][0]:.3f}  C/E {eresp[E][1]:.3f}")
lE = np.log(sorted(eresp)); eS = [eresp[E][0] for E in sorted(eresp)]; eC = [eresp[E][1] for E in sorted(eresp)]
e_at = lambda E: (np.interp(np.log(E), lE, eS), np.interp(np.log(E), lE, eC))    # flat beyond the measured range

d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json"))
dS = dict(zip(d11["energy_GeV"], d11["S_over_E"])); dQ = dict(zip(d11["energy_GeV"], d11["QS_over_E"]))
rows = {}
print("\npion response at positron depth, divided by our positron response at the same energy (pi/e), vs DREAM")
print(f"{'E':>4s} | {'S/E':>6s} {'S pi/e':>7s} {'DREAM':>6s} | {'C/E':>6s} {'C pi/e':>7s} {'DREAM':>6s} | {'Q/S/E':>6s} {'Q/S pi/e':>9s} {'DREAM':>6s}")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    d = np.load(f"{SCR}/hodo_{run}.npz")
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"] >= 1) & (d["TTU_Hodo_nHitX"] <= 3) & (d["TTU_Hodo_nHitY"] >= 1) & (d["TTU_Hodo_nHitY"] <= 3) & (d["FERS_energyMix_Sci_calib_sum"] > 0.3 * E)
    SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py; ix, iy = px - np.median(px), py - np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix * ix, iy * iy, ix * iy])

    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX) < 3) & (np.abs(dY) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
        for _ in range(5):
            c_, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            r = yy - A @ c_; s_ = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s_)
        return c_
    cS, cC = fit(p["S"]), fit(p["C"])
    wS, wC = np.exp(-(A[:, 1:] @ cS[1:])), np.exp(-(A[:, 1:] @ cC[1:]))
    sE, cE = e_at(E)
    S, C, So = p["S"] * wS, p["C"] * wC, SO * wS
    q_raw = qs(S, C, So).mean() / E
    q_pe = qs(S / sE, C / cE, So / sE).mean() / E           # S and C on the measured electron scale before the Q/S step
    rows[E] = dict(S=S.mean() / E, S_pe=S.mean() / E / sE, C=C.mean() / E, C_pe=C.mean() / E / cE, Q=q_raw, Q_pe=q_pe)
    r = rows[E]
    print(f"{E:4.0f} | {r['S']:6.3f} {r['S_pe']:7.3f} {dS.get(int(E), np.nan):6.3f} | {r['C']:6.3f} {r['C_pe']:7.3f} "
          f"{(0.640 if E == 100 else np.nan):6.3f} | {r['Q']:6.3f} {r['Q_pe']:9.3f} {dQ.get(int(E), np.nan):6.3f}")

# ---------------- item 2 ----------------
OFF = [1409, 1389, 1390, 1398, 1399, 1404, 1405, 1406, 1408]
print("\nitem 2: 80 GeV positrons at different table positions (core peak; 1409 is the central position)")
pos = {}
for run in OFF:
    if not os.path.exists(f"{SCR}/hodo_{run}.npz"):
        print(f"   run {run}: no dump"); continue
    e = load(run, "electron"); s, c = e["S"] / e["E"], e["C"] / e["E"]
    k = (s > 0.6) & (s < 1.5)
    if k.sum() < 300:
        print(f"   run {run}: {k.sum()} events"); continue
    pos[run] = dict(tx=num(rl[str(run)]["Table position X"]), ty=num(rl[str(run)]["Table position Y"]), S=core(s[k])[0], C=core(c[k])[0],
                    ex=float(np.median(e["ex"][k])), ey=float(np.median(e["ey"][k])), N=int(k.sum()))
ref = pos.get(1409)
for run, v in pos.items():
    rel = f"  relative to centre: S {v['S'] / ref['S']:.3f}, C {v['C'] / ref['C']:.3f}" if ref else ""
    print(f"   run {run}: table ({v['tx']:+5.0f}, {v['ty']:5.0f}) mm, shower centre ({v['ex']:+5.2f}, {v['ey']:+5.2f}) cm, N={v['N']:5d}: "
          f"S/E {v['S']:.3f}  C/E {v['C']:.3f}{rel}")

# electron scale as seen by pions: weight the positron response map by where pion (vs central-positron) energy goes
if ref and os.path.exists(f"{SCR}/chan_1432.npz") and os.path.exists(f"{SCR}/chan_1409.npz"):
    def profile(run, part):
        d = np.load(f"{SCR}/chan_{run}.npz"); meta = d["meta"]
        sel = np.ones(len(d["vals"]), bool)
        for f in SEL[part]: sel &= d[f].astype(bool)
        v = d["vals"][sel].mean(0)
        x, y = meta[:, 4] * 1.2, meta[:, 5] * 1.6                       # tower centres, same scale as the shower centre
        return x, y, meta[:, 2], np.clip(v, 0, None)
    # response at each channel position: from the positron run whose shower centre is nearest, after shifting the
    # central position by the table offset (table X -> detector X, table Y -> detector Y; signs fixed by the centres)
    centres = {run: (v["ex"], v["ey"]) for run, v in pos.items()}
    def resp_at(x, y, key):
        out = np.empty(len(x))
        for i, (a, b) in enumerate(zip(x, y)):
            run = min(centres, key=lambda r: (centres[r][0] - a) ** 2 + (centres[r][1] - b) ** 2)
            out[i] = pos[run][key] / ref[key]
        return out
    for key, isc in (("S", 0), ("C", 1)):
        xp, yp, cerp, vp = profile(1432, "pion"); xe, ye, cere, ve = profile(1409, "electron")
        mp, me = cerp == isc, cere == isc
        Rp = resp_at(xp[mp], yp[mp], key); Re = resp_at(xe[me], ye[me], key)
        sp = (vp[mp] * Rp).sum() / vp[mp].sum(); se = (ve[me] * Re).sum() / ve[me].sum()
        far = np.hypot(xp[mp] - centres[1409][0], yp[mp] - centres[1409][1]) > 5.0
        print(f"   {key}: energy-weighted relative response seen by 100 GeV pions {sp:.3f}, by central 80 GeV positrons {se:.3f}"
              f" -> pions vs positrons {sp / se:.3f};  pion {key} fraction beyond 5 cm of the beam axis: {vp[mp][far].sum() / vp[mp].sum():.2f}")
json.dump({"eresp": {str(k): v for k, v in eresp.items()}, "rows": {str(k): v for k, v in rows.items()}, "pos": {str(k): v for k, v in pos.items()}},
          open(f"{SCR}/escale_items.json", "w"), indent=1)
