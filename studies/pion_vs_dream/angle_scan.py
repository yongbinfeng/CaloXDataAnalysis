"""Pion response and resolution vs an assumed constant total tilt (direction fixed by the measured dY/dX = 0.71)."""
import json
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final_res.py")).read(); exec(src[:src.index("d32 = json.load")])     # pion() with ARR, globals T, R_TILT
d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json")); dres = {int(e): q for e, q in zip(d32["energy_GeV"], d32["QS_res_pct"])}
def coreres(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return 100 * s / m
out = {}
ANG = (2.0, 2.5, 3.0, 3.5, 4.0)
print("total tilt -> per energy: S/e, C/e, C/S, Q/S/e, Q/S rms/core  (DREAM S/E, Q/S/E, Q/S res)")
for run in (1438, 1437, 1429, 1432, 1433, 1434):
    for th in ANG:
        T = np.tan(np.radians(th)) / np.sqrt(1 + R_TILT ** 2)        # X component of the tilt, so that the total is th
        ARR.clear(); E, v = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], 1)
        S, C, So = ARR["final"]; Q = qsfun(S, C, So)
        out.setdefault(str(int(E)), {})[str(th)] = dict(S=v["S"], C=v["C"], CS=v["CS"], Q=v["Q"], Qrms=100*Q.std()/Q.mean(), Qcore=coreres(Q))
        r = out[str(int(E))][str(th)]
        print(f"{E:4.0f} GeV  {th:.1f} deg: S/e {r['S']:.3f}  C/e {r['C']:.3f}  C/S {r['CS']:.3f}  Q/S/e {r['Q']:.3f}  Q/S {r['Qrms']:.1f}/{r['Qcore']:.1f}"
              + (f"   | DREAM {dS.get(int(E), np.nan):.3f} / {dQ.get(int(E), np.nan):.3f} / {dres[int(E)]:.1f}" if th == ANG[0] and int(E) in dres else ""))
json.dump(out, open(f"{SCR}/angle_scan.json", "w"), indent=1)
