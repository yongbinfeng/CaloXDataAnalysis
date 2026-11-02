"""DREAM-like cut on the depth of the light production (their Section 4.2: centre of gravity in the first metre, i.e. ~70 cm
beyond the 30 cm depth of their 40 GeV calibration electrons). Same correction as final_analysis.py (2.5 deg)."""
import json
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final_res.py")).read()
src = src.replace('ARR["raw"] = (S/sE, C/cE, So/sE)', 'ARR["raw"] = (S/sE, C/cE, So/sE); ARR["z"] = (dX + R_TILT*dY) / ((1 + R_TILT**2) * T)')
exec(src[:src.index("d32 = json.load")])
T = np.tan(np.radians(2.5)) / np.sqrt(1 + R_TILT ** 2)
def coreres(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return 100 * s / m
CUTS = (None, 70, 50, 40)
out = {}
print("depth z relative to the positrons [cm]; cut z < 70 cm = DREAM's 'first metre'")
for run in (1438, 1437, 1429, 1432, 1433, 1434):
    ARR.clear(); E, _ = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], 1)
    S, C, So = ARR["final"]; z = ARR["z"]
    print(f"{E:4.0f} GeV: z median {np.median(z):5.1f} cm, 90% {np.percentile(z, 90):5.1f}, 99% {np.percentile(z, 99):5.1f}")
    for cut in CUTS:
        k = np.ones(len(z), bool) if cut is None else z < cut
        s, c, so = S[k], C[k], So[k]; dr = (s - 0.29*c)/0.71; qs = qsfun(s, c, so)
        r = {n: (100*v.std()/v.mean(), coreres(v), v.mean()/E) for n, v in (("S", s), ("C", c), ("DR", dr), ("QS", qs))}
        out.setdefault(str(int(E)), {})[str(cut)] = dict(kept=float(k.mean()), **{n: list(v) for n, v in r.items()})
        print(f"     z < {str(cut):>4s}: kept {100*k.mean():5.1f}% | rms S {r['S'][0]:5.1f} C {r['C'][0]:5.1f} DR {r['DR'][0]:5.1f} Q/S {r['QS'][0]:5.1f}"
              f" | core Q/S {r['QS'][1]:4.1f} | response S {r['S'][2]:.3f} Q/S {r['QS'][2]:.3f}")
json.dump(out, open(f"{SCR}/depth_cut.json", "w"), indent=1)
