"""Does DR / Q/S improve if the Cherenkov signal is NOT corrected for attenuation (S corrected as in final_analysis.py)?
C variants: fully corrected (attenuation + impact position), impact position only (no attenuation), raw.
Same events as final_analysis.py (2.5 deg tilt, depth cut z < 70 cm), electron scale."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json
import numpy as np

src = open(code("pion_final_res.py")).read()
src = src.replace('ARR["raw"] = (S/sE, C/cE, So/sE)', 'ARR["raw"] = (S/sE, C/cE, So/sE); ARR["z"] = (dX + R_TILT*dY) / ((1 + R_TILT**2) * T)')
exec(src[:src.index("d32 = json.load")])                   # pion(): S, C, S_outer per event (ARR), L = attenuation lengths
T = np.tan(np.radians(2.5)) / np.sqrt(1 + R_TILT ** 2)
ZCUT, NOATT = 70.0, 1e9                                     # lambda = 1e9 m: no attenuation correction


def coreres(x):
    m, s = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2 * s; m, s = x[k].mean(), x[k].std()
    return 100 * s / m


def best_chi(s, c):
    grid = np.arange(-0.2, 0.81, 0.01)
    r = [((s - x * c) / (1 - x)).std() / ((s - x * c) / (1 - x)).mean() for x in grid]
    i = int(np.argmin(r)); return grid[i], 100 * r[i]


out = {}
print("S fully corrected in all cases. DR: chi = 0.29; best chi: minimises rms/mean of (S - chi C)/(1 - chi). [%] = rms/mean, core")
for run in (1439, 1438, 1437, 1429, 1432, 1433, 1434):
    ARR.clear(); E, _ = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], 1)
    keep = ARR["z"] < ZCUT; S, Cfull, So = (a[keep] for a in ARR["final"]); Craw = ARR["raw"][1][keep]
    ARR.clear(); pion(run, L["Sci"][0], NOATT, NOATT, 1)
    Cimp = ARR["final"][1][keep]                            # same events (same selection order), C without attenuation correction
    rows = {}
    for lab, C in (("C corrected", Cfull), ("C impact only", Cimp), ("C raw", Craw)):
        dr = (S - 0.29 * C) / 0.71; qs = qsfun(S, C, So); bc, br = best_chi(S, C)
        rows[lab] = dict(C=[100 * C.std() / C.mean(), coreres(C)], DR=[100 * dr.std() / dr.mean(), coreres(dr), dr.mean() / E],
                         QS=[100 * qs.std() / qs.mean(), coreres(qs), qs.mean() / E], chi=[bc, br])
    out[str(int(E))] = rows
    print(f"{E:4.0f} GeV")
    for lab, r in rows.items():
        print(f"   {lab:14s}: C {r['C'][0]:5.1f} | DR {r['DR'][0]:5.1f} / {r['DR'][1]:4.1f} (mean {r['DR'][2]:.3f}) | best chi {r['chi'][0]:.2f} -> "
              f"{r['chi'][1]:5.1f} | Q/S {r['QS'][0]:5.1f} / {r['QS'][1]:4.1f} (mean {r['QS'][2]:.3f})")
json.dump(out, open(f"{DATA_DIR}/c_attenuation_check.json", "w"), indent=1, default=float)
