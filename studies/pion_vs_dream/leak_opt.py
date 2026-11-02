"""Optimise DREAM's lateral-leakage correction S_c = S (1 + x f_out) for the CaloX pions: the factor x and the region used
for f_out, and how the outer-ring sum is treated by the impact-position correction. Everything else as final_analysis.py
(2.5 deg tilt, measured lambda per fibre type, positron-depth reference, depth cut, own electron scale)."""
import json, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final.py")).read(); exec(src[:src.index("def pion(")])         # e_at, cx, cy, x0, y0, SEL, qsfun, ATT
LAM = {t: ATT[t]["lam"] * 100 for t in ("Sci", "CerQuartz", "CerPlastic")}           # cm
R_TILT, TILT, ZCUT = 0.71, 2.5, 70.0
TX = np.tan(np.radians(TILT)) / np.sqrt(1 + R_TILT ** 2)
AXIS = (2.32, 1.23)                                                                   # pion shower axis in the module frame [cm]
RINGS = {**{f"square>={a}": ("sq", a) for a in (6.0, 8.0, 9.6, 11.0, 12.5)}, **{f"axis r>={a}": ("ax", a) for a in (8.0, 10.0, 12.0, 14.0)}}
XS = np.round(np.arange(0.0, 3.01, 0.25), 2)

def coreres(v):
    m, s = np.median(v), 1.4826*np.median(np.abs(v - np.median(v)))
    for _ in range(5):
        k = np.abs(v - m) < 2*s; m, s = v[k].mean(), v[k].std()
    return 100 * s / m

def prepare(run):
    d = np.load(f"{SCR}/hodoc_{run}.npz"); E = float(d["energy"]); sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"]>=1)&(d["TTU_Hodo_nHitX"]<=3)&(d["TTU_Hodo_nHitY"]>=1)&(d["TTU_Hodo_nHitY"]<=3)&(d["FERS_energyMix_Sci_calib_sum"]>0.3*E)
    g = lambda k: d[k][sel].astype(float)
    S, Cq, Cp = g("FERS_energyMix_Sci_calib_sum"), g("FERS_energyMix_CerQuartz_calib_sum"), g("FERS_energyMix_CerPlastic_calib_sum")
    So_fw = g("FERS_energyMix_Sci_calib_sum_OuterRing")                                     # the framework's outer ring
    V = d["sci_vals"][sel].astype(float); meta = d["sci_meta"]; x, y = meta[:, 2] * 1.2, meta[:, 3] * 1.6
    rings = {}
    for name, (kind, a) in RINGS.items():
        m = (np.maximum(np.abs(x), np.abs(y)) >= a) if kind == "sq" else (np.hypot(x - AXIS[0], y - AXIS[1]) >= a)
        rings[name] = V[:, m].sum(1)
    px = cx[0]*g("hodoX_c")+cx[1]*g("hodoY_c")+cx[2]; py = cy[0]*g("hodoX_c")+cy[1]*g("hodoY_c")+cy[2]
    dX, dY = g("FERS_energyMix_Sci_calib_weighted_center_X")-px-x0, g("FERS_energyMix_Sci_calib_weighted_center_Y")-py-y0
    z = (dX + R_TILT * dY) / ((1 + R_TILT ** 2) * TX)
    ix, iy = px - np.median(px), py - np.median(py)
    A = np.column_stack([np.ones(len(z)), ix, iy, ix*ix, iy*iy, ix*iy])
    ok = (np.abs(dX)<3)&(np.abs(dY)<3)&(np.abs(ix)<3)&(np.abs(iy)<3)&(S>0)&(Cq+Cp>0)
    att = {t: np.exp(-z / LAM[t]) for t in LAM}
    def impact_factor(sig):
        y_ = np.log(np.maximum(sig, 1e-3)); m = ok.copy()
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], y_[m], rcond=None)
            r = y_ - A @ c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m & (np.abs(r) < 2.5*s)
        return np.exp(-(A[:, 1:] @ c[1:]))
    sE, cE = e_at(E); keep = z < ZCUT
    fS = impact_factor(S * att["Sci"])
    base = dict(E=E, keep=keep,
                S_imp=S*att["Sci"]*fS/sE, C_imp=(Cq*att["CerQuartz"]*impact_factor(Cq*att["CerQuartz"]) + Cp*att["CerPlastic"]*impact_factor(Cp*att["CerPlastic"]))/cE,
                S_att=S*att["Sci"]/sE, C_att=(Cq*att["CerQuartz"] + Cp*att["CerPlastic"])/cE)
    # ring treatments: A own impact fit (current); B same factor as S; C attenuation only
    base["rings"] = {name: dict(A=R*att["Sci"]*impact_factor(np.maximum(R, 1e-3)*att["Sci"])/sE, B=R*att["Sci"]*fS/sE, C=R*att["Sci"]/sE)
                     for name, R in rings.items()}
    base["fw_check"] = float(np.abs(rings["square>=9.6"] - So_fw).max())
    return base

def evaluate(b, ring, treat, x):
    k = b["keep"]
    S, C = (b["S_att"], b["C_att"]) if treat == "C" else (b["S_imp"], b["C_imp"])
    R = b["rings"][ring][treat]
    if x == "eq10":       # DREAM Eq. (10): S_corr = S [1 + 0.2 (1 - exp(-y f_out / sqrt(S)))], y = 60, S in GeV
        Sk, fout = S[k], R[k] / np.maximum(S[k], 1e-3)
        q = qsfun(Sk, C[k], Sk * 0.2 * (1 - np.exp(-60 * fout / np.sqrt(np.maximum(Sk, 1e-3)))))
    else:
        q = qsfun(S[k], C[k], x * R[k])
    return 100 * q.std() / q.mean(), coreres(q), q.mean() / b["E"]

RUNS = [int(r) for r in sys.argv[1:]] or [1438, 1437, 1429, 1432, 1433, 1434]
res = {}
for run in RUNS:
    b = prepare(run); E = int(b["E"])
    print(f"run {run} ({E} GeV): framework outer ring == square>=9.6 to {b['fw_check']:.1e} GeV")
    for ring in RINGS:
        for treat in "ABC":
            for x in list(XS) + ["eq10"]:
                res.setdefault(str(E), {}).setdefault(ring, {}).setdefault(treat, {})[str(x)] = evaluate(b, ring, treat, x)
    cur = res[str(E)]["square>=9.6"]["A"]["1.0"]
    print(f"   current (square>=9.6, A, x=1): Q/S rms {cur[0]:.2f}  core {cur[1]:.2f}  mean {cur[2]:.3f}")
    best = min(((ring, t, x, v) for ring, d1 in res[str(E)].items() for t, d2 in d1.items() for x, v in d2.items()), key=lambda r: r[3][0])
    print(f"   best rms: {best[0]}, {best[1]}, x={best[2]}: rms {best[3][0]:.2f} core {best[3][1]:.2f} mean {best[3][2]:.3f}")
json.dump(res, open(f"{SCR}/leak_opt.json", "w"))
