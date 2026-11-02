"""Angle-free depth correction: the slope of ln S vs displacement (dX, dY) is measured in the pion data at fixed C/S
(fixed em fraction -> attenuation only); C uses that slope scaled by lambda_S / lambda_quartz and lambda_S / lambda_plastic
from the 90-degree positron runs. No tilt angle enters. Implied effective angle = atan(|d| / depth), depth from lambda_S."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final.py")).read(); exec(src[:src.index("def pion(")])
d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
dres = {int(e): (s, c, q) for e, s, c, q in zip(d32["energy_GeV"], d32["S_res_pct"], d32["C_res_pct"], d32["QS_res_pct"])}
L = {t: (ATT[t]["lam"], ATT[t]["elam"]) for t in ("Sci", "CerQuartz", "CerPlastic")}
def coreres(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return 100 * s / m

def run_af(run, lam):
    d = np.load(f"{SCR}/hodoq_{run}.npz"); E = float(d["energy"]); sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"]>=1)&(d["TTU_Hodo_nHitX"]<=3)&(d["TTU_Hodo_nHitY"]>=1)&(d["TTU_Hodo_nHitY"]<=3)&(d["FERS_energyMix_Sci_calib_sum"]>0.3*E)
    g = lambda k: d[k][sel].astype(float)
    S, Cq, Cp, So = g("FERS_energyMix_Sci_calib_sum"), g("FERS_energyMix_CerQuartz_calib_sum"), g("FERS_energyMix_CerPlastic_calib_sum"), g("FERS_energyMix_Sci_calib_sum_OuterRing")
    C = Cq + Cp
    px = cx[0]*g("hodoX_c")+cx[1]*g("hodoY_c")+cx[2]; py = cy[0]*g("hodoX_c")+cy[1]*g("hodoY_c")+cy[2]
    dX, dY = g("FERS_energyMix_Sci_calib_weighted_center_X")-px-x0, g("FERS_energyMix_Sci_calib_weighted_center_Y")-py-y0
    ix, iy = px-np.median(px), py-np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix*ix, iy*iy, ix*iy])
    ok = (np.abs(dX)<3)&(np.abs(dY)<3)&(np.abs(ix)<3)&(np.abs(iy)<3)&(S>0)&(C>0)
    yS = np.log(np.maximum(S, 1e-3))
    def rfit(y, m):
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            r = y - A @ c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m & (np.abs(r) < 2.5*s)
        return c
    # slope vector at fixed C/S: fit within quantile bins of C/S, signal-weighted average; spread over bins -> error
    q = C / S; nb = int(np.clip(ok.sum() // 400, 3, 8)); edges = np.quantile(q[ok], np.linspace(0, 1, nb + 1))
    ks, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (q >= lo) & (q < hi)
        if m.sum() > 60: ks.append(rfit(yS, m)[1:3]); ws.append(m.sum())
    ks = np.array(ks); k = np.average(ks, axis=0, weights=ws); ek = ks.std(axis=0) / np.sqrt(len(ks))
    def corr(sig, kk):
        y = np.log(np.maximum(sig, 1e-3)) - A[:, 1:3] @ kk; m = ok.copy()
        idx = [0, 3, 4, 5, 6, 7]
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m][:, idx], y[m], rcond=None)
            r = y - A[:, idx] @ c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m & (np.abs(r) < 2.5*s)
        return sig * np.exp(-(A[:, 1:3] @ kk) - (A[:, 3:] @ c[1:]))
    kq, kp = k * lam["Sci"] / lam["CerQuartz"], k * lam["Sci"] / lam["CerPlastic"]
    Sx, Qx, Px, Sox = corr(S, k), corr(Cq, kq), corr(Cp, kp), corr(So, k)
    Cx = Qx + Px; sE, cE = e_at(E); Se, Ce, Soe = Sx / sE, Cx / cE, Sox / sE
    QS = qsfun(Se, Ce, Soe); DR = (Se - 0.29 * Ce) / 0.71
    dm = np.array([np.median(dX[ok]), np.median(dY[ok])])
    depth = lam["Sci"] * 100 * (k @ dm)                      # cm deeper than the positrons, from the attenuation of S
    ang = np.degrees(np.arctan(np.hypot(*dm) / depth)) if depth > 0 else np.nan
    res = {key: (100 * v.std() / v.mean(), coreres(v)) for key, v in (("S", Se), ("C", Ce), ("DR", DR), ("QS", QS))}
    return E, dict(N=int(sel.sum()), k=k.tolist(), ek=ek.tolist(), dmed=dm.tolist(), depth=depth, angle=ang,
                   S=Se.mean() / E, C=Ce.mean() / E, Q=QS.mean() / E, CS=Cx.mean() / Sx.mean(), res=res)

out = {}
print("slope of ln S vs displacement at fixed C/S (per cm), pion median displacement vs positrons, implied depth and angle (lambda_S = 3.07 m)")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    E, v = run_af(run, {t: L[t][0] for t in L})
    # systematic: lambda ratio uncertainties and the hottest-channel lengths
    _, v_core = run_af(run, {"Sci": 2.5, "CerQuartz": 3.3, "CerPlastic": 5.2})
    v["core"] = {k2: v_core[k2] for k2 in ("S", "C", "Q", "CS")}
    out[str(int(E))] = v; dd = dres.get(int(E))
    print(f"{E:4.0f} N {v['N']:5d} | k = ({v['k'][0]:.3f}+-{v['ek'][0]:.3f}, {v['k'][1]:.3f}+-{v['ek'][1]:.3f}) | d = ({v['dmed'][0]:.2f}, {v['dmed'][1]:.2f}) cm"
          f" -> {v['depth']:5.1f} cm deeper, angle {v['angle']:.1f} deg | S/e {v['S']:.3f} C/e {v['C']:.3f} C/S {v['CS']:.3f} Q/S/e {v['Q']:.3f}"
          f" (core-ch. lambda: {v['core']['S']:.3f} {v['core']['C']:.3f} {v['core']['Q']:.3f}) | res rms/core QS {v['res']['QS'][0]:.1f}/{v['res']['QS'][1]:.1f}"
          f" S {v['res']['S'][0]:.1f} C {v['res']['C'][0]:.1f}" + (f" | DREAM S {dd[0]} C {dd[1]} QS {dd[2]}" if dd else ""))
json.dump(out, open(f"{SCR}/pion_anglefree.json", "w"), indent=1, default=float)
