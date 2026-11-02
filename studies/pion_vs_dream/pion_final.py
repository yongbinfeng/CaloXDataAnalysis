"""Pion response vs DREAM after (i) normalising to the depth of the positrons taken just before the pion scan,
using the attenuation lengths measured on the 2025 90-degree positron runs (Sci, quartz and plastic separately),
and (ii) dividing by our own positron response at the same energy (test 1)."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("edepth_test.py")).read(); exec(src[:src.index("d11 = json.load")])    # load (electrons), fit_map, cx, cy, SEL
d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json")); dS = dict(zip(d11["energy_GeV"], d11["S_over_E"])); dQ = dict(zip(d11["energy_GeV"], d11["QS_over_E"]))
ATT = json.load(open(f"{SCR}/att90_Y1064.json"))["FERS"]
T = np.tan(np.radians(3.3))
# The pion displacement relative to positrons has dY/dX = 0.71 at all energies and dX, dY correlate event by event:
# a vertical tilt of about 2.3 deg besides the logged 3.3 deg. Option 1 projects (dX, dY) on that direction,
# depth = (dX + r dY) / ((1 + r^2) tan 3.3); option 0 uses dX alone.
R_TILT = 0.71
def fem_from(q, ehS=1.3, ehC=5.0):
    a, b = 1/ehC, 1/ehS; return (a - q*b)/(q*(1-b) - (1-a))
def qsfun(S, C, So):
    Sc = S + So; f = np.clip(fem_from(C/Sc), 0, 1); return Sc/(f + (1-f)/1.3)
def core(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return m, s
def off(e):
    px = cx[0]*e["hx"]+cx[1]*e["hy"]+cx[2]; py = cy[0]*e["hx"]+cy[1]*e["hy"]+cy[2]; return e["ex"]-px, e["ey"]-py

# positron reference: depth (displacement) and response, from the positron runs of the same period / energies
EARLY = (1409, 1411, 1416, 1423)
x0 = np.mean([np.median(off(load(r, "electron"))[0]) for r in EARLY]); y0 = np.mean([np.median(off(load(r, "electron"))[1]) for r in EARLY])
POS = {20: 1423, 30: 1416, 40: 1514, 60: 1526, 80: 1409, 100: 1527, 120: 1411}
esc = {}
for E, r in POS.items():
    e = load(r, "electron"); s, c = e["S"]/E, e["C"]/E; k = (s > 0.7) & (s < 1.4); esc[E] = (core(s[k])[0], core(c[k])[0])
lE = np.log(sorted(esc)); e_at = lambda E: (np.interp(np.log(E), lE, [esc[x][0] for x in sorted(esc)]), np.interp(np.log(E), lE, [esc[x][1] for x in sorted(esc)]))

def pion(run, lamS, lamQ, lamP, yd):
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
    def impact(sig, k):
        # impact-position terms, fitted with the depth term fixed to the measured attenuation
        y = np.log(np.maximum(sig, 1e-3)) - (A[:, 1:3] @ k); m = ok.copy()
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m][:, [0, 3, 4, 5, 6, 7]], y[m], rcond=None)
            r = y - A[:, [0, 3, 4, 5, 6, 7]] @ c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m & (np.abs(r) < 2.5*s)
        return sig*np.exp(-(A[:, 1:3] @ k) - (A[:, 3:] @ c[1:]))
    r = R_TILT if yd else 0.0
    kv = lambda lam: np.array([1, r]) / (1 + r * r) / (T * lam * 100)     # per cm of displacement (dX, dY)
    Sx, Qx, Px, Sox = impact(S, kv(lamS)), impact(Cq, kv(lamQ)), impact(Cp, kv(lamP)), impact(So, kv(lamS))
    Cx = Qx + Px; sE, cE = e_at(E)
    return E, dict(S=Sx.mean()/E/sE, C=Cx.mean()/E/cE, Q=qsfun(Sx/sE, Cx/cE, Sox/sE).mean()/E, CS=Cx.mean()/Sx.mean(),
                   Sraw=S.mean()/E, Craw=C.mean()/E, fq=Cq.sum()/C.sum(), N=int(sel.sum()))

RUNS = [r for r in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434) if os.path.exists(f"{SCR}/hodoq_{r}.npz")]
L = {t: (ATT[t]["lam"], ATT[t]["elam"]) for t in ("Sci", "CerQuartz", "CerPlastic")}
out = {"esc": {str(k): v for k, v in esc.items()}, "ref0": [x0, y0], "lam": L, "rows": {}}
print(f"lambda S {L['Sci'][0]:.2f}+-{L['Sci'][1]:.2f}, quartz {L['CerQuartz'][0]:.2f}+-{L['CerQuartz'][1]:.2f}, plastic {L['CerPlastic'][0]:.2f}+-{L['CerPlastic'][1]:.2f} m; "
      f"positron depth reference ({x0:+.2f}, {y0:+.2f}) cm from runs {EARLY}")
for yd in (0, 1):
    print(f"\n depth from {'dX and dY (projected, r = 0.71)' if yd else 'dX only'}:  values on our electron scale; +- from the attenuation-length uncertainties")
    print(f"{'E':>4s} {'N':>6s} {'f_quartz':>8s} | {'S/e':>13s} {'DREAM':>6s} | {'C/e':>13s} {'DREAM':>6s} | {'C/S':>5s} | {'Q/S/e':>13s} {'DREAM':>6s}")
    for run in RUNS:
        E, nom = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], yd)
        var = []
        for t in ("Sci", "CerQuartz", "CerPlastic"):
            for sgn in (+1, -1):
                lam = {k: v[0] for k, v in L.items()}; lam[t] = max(lam[t] + sgn * L[t][1], 0.5)
                var.append(pion(run, lam["Sci"], lam["CerQuartz"], lam["CerPlastic"], yd)[1])
        err = {k: np.sqrt(sum(((v[k] - nom[k]) ** 2) for v in var) / 2) for k in ("S", "C", "Q")}
        out["rows"].setdefault(str(yd), {})[str(int(E))] = dict(nom, eS=err["S"], eC=err["C"], eQ=err["Q"], run=run)
        print(f"{E:4.0f} {nom['N']:6d} {nom['fq']:8.2f} | {nom['S']:6.3f}+-{err['S']:.3f} {dS.get(int(E), np.nan):6.3f} | {nom['C']:6.3f}+-{err['C']:.3f} "
              f"{(0.640 if E == 100 else np.nan):6.3f} | {nom['CS']:5.3f} | {nom['Q']:6.3f}+-{err['Q']:.3f} {dQ.get(int(E), np.nan):6.3f}")
json.dump(out, open(f"{SCR}/pion_final.json", "w"), indent=1, default=float)
