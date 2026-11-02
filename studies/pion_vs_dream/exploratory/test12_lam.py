import sys
YD = float(sys.argv[1])
"""Test 1: pion response (at early-period positron depth, pure-attenuation slope) on our own measured electron scale, vs DREAM.
Test 2: 80 GeV positron response at off-centre table positions relative to the central run."""
import json, os
import numpy as np
src = open("edepth_test.py").read(); exec(src[:src.index("d11 = json.load")])
rl = json.load(open("/data1/yofeng/CaloX/CaloXDataAnalysis/data/Runlist.json"))
num = lambda s: float(str(s).replace("mm", "").strip())
d11 = json.load(open("dream_fig11.json")); dS = dict(zip(d11["energy_GeV"], d11["S_over_E"])); dQ = dict(zip(d11["energy_GeV"], d11["QS_over_E"]))
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

# ---- test 1: electron scale at the central position
POS = {20: 1423, 30: 1416, 40: 1514, 60: 1526, 80: 1409, 100: 1527, 120: 1411}
esc = {}
print("TEST 1. positron core peak, central position:")
for E, r in POS.items():
    e = load(r, "electron"); s, c = e["S"]/E, e["C"]/E; k = (s > 0.7) & (s < 1.4)
    esc[E] = (core(s[k])[0], core(c[k])[0])
    print(f"   {E:4d} GeV run {r} ({'early' if r < 1500 else 'late '}): S/E {esc[E][0]:.3f}  C/E {esc[E][1]:.3f}  C/S {esc[E][1]/esc[E][0]:.3f}")
lE = np.log(sorted(esc)); e_at = lambda E: (np.interp(np.log(E), lE, [esc[x][0] for x in sorted(esc)]), np.interp(np.log(E), lE, [esc[x][1] for x in sorted(esc)]))
x0, y0 = np.mean([np.median(off(load(r, "electron"))[0]) for r in (1409, 1411, 1416, 1423)]), np.mean([np.median(off(load(r, "electron"))[1]) for r in (1409, 1411, 1416, 1423)])
print(f"\n   [YD={YD}] pions at early positron depth (dX0 {x0:+.2f}, dY0 {y0:+.2f}), pure-attenuation slope; '/e' = divided by our e+ response at that energy")
print(f"{'E':>4s} | {'S/E':>6s} {'S/e':>6s} {'DREAM':>6s} | {'C/E':>6s} {'C/e':>6s} {'DREAM':>6s} | {'Q/S':>6s} {'Q/S /e':>6s} {'DREAM':>6s}")
rows = {}
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    d = np.load(f"hodo_{run}.npz"); sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"]>=1)&(d["TTU_Hodo_nHitX"]<=3)&(d["TTU_Hodo_nHitY"]>=1)&(d["TTU_Hodo_nHitY"]<=3)&(d["FERS_energyMix_Sci_calib_sum"]>0.3*E)
    So = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    px = cx[0]*p["hx"]+cx[1]*p["hy"]+cx[2]; py = cy[0]*p["hx"]+cy[1]*p["hy"]+cy[2]
    dX, dY = p["ex"]-px-x0, p["ey"]-py-y0; ix, iy = px-np.median(px), py-np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix*ix, iy*iy, ix*iy])
    ok = (np.abs(dX)<3)&(np.abs(dY)<3)&(np.abs(ix)<3)&(np.abs(iy)<3)&(p["S"]>0)&(p["C"]>0)
    yS, yC = np.log(np.maximum(p["S"],1e-3)), np.log(np.maximum(p["C"],1e-3))
    def fit(y, m):
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            r_ = y-A@c; s = 1.4826*np.median(np.abs(r_[m]-np.median(r_[m]))); m = m&(np.abs(r_)<2.5*s)
        return c
    cS, cC = fit(yS, ok), fit(yC, ok)
    q = p["C"]/p["S"]; edges = np.quantile(q[ok], np.linspace(0, 1, 9)); acc, w = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok&(q>=lo)&(q<hi)
        if m.sum() > 100: acc.append(fit(yS, m)[1:3]); w.append(m.sum())
    att = np.average(np.array(acc), axis=0, weights=w) if w else cS[1:3]
    T = np.tan(np.radians(3.3)); kS, kC = 1/(T*307.0), 1/(T*650.0)       # per cm of dX: lambda_S 3.07 m, lambda_C ~6.5 m (2025 FERS, Y~1064 row)
    attS = np.array([kS, kS*YD]); attC = np.array([kC, kC*YD])          # YD = 1: dY treated as depth too; 0: dX only
    corr = lambda sig, c, a: sig*np.exp(-(A[:, 1:3]@a) - (A[:, 3:]@c[3:]))
    S, C, So_ = corr(p["S"], cS, attS), corr(p["C"], cC, attC), corr(So, cS, attS)
    sE, cE = e_at(E)
    r = dict(S=S.mean()/E, Se=S.mean()/E/sE, C=C.mean()/E, Ce=C.mean()/E/cE, Q=qsfun(S, C, So_).mean()/E, Qe=qsfun(S/sE, C/cE, So_/sE).mean()/E)
    rows[E] = r
    print(f"{E:4.0f} | {r['S']:6.3f} {r['Se']:6.3f} {dS.get(int(E), np.nan):6.3f} | {r['C']:6.3f} {r['Ce']:6.3f} {(0.640 if E == 100 else np.nan):6.3f} | "
          f"{r['Q']:6.3f} {r['Qe']:6.3f} {dQ.get(int(E), np.nan):6.3f}")

