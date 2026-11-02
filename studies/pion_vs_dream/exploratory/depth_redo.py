"""Pion S, C, Q/S at positron depth using the pure-attenuation slope (fit at fixed C/S), for different d = 0 references."""
import json, os, sys
import numpy as np
src = open("edepth_test.py").read(); exec(src[:src.index("d11 = json.load")])
d11 = json.load(open("dream_fig11.json")); dS = dict(zip(d11["energy_GeV"], d11["S_over_E"])); dQ = dict(zip(d11["energy_GeV"], d11["QS_over_E"]))
def fem_from(q, ehS=1.3, ehC=5.0):
    a, b = 1/ehC, 1/ehS; return (a - q*b)/(q*(1-b) - (1-a))
def qsfun(S, C, So):
    Sc = S + So; f = np.clip(fem_from(C/Sc), 0, 1); return Sc/(f + (1-f)/1.3)
refs = {}
for r in (1527, 1526, 1409, 1411, 1416, 1423, 1514):
    if os.path.exists(f"hodo_{r}.npz"):
        e = load(r, "electron"); px = cx[0]*e["hx"]+cx[1]*e["hy"]+cx[2]; py = cy[0]*e["hx"]+cy[1]*e["hy"]+cy[2]
        refs[r] = (np.median(e["ex"]-px), np.median(e["ey"]-py))
print("positron d=0 references (dX, dY in cm on the 1527 map):", {r: (round(a, 2), round(b, 2)) for r, (a, b) in refs.items()})
early = [r for r in (1409, 1411, 1416, 1423) if r in refs]
use = {"late (1527)": refs[1527], f"early ({','.join(map(str, early))})": tuple(np.mean([refs[r] for r in early], axis=0))}
for lab, (x0, y0) in use.items():
    print(f"\n=== reference {lab}: dX0 {x0:+.2f}, dY0 {y0:+.2f} cm.  Columns: old = in-sample slope (current page), new = fixed-C/S slope, depth-only part")
    print(f"{'E':>4s} | {'S/E old':>7s} {'new':>6s} {'DREAM':>6s} | {'C/E old':>7s} {'new':>6s} | {'C/S new':>7s} | {'Q/S old':>7s} {'new':>6s} {'DREAM':>6s} | {'S boost':>7s} {'C boost':>7s}")
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
        cS_old, cC_old = fit(yS, ok), fit(yC, ok)
        # pure-attenuation depth slope: fit within C/S bins (f_em fixed), average; impact terms from the global fit
        q = p["C"]/p["S"]; edges = np.quantile(q[ok], np.linspace(0, 1, 9)); acc, w = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = ok&(q>=lo)&(q<hi)
            if m.sum() > 100: acc.append(fit(yS, m)[1:3]); w.append(m.sum())
        att = np.average(np.array(acc), axis=0, weights=w) if w else cS_old[1:3]
        if not w: print(f'   {E:.0f} GeV: N={ok.sum()} too few per C/S bin, in-sample slope kept')
        def corr(sig, c, dep):
            return sig*np.exp(-(A[:, 1:3]@dep) - (A[:, 3:]@c[3:]))
        S_o, C_o, So_o = corr(p["S"], cS_old, cS_old[1:3]), corr(p["C"], cC_old, cC_old[1:3]), corr(So, cS_old, cS_old[1:3])
        S_n, C_n, So_n = corr(p["S"], cS_old, att), corr(p["C"], cC_old, att), corr(So, cS_old, att)
        bS = np.mean(p["S"]*np.exp(-(A[:,3:]@cS_old[3:])))/S_n.mean(); bC = np.mean(p["C"]*np.exp(-(A[:,3:]@cC_old[3:])))/C_n.mean()
        print(f"{E:4.0f} | {S_o.mean()/E:7.3f} {S_n.mean()/E:6.3f} {dS.get(int(E), np.nan):6.3f} | {C_o.mean()/E:7.3f} {C_n.mean()/E:6.3f} | "
              f"{C_n.mean()/S_n.mean():7.3f} | {qsfun(S_o, C_o, So_o).mean()/E:7.3f} {qsfun(S_n, C_n, So_n).mean()/E:6.3f} {dQ.get(int(E), np.nan):6.3f} | {bS:7.3f} {bC:7.3f}")
