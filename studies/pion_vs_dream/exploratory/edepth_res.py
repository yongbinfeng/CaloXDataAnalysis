import numpy as np
src = open("edepth_test.py").read()
exec(src[:src.index("cx, cy = fit_map")])
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6: break
        m, s = m2, s2
    return m, s
def qs(S, C, SO):
    Sc = S + SO; b, a = 1 / 1.3, 0.2
    f = np.clip((a - C / Sc * b) / (C / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    return Sc / (f + (1 - f) * b)
cx, cy = fit_map(load(1527, "electron"))
print("resolution sigma_rms/mean / core [%] and mean / E; correction normalised to the pion median (as in the plots) vs to positron depth (d = 0)")
print(f"{'E':>4s} | {'S':>16s} | {'DR, median':>22s} {'DR, e+ depth':>22s} | {'Q/S, median':>22s} {'Q/S, e+ depth':>22s}")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    d = np.load(f"{SCR}/hodo_{run}.npz")
    # outer ring for the same events (load() does not return it): redo the selection
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL["pion"]: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"] >= 1) & (d["TTU_Hodo_nHitX"] <= 3) & (d["TTU_Hodo_nHitY"] >= 1) & (d["TTU_Hodo_nHitY"] <= 3) & (d["FERS_energyMix_Sci_calib_sum"] > 0.3 * E)
    SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py; ix, iy = px - np.median(px), py - np.median(py)
    def build(vx, vy):
        return np.column_stack([np.ones(len(vx)), vx, vy, ix, iy, ix * ix, iy * iy, ix * iy])
    A0 = build(dX, dY)
    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX) < 3) & (np.abs(dY) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A0[m], yy[m], rcond=None)
            r = yy - A0 @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    out = {}
    for lab, A in (("median", build(dX - np.median(dX), dY - np.median(dY))), ("edepth", A0)):
        wS, wC = np.exp(-(A[:, 1:] @ cS[1:])), np.exp(-(A[:, 1:] @ cC[1:]))
        S, C, So = p["S"] * wS, p["C"] * wC, SO * wS
        for k, x in (("S", S), ("DR", (S - 0.29 * C) / 0.71), ("QS", qs(S, C, So))):
            m_, s_ = core(x); out[(lab, k)] = (x.std() / x.mean() * 100, s_ / m_ * 100, x.mean() / E)
    f = lambda t: f"{t[0]:5.1f}/{t[1]:5.1f} ({t[2]:.3f})"
    print(f"{E:4.0f} | {f(out[('median', 'S')])[:11]:>16s} | {f(out[('median', 'DR')]):>22s} {f(out[('edepth', 'DR')]):>22s} | "
          f"{f(out[('median', 'QS')]):>22s} {f(out[('edepth', 'QS')]):>22s}")
