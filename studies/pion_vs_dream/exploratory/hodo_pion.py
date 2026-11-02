import numpy as np
exec(open("hodo_ana.py").read().split("e = get(1527")[0])
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6: break
        m, s = m2, s2
    return m, s
cx, cy = np.load("hodo_map.npy")
p = get(1432, "pion"); E = p["E"]
k = p["S"] > 0.3 * E
for key in list(p): 
    if isinstance(p[key], np.ndarray): p[key] = p[key][k]
px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
dX, dY = p["ex"] - px, p["ey"] - py                         # shower centre minus mapped impact point
mad = lambda v: 1.4826 * np.median(np.abs(v - np.median(v)))
print(f"pions 100 GeV, {len(dX)} events (pion sel., 1-3 hodo fibres per axis, S>0.3E; no hole veto needed)")
print(f"   spread of the shower centre: X {mad(p['ex']) * 10:.1f} mm, Y {mad(p['ey']) * 10:.1f} mm  ->  displacement from impact point: X {mad(dX) * 10:.1f} mm, Y {mad(dY) * 10:.1f} mm")
print(f"   median displacement: X {np.median(dX) * 10:+.1f} mm, Y {np.median(dY) * 10:+.1f} mm  (positrons: 0 by construction)")
S = p["S"] / E
for name, v in (("displacement X", dX), ("displacement Y", dY), ("shower centre X (no hodoscope)", p["ex"] - np.median(p["ex"]))):
    edges = np.percentile(v, [0, 10, 25, 50, 75, 90, 100])
    print(f"   S/E vs {name}: " + " ".join(f"[{a * 10:+.0f},{b * 10:+.0f}) mm {np.median(S[(v >= a) & (v < b)]):.3f}" for a, b in zip(edges[:-1], edges[1:])))
# position correction using the displacement, vs using the shower-centre offset alone
def corrected(var_x, var_y):
    def fit(sig):
        y = np.log(np.maximum(sig, 1e-3)); m = (np.abs(var_x) < 3) & (np.abs(var_y) < 3)
        A = np.column_stack([np.ones(len(y)), var_x, var_y])
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            r = y - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    wS, wC = np.exp(-(cS[1] * var_x + cS[2] * var_y)), np.exp(-(cC[1] * var_x + cC[2] * var_y))
    return p["S"] * wS, p["C"] * wC, p["SO"] * wS, cS, cC
def views(S_, C_, SO_):
    Sc = S_ + SO_; b, a = 1 / 1.3, 0.2
    f = np.clip((a - C_ / Sc * b) / (C_ / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    return {"S": S_ / E, "C": C_ / E, "DR": (S_ - 0.29 * C_) / 0.71 / E, "Q/S": Sc / (f + (1 - f) * b) / E}
ex0, ey0 = p["ex"] - np.median(p["ex"]), p["ey"] - np.median(p["ey"])
rows = [("no correction", views(p["S"], p["C"], p["SO"]), None)]
for lab, vx, vy in (("shower-centre offset", ex0, ey0), ("displacement from hodoscope", dX - np.median(dX), dY - np.median(dY))):
    S_, C_, SO_, cS, cC = corrected(vx, vy); rows.append((lab, views(S_, C_, SO_), (cS, cC)))
print("\n   resolution [%] sigma_rms/mean / Gaussian core, and S high tail (>2.5 sigma):")
for lab, V, c in rows:
    m0, s0 = core(V["S"]); tail = np.mean((V["S"] - m0) / s0 > 2.5) * 100
    sl = f"  (S slopes {c[0][1] * 100:+.1f}%/cm x, {c[0][2] * 100:+.1f}%/cm y; C {c[1][1] * 100:+.1f}, {c[1][2] * 100:+.1f})" if c else ""
    print(f"   {lab:30s} tail {tail:4.1f}% | " + " | ".join(f"{k} {v.std() / v.mean() * 100:5.1f}/{(lambda mm, ss: ss / mm * 100)(*core(v)):5.1f}" for k, v in V.items()) + sl)
