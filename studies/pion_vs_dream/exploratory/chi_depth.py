import numpy as np
src = open("edepth_res.py").read()
exec(src[:src.index('print("resolution sigma_rms')])
def fem_from(qs, ehS=1.3, ehC=5.0):
    a, b = 1 / ehC, 1 / ehS
    return (a - qs * b) / (qs * (1 - b) - (1 - a))
print(f"{'E':>4s} {'version':>10s} | {'<S>/E':>6s} {'<C>/E':>6s} {'chi centroid':>12s} | {'p1/p0 (Fig.27)':>14s} {'-> e/h_S':>8s} {'chi':>5s} | {'best chi (DR rms)':>17s}")
store = {}
for run in (1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    px = cx[0] * p["hx"] + cx[1] * p["hy"] + cx[2]; py = cy[0] * p["hx"] + cy[1] * p["hy"] + cy[2]
    dX, dY = p["ex"] - px, p["ey"] - py; ix, iy = px - np.median(px), py - np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix * ix, iy * iy, ix * iy])
    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX) < 3) & (np.abs(dY) < 3) & (np.abs(ix) < 3) & (np.abs(iy) < 3)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            r = yy - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    # impact-only correction (no depth) vs impact + depth normalised to positron depth
    Ai = A.copy(); Ai[:, 1:3] = 0
    versions = {"impact": (p["S"] * np.exp(-(Ai[:, 1:] @ cS[1:])), p["C"] * np.exp(-(Ai[:, 1:] @ cC[1:]))),
                "e+ depth": (p["S"] * np.exp(-(A[:, 1:] @ cS[1:])), p["C"] * np.exp(-(A[:, 1:] @ cC[1:])))}
    for lab, (S, C) in versions.items():
        s, c = S / E, C / E
        chi_c = (1 - s.mean()) / (1 - c.mean())
        f = fem_from(c / s); ok = (f > 0) & (f < 1)
        p1, p0 = np.polyfit(f[ok], s[ok], 1); eh = 1 + p1 / p0
        chis = np.arange(-0.4, 0.81, 0.02)
        res = [((s - x * c) / (1 - x)).std() / ((s - x * c) / (1 - x)).mean() for x in chis]
        print(f"{E:4.0f} {lab:>10s} | {s.mean():6.3f} {c.mean():6.3f} {chi_c:12.3f} | {p1 / p0:14.3f} {eh:8.2f} {(1 - 1 / eh) / (1 - 0.2):5.2f} | {chis[int(np.argmin(res))]:+17.2f}")
        store[(E, lab)] = (s, c)
np.savez("chi_depth_store.npz", **{f"{int(k[0])}_{k[1].replace(' ', '').replace('+', 'p')}_{w}": v[i] for k, v in store.items() for i, w in enumerate("sc")})
