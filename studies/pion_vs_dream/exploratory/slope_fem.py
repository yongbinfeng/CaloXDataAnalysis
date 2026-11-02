import numpy as np
src = open("edepth_test.py").read(); exec(src[:src.index("d11 = json.load")])
T = np.tan(np.radians(3.3))
print("d lnS/d dX (per cm): all events vs within fixed C/S bins (f_em fixed); pure attenuation expectation = 1/(tan3.3 * lambda)")
print(f"   lambda 1.8 m -> {1/(T*180):.3f}/cm,  2.5 m -> {1/(T*250):.3f}/cm")
for run in (1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    px = cx[0]*p["hx"]+cx[1]*p["hy"]+cx[2]; py = cy[0]*p["hx"]+cy[1]*p["hy"]+cy[2]
    dX, dY = p["ex"]-px, p["ey"]-py; ix, iy = px-np.median(px), py-np.median(py)
    ok = (np.abs(dX)<3)&(np.abs(dY)<3)&(np.abs(ix)<3)&(np.abs(iy)<3)&(p["S"]>0)&(p["C"]>0)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix*ix, iy*iy, ix*iy])
    yS, yC = np.log(np.maximum(p["S"],1e-3)), np.log(np.maximum(p["C"],1e-3))
    def slope(m, y):
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); return c[1], c[2]
    a = slope(ok, yS); ac = slope(ok, yC)
    q = p["C"]/p["S"]; edges = np.quantile(q[ok], np.linspace(0,1,9))
    bs, bc, w = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok&(q>=lo)&(q<hi)
        if m.sum() > 200:
            bs.append(slope(m, yS)); bc.append(slope(m, yC)); w.append(m.sum())
    bs, bc = np.average(np.array(bs), axis=0, weights=w), np.average(np.array(bc), axis=0, weights=w)
    print(f"{E:4.0f} GeV  S: all dX {a[0]:+.3f} dY {a[1]:+.3f} | fixed C/S dX {bs[0]:+.3f} dY {bs[1]:+.3f}    "
          f"C: all dX {ac[0]:+.3f} dY {ac[1]:+.3f} | fixed C/S dX {bc[0]:+.3f} dY {bc[1]:+.3f}")
