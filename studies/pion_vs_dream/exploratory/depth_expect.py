import numpy as np
src = open("edepth_test.py").read(); exec(src[:src.index("d11 = json.load")])
T = np.tan(np.radians(3.3))
print("mean displacement of pions vs 100 GeV positrons (cm), implied extra depth dz = dX/tan(3.3 deg), and S boosts")
print(f"{'E':>4s} {'<dX>':>6s} {'<dY>':>6s} {'dz cm':>6s} | {'fit slope dlnS/dX':>17s} {'applied S':>9s} {'applied C':>9s} | {'exp l=1.8m':>10s} {'l=2.5m':>7s}")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    p = load(run, "pion"); E = p["E"]
    px = cx[0]*p["hx"]+cx[1]*p["hy"]+cx[2]; py = cy[0]*p["hx"]+cy[1]*p["hy"]+cy[2]
    dX, dY = p["ex"]-px, p["ey"]-py; ix, iy = px-np.median(px), py-np.median(py)
    A = np.column_stack([np.ones(len(dX)), dX, dY, ix, iy, ix*ix, iy*iy, ix*iy])
    def fit(sig):
        yy = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dX)<3)&(np.abs(dY)<3)&(np.abs(ix)<3)&(np.abs(iy)<3)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(A[m], yy[m], rcond=None)
            r = yy-A@c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m&(np.abs(r)<2.5*s)
        return c
    cS, cC = fit(p["S"]), fit(p["C"])
    # applied depth-only factor (displacement terms only)
    aS = p["S"].mean()/(p["S"]*np.exp(-(A[:,1:3]@cS[1:3]))).mean(); aC = p["C"].mean()/(p["C"]*np.exp(-(A[:,1:3]@cC[1:3]))).mean()
    mx, my = np.median(dX), np.median(dY); dz = abs(mx)/T
    print(f"{E:4.0f} {mx:6.2f} {my:6.2f} {dz:6.1f} | {cS[1]:+8.3f} {cS[2]:+8.3f} {aS:9.3f} {aC:9.3f} | {np.exp(dz/180):10.3f} {np.exp(dz/250):7.3f}")
# positrons at other energies: displacement shift with energy, as a check of the depth scale (shower max moves ~ X0 ln E)
import os
for run in (1423, 1416, 1514, 1526, 1409, 1527, 1411):
    if os.path.exists(f"hodo_{run}.npz"):
        e = load(run, "electron"); px = cx[0]*e["hx"]+cx[1]*e["hy"]+cx[2]; py = cy[0]*e["hx"]+cy[1]*e["hy"]+cy[2]
        print(f"e+ run {run} {e['E']:4.0f} GeV: median dX {np.median(e['ex']-px):+.3f}  dY {np.median(e['ey']-py):+.3f} cm")
