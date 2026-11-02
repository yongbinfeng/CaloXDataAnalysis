import numpy as np
SEL = {"electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"],
       "pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]}
def get(run, part, hodo="c"):
    d = np.load(f"hodo_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL[part]: sel &= d[f].astype(bool)
    nx, ny = d["TTU_Hodo_nHitX"], d["TTU_Hodo_nHitY"]
    print(f"run {run} ({part}): {sel.sum()} selected; hodoscope hits per axis X: " +
          " ".join(f"{k}:{np.mean(nx[sel] == k) * 100:.0f}%" for k in range(0, 5)) + f" >4:{np.mean(nx[sel] > 4) * 100:.0f}% | Y: " +
          " ".join(f"{k}:{np.mean(ny[sel] == k) * 100:.0f}%" for k in range(0, 5)) + f" >4:{np.mean(ny[sel] > 4) * 100:.0f}%")
    good = sel & (nx >= 1) & (nx <= 3) & (ny >= 1) & (ny <= 3)
    hx, hy = (d["hodoX_c"], d["hodoY_c"]) if hodo == "c" else (d["TTU_Hodo_X"], d["TTU_Hodo_Y"])
    return dict(E=E, hx=hx[good], hy=hy[good], ex=d["FERS_energyMix_Sci_calib_weighted_center_X"][good],
                ey=d["FERS_energyMix_Sci_calib_weighted_center_Y"][good], cx=d["FERS_energyMix_Cer_calib_weighted_center_X"][good],
                cy=d["FERS_energyMix_Cer_calib_weighted_center_Y"][good], S=d["FERS_energyMix_Sci_calib_sum"][good],
                C=d["FERS_energyMix_Cer_calib_sum"][good], SO=d["FERS_energyMix_Sci_calib_sum_OuterRing"][good])
def fit_map(h, ex, ey):
    """EWC = A @ hodo + b, robust (iterative 3-sigma clipping)."""
    X = np.column_stack([h["hx"], h["hy"], np.ones(len(h["hx"]))]); m = np.ones(len(X), bool)
    for _ in range(6):
        cx, *_ = np.linalg.lstsq(X[m], ex[m], rcond=None); cy, *_ = np.linalg.lstsq(X[m], ey[m], rcond=None)
        rx, ry = ex - X @ cx, ey - X @ cy
        sx, sy = 1.4826 * np.median(np.abs(rx[m])), 1.4826 * np.median(np.abs(ry[m]))
        m = (np.abs(rx) < 3 * sx) & (np.abs(ry) < 3 * sy)
    return cx, cy, sx, sy, m.mean()
e = get(1527, "electron")
print(f"positrons: {len(e['hx'])} events with 1-3 hit fibres per axis; hodo X range {np.percentile(e['hx'], [2, 98]).round(2)} cm, Y {np.percentile(e['hy'], [2, 98]).round(2)} cm")
print(f"   raw correlations: corr(EWC_x, hodo_x) {np.corrcoef(e['ex'], e['hx'])[0,1]:+.2f}, corr(EWC_x, hodo_y) {np.corrcoef(e['ex'], e['hy'])[0,1]:+.2f}, "
      f"corr(EWC_y, hodo_y) {np.corrcoef(e['ey'], e['hy'])[0,1]:+.2f}, corr(EWC_y, hodo_x) {np.corrcoef(e['ey'], e['hx'])[0,1]:+.2f}")
cx, cy, sx, sy, fin = fit_map(e, e["ex"], e["ey"])
A = np.array([[cx[0], cx[1]], [cy[0], cy[1]]])
print(f"   map EWC = A hodo + b:  A = [[{A[0,0]:+.3f}, {A[0,1]:+.3f}], [{A[1,0]:+.3f}, {A[1,1]:+.3f}]],  b = ({cx[2]:+.2f}, {cy[2]:+.2f}) cm;"
      f"  residual width X {sx * 10:.2f} mm, Y {sy * 10:.2f} mm ({fin * 100:.0f}% inside 3 sigma)")
sxv, syv = np.std(e["ex"]), np.std(e["ey"])
print(f"   spread of the positron shower centre without the hodoscope: X {1.4826 * np.median(np.abs(e['ex'] - np.median(e['ex']))) * 10:.2f} mm, "
      f"Y {1.4826 * np.median(np.abs(e['ey'] - np.median(e['ey']))) * 10:.2f} mm")
np.save("hodo_map.npy", np.array([cx, cy]))
