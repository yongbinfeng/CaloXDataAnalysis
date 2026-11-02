import numpy as np
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"]
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6: break
        m, s = m2, s2
    return m, s
for run in (1429, 1432, 1433):
    d = np.load(f"tail_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION: sel &= d[f].astype(bool)
    S = d["FERS_energyMix_Sci_calib_sum"][sel]; C = d["FERS_energyMix_Cer_calib_sum"][sel]; SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    X = d["FERS_energyMix_Sci_calib_weighted_center_X"][sel]; Y = d["FERS_energyMix_Sci_calib_weighted_center_Y"][sel]
    shw = S > 0.3 * E
    S, C, SO, X, Y = S[shw], C[shw], SO[shw], X[shw], Y[shw]
    Sc = S + SO; b = 1 / 1.3; a = 0.2
    f = np.clip((a - C / Sc * b) / (C / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    V = {"S": S / E, "C": C / E, "DR": (S - 0.29 * C) / 0.71 / E, "Q/S": Sc / (f + (1 - f) * b) / E}
    x0, y0 = np.median(X), np.median(Y)
    print(f"== run {run}, {E:.0f} GeV pions (pion sel. + hole veto + S>0.3E): N={len(S)}, median shower centre X={x0:.2f}, Y={y0:.2f} cm")
    print("   S/E vs shower-centre X (cm from median):", " ".join(f"[{lo:+.1f},{hi:+.1f}) {np.median(V['S'][(X - x0 >= lo) & (X - x0 < hi)]):.3f}/{np.sum((X - x0 >= lo) & (X - x0 < hi))}"
                                                             for lo, hi in ((-9, -1.5), (-1.5, -0.5), (-0.5, 0.5), (0.5, 1.5), (1.5, 9))))
    print("   S/E vs shower-centre Y (cm from median):", " ".join(f"[{lo:+.1f},{hi:+.1f}) {np.median(V['S'][(Y - y0 >= lo) & (Y - y0 < hi)]):.3f}/{np.sum((Y - y0 >= lo) & (Y - y0 < hi))}"
                                                             for lo, hi in ((-9, -1.5), (-1.5, -0.5), (-0.5, 0.5), (0.5, 1.5), (1.5, 9))))
    m0, s0 = core(V["S"])
    print(f"   {'cut on shower-centre distance r':34s} {'kept':>6s} {'S tail>2.5sig':>13s} | " + " | ".join(f"{k:>4s} rms/mean core" for k in V))
    for rcut in (99, 3.0, 2.0, 1.5, 1.0, 0.7):
        k = np.hypot(X - x0, Y - y0) < rcut
        tail = np.mean((V["S"][k] - m0) / s0 > 2.5) * 100
        cells = []
        for key, x in V.items():
            xx = x[k]; cm, cs = core(xx)
            cells.append(f"{xx.std() / xx.mean() * 100:5.1f} {cs / cm * 100:5.1f}")
        print(f"   r < {rcut:<4g} cm {'':23s} {k.mean() * 100:5.1f}% {tail:12.1f}% | " + " | ".join(f"{c:>17s}" for c in cells))
