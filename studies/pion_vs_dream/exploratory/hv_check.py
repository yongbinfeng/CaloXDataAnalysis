import numpy as np
BASE = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6: break
        m, s = m2, s2
    return m, s
def fit(sig, vx, vy):
    y = np.log(np.maximum(sig, 1e-3)); m = (np.abs(vx) < 3) & (np.abs(vy) < 3)
    A = np.column_stack([np.ones(len(y)), vx, vy])
    for _ in range(5):
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
        r = y - A @ c; s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = m & (np.abs(r) < 2.5 * s)
    return c
print(f"{'E':>4s} {'sample':10s} {'N':>6s} | {'S raw':>12s} {'S corr':>12s} | {'Q/S raw':>12s} {'Q/S corr':>12s}   [rms/mean / core, %]")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    d = np.load(f"tail_{run}.npz"); E = float(d["energy"])
    for hv in (True, False):
        sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
        for f in BASE: sel &= d[f].astype(bool)
        if hv: sel &= d["is_HoleVeto_vetoed"].astype(bool)
        S = d["FERS_energyMix_Sci_calib_sum"][sel]; C = d["FERS_energyMix_Cer_calib_sum"][sel]; SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
        X = d["FERS_energyMix_Sci_calib_weighted_center_X"][sel]; Y = d["FERS_energyMix_Sci_calib_weighted_center_Y"][sel]
        k = S > 0.3 * E; S, C, SO, X, Y = S[k], C[k], SO[k], X[k], Y[k]
        if len(S) < 100: print(f"{E:4.0f} {'hole veto' if hv else 'no veto':10s} {len(S):6d} | too few"); continue
        vx, vy = X - np.median(X), Y - np.median(Y)
        cS, cC = fit(S, vx, vy), fit(C, vx, vy); wS, wC = np.exp(-(cS[1] * vx + cS[2] * vy)), np.exp(-(cC[1] * vx + cC[2] * vy))
        def qs(S_, C_, SO_):
            Sc = S_ + SO_; b, a = 1 / 1.3, 0.2; f = np.clip((a - C_ / Sc * b) / (C_ / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
            return Sc / (f + (1 - f) * b)
        cells = []
        for x in (S, S * wS, qs(S, C, SO), qs(S * wS, C * wC, SO * wS)):
            m_, s_ = core(x); cells.append(f"{x.std() / x.mean() * 100:5.1f}/{s_ / m_ * 100:5.1f}")
        print(f"{E:4.0f} {'hole veto' if hv else 'no veto':10s} {len(S):6d} | {cells[0]:>12s} {cells[1]:>12s} | {cells[2]:>12s} {cells[3]:>12s}")
