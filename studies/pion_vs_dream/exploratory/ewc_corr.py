import numpy as np
exec(open("ewc_cut.py").read().split("for run in")[0])
for run in (1429, 1432, 1433):
    d = np.load(f"tail_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION: sel &= d[f].astype(bool)
    S = d["FERS_energyMix_Sci_calib_sum"][sel]; C = d["FERS_energyMix_Cer_calib_sum"][sel]; SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    X = d["FERS_energyMix_Sci_calib_weighted_center_X"][sel]; Y = d["FERS_energyMix_Sci_calib_weighted_center_Y"][sel]
    shw = S > 0.3 * E
    S, C, SO, X, Y = S[shw], C[shw], SO[shw], X[shw], Y[shw]
    dx, dy = X - np.median(X), Y - np.median(Y)
    # robust linear fit of ln(signal) vs (dx, dy): least squares on events within 3 cm, iterated with 2.5-sigma clipping
    def fit(sig):
        y = np.log(np.maximum(sig, 1e-3)); m = (np.abs(dx) < 3) & (np.abs(dy) < 3)
        A = np.column_stack([np.ones(len(y)), dx, dy])
        for _ in range(5):
            coef, *_ = np.linalg.lstsq(A[m], y[m], rcond=None)
            res = y - A @ coef; s = 1.4826 * np.median(np.abs(res[m] - np.median(res[m])))
            m = m & (np.abs(res) < 2.5 * s)
        return coef
    cS, cC, cO = fit(S), fit(C), fit(SO + 1e-3)
    Sx = S * np.exp(-(cS[1] * dx + cS[2] * dy)); Cx = C * np.exp(-(cC[1] * dx + cC[2] * dy)); SOx = SO * np.exp(-(cS[1] * dx + cS[2] * dy))
    def views(S_, C_, SO_):
        Sc = S_ + SO_; b = 1 / 1.3; a = 0.2
        f = np.clip((a - C_ / Sc * b) / (C_ / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
        return {"S": S_ / E, "C": C_ / E, "DR": (S_ - 0.29 * C_) / 0.71 / E, "Q/S": Sc / (f + (1 - f) * b) / E}
    print(f"== run {run}, {E:.0f} GeV: fitted slopes per cm of centre offset: S x {cS[1] * 100:+.1f}%, y {cS[2] * 100:+.1f}% | C x {cC[1] * 100:+.1f}%, y {cC[2] * 100:+.1f}%")
    for lab, V in (("no correction", views(S, C, SO)), ("position-corrected", views(Sx, Cx, SOx))):
        m0, s0 = core(V["S"]); tail = np.mean((V["S"] - m0) / s0 > 2.5) * 100
        cells = []
        for key, x in V.items():
            cm, cs = core(x); cells.append(f"{key} {x.std() / x.mean() * 100:5.1f} / {cs / cm * 100:5.1f}")
        print(f"   {lab:20s} all {len(S)} events, S tail {tail:4.1f}% | " + " | ".join(cells))
