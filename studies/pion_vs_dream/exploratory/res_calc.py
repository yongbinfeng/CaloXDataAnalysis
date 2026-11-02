import json, numpy as np
RUNS = [1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434]
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-5 and abs(s2 - s) < 1e-5: break
        m, s = m2, s2
    return s / m
def fem_from(qs, ehS=1.3, ehC=5.0):
    a, b = 1 / ehC, 1 / ehS
    return (a - qs * b) / (qs * (1 - b) - (1 - a))
out = {}
for run in RUNS:
    d = np.load(f"resp_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION: sel &= d[f].astype(bool)
    S, C = d["FERS_energyMix_Sci_calib_sum"][sel], d["FERS_energyMix_Cer_calib_sum"][sel]
    SO = d["FERS_energyMix_Sci_calib_sum_OuterRing"][sel]
    Sc = S + SO                                   # = S (1 + f_out), DREAM Eq. 6 with x = 1
    f = np.clip(fem_from(C / Sc), -0.5, 1.5)
    Sq = Sc / (f + (1 - f) / 1.3)                 # DREAM Eq. 7 / 9
    shw = S > 0.3 * E                             # drop non-showering / very late pions
    r = {"E": E, "N": int(sel.sum()), "frac_cut": float(1 - shw.mean())}
    for k, x in (("S", S), ("C", C), ("QS", Sq)):
        r[f"{k}_rms_all"] = float(x.std() / x.mean() * 100)
        r[f"{k}_rms"] = float(x[shw].std() / x[shw].mean() * 100)
        r[f"{k}_core"] = float(core(x[shw]) * 100)
    out[run] = r
    print(f"{run} {E:5.0f} GeV N={r['N']:6d} cut {r['frac_cut']*100:4.1f}% | rms/mean all: S {r['S_rms_all']:5.1f} C {r['C_rms_all']:5.1f} QS {r['QS_rms_all']:5.1f}"
          f" | after cut: S {r['S_rms']:5.1f} C {r['C_rms']:5.1f} QS {r['QS_rms']:5.1f} | core: S {r['S_core']:5.1f} C {r['C_core']:5.1f} QS {r['QS_core']:5.1f}")
json.dump(out, open("calox_res.json", "w"), indent=1)
