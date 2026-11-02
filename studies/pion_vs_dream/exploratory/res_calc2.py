"""Resolutions (pion selection + hole veto, S > 0.3 E) for S, C, raw-S DR and Q/S; also saves per-event values
(pion selection + hole veto, no S cut) for the distribution plots."""
import json, numpy as np
RUNS = [1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434]
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"]
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-5 and abs(s2 - s) < 1e-5: break
        m, s = m2, s2
    return s / m
out, events = {}, {}
for run in RUNS:
    d = np.load(f"resp_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in PION: sel &= d[f].astype(bool)
    S, C, SO = (d[k][sel] for k in ("FERS_energyMix_Sci_calib_sum", "FERS_energyMix_Cer_calib_sum", "FERS_energyMix_Sci_calib_sum_OuterRing"))
    Sc = S + SO
    a, b = 1 / 5.0, 1 / 1.3
    f = np.clip((a - C / Sc * b) / (C / Sc * (1 - b) - (1 - a)), -0.5, 1.5)
    vals = {"S": S, "C": C, "DR": (S - 0.29 * C) / 0.71, "QS": Sc / (f + (1 - f) * b)}
    events[run] = {k: v / E for k, v in vals.items()}
    shw = S > 0.3 * E
    r = {"E": E, "N": int(sel.sum())}
    for k, x in vals.items():
        r[f"{k}_rms"] = float(x[shw].std() / x[shw].mean() * 100); r[f"{k}_core"] = float(core(x[shw]) * 100)
        r[f"{k}_mean"] = float(x[shw].mean() / E)
    out[run] = r
    print(f"{run} {E:5.0f} GeV N={r['N']:5d} | rms/mean: S {r['S_rms']:5.1f} C {r['C_rms']:5.1f} DR {r['DR_rms']:5.1f} QS {r['QS_rms']:5.1f}"
          f" | core: S {r['S_core']:5.1f} C {r['C_core']:5.1f} DR {r['DR_core']:5.1f} QS {r['QS_core']:5.1f}")
json.dump(out, open("calox_res_holeveto.json", "w"), indent=1)
np.savez("calox_events_holeveto.npz", **{f"{run}_{k}": v for run, d in events.items() for k, v in d.items()})
