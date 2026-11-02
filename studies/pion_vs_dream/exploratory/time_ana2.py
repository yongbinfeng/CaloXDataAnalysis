import numpy as np
exec(open("time_ana.py").read().split("for run, part in")[0])
for run, part in ((1432, "pion"), (1527, "electron")):
    d = np.load(f"time_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in SEL[part]: sel &= d[f].astype(bool)
    S = d["FERS_energyMix_Sci_calib_sum"][sel] / E
    m, s = core(S); z = (S - m) / s
    print(f"== run {run} ({part}): per-pulse time minus channel median [TS], in bins of pulse amplitude [ADC]; core |z|<1 vs high tail z>2.5")
    for key, metak in (("Sci", "sci_meta"), ("Cer", "cer_meta")):
        meta = d[metak]; T = d[f"{key}_t"][sel]; A = d[f"{key}_a"][sel]
        for typ in sorted(set(meta[:, 2])):
            idx = np.nonzero((meta[:, 2] == typ) & (meta[:, 3] == "0"))[0]
            t, a = T[:, idx], A[:, idx]
            good = (a > AMIN) & (a < AMAX)
            ref = np.array([np.median(t[(np.abs(z) < 1) & good[:, j], j]) if ((np.abs(z) < 1) & good[:, j]).sum() > 30 else np.nan for j in range(len(idx))])
            dt = t - ref
            line = []
            for lo, hi in ((100, 200), (200, 400), (400, 800), (800, 1500)):
                pb = good & (a >= lo) & (a < hi) & np.isfinite(dt) & (np.abs(dt) < 15)
                c = dt[pb & (np.abs(z) < 1)[:, None]]; tl = dt[pb & (z > 2.5)[:, None]]
                line.append(f"{lo}-{hi}: core {np.median(c):+5.2f} ({len(c)}), tail {np.median(tl):+5.2f} ({len(tl)})" if len(tl) >= 15 and len(c) >= 15 else f"{lo}-{hi}: -")
            print(f"   {typ:10s} 3mm | " + " | ".join(line))
