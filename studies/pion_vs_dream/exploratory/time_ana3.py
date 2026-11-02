import numpy as np
exec(open("time_ana.py").read().split("for run, part in")[0])
MCP_MIN = 50.0
for run, part in ((1432, "pion"), (1527, "electron")):
    d = np.load(f"time_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in SEL[part]: sel &= d[f].astype(bool)
    nall = sel.sum()
    sel &= d["MCP_DS_1_peak_value"] > MCP_MIN
    S = d["FERS_energyMix_Sci_calib_sum"][sel] / E
    m, s = core(S); z = (S - m) / s
    print(f"== run {run} ({part}): {nall} selected, {sel.sum()} with an MCP pulse > {MCP_MIN:.0f} ADC; core S/E {m:.3f}, sigma {s / m * 100:.1f}%; "
          f"high tail {np.sum(z > 2.5)} events ({np.mean(z > 2.5) * 100:.1f}%)")
    for key, metak in (("Sci", "sci_meta"), ("Cer", "cer_meta")):
        meta = d[metak]; T = d[f"{key}_t"][sel]; A = d[f"{key}_a"][sel]
        for typ in sorted(set(meta[:, 2])):
            idx = np.nonzero((meta[:, 2] == typ) & (meta[:, 3] == "0"))[0]
            t, a = T[:, idx], A[:, idx]
            good = (a > AMIN) & (a < 1400)                 # below the DRS clipping (~1450-1530 ADC)
            ref = np.array([np.median(t[(np.abs(z) < 1) & good[:, j], j]) if ((np.abs(z) < 1) & good[:, j]).sum() > 20 else np.nan
                            for j in range(len(idx))])
            dt = t - ref
            print(f"   {typ:10s} 3mm: median valid time {np.nanmedian(t[good]):.1f} TS")
            for lo, hi in ((100, 300), (300, 700), (700, 1400)):
                pb = good & (a >= lo) & (a < hi) & np.isfinite(dt) & (np.abs(dt) < 20)
                c = dt[pb & (np.abs(z) < 1)[:, None]]; tl = dt[pb & (z > 2.5)[:, None]]
                if len(c) >= 15 and len(tl) >= 10:
                    err = 1.2533 * 0.5 * np.subtract(*np.percentile(tl, [84, 16])) / np.sqrt(len(tl))
                    print(f"      amplitude {lo:4d}-{hi:<4d}: core {np.median(c):+5.2f} ({len(c):5d} pulses)   tail {np.median(tl):+5.2f} +- {err:.2f} ({len(tl):4d} pulses)"
                          f"   tail - core {np.median(tl) - np.median(c):+5.2f} TS")
                else:
                    print(f"      amplitude {lo:4d}-{hi:<4d}: too few tail pulses ({len(tl)})")
