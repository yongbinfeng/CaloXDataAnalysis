import numpy as np
SEL = {"pion": ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"],
       "electron": ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"]}
AMIN, AMAX = 100.0, 2000.0
def core(x):
    m, s = np.median(x), 0.5 * np.subtract(*np.percentile(x, [84, 16]))
    for _ in range(50):
        w = np.abs(x - m) < 2 * s
        m2, s2 = x[w].mean(), x[w].std() / 0.88
        if abs(m2 - m) < 1e-6: break
        m, s = m2, s2
    return m, s
for run, part in ((1432, "pion"), (1527, "electron")):
    d = np.load(f"time_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["FERS_energyMix_Sci_calib_sum"]), bool)
    for f in SEL[part]: sel &= d[f].astype(bool)
    S = d["FERS_energyMix_Sci_calib_sum"][sel] / E; C = d["FERS_energyMix_Cer_calib_sum"][sel] / E
    m, s = core(S); z = (S - m) / s
    groups = {}
    for key, metak in (("Sci", "sci_meta"), ("Cer", "cer_meta")):
        meta = d[metak]; T = d[f"{key}_t"][sel]; A = d[f"{key}_a"][sel]
        for typ in sorted(set(meta[:, 2])):
            for six in sorted(set(meta[:, 3])):
                idx = np.nonzero((meta[:, 2] == typ) & (meta[:, 3] == six))[0]
                if len(idx) == 0: continue
                t, a = T[:, idx], A[:, idx]
                good = (a > AMIN) & (a < AMAX)
                ref = np.array([np.median(t[(np.abs(z) < 1) & good[:, j], j]) if ((np.abs(z) < 1) & good[:, j]).sum() > 30 else np.nan for j in range(len(idx))])
                ok = np.isfinite(ref)
                dt = t[:, ok] - ref[ok]; w = a[:, ok] * good[:, ok] * (np.abs(dt) < 15)
                tw = (w * dt).sum(1) / np.maximum(w.sum(1), 1e-6); has = w.sum(1) > 0
                groups[f"{typ} {'6mm' if six == '1' else '3mm'} ({ok.sum()} ch)"] = (tw, has)
    print(f"== run {run} ({part}, {E:.0f} GeV): N={sel.sum()}, high tail z>2.5: {np.mean(z > 2.5) * 100:.1f}%")
    print(f"   {'channel type':28s} {'core |z|<1':>14s} {'high tail z>2.5':>16s} {'low tail z<-2.5':>16s}   [amplitude-weighted time minus channel median, TS (0.2 ns)]")
    for name, (tw, has) in groups.items():
        cells = []
        for msk in (np.abs(z) < 1, z > 2.5, z < -2.5):
            v = tw[msk & has]
            cells.append(f"{np.median(v):+6.2f} ({len(v):4d})" if len(v) >= 10 else "      -       ")
        print(f"   {name:28s} " + " ".join(f"{c:>16s}" for c in cells))
