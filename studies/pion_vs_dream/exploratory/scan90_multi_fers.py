"""tb2025 90/93.3-degree positrons, FERS: fit ln(signal) vs table X over all runs at the same Y (not pairs)."""
import json, os
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
rl = json.load(open("/data1/yofeng/CaloX/CaloXDataAnalysis/data/Runlist.json"))
ELECTRON = ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
num = lambda s: float(str(s).replace("mm", "").strip())
# corrected runlist entries 15120 (run 1512 at Y = 936) and 15130 (run 1513 at X = -400.3)
GROUPS = {"Y~1064": [1500, 1501, 1507, 1511, 1513], "Y=1128": [1499, 1506, 1510], "Y=1000": [1502, 1508], "Y=936": [1503, 1504, 1509, 1512]}
X1513 = {"corrected -400.3": -400.3, "original -54.5": -54.5}
TYPES = (("Sci", lambda m: m[:, 2] == 0), ("C quartz", lambda m: (m[:, 2] == 1) & (m[:, 3] == 1)), ("C plastic", lambda m: (m[:, 2] == 1) & (m[:, 3] == 0)))
data = {}
for runs in GROUPS.values():
    for r in runs:
        p = f"{SCR}/scan90ch_{r}.npz"
        if not os.path.exists(p): continue
        d = np.load(p); sel = np.ones(len(d["vals"]), bool)
        for f in ELECTRON: sel &= d[f].astype(bool)
        v = d["vals"][sel]; S = v[:, d["meta"][:, 2] == 0].sum(1)
        core = S > 0.5 * np.median(S)                      # drop the low tail (non-showering / leaking events)
        data[r] = (d["meta"], v[core].mean(0), int(core.sum()), float(np.median(S)))
def lam(s): return abs(1 / s) if s else np.inf
for gname, runs in GROUPS.items():
    rs = [r for r in runs if r in data]
    for xlab, x13 in (X1513.items() if 1513 in rs else [("", None)]):
        xs = np.array([x13 if r == 1513 else num(rl[str(r)]["Table position X"]) for r in rs]) / 1000.0
        if len(set(xs)) < 2: continue
        meta = data[rs[0]][0]
        print(f"== {gname} {('(1513 at ' + xlab + ')') if xlab else ''}: runs " +
              ", ".join(f"{r}@{1000*x:.0f}mm/{rl[str(r)]['Table Degree']}deg" for r, x in zip(rs, xs)))
        for tl, tm in TYPES:
            m = tm(meta); means = np.array([data[r][1] for r in rs])
            ok = m & (means.min(0) > 0.3)
            if ok.sum() == 0: print(f'   {tl:9s} no channels above 0.3 GeV in all runs'); continue
            y = np.log(means[:, ok].sum(1)); p = np.polyfit(xs, y, 1); res = y - np.polyval(p, xs)
            per = [np.polyfit(xs, np.log(means[:, c]), 1)[0] for c in np.where(ok)[0]]
            w = means[:, ok].mean(0)
            print(f"   {tl:9s} n {ok.sum():3d} | summed: slope {p[0]:+.3f}/m (lambda {lam(p[0]):5.2f} m), residuals % " +
                  " ".join(f"{100*e:+.1f}" for e in res) +
                  f" | per-channel median {np.median(per):+.3f}/m, weighted {np.average(per, weights=w):+.3f}/m")
