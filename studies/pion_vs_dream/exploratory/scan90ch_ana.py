"""Per-channel attenuation from the tb2025 90-degree positron runs: channel mean vs table X, same Y and angle."""
import json
import numpy as np

import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
rl = json.load(open("/data1/yofeng/CaloX/CaloXDataAnalysis/data/Runlist.json"))
ELECTRON = ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
PAIRS = [("90 deg, Y=1128", 1499, 1506), ("90 deg, Y~1064", 1501, 1507), ("93.3 deg, Y=1000", 1508, 1512),
         ("93.3 deg, Y=1064 (X of 1513 uncertain)", 1511, 1513)]
MIN_MEAN = 0.3          # GeV: channels with a real share of the shower in both runs


def num(s):
    return float(str(s).replace("mm", "").strip())


def load(run):
    d = np.load(f"{SCR}/scan90ch_{run}.npz")
    sel = np.ones(len(d["vals"]), bool)
    for f in ELECTRON:
        sel &= d[f].astype(bool)
    v = d["vals"][sel]
    return d["meta"], v.mean(0), v.std(0) / np.sqrt(len(v)), int(sel.sum()), len(d["vals"])


summary = {}
for name, ra, rb in PAIRS:
    try:
        meta, ma, ea, na, ta = load(ra)
        _, mb, eb, nb, tb = load(rb)
    except FileNotFoundError as e:
        print(f"{name}: missing {e}")
        continue
    xa, xb = num(rl[str(ra)]["Table position X"]), num(rl[str(rb)]["Table position X"])
    dx = xb - xa
    print(f"== {name}: run {ra} (X={xa:g}, {na}/{ta} positrons) vs run {rb} (X={xb:g}, {nb}/{tb}); dX = {dx:+g} mm")
    for label, mask in (("Sci", meta[:, 2] == 0), ("C quartz", (meta[:, 2] == 1) & (meta[:, 3] == 1)),
                        ("C plastic", (meta[:, 2] == 1) & (meta[:, 3] == 0))):
        ok = mask & (ma > MIN_MEAN) & (mb > MIN_MEAN)
        if ok.sum() < 3:
            print(f"   {label:9s}: only {ok.sum()} channels above {MIN_MEAN} GeV")
            continue
        r = np.log(mb[ok] / ma[ok])                          # per-channel log ratio
        w = (ma[ok] + mb[ok])                                # weight by signal size
        slope = np.average(r, weights=w) / dx * 100          # per 100 mm
        q = np.percentile(r / dx * 100, [16, 50, 84])
        tot = np.log(mb[ok].sum() / ma[ok].sum()) / dx * 100
        lam = lambda s: abs(100.0 / s) / 1000 if s else np.inf
        print(f"   {label:9s}: {ok.sum():3d} channels | change per 100 mm: weighted {np.expm1(slope) * 100:+5.1f}% "
              f"(|lambda| {lam(slope):5.2f} m), channel median {np.expm1(q[1]) * 100:+5.1f}% "
              f"[16-84%: {np.expm1(q[0]) * 100:+5.1f} .. {np.expm1(q[2]) * 100:+5.1f}], summed {np.expm1(tot) * 100:+5.1f}%")
        summary.setdefault(label, []).append((name, float(slope), int(ok.sum())))
json.dump(summary, open(f"{SCR}/scan90ch_summary.json", "w"), indent=1)
