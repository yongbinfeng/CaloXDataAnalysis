"""Per-channel attenuation length from the 2026 90-degree DRS position scan (40 GeV e+, runs 1897-1931)."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json, os, sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
build_drs_boards = P.build_drs_boards
rl = json.load(open("data/Runlist.json"))
runs = [r for r in P.SCAN_RUNS if os.path.exists(f"results/root/Run{r}/drs_energy_stats.json")]
X = {r: float(rl[str(r)]["Table position X"].replace("mm", "")) for r in runs}
st = {r: json.load(open(f"results/root/Run{r}/drs_energy_stats.json")) for r in runs}
order = P._build_channel_order(build_drs_boards(run_number=P.REFERENCE_RUN))
res = {}
for ch, b, g, c, var, tx, ty, is6 in order:
    pts = [(X[r], *st[r][ch]) for r in runs if ch in st[r]]
    if len(pts) < 15: continue
    x, m, n, sat = (np.array([p[i] for p in pts]) for i in range(4))
    if sat.max() > 0.01 or np.median(m) < 100: continue
    p1, p0 = np.polyfit(x, np.log(m), 1)
    resid = np.log(m) - (p0 + p1 * x)
    res[ch] = dict(var=var, tower=(tx, ty), lam_m=(1 / p1) / 1000 if p1 != 0 else np.inf, slope_per_m=p1 * 1000, med=float(np.median(m)),
                   rms_resid=float(resid.std()), is6mm=bool(is6))
print(f"{len(runs)} runs, X {min(X.values()):.0f}..{max(X.values()):.0f} mm; channels kept {len(res)} of {len(order)} (sat<=1%, median mean >=100 ADC)")
for var in P._TYPE_ORDER:
    s = np.array([v["slope_per_m"] for v in res.values() if v["var"] == var]); w = np.array([v["med"] for v in res.values() if v["var"] == var])
    if not len(s): continue
    sw = np.average(s, weights=w)
    q = np.percentile(s, [16, 50, 84])
    print(f"  {var:10s}: n {len(s):3d} | d lnQ/dX per m: median {q[1]:+.3f} [16-84%: {q[0]:+.3f} .. {q[2]:+.3f}], signal-weighted {sw:+.3f}"
          f"  ->  lambda median {abs(1/q[1]):.2f} m, weighted {abs(1/sw):.2f} m  | fit residual rms median {np.median([v['rms_resid'] for v in res.values() if v['var']==var]):.3f}")
json.dump(res, open(f"{DATA_DIR}/scan90_att.json", "w"), indent=1, default=float)
