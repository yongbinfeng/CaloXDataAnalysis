"""2025 90-degree positron runs: DRS per-channel attenuation from run pairs at the same Y and angle (dX = 50 mm),
for comparison with the FERS result on the same pairs and with the 2026 DRS scan."""
import json, os, sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
rl = json.load(open("data/Runlist.json"))
xof = lambda r: float(rl[str(r)]["Table position X"].replace("mm", ""))
PAIRS = [(1499, 1506), (1501, 1507), (1508, 1512)]
FERS = {(1499, 1506): (1.84, 1.93, 0.7), (1501, 1507): (2.5, 5.4, None), (1508, 1512): (1.7, 2.5, None)}   # lambda_S, quartz, plastic (m)
boards = P.build_drs_boards(run_number=1499); order = P._build_channel_order(boards)
for a, b in PAIRS:
    fa, fb = (f"results/root/Run{r}/drs_energy_stats.json" for r in (a, b))
    if not (os.path.exists(fa) and os.path.exists(fb)):
        print(f"pair {a}/{b}: stats missing"); continue
    sa, sb = json.load(open(fa)), json.load(open(fb)); dx = (xof(a) - xof(b)) / 1000.0   # m; a is at larger X
    mpv = P._collect_mpv_data(boards, {a: xof(a), b: xof(b)})
    print(f"pair {a} (X={xof(a):.0f}) vs {b} (X={xof(b):.0f}), Y={rl[str(a)]['Table position Y']}, angle {rl[str(a)]['Table Degree']}:"
          f" N(e+) {sa[order[0][0]][1]:.0f} / {sb[order[0][0]][1]:.0f}")
    for var in P._TYPE_ORDER:
        sl, w, dt = [], [], []
        for ch, _b, _g, _c, v, _tx, _ty, _is6 in order:
            if v != var or ch not in sa or ch not in sb: continue
            ma, mb = sa[ch][0], sb[ch][0]
            if max(sa[ch][2], sb[ch][2]) > 0.01 or min(ma, mb) < 100: continue
            sl.append(np.log(ma / mb) / dx); w.append(np.sqrt(ma * mb))
            pts = dict((x, m) for x, m, _e in mpv.get(ch, []))
            if len(pts) == 2: dt.append((pts[xof(a)] - pts[xof(b)]) * 0.2 / dx)
        if not sl: print(f"   {var:10s}: no usable channels"); continue
        sl = np.array(sl); sw = np.average(sl, weights=w); q = np.percentile(sl, [16, 50, 84])
        f = FERS[(a, b)][{"Sci": 0, "CerQuartz": 1, "CerPlastic": 2}[var]]
        print(f"   {var:10s}: n {len(sl):3d} | dlnQ/dX median {q[1]:+.3f}/m [16-84%: {q[0]:+.3f}..{q[2]:+.3f}], weighted {sw:+.3f}/m -> lambda {abs(1/sw):5.2f} m"
              f" | FERS same pair: {f if f else 'n/a'} m | time slope median {np.median(dt) if dt else np.nan:+.2f} ns/m")
