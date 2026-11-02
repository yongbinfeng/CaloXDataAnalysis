"""Type-summed DRS mean signal vs table X: 2026 scan (fit over 33 cm) and 2025 run pairs (50 mm)."""
import json, os, sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
rl = json.load(open("data/Runlist.json"))
xof = lambda r: float(rl[str(r)]["Table position X"].replace("mm", ""))
def typesums(run, order):
    st = json.load(open(f"results/root/Run{run}/drs_energy_stats.json")); out = {}
    for ch, _b, _g, _c, var, _tx, _ty, _6 in order:
        if ch in st: out[var] = out.get(var, 0.0) + st[ch][0]
    return out
o26 = P._build_channel_order(P.build_drs_boards(run_number=1897))
runs = [r for r in P.SCAN_RUNS if os.path.exists(f"results/root/Run{r}/drs_energy_stats.json")]
x = np.array([xof(r) for r in runs]); T = [typesums(r, o26) for r in runs]
print("2026 scan, 34 runs over 33 cm:")
for var in P._TYPE_ORDER:
    y = np.log([t[var] for t in T]); p = np.polyfit(x, y, 1); res = y - np.polyval(p, x)
    print(f"   {var:10s}: slope {p[0]*1000:+.3f}/m -> lambda {abs(1/(p[0]*1000)):5.1f} m, residual rms {res.std()*100:.2f}%, "
          f"largest run-to-run jump between neighbours {np.max(np.abs(np.diff(res)))*100:.2f}%")
o25 = P._build_channel_order(P.build_drs_boards(run_number=1499))
print("2025 pairs (dX = 50 mm):")
for a, b in ((1499, 1506), (1501, 1507), (1508, 1512)):
    if not all(os.path.exists(f"results/root/Run{r}/drs_energy_stats.json") for r in (a, b)): print(f"   {a}/{b}: pending"); continue
    ta, tb = typesums(a, o25), typesums(b, o25); dx = (xof(a) - xof(b)) / 1000
    print(f"   {a}/{b}: " + ", ".join(f"{v} {100*(ta[v]/tb[v]-1):+.1f}% (slope {np.log(ta[v]/tb[v])/dx:+.2f}/m)" for v in P._TYPE_ORDER))
