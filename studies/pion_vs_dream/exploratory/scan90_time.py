import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json, sys
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
rl = json.load(open("data/Runlist.json"))
run_x = {r: float(rl[str(r)]["Table position X"].replace("mm", "")) for r in P.SCAN_RUNS}
boards = P.build_drs_boards(run_number=P.REFERENCE_RUN); order = P._build_channel_order(boards)
mpv = P._collect_mpv_data(boards, run_x)
att = json.load(open(f"{DATA_DIR}/scan90_att.json"))
for var in P._TYPE_ORDER:
    sl, sa = [], []
    for ch, b, g, c, v, tx, ty, is6 in order:
        if v != var or ch not in mpv or ch not in att: continue
        pts = np.array(mpv[ch]); ok = (pts[:, 1] > 0) & (pts[:, 1] < 2000)
        if ok.sum() < 15: continue
        p1 = np.polyfit(pts[ok, 0], pts[ok, 1], 1)[0]          # TS per mm
        sl.append(p1 * 0.2 * 1000); sa.append(att[ch]["slope_per_m"])   # ns per m
    sl = np.array(sl)
    print(f"{var:10s}: n {len(sl)}  arrival-time slope median {np.median(sl):+.2f} ns/m [16-84%: {np.percentile(sl,16):+.2f} .. {np.percentile(sl,84):+.2f}]"
          f"  (1/v for v~0.6c: 5.6 ns/m);  corr(time slope, ln-signal slope) {np.corrcoef(sl, sa)[0,1]:+.2f}")
