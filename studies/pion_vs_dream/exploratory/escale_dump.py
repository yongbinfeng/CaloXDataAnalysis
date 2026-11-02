"""Per-tower S and C (Mix, calibrated) and gain-switch flags, for the electron-scale check."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import sys
run, particle, last = int(sys.argv[1]), sys.argv[2], sys.argv[3]
sys.argv = ["x", "--run", str(run), "--last-event", last]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import ROOT, numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import get_fers_saturation_value
args = get_args()
m = CaloXAnalysisManager(args).prepare().calibrate_fers()
m.apply_beam_pid_selection(flag_only=True, particle=particle)
rdf = m.get_rdf()
sat = get_fers_saturation_value()
towers = {}
for b in m.fersboards.values():
    for ch in b:
        key = (round(ch.i_tower_x, 3), round(ch.i_tower_y, 3))
        towers.setdefault(key, {})["C" if ch.isCer else "S"] = ch
cols, meta = [], []
for i, (key, d) in enumerate(sorted(towers.items())):
    if "S" not in d or "C" not in d:
        continue
    for t in ("S", "C"):
        ch = d[t]
        cal = ch.get_channel_name(gain="Mix", pdsub=True, calib=True)
        hg = ch.get_channel_name(gain="HG", pdsub=True)
        rdf = rdf.Define(f"t{i}_{t}", f"(float){cal}").Define(f"t{i}_{t}_lg", f"{hg} >= {sat}")
        cols += [f"t{i}_{t}", f"t{i}_{t}_lg"]
    meta.append((i, key[0], key[1]))
pid = [c for c in map(str, rdf.GetColumnNames()) if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy(cols + pid)
np.savez(f"{DATA_DIR}/escale_{run}.npz",
         meta=np.array(meta), pid=np.array(pid), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]), "events", len(meta), "towers")
