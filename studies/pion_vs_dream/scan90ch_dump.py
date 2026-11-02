"""90-degree positron runs (tb2025): per-channel calibrated FERS values plus the existing electron PID flags."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import sys
run, last = int(sys.argv[1]), sys.argv[2]
label = sys.argv[3] if len(sys.argv) > 3 else str(run)   # output name (e.g. 15130 for the second file of 1513)
sys.argv = ["x", "--run", str(run), "--last-event", last]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
args = get_args()
m = CaloXAnalysisManager(args).prepare().calibrate_fers()
m.apply_beam_pid_selection(flag_only=True, particle="electron")
rdf = m.get_rdf()
names, meta = [], []
for b in m.fersboards.values():
    for ch in b:
        names.append(ch.get_channel_name(gain="Mix", pdsub=True, calib=True))
        meta.append((b.board_no, ch.channel_no, int(ch.isCer), int(ch.isQuartz), ch.i_tower_x, ch.i_tower_y))
rdf = rdf.Define("chvals", "ROOT::RVec<float>{" + ", ".join(f"(float){n}" for n in names) + "}")
pid = [c for c in map(str, rdf.GetColumnNames()) if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy(["chvals"] + pid)
vals = np.stack([np.asarray(x, dtype=np.float32) for x in d.pop("chvals")])
np.savez(f"{DATA_DIR}/scan90ch_{label}.npz",
         energy=np.float32(m.beam_energy), vals=vals, meta=np.array(meta, dtype=float), **{k: np.asarray(v) for k, v in d.items()})
print("done", label, vals.shape, pid)
