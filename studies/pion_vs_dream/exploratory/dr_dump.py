import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import sys
sys.argv = ["x", "--run", "1465"]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import ROOT, numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
from configs.selection_config import get_particle_selection
args = get_args()
m = (CaloXAnalysisManager(args).prepare().calibrate_fers()
     .define_physics_variables(gain="Mix", calib=True, pdsub=True)
     .apply_hole_veto(flag_only=True))
m.apply_beam_pid_selection(flag_only=True, particle="pion")
m.apply_beam_pid_selection(flag_only=True, particle="muon")
sel = m._get_or_create_sel_mgr()
for det in get_particle_selection("muon", 1465):
    pass
rdf = m.get_rdf()
cols = [str(c) for c in rdf.GetColumnNames()]
print([c for c in cols if "energyMix" in c and "Board" not in c][:40])
want = [c for c in cols if (c.startswith("FERS_energyMix_") and "calib_sum" in c and "Board" not in c)
        or c.startswith("is") and ("Fired" in c or "_vetoed" in c or "_fired" in c)
        or ("_firingVal" in c) or ("EWC" in c.upper() and "Mix" in c) or c in ("COverS", "fEM")]
want = sorted(set(want))
print({w: rdf.GetColumnType(w) for w in want if "energy" in w or "EWC" in w.upper()})
print("columns:", want)
d = rdf.AsNumpy(want)
np.savez(f"{DATA_DIR}/dr_1465.npz", **{k: np.asarray(v) for k, v in d.items()})
print("events", len(d[want[0]]))
