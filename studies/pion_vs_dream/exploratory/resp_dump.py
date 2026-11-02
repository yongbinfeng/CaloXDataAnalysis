"""Per-event calibrated S and C sums and the pion-selection flags, for the response-vs-energy plot."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import sys
run, last = int(sys.argv[1]), sys.argv[2]
sys.argv = ["x", "--run", str(run), "--last-event", last]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
args = get_args()
m = (CaloXAnalysisManager(args).prepare().calibrate_fers()
     .define_physics_variables(gain="Mix", calib=True, pdsub=True)
     .apply_hole_veto(flag_only=True))
m.apply_beam_pid_selection(flag_only=True, particle="pion")
rdf = m.get_rdf()
cols = ["FERS_energyMix_Sci_calib_sum", "FERS_energyMix_Cer_calib_sum",
        "FERS_energyMix_Sci_calib_sum_LeakCorr", "FERS_energyMix_Cer_calib_sum_LeakCorr",
        "FERS_energyMix_Sci_calib_sum_OuterRing"]
cols += [c for c in map(str, rdf.GetColumnNames()) if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy(cols)
np.savez(f"{DATA_DIR}/resp_{run}.npz",
         energy=np.float32(m.beam_energy), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, m.beam_energy, len(d[cols[0]]))
