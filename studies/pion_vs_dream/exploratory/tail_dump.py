"""Per-event variables to compare the high-S tail with the core (pion runs)."""
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
from configs.plot_config import get_fers_saturation_value
args = get_args()
m = (CaloXAnalysisManager(args).prepare().calibrate_fers()
     .define_physics_variables(gain="Mix", calib=True, pdsub=True)
     .apply_hole_veto(flag_only=True))
m.apply_beam_pid_selection(flag_only=True, particle="pion")
rdf = m.get_rdf()
allcols = [str(c) for c in rdf.GetColumnNames()]
if run == 1432:
    print("time/event-like columns:", [c for c in allcols if any(k in c.lower() for k in ("time", "evt", "event", "trig", "stamp"))
                                       and "FERS_Board" not in c and "DRS" not in c][:30])
sci, hg = [], []
for b in m.fersboards.values():
    for ch in b:
        if not ch.isCer:
            sci.append(ch.get_channel_name(gain="Mix", pdsub=True, calib=True)); hg.append(ch.get_channel_name(gain="HG", pdsub=True))
sat = get_fers_saturation_value()
rdf = (rdf.Define("Svec", "ROOT::RVec<double>{" + ", ".join(sci) + "}")
          .Define("Shg", "ROOT::RVec<double>{" + ", ".join(hg) + "}")
          .Define("S_max", "ROOT::VecOps::Max(Svec)")
          .Define("S_nhit", f"(int)ROOT::VecOps::Sum(Svec > 0.1)")
          .Define("S_lgfrac", f"ROOT::VecOps::Sum(Svec[Shg >= {sat}]) / (ROOT::VecOps::Sum(Svec) + 1e-9)"))
cols = ["FERS_energyMix_Sci_calib_sum", "FERS_energyMix_Cer_calib_sum", "FERS_energyMix_Sci_calib_sum_OuterRing",
        "FERS_energyMix_Sci_calib_weighted_center_X", "FERS_energyMix_Sci_calib_weighted_center_Y",
        "S_max", "S_nhit", "S_lgfrac"]
cols += [c for c in allcols if c.endswith("_firingVal") or (c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed")))]
tcols = [c for c in allcols if c in ("event_n", "EventNumber", "event_number", "TriggerTime", "trigger_time", "timestamp", "Timestamp")]
d = rdf.AsNumpy(cols + tcols)
np.savez(f"{DATA_DIR}/tail_{run}.npz",
         energy=np.float32(m.beam_energy), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]), "time columns:", tcols)
