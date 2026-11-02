"""90-degree positron scan: per-event FERS S, C(quartz), C(plastic) sums and selection flags."""
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
m = CaloXAnalysisManager(args).prepare(do_drs=False).calibrate_fers().define_physics_variables(gain="Mix", calib=True, pdsub=True)
# no particle selection: it needs the service DRS; showering positrons are selected offline from the S peak
rdf = m.get_rdf()
S, cq, cp = [], [], []
for b in m.fersboards.values():
    for ch in b:
        n = ch.get_channel_name(gain="Mix", pdsub=True, calib=True)
        (S if not ch.isCer else (cq if ch.isQuartz else cp)).append(n)
rdf = (rdf.Define("FERS_S", "(" + " + ".join(S) + ")").Define("FERS_Cq", "(" + " + ".join(cq) + ")")
          .Define("FERS_Cp", "(" + " + ".join(cp) + ")"))
cols = ["FERS_S", "FERS_Cq", "FERS_Cp"] + [c for c in map(str, rdf.GetColumnNames())
                                           if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy(cols)
np.savez(f"{DATA_DIR}/scan90_{run}.npz",
         energy=np.float32(m.beam_energy), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]), [c for c in cols if c.startswith("is_")])
