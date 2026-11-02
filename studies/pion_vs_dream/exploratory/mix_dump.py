"""Raw pedestal-subtracted HG and LG per channel, to test the HG->LG conversion used by the Mix gain."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import sys, json
run, last = int(sys.argv[1]), sys.argv[2]
sys.argv = ["x", "--run", str(run), "--last-event", last]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
args = get_args()
m = CaloXAnalysisManager(args).prepare().calibrate_fers()
rdf = m.get_rdf()
cols, meta = [], []
for b in m.fersboards.values():
    for ch in b:
        hg, lg = ch.get_channel_name(gain="HG", pdsub=True), ch.get_channel_name(gain="LG", pdsub=True)
        key = f"FERS_Board{b.board_no}_{ch.channel_no}"
        i = len(meta)
        rdf = rdf.Define(f"c{i}_hg", f"(float){hg}").Define(f"c{i}_lg", f"(float){lg}")
        cols += [f"c{i}_hg", f"c{i}_lg"]; meta.append((key, int(ch.isCer)))
d = rdf.AsNumpy(cols)
np.savez(f"{DATA_DIR}/mix_{run}.npz",
         keys=np.array([k for k, _ in meta]), iscer=np.array([c for _, c in meta]), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]), "events", len(meta), "channels")
