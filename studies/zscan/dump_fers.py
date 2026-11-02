"""Per-event FERS dump for the z-scan analysis (datasets with calorimeter FERS): calibrated Mix value of every channel, the raw LG
(saturation check), the hodoscope, the beam-PID flags and the MCP pulse heights.
Usage (repository root, or the repository mirror for a second data file): python3 dump_fers.py RUN [LAST_EVENT] [LABEL]
Output: $ZSCAN_DATA/fers/fers_<LABEL>.npz"""
import os, sys
run = int(sys.argv[1]); last = sys.argv[2] if len(sys.argv) > 2 else "40000"
label = sys.argv[3] if len(sys.argv) > 3 else str(run)          # 15130: second file of 1513, analysed as 1513 in the mirror
OUT = os.path.join(os.environ.get("ZSCAN_DATA", "/data1/yofeng/CaloX/zscan_dumps"), "fers")
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
mix, lg, meta = [], [], []
for b in m.fersboards.values():
    for ch in b:
        mix.append(ch.get_channel_name(gain="Mix", pdsub=True, calib=True)); lg.append(ch.get_channel_name(gain="LG"))
        meta.append((b.board_no, ch.channel_no, int(ch.isCer), int(ch.isQuartz), ch.i_tower_x, ch.i_tower_y, int(ch.is6mm)))
vec = lambda ns: "ROOT::RVec<float>{" + ", ".join(f"(float){n}" for n in ns) + "}"
rdf = rdf.Define("v_mix", vec(mix)).Define("v_lg", vec(lg))
cols = set(map(str, rdf.GetColumnNames()))
keep = [c for c in ("TTU_Hodo_X", "TTU_Hodo_Y", "TTU_Hodo_nHitX", "TTU_Hodo_nHitY") if c in cols]
keep += sorted(c for c in cols if c.startswith("MCP_") and c.endswith("_peak_value"))
keep += sorted(c for c in cols if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed")))
d = rdf.AsNumpy(["v_mix", "v_lg"] + keep)
arr = {k: np.stack([np.asarray(x, dtype=np.float32) for x in d.pop(f"v_{k}")]) for k in ("mix", "lg")}
os.makedirs(OUT, exist_ok=True)
np.savez(f"{OUT}/fers_{label}.npz", meta=np.array(meta, dtype=np.float32), names=np.array(mix), energy=np.float32(m.beam_energy),
         vals=arr["mix"], lg=arr["lg"].astype(np.int16), **{k: np.asarray(v) for k, v in d.items()})
print("done", label, arr["mix"].shape)
