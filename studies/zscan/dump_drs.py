"""Per-event DRS dump for the z-scan analysis (and the integral checks in checks/), per DRS channel and event:
  E, Es, A, T        default chain: constant baseline (median of samples 0-200), peak searched over the whole 1024-sample record, long
                     integral (50 samples, -5..+45 from the CFD crossing) and short integral (15 samples, -4..+11), peak amplitude, CFD time
                     relative to the MCP [TS]                                       (used by the analysis)
  El, Esl            linear baseline (line through the medians of samples 0-120 and 800-1000), whole-record search      (checks)
  Eli, Esli, Ali, Tli linear baseline, peak searched only in raw samples 380-580 (in-time search); Tli in the MCP frame  (checks)
plus the beam-PID flags, the hodoscope and the MCP pulse heights / CFD times. Dumps written with implicit multithreading have their
own event order: columns of two dumps must never be combined row by row.
Usage (repository root, or the repository mirror for a second data file): python3 dump_drs.py RUN [LAST_EVENT] [LABEL] [DATASET]
Output: $ZSCAN_DATA/drs/drs_<LABEL>.npz; for a dataset with its own channel map (datasets.py, e.g. testfibres) only those channels, in
the dataset's dump directory"""
import os, sys
run = int(sys.argv[1]); last = sys.argv[2] if len(sys.argv) > 2 else "40000"
label = sys.argv[3] if len(sys.argv) > 3 else str(run)          # e.g. 15130: second data file of 1513, dumped as 1513 from the mirror
dataset = sys.argv[4] if len(sys.argv) > 4 else None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from datasets import DATASETS, drs_dump
OUT = os.path.dirname(drs_dump(dataset, label)) if dataset else os.path.join(os.environ.get("ZSCAN_DATA", "/data1/yofeng/CaloX/zscan_dumps"), "drs")
cmap = DATASETS[dataset].get("channel_map") if dataset else None
sys.argv = ["x", "--run", str(run), "--last-event", last]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import ROOT
import numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager

SEARCH = (380, 580)
ROOT.gInterpreter.Declare(r"""
// subtract a straight line through the medians of samples [0,120) and [800,1000) (two-anchor linear baseline)
ROOT::RVec<float> att_linsub(const ROOT::RVec<float>& w) {
    auto med = [&](int lo, int hi) { std::vector<float> v(w.begin() + lo, w.begin() + std::min<int>(hi, w.size()));
                                     std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end()); return v[v.size() / 2]; };
    const float a = med(0, 120), b = med(800, 1000);
    ROOT::RVec<float> o(w.size());
    for (size_t i = 0; i < w.size(); ++i) o[i] = w[i] - (a + (b - a) * (float(i) - 60.f) / 840.f);
    return o;
}""")
args = get_args()
m = CaloXAnalysisManager(args).prepare()
m.apply_beam_pid_selection(flag_only=True, particle="electron")
rdf = m.get_rdf()
names = [ch.get_channel_name(blsub=False) for b in m.drsboards.values() for ch in b if not ch.is_reference]
cols = set(map(str, rdf.GetColumnNames()))
names = [n for n in names if f"{n}_blsub" in cols]
if cmap:
    from channels.maps.paths import load_channel_map
    keep_ch = set(load_channel_map(cmap).values()); names = [n for n in names if n in keep_ch]
for n in names:
    rdf = (rdf.Define(f"{n}_bll", f"att_linsub({n}_blsub)")
              .Define(f"{n}_CFDl", f"compute_cfd_integral({n}_bll)")
              .Define(f"{n}_CFDli", f"compute_cfd_integral({n}_bll, true, 5.0f, {SEARCH[0]}, {SEARCH[1]})")
              # MCP frame for the in-time-search time: same per-event offset as the default chain
              .Define(f"{n}_Tli", f"{n}_CFDli.is_valid && {n}_CFD.is_valid ? {n}_CFDli.time_slice + ({n}_TS_cfd_mcp - {n}_TS_cfd) : -99999.f"))
def vec(expr):
    return "ROOT::RVec<float>{" + ", ".join("(float)" + expr.format(n=n) for n in names) + "}"
VARS = {"E": "{n}_energy", "Es": "{n}_energy_short", "A": "{n}_peak_value", "T": "{n}_TS_cfd_mcp",
        "El": "{n}_CFDl.energy", "Esl": "{n}_CFDl.energy_short",
        "Eli": "{n}_CFDli.energy", "Esli": "{n}_CFDli.energy_short", "Ali": "{n}_CFDli.peak_value", "Tli": "{n}_Tli"}
for k, e in VARS.items():
    rdf = rdf.Define(f"v_{k}", vec(e))
keep = [c for c in ("TTU_Hodo_X", "TTU_Hodo_Y", "TTU_Hodo_nHitX", "TTU_Hodo_nHitY") if c in cols]
# MCP pulse heights (in their search windows): the time reference needs a real pulse (see mcp_check.py)
keep += [c for c in cols if c.startswith("MCP_") and (c.endswith("_peak_value") or c.endswith("_TS_cfd"))]
pid = [c for c in cols if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy([f"v_{k}" for k in VARS] + keep + pid)
arr = {k: np.stack([np.asarray(x, dtype=np.float32) for x in d.pop(f"v_{k}")]) for k in VARS}
os.makedirs(OUT, exist_ok=True)
np.savez(f"{OUT}/drs_{label}.npz", channels=np.array(names), energy=np.float32(m.beam_energy), **arr, **{k: np.asarray(v) for k, v in d.items()})
print("done", label, arr["E"].shape)
