"""DRS timing (CFD, wrt the MCP) of the calorimeter towers per event, plus FERS sums, for tail-vs-core comparisons."""
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
import numpy as np
from utils.parser import get_args
from utils.utils import get_channel_var
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import get_fers_saturation_value
from channels.maps.services import get_mcp_reference
args = get_args()
m = (CaloXAnalysisManager(args).prepare().calibrate_fers()
     .define_physics_variables(gain="Mix", calib=True, pdsub=True)
     .apply_hole_veto(flag_only=True))
m.apply_beam_pid_selection(flag_only=True, particle=particle)
rdf = m.get_rdf()
cols_all = set(map(str, rdf.GetColumnNames()))
ref = get_mcp_reference(run)
grp = {"Sci": [], "Cer": []}; meta = []
for b in m.drsboards.values():
    for ch in b:
        if ch.is_reference: continue
        n = ch.get_channel_name(blsub=False); v = get_channel_var(ch)
        key = "Sci" if v == "Sci" else ("Cer" if v.startswith("Cer") else None)
        if key and f"{n}_TS_cfd_mcp" in cols_all and f"{n}_peak_value" in cols_all:
            grp[key].append(n); meta.append((key, n, v, int(ch.is6mm)))
print("DRS channels with timing:", {k: len(v) for k, v in grp.items()}, "MCP reference:", ref)
for key, chans in grp.items():
    t = "ROOT::RVec<float>{" + ", ".join(f"(float){c}_TS_cfd_mcp" for c in chans) + "}"
    a = "ROOT::RVec<float>{" + ", ".join(f"(float){c}_peak_value" for c in chans) + "}"
    rdf = (rdf.Define(f"{key}_t", t).Define(f"{key}_a", a)
              .Define(f"{key}_imax", f"(int)ROOT::VecOps::ArgMax({key}_a)")
              .Define(f"{key}_tmax", f"{key}_t[{key}_imax]").Define(f"{key}_amax", f"{key}_a[{key}_imax]")
              .Define(f"{key}_w", f"{key}_a * ({key}_a > 50.f)")
              .Define(f"{key}_tw", f"ROOT::VecOps::Sum({key}_w * {key}_t) / (ROOT::VecOps::Sum({key}_w) + 1e-6)")
              .Define(f"{key}_nsat", f"(int)ROOT::VecOps::Sum({key}_a > 2000.f)"))
# FERS Cherenkov sums per fibre type (quartz / plastic), calibrated Mix gain
cq, cp = [], []
for b in m.fersboards.values():
    for ch in b:
        if ch.isCer:
            (cq if ch.isQuartz else cp).append(ch.get_channel_name(gain="Mix", pdsub=True, calib=True))
rdf = rdf.Define("FERS_Cq", "(" + " + ".join(cq) + ")").Define("FERS_Cp", "(" + " + ".join(cp) + ")")
cols_arr = ["Sci_t", "Sci_a", "Cer_t", "Cer_a"]
cols = ["FERS_Cq", "FERS_Cp", "FERS_energyMix_Sci_calib_sum", "FERS_energyMix_Cer_calib_sum", "FERS_energyMix_Sci_calib_sum_OuterRing",
        "Sci_tmax", "Sci_amax", "Sci_tw", "Sci_nsat", "Cer_tmax", "Cer_amax", "Cer_tw", "Cer_nsat"]
cols += [c for c in map(str, rdf.GetColumnNames()) if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
if ref and f"{ref}_peak_value" in cols_all: cols += [f"{ref}_peak_value"]
d = rdf.AsNumpy(cols + cols_arr)
arrs = {k: np.stack([np.asarray(x, dtype=np.float32) for x in d.pop(k)]) for k in cols_arr}
np.savez(f"{DATA_DIR}/time_{run}.npz",
         energy=np.float32(m.beam_energy), ref=np.array(str(ref)),
         sci_meta=np.array([x for x in meta if x[0] == "Sci"]), cer_meta=np.array([x for x in meta if x[0] == "Cer"]),
         **arrs, **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]))
