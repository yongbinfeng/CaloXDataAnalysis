"""Why ~25% of the events have all calorimeter times shifted against the MCP reference (MCP_DS_1).
Per event: the reference MCP's CFD time / peak in its search window (framework), its largest pulse in the whole record, and the
event's calorimeter time = median TS_cfd_ref of bright scintillating channels. Shift = calo - MCP - median(calo - MCP).
Usage (repository root): python3 studies/zscan/checks/mcp_check.py RUN DATASET; output $ZSCAN_DATA/checks/mcp_check_<run>.npz"""
import json, os, sys
run = int(sys.argv[1]) if len(sys.argv) > 1 else 1910; dataset = sys.argv[2] if len(sys.argv) > 2 else "2026A"
sys.argv = ["x", "--run", str(run), "--last-event", "30000"]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=1, batch_mode=True, load_functions=True)      # single thread: event order kept
import ROOT
import numpy as np
from core.analysis_manager import CaloXAnalysisManager
from utils.parser import get_args
from channels.maps.services import get_mcp_channels, get_mcp_reference
DATA = os.environ.get("ZSCAN_DATA", "/data1/yofeng/CaloX/zscan_dumps")
J = json.load(open(f"{DATA}/zscan_{dataset}.json")); C = J["selections"][J["main_selection"]]["drs"]["channels"]
m = CaloXAnalysisManager(get_args()).prepare()
m.apply_beam_pid_selection(flag_only=True, particle="electron")
rdf = m.get_rdf()
ref = get_mcp_reference(run); mcps = get_mcp_channels(run)
bright = sorted([(np.median(p["lpe"]), c) for c, p in C.items() if p["type"] == "Sci"])[-12:]
bright = [c for _, c in bright]
for det, ch in mcps.items():
    rdf = rdf.Define(f"{det}_full", f"compute_cfd_integral({ch}_blsub, true, 10.0f, 0, 1024)")
    rdf = rdf.Define(f"{det}_full_TS", f"{det}_full.time_slice").Define(f"{det}_full_A", f"{det}_full.peak_value")
    rdf = rdf.Define(f"{det}_wf", f"ROOT::RVec<float>({ch}_blsub.begin() + 300, {ch}_blsub.begin() + 750)")
cols = [f"{c}_TS_cfd_ref" for c in bright] + [f"{c}_peak_value" for c in bright]
for det in mcps:
    cols += [f"{det}_TS_cfd", f"{det}_TS_cfd_ref", f"{det}_peak_value", f"{det}_full_TS", f"{det}_full_A", f"{det}_wf"]
cols += [c for c in ("is_Cer474_fired", "is_Cer519_fired", "is_TTUMuonVeto_vetoed", "TTU_Hodo_X", "TTU_Hodo_Y") if c in set(map(str, rdf.GetColumnNames()))]
d = rdf.AsNumpy(cols)
out = {k: (np.stack([np.asarray(x, dtype=np.float32) for x in v]) if v.dtype == object else np.asarray(v)) for k, v in d.items()}
np.savez_compressed(f"{DATA}/checks/mcp_check_{run}.npz", ref=ref, mcps=np.array(list(mcps)), bright=np.array(bright), **out)
print("saved", ref, list(mcps))
