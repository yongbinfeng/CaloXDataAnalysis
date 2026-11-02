"""Per-event dump for the electronics part of the timing study: for every DRS group the reference (Channel8) time as used by the
framework (interpolated 50% crossing), the integer sample (the former definition) and a line fit to the 20-80% edge, the
StartIndexCell; for every MCP copy its raw and reference-corrected CFD time and amplitude; for every physics channel the raw CFD time
and amplitude; the beam-PID flags. No event selection.
Usage (repository root): python3 studies/timing/dump_refs.py RUN [NEVENTS]      Output: $TIMING_DATA/refs/refs_<RUN>.npz"""
import os, sys
run = int(sys.argv[1]); nev = sys.argv[2] if len(sys.argv) > 2 else None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import DATA, REF_EVENTS
sys.argv = ["x", "--run", str(run), "--last-event", str(nev or REF_EVENTS)]
sys.path.insert(0, ".")
from utils.root_setup import setup_root
setup_root(n_threads=1, batch_mode=True, load_functions=True)
import ROOT
import numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
from channels.maps.services import get_mcp_channels

ROOT.gInterpreter.Declare(r"""
// alternative reference times: [0] integer first sample above 50% (former framework definition), [1] line fit to the samples between
// 20% and 80% of the amplitude, crossing of 50%; -9999 below 500 ADC
ROOT::RVec<float> tim_ref_alt(const ROOT::RVec<float>& w) {
    ROOT::RVec<float> o(2, -9999.f);
    if (w.empty()) return o;
    int ip = std::max_element(w.begin(), w.end()) - w.begin(); float A = w[ip];
    if (A < 500.f) return o;
    int i = ip; while (i > 0 && w[i] > 0.5f * A) i--;
    o[0] = float(i + 1);
    double sx = 0, sy = 0, sxx = 0, sxy = 0; int n = 0;
    for (int k = i - 10; k <= ip; ++k) {
        if (k < 0 || w[k] <= 0.2f * A || w[k] >= 0.8f * A) continue;
        sx += k; sy += w[k]; sxx += double(k) * k; sxy += double(k) * w[k]; n++;
    }
    if (n >= 2) { double b = (n * sxy - sx * sy) / (n * sxx - sx * sx), a = (sy - b * sx) / n; o[1] = float((0.5 * A - a) / b); }
    return o;
}""")
m = CaloXAnalysisManager(get_args()).prepare()
rdf = m.get_rdf(); cols = set(map(str, rdf.GetColumnNames()))
groups = sorted(c[:-len("_Channel8_blsub")] for c in cols if c.endswith("_Channel8_blsub"))
mcps = {k: v for k, v in get_mcp_channels(run).items() if f"{k}_TS_cfd_ref" in cols}
phys = sorted(c[:-6] for c in cols if c.startswith("DRS_") and c.endswith("_blsub") and not c.endswith("Channel8_blsub") and f"{c[:-6]}_TS_cfd" in cols
              and c[:-6] not in mcps.values())
vec = lambda names: "ROOT::RVec<float>{" + ", ".join(f"(float)({n})" for n in names) + "}"
rdf = rdf.Define("ref_fw", vec([f"{g}_ref_TS" for g in groups]))
for g in groups: rdf = rdf.Define(f"{g}_alt", f"tim_ref_alt({g}_Channel8_blsub)")
rdf = (rdf.Define("ref_int", vec([f"{g}_alt[0]" for g in groups])).Define("ref_fit", vec([f"{g}_alt[1]" for g in groups]))
          .Define("start_cell", vec([f"{g}_StartIndexCell" for g in groups]))
          .Define("mcp_T", vec([f"{k}_TS_cfd" for k in mcps])).Define("mcp_Tref", vec([f"{k}_TS_cfd_ref" for k in mcps]))
          .Define("mcp_A", vec([f"{k}_peak_value" for k in mcps]))
          .Define("phys_T", vec([f"{p}_TS_cfd" for p in phys])).Define("phys_A", vec([f"{p}_peak_value" for p in phys])))
flags = sorted(c for c in cols if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed")))
arr = ["ref_fw", "ref_int", "ref_fit", "start_cell", "mcp_T", "mcp_Tref", "mcp_A", "phys_T", "phys_A"]
d = rdf.AsNumpy(arr + flags)
os.makedirs(f"{DATA}/refs", exist_ok=True)
np.savez(f"{DATA}/refs/refs_{run}.npz", groups=np.array(groups), mcps=np.array([f"{k}={v}" for k, v in mcps.items()]), phys=np.array(phys),
         **{k: np.stack([np.asarray(x, dtype=np.float32) for x in d[k]]) for k in arr}, **{k: np.asarray(d[k]) for k in flags})
print("done", run, len(groups), "groups", len(phys), "channels", list(mcps))
