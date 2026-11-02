"""Fibre hodoscope (TTU, FERS boards 0/1) positions per event, with FERS shower centres and sums,
plus the Cherenkov sum split into quartz and plastic fibres."""
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
import ROOT, numpy as np
from utils.parser import get_args
from core.analysis_manager import CaloXAnalysisManager
ROOT.gInterpreter.Declare(r"""
// amplitude-weighted fibre index of channels above threshold (fibre index from the channel via get_x/y_index)
float hodo_centroid(const ROOT::RVec<float>& hg, float thr, bool isx) {
    double sw = 0, sx = 0;
    for (int i = 0; i < (int)hg.size() && i < 64; ++i) {
        if (hg[i] <= thr) continue;
        int idx = isx ? get_x_index(i) : get_y_index(i);
        if (idx < 0) continue;
        sw += hg[i] - thr; sx += (hg[i] - thr) * idx;
    }
    return sw > 0 ? 0.064f * sx / sw : -1.f;
}""")
args = get_args()
m = (CaloXAnalysisManager(args).prepare().calibrate_fers()
     .define_physics_variables(gain="Mix", calib=True, pdsub=True)
     .apply_hole_veto(flag_only=True))
m.apply_beam_pid_selection(flag_only=True, particle=particle)
rdf = (m.get_rdf()
       .Define("hodoX_c", "hodo_centroid(ROOT::RVec<float>(FERS_Board1_energyHG.begin(), FERS_Board1_energyHG.end()), 4000.f, true)")
       .Define("hodoY_c", "hodo_centroid(ROOT::RVec<float>(FERS_Board0_energyHG.begin(), FERS_Board0_energyHG.end()), 4000.f, false)"))
for name, quartz in (("FERS_energyMix_CerQuartz_calib_sum", True), ("FERS_energyMix_CerPlastic_calib_sum", False)):
    chans = [ch.get_channel_name(gain="Mix", pdsub=True, calib=True) for b in m.fersboards.values() for ch in b
             if ch.isCer and bool(ch.isQuartz) == quartz]
    rdf = rdf.Define(name, " + ".join(f"(double){c}" for c in chans))
cols = ["TTU_Hodo_X", "TTU_Hodo_Y", "TTU_Hodo_nHitX", "TTU_Hodo_nHitY", "hodoX_c", "hodoY_c",
        "FERS_energyMix_Sci_calib_weighted_center_X", "FERS_energyMix_Sci_calib_weighted_center_Y",
        "FERS_energyMix_Cer_calib_weighted_center_X", "FERS_energyMix_Cer_calib_weighted_center_Y",
        "FERS_energyMix_Sci_calib_sum", "FERS_energyMix_Cer_calib_sum", "FERS_energyMix_Sci_calib_sum_OuterRing",
        "FERS_energyMix_CerQuartz_calib_sum", "FERS_energyMix_CerPlastic_calib_sum"]
cols += [c for c in map(str, rdf.GetColumnNames()) if c.startswith("is_") and (c.endswith("_fired") or c.endswith("_vetoed"))]
d = rdf.AsNumpy(cols)
np.savez(f"{DATA_DIR}/hodoq_{run}.npz",
         energy=np.float32(m.beam_energy), **{k: np.asarray(v) for k, v in d.items()})
print("done", run, len(d[cols[0]]))
