"""
FERS physics analysis.

Books all lazy FERS histograms (Cer vs Sci, dual readout, energy weighted
centre, EWC vs hodoscope, shower shape) using the calibrated Mix-gain
physics variables, triggers ROOT.RDF.RunGraphs once, saves ROOT files,
then generates HTML plot pages.

Usage
-----
  python scripts/make_fers.py                            # default FERS_SEQUENCES, inclusive
  python scripts/make_fers.py --particle pion            # pion-selected sample
  python scripts/make_fers.py --sequences fers_cer_vs_sci fers_dr
  python scripts/make_fers.py --particle pion --sequences fers_energy_sum fers_dr

Available sequences (see analysis/registry.py  FERS_SEQUENCES):
  fers_energy_sum    — 1D total (all-boards) energy sum distributions (all gains)
  fers_energy_sum_board — 1D per-board energy sum distributions (opt-in; not run by default)
  fers_lg_vs_mix     — LG vs Mix 2D per-channel correlations (opt-in; not run by default)
  fers_stats         — per-channel mean / max / saturation maps
  fers_cer_vs_sci    — Cer vs Sci 2D correlations (all gains)
  fers_dr            — dual-readout histograms (Mix, calibrated)
  fers_ewc           — energy-weighted centre distributions
  fers_ewc_vs_hodo   — EWC vs TTU hodoscope correlations
  fers_shower_shape  — per-channel shower-shape profiles
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase, run_plot_phase
from analysis.registry import FERS_SEQUENCES, default_fers_sequences, select_sequences
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame: FERS calibration + physics variables for all gains
# ------------------------------------------------------------------
manager = (
    CaloXAnalysisManager(args)
    .prepare()
    .calibrate_fers()
    .apply_hole_veto(flag_only=True)
    #.apply_beam_pid_selection(particle=args.particle or None)
    .define_physics_variables(gain="HG",  calib=False, pdsub=True)
    .define_physics_variables(gain="LG",  calib=False, pdsub=True)
    .define_physics_variables(gain="Mix", calib=True,  pdsub=True)
)

particle = manager.get_particle_type(args.particle)
if particle:
    manager.set_output_subdir(particle)

if args.sequences:
    sequences = select_sequences(args.sequences, pool=FERS_SEQUENCES)
else:
    sequences = default_fers_sequences()

print("\nRunning FERS sequences:")
for seq in sequences:
    print(f"  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
