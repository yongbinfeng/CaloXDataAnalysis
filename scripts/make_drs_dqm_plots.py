"""
CaloX DRS DQM plot generation.

Reads the ROOT / JSON files produced by make_drs_dqm_hists.py and generates
HTML plot pages.  For FERS-only runs, use make_fers_dqm_plots.py instead.

Usage
-----
  python scripts/make_drs_dqm_plots.py                          # default sequences
  python scripts/make_drs_dqm_plots.py --sequences drs_stats drs_waveforms
  python scripts/make_drs_dqm_plots.py --sequences drs_peak_ts

Available sequences (see analysis/registry.py  DRS_DQM_SEQUENCES):
  drs_mapping            — channel layout map (plot-only)
  drs_waveforms          — per-channel DRS waveform 2D (vs TS)
  drs_profiles           — per-channel DRS waveform profiles
  drs_prof_combined      — combined waveform profiles (book-only)
  drs_stats              — per-channel stats maps
  drs_finebins_combined  — fine-binned combined profiles (off by default)
  drs_peak_ts            — peak time-slice distributions (off by default)
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_plot_phase
from analysis.registry import DRS_DQM_SEQUENCES, default_drs_dqm_sequences, select_sequences
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=1, batch_mode=True, load_functions=False)

args = get_args()

# ------------------------------------------------------------------
# Build manager without loading data (plot phase only reads saved files)
# ------------------------------------------------------------------
manager = CaloXAnalysisManager(args, load_data=False)

# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences, pool=DRS_DQM_SEQUENCES)
else:
    sequences = default_drs_dqm_sequences()

print("\nRunning DRS plot sequences:")
for seq in sequences:
    has_plot = "plot" if seq.make_plots else "    "
    print(f"  [{has_plot}]  {seq.name}")
print()


def main():
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
