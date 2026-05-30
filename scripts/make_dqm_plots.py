"""
CaloX plot-generation entry point.

Reads the ROOT / JSON files produced by make_dqm_hists.py and generates
HTML plot pages.  Sequences are selected the same way as the hist script.

Usage
-----
  python scripts/make_dqm_plots.py                          # default sequences
  python scripts/make_dqm_plots.py --sequences monitor_conditions drs_stats
  python scripts/make_dqm_plots.py --sequences fers_channels fers_stats fers_max

Available sequences (see analysis/registry.py):
  monitor_conditions  fers_esum_vs_event  fers_mapping  fers_channels  fers_stats
  fers_max  fers_2d  fers_track  drs_waveforms  drs_stats  drs_peak_ts
  drs_sum_vs_fers  drs_peak_vs_fers
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_plot_phase
from analysis.registry import default_sequences, select_sequences
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
    sequences = select_sequences(args.sequences)
else:
    sequences = default_sequences()

print("\nRunning plot sequences:")
for seq in sequences:
    has_plot = "plot" if seq.make_plots else "    "
    print(f"  [{has_plot}]  {seq.name}")
print()


def main():
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
