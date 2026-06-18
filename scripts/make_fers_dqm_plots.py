"""
CaloX FERS-only DQM plot generation.

Reads the ROOT / JSON files produced by make_fers_dqm_hists.py and generates
HTML plot pages.  No DRS data required.

Usage
-----
  python scripts/make_fers_dqm_plots.py                          # default sequences
  python scripts/make_fers_dqm_plots.py --sequences monitor_conditions fers_stats
  python scripts/make_fers_dqm_plots.py --sequences fers_channels fers_stats fers_max

Available sequences (see analysis/registry.py):
  monitor_conditions  fers_esum_vs_event  fers_mapping  fers_channels
  fers_stats  fers_max  fers_2d  fers_track
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_plot_phase
from analysis.registry import FERS_DQM_SEQUENCES, default_fers_dqm_sequences, select_sequences
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
    sequences = select_sequences(args.sequences, pool=FERS_DQM_SEQUENCES)
else:
    sequences = default_fers_dqm_sequences()

print("\nRunning FERS plot sequences:")
for seq in sequences:
    has_plot = "plot" if seq.make_plots else "    "
    print(f"  [{has_plot}]  {seq.name}")
print()


def main():
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
