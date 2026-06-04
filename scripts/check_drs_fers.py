"""
DRS integral vs FERS energy correlations.

Books and plots 2D histograms of DRS sum/peak vs FERS HG and LG per channel,
for both Cer and Sci. Requires both DRS (with integral variables) and FERS to
be processed.

Usage
-----
  python scripts/check_drs_fers.py
  python scripts/check_drs_fers.py --sequences drs_sum_vs_fers
  python scripts/check_drs_fers.py --sequences drs_sum_vs_fers drs_peak_vs_fers

Available sequences (see analysis/registry.py):
  drs_sum_vs_fers   — DRS integral sum vs FERS HG/LG 2D correlations
  drs_peak_vs_fers  — DRS peak value vs FERS HG/LG 2D correlations
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase, run_plot_phase
from analysis.registry import DRS_FERS_SEQUENCES, select_sequences
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame: DRS + FERS processing
# ------------------------------------------------------------------
manager = (
    CaloXAnalysisManager(args)
    .prepare()
    .calibrate_fers()
)

# Define DRS integral sum and peak aliases expected by book_drs_sum_vs_fers
# and book_drs_peak_vs_fers (get_channel_sum_name / get_channel_peak_name).
# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences, pool=DRS_FERS_SEQUENCES)
else:
    sequences = list(DRS_FERS_SEQUENCES)

print("\nRunning DRS vs FERS sequences:")
for seq in sequences:
    print(f"  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
