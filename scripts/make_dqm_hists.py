"""
CaloX DRS DQM histogram-booking entry point.

Builds the RDataFrame via CaloXAnalysisManager with full DRS processing, then
books DRS monitoring histograms.  For FERS-only runs (no DRS data), use
make_fers_dqm_hists.py instead.

Usage
-----
  python scripts/make_dqm_hists.py                          # default sequences
  python scripts/make_dqm_hists.py --sequences drs_stats drs_waveforms
  python scripts/make_dqm_hists.py --sequences drs_peak_ts

Available sequences (see analysis/registry.py):
  drs_mapping  drs_waveforms  drs_profiles  drs_stats  drs_peak_ts
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase
from analysis.registry import DRS_DQM_SEQUENCES, default_drs_dqm_sequences, select_sequences
from utils.root_setup import setup_root
from utils.parser import get_args
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root()

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame with DRS processing
# ------------------------------------------------------------------
manager = CaloXAnalysisManager(args).prepare()

# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences, pool=DRS_DQM_SEQUENCES)
else:
    sequences = default_drs_dqm_sequences()

print("\nRunning DRS histogram sequences:")
for seq in sequences:
    has_book = "hist" if seq.book_hists else "    "
    print(f"  [{has_book}]  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)


if __name__ == "__main__":
    main()
