"""
CaloX FERS-only DQM histogram booking.

Books FERS monitoring histograms without requiring DRS data.  Use this script
for FERS-only runs where no DRS detector data is present.

Usage
-----
  python scripts/make_fers_dqm_hists.py                          # default sequences
  python scripts/make_fers_dqm_hists.py --sequences monitor_conditions fers_stats
  python scripts/make_fers_dqm_hists.py --sequences fers_channels fers_stats fers_max

Available sequences (see analysis/registry.py):
  monitor_conditions  fers_esum_vs_event  fers_mapping  fers_channels
  fers_stats  fers_max  fers_2d  fers_track
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase
from analysis.registry import FERS_DQM_SEQUENCES, default_fers_dqm_sequences, select_sequences
from utils.root_setup import setup_root
from utils.parser import get_args
from utils.timing import auto_timer
from variables.fers import get_fers_energy_max, get_fers_energy_sum

auto_timer("Total Execution Time")

setup_root()

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame: FERS only, skip DRS processing entirely
# ------------------------------------------------------------------
manager = CaloXAnalysisManager(args).prepare(do_drs=False)

manager.rdf = get_fers_energy_max(manager.rdf, manager.fersboards, gain="HG")
manager.rdf = get_fers_energy_max(manager.rdf, manager.fersboards, gain="LG")
manager.rdf = get_fers_energy_sum(manager.rdf, manager.fersboards, gain="HG")
manager.rdf = get_fers_energy_sum(manager.rdf, manager.fersboards, gain="LG")

# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences, pool=FERS_DQM_SEQUENCES)
else:
    sequences = default_fers_dqm_sequences()

print("\nRunning FERS histogram sequences:")
for seq in sequences:
    has_book = "hist" if seq.book_hists else "    "
    print(f"  [{has_book}]  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)


if __name__ == "__main__":
    main()
