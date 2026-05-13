"""
CaloX histogram-booking entry point.

Builds the RDataFrame via CaloXAnalysisManager, selects which sequences to run
(via --sequences or the registry defaults), books all lazy ROOT objects, triggers
ROOT.RDF.RunGraphs once, then saves every output file.

Usage
-----
  python scripts/make_dqm_hists.py                          # default sequences
  python scripts/make_dqm_hists.py --sequences monitor_conditions drs_stats
  python scripts/make_dqm_hists.py --sequences fers_channels fers_stats fers_max

Available sequences (see analysis/registry.py):
  monitor_conditions  fers_energy_sum  fers_channels  fers_stats  fers_max
  fers_2d  fers_track  drs_waveforms  drs_stats  drs_peak_ts
  drs_sum_vs_fers  drs_peak_vs_fers
"""

from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase
from analysis.registry import default_sequences, select_sequences
from utils.root_setup import setup_root
from utils.parser import get_args
from utils.timing import auto_timer
from variables.fers import get_fers_energy_max, get_fers_energy_sum

auto_timer("Total Execution Time")

setup_root()

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame and derive FERS columns
# ------------------------------------------------------------------
manager = CaloXAnalysisManager(args).prepare()

manager.rdf = get_fers_energy_max(manager.rdf, manager.fersboards, gain="HG")
manager.rdf = get_fers_energy_max(manager.rdf, manager.fersboards, gain="LG")
manager.rdf = get_fers_energy_sum(manager.rdf, manager.fersboards, gain="HG")
manager.rdf = get_fers_energy_sum(manager.rdf, manager.fersboards, gain="LG")

# MCP filter
channels_mcp = get_mcp_channels(manager.run_number)
ref_det  = list(channels_mcp.keys())[0]
ref_det2 = list(channels_mcp.keys())[4]
manager.rdf = manager.rdf.Filter(
    f"{ref_det}_integral_to_peak > 4 && {ref_det}_peak_value > 20")
manager.rdf = manager.rdf.Filter(
    f"{ref_det2}_integral_to_peak > 4 && {ref_det2}_peak_value > 20")

# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences)
else:
    sequences = default_sequences()

print("\nRunning histogram sequences:")
for seq in sequences:
    has_book = "hist" if seq.book_hists else "    "
    print(f"  [{has_book}]  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)


if __name__ == "__main__":
    main()
