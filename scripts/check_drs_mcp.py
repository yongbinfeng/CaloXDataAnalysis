"""
DRS analysis under MCP timing selection.

Books all lazy DRS histograms (waveforms, pulse stats, peak TS) using only
events that pass the MCP clean-pulse filter, triggers ROOT.RDF.RunGraphs once,
saves output files, then generates HTML plot pages.

Usage
-----
  python scripts/check_drs_mcp.py                   # all DRS_MCP_SEQUENCES

Available sequences (see analysis/registry.py  DRS_MCP_SEQUENCES):
  drs_waveforms  — DRS waveform 2D histograms vs TS
  drs_stats      — pulse peak, energy, and timing distributions
  drs_peak_ts    — peak time-slice distributions
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase, run_plot_phase
from analysis.registry import DRS_MCP_SEQUENCES, default_drs_mcp_sequences
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
args.jsroot = True  # default to JSROOT for this script; pass --no-jsroot is not supported, edit this line to disable

# ------------------------------------------------------------------
# Build RDataFrame: DRS baseline subtraction + MCP clean-pulse +
# timing-difference selection applied globally to all sequences
# ------------------------------------------------------------------
manager = (
    CaloXAnalysisManager(args)
    .prepare(do_fers=False)
    #.apply_hole_veto(flag_only=True)
    .apply_mcp_diff_selection()
)

sequences = default_drs_mcp_sequences()

print("\nRunning DRS MCP sequences:")
for seq in sequences:
    print(f"  {seq.name}")
print()


def main():
    run_hist_phase(sequences, manager)
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
