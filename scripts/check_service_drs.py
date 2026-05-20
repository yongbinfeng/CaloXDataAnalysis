"""
Service DRS analysis entry point.

Books all lazy histograms, triggers ROOT.RDF.RunGraphs once, saves output
files, then generates HTML plot pages — all in a single script run.

Usage
-----
  python scripts/check_service_drs.py                        # default sequences
  python scripts/check_service_drs.py --sequences service_drs_pid service_drs_mcp

Available sequences (see analysis/registry.py  SERVICE_DRS_SEQUENCES):
  service_drs_pid        — PID detector pulse analysis (HoleVeto, PSD, Cer, …)
  service_drs_mcp        — MCP pulse analysis
  service_drs_mcp_timing — MCP CFD timing differences
  service_drs_hodo       — hodoscope peak positions (off by default)
"""

from core.analysis_manager import CaloXAnalysisManager
from core.sequence import run_hist_phase, run_plot_phase
from analysis.registry import (
    SERVICE_DRS_SEQUENCES,
    default_service_drs_sequences,
    select_sequences,
)
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()

# ------------------------------------------------------------------
# Build RDataFrame: DRS baseline subtraction + MCP clean-pulse
# selection applied globally to all sequences
# ------------------------------------------------------------------
manager = (
    CaloXAnalysisManager(args)
    .prepare()
    .apply_hole_veto(flag_only=True)
    .apply_mcp_selection()
)

# ------------------------------------------------------------------
# Select sequences
# ------------------------------------------------------------------
if args.sequences:
    sequences = select_sequences(args.sequences, pool=SERVICE_DRS_SEQUENCES)
else:
    sequences = default_service_drs_sequences()

print("\nRunning service DRS sequences:")
for seq in sequences:
    print(f"  {seq.name}")
print()


def main():
    # Hist phase: book → RunGraphs → save
    run_hist_phase(sequences, manager)
    # Plot phase: read saved files → HTML
    run_plot_phase(sequences, manager)


if __name__ == "__main__":
    main()
