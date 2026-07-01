"""
DRS analysis under MCP timing selection.

Books all lazy DRS histograms (waveform profiles, pulse stats, CFD timing) using
only events that pass the MCP clean-pulse filter, triggers ROOT.RDF.RunGraphs
once, saves output files, then generates HTML plot pages.

Usage
-----
  python scripts/check_drs_mcp.py                              # default DRS_MCP_SEQUENCES
  python scripts/check_drs_mcp.py --sequences drs_stats drs_cfd_mpv
  python scripts/check_drs_mcp.py --channels data/channel_maps/testingfibers.json

  --channels FILE  JSON file mapping labels to channel names; only those channels
                   are histogrammed (reference channels are always kept).

Available sequences (see analysis/registry.py  DRS_MCP_SEQUENCES):
  drs_profiles               — per-channel DRS waveform profiles vs TS
  drs_prof_combined          — combined waveform profiles per fiber size/type (book-only)
  drs_prof_corr_combined     — combined waveform profiles vs time [ns] (book-only, off by default)
  drs_stats                  — peak/energy/timing/noise distributions + sum/saturation/noise board maps
  drs_finebins_combined      — fine-binned CFD-TS combined distributions (book-only, off by default)
  drs_finebins_corr_combined — fine-binned CFD-time [ns] combined distributions (book-only, off by default)
  drs_cfd_mpv                — CFD MPV board maps (plot-only; needs drs_stats finebins)
"""

import json

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
manager = CaloXAnalysisManager(args)

if args.channels:
    with open(args.channels) as f:
        ch_map = json.load(f)
    allowed = set(ch_map.values())
    for board in manager.drsboards.values():
        board.channels = [
            ch for ch in board.channels
            if ch.is_reference or ch.get_channel_name(blsub=False) in allowed
        ]
    print(f"Channel filter: {len(allowed)} channels from {args.channels}")

manager.prepare(do_fers=False)
#manager.apply_hole_veto(flag_only=True)
#manager.apply_mcp_diff_selection()

if args.mcp_clean:
    manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
if args.particle:
    manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

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
