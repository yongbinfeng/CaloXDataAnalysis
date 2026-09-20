"""
In-time vs out-of-time energy and amplitude distributions for the test fibres
(standalone entry point; see analysis/intime_analysis.py for what is produced).

Usage
-----
  python scripts/compare/plot_drs_intime_energy_testfibers.py --run 2030 \
      --channels channels/data/testingfibers.json --mcp-clean

The same analysis runs as part of scripts/check_drs_mcp.py with --intime,
right after the DRS stats it depends on have been written.
"""

import json
import os

from analysis.intime_analysis import run_intime_analysis
from core.analysis_manager import CaloXAnalysisManager
from core.plot_manager import PlotManager
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer
auto_timer("Total Execution Time")


def main():
    setup_root(n_threads=10, batch_mode=True, load_functions=True)
    args = get_args()
    if not args.channels:
        raise SystemExit("--channels FILE is required (label -> channel map)")

    stats_path = f"results/root/Run{args.run}/drs_stats.root"
    if not os.path.exists(stats_path):
        raise SystemExit(f"{stats_path} not found; run check_drs_mcp.py for this run first "
                         f"(or use check_drs_mcp.py --intime to do both)")

    manager = CaloXAnalysisManager(args)
    ch_map = json.load(open(args.channels))
    allowed = set(ch_map.values())
    for board in manager.drsboards.values():
        board.channels = [c for c in board.channels
                          if c.is_reference or c.get_channel_name(blsub=False) in allowed]
    manager.prepare(do_fers=False)
    if args.mcp_clean:
        manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
    if args.particle:
        manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

    run_intime_analysis(manager, args)
    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
