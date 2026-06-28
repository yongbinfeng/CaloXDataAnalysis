"""
DRS SumMean (per-channel mean DRS energy) position scan.

For each DRS channel, plots DRS_SumMean -- the per-channel mean DRS CFD energy
(vals[0] of drs_energy_stats.json, same quantity shown on the DRS_Sum_Map Mean
board map) -- as a function of the table position X [mm] across a set of scan
runs, one graph per channel.

Runs are configurable via SCAN_RUNS below; their table position X is read from
data/Runlist.json. Runs whose JSON or position is missing are skipped. Output:
one HTML page per fiber type under results/html/PositionScan.

Usage
-----
  python scripts/compare/plot_drs_summean_scan.py
"""

import os
import json
import ROOT

from channels.channel_map import build_drs_boards
from core.plot_manager import PlotManager
from plotting.my_function import DrawHistos
from plotting.calox_plot_helper import create_pave_text
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var

auto_timer("Total Execution Time")
setup_root(batch_mode=True, load_functions=False)

SCAN_RUNS = []
for run_number in range(1897, 1928):
    SCAN_RUNS.append(run_number)
REFERENCE_RUN = 1897   # used to build the DRS board/channel list

PLOTDIR = "results/plots/PositionScan"
HTMLDIR = "results/html/PositionScan"


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def _collect_summean(run_x):
    """Read DRS_SumMean (mean energy, vals[0]) per channel from each run's JSON.

    run_x maps run -> table position X [mm]. Returns {channel: [(x, mean), ...]}.
    """
    data = {}
    for run, x_val in sorted(run_x.items(), key=lambda kv: kv[1]):
        path = f"results/root/Run{run}/drs_energy_stats.json"
        if not os.path.exists(path):
            print(f"Warning: {path} not found, skipping run {run}")
            continue
        with open(path) as f:
            stats = json.load(f)
        for ch, vals in stats.items():
            data.setdefault(ch, []).append((x_val, vals[0]))
    return data


def _build_channel_order(drsboards):
    """Return list of (ch, board_no, group_no, chan_no, var, tower_x, tower_y)."""
    order = []
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            order.append((
                chan.get_channel_name(blsub=False),
                board.board_no, chan.group_no, chan.channel_no,
                get_channel_var(chan), chan.i_tower_x, chan.i_tower_y))
    return order


# ---------------------------------------------------------------------------
# Per-channel scan plots
# ---------------------------------------------------------------------------

def _plot_channels(channel_order, data, pm, xmin, xmax):
    """Draw one SumMean-vs-run TGraph per channel. Returns plots grouped by type."""
    plots_by_type = {}
    for ch, b, g, c, var, tx, ty in channel_order:
        pts = sorted(data.get(ch, []))
        if not pts:
            continue

        gr = ROOT.TGraph(len(pts))
        gr.SetName(f"gr_summean_{ch}")
        gr.SetTitle("")
        for i, (x, mean) in enumerate(pts):
            gr.SetPoint(i, x, mean)
        gr.SetLineWidth(2)

        ys = [m for _, m in pts]
        ymid = (max(ys) + min(ys)) / 2
        yhalf = max(50.0, (max(ys) - min(ys)) / 2 * 1.4 + 20)

        pave = create_pave_text(0.18, 0.77, 0.60, 0.87)
        pave.AddText(f"B: {b}, G: {g}, C: {c}")
        pave.AddText(f"Tower: ({tx}, {ty})")

        plot_name = f"drs_summean_scan_{ch}_{var}"
        DrawHistos(
            [gr], [],
            xmin, xmax, "Table position X [mm]",
            ymid - yhalf, ymid + yhalf, "Mean DRS Sum [ADC]",
            plot_name,
            outdir=pm.get_output_dir(),
            drawoptions="EPL",
            markerstyles=[20],
            mycolors=[ROOT.kBlue + 1],
            dology=False,
            usePDF=False,
            extra_text=var,
            extraToDraw=[pave],
        )
        plots_by_type.setdefault(var, []).append(plot_name)
    return plots_by_type


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    with open("data/Runlist.json") as f:
        runlist = json.load(f)

    def _parse_x(run):
        raw = runlist[str(run)]["Table position X"]
        return float(raw.replace("mm", "").strip())

    run_x = {r: _parse_x(r) for r in SCAN_RUNS if str(r) in runlist}
    if not run_x:
        print("No table positions found for the scan runs in data/Runlist.json.")
        return

    drsboards = build_drs_boards(run_number=REFERENCE_RUN)
    data = _collect_summean(run_x)
    channel_order = _build_channel_order(drsboards)

    if not any(data.values()):
        print("No drs_energy_stats.json data found for any scan run.")
        return
    xs_all = list(run_x.values())
    xmin, xmax = min(xs_all) - 20, max(xs_all) + 20

    pm = PlotManager(
        rootdir=f"results/root/Run{REFERENCE_RUN}",
        plotdir=PLOTDIR,
        htmldir=HTMLDIR,
    )
    pm.set_output_dir("DRS_SumMean_Scan")
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    plots_by_type = _plot_channels(channel_order, data, pm, xmin, xmax)

    run_list = ", ".join(str(r) for r in sorted(run_x))
    for type_name, plot_names in plots_by_type.items():
        pm.reset_plots()
        for name in plot_names:
            pm.add_plot(name)
        pm.generate_html(
            f"DRS_SumMean_Scan_{type_name}.html",
            plots_per_row=4,
            title=f"DRS SumMean Scan — {type_name}",
            intro_text=(
                "Per-channel mean DRS CFD energy (DRS_SumMean) vs table "
                f"position X [mm] (runs {run_list})."
            ),
        )

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
