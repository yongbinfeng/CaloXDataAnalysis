#!/usr/bin/env python3
"""Per-run 1D overlay of mean fitted TOA, one histogram per 3mm column.

The tb2026 3mm DRS channels sit in exactly four tower-X columns (any
material - Sci, CerQuartz, CerPlastic - shares the same four columns). For
each scan run this fills one TH1D per column with the mean fitted TOA (the
same per-run, per-channel timing centre used in make_toa_position_maps.py /
plot_velocityfits.py) of every 3mm channel in that column, then overlays the
four histograms - one colour per column, left to right - on a single page.
Pages are ordered by table position, one per scan run, mirroring
make_toa_position_maps.py.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from channels.channel_map import build_drs_boards
from utils.root_setup import setup_root
from utils.visualization import _X_SCALE

# Left-to-right column colours, reused from compare_timing_runs.py's
# run-colour convention so plots across the repo read consistently.
COLUMN_COLORS = [ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Per-run 1D overlay of mean fitted TOA, one histogram per 3mm "
            "tower-X column."
        )
    )
    parser.add_argument(
        "--input",
        default=(
            "results/timing_velocity_smoke/runs_1897_1943_reference_mpv/"
            "timing_summary.json"
        ),
        help="timing_summary.json holding per-run, per-channel fitted TOA.",
    )
    parser.add_argument(
        "--outdir",
        default="results/timing_velocity_smoke/TimingHeatmaps3mm",
    )
    parser.add_argument(
        "--x-min", type=float, default=None,
        help="Lowest table position to include, in mm. Default: no cut (all runs).",
    )
    parser.add_argument(
        "--x-max", type=float, default=None,
        help="Highest table position to include, in mm. Default: no cut (all runs).",
    )
    parser.add_argument(
        "--map-run", type=int, default=1897,
        help=(
            "Run whose channel map supplies the tower-X column each 3mm "
            "channel belongs to. Must be a tb2026 run (>= 1896). Default: 1897."
        ),
    )
    parser.add_argument(
        "--toa-bins", type=int, default=50,
        help="Bins for each column's TOA histogram (default: 50).",
    )
    parser.add_argument(
        "--toa-min", type=float, default=None,
        help="Override the common TOA-axis minimum (ns).",
    )
    parser.add_argument(
        "--toa-max", type=float, default=None,
        help="Override the common TOA-axis maximum (ns).",
    )
    parser.add_argument(
        "--range-percentile", type=float, default=1.0,
        help=(
            "Percentile trimmed from each end when fixing the shared TOA "
            "axis range across all pages. Default: 1."
        ),
    )
    return parser.parse_args()


def load_points(path: Path) -> Tuple[Dict[int, float], Dict[int, Dict[str, float]]]:
    """Return (position by run, mean TOA [ns] by run and channel).

    Same source as make_toa_position_maps.py: timing_summary.json's
    "run_positions_mm" and "timing_fits".
    """
    with path.open() as handle:
        source = json.load(handle)
    positions = {int(run): float(x) for run, x in source["run_positions_mm"].items()}
    values: Dict[int, Dict[str, float]] = {}
    for branch, per_run in source["timing_fits"].items():
        for run_str, entry in per_run.items():
            if not entry.get("fit_ok"):
                continue
            mean = entry.get("mean_ns")
            if mean is None or not math.isfinite(mean):
                continue
            run = int(run_str)
            if run not in positions:
                continue
            values.setdefault(run, {})[branch] = float(mean)
    return positions, values


def build_columns(map_run: int) -> List[Tuple[str, set]]:
    """Group every 3mm channel's branch name by tower-X column.

    Returns [(label, branches), ...] ordered left to right (ascending
    tower-X, matching the heatmaps' X axis), pooling all materials (Sci,
    CerQuartz, CerPlastic) that share a column.
    """
    drsboards = build_drs_boards(run_number=map_run)
    by_column: Dict[float, set] = {}
    for _, board in drsboards.items():
        for channel in board:
            if channel.is_reference or channel.is6mm:
                continue
            branch = channel.get_channel_name(blsub=False)
            by_column.setdefault(channel.i_tower_x, set()).add(branch)
    columns = []
    for i_tower_x in sorted(by_column):
        x_cm = i_tower_x * _X_SCALE
        label = "X = {:+.1f} cm".format(x_cm)
        columns.append((label, by_column[i_tower_x]))
    return columns


def shared_toa_range(
    values: Dict[int, Dict[str, float]],
    runs: List[int],
    columns: List[Tuple[str, set]],
    args: argparse.Namespace,
) -> Tuple[float, float]:
    if args.toa_min is not None and args.toa_max is not None:
        return args.toa_min, args.toa_max
    all_branches = set().union(*(branches for _, branches in columns))
    pool = sorted(
        value
        for run in runs
        for branch, value in values.get(run, {}).items()
        if branch in all_branches
    )
    if not pool:
        return 0.0, 1.0
    frac = max(0.0, min(0.45, args.range_percentile / 100.0))
    lo = pool[int(frac * (len(pool) - 1))]
    hi = pool[int((1.0 - frac) * (len(pool) - 1))]
    if hi <= lo:
        lo, hi = pool[0], max(pool[-1], pool[0] + 1.0)
    return (args.toa_min if args.toa_min is not None else lo,
            args.toa_max if args.toa_max is not None else hi)


def draw_pad_label(text: str, subtext: str) -> List[Any]:
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.040)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.12, 0.945, text)
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.88, 0.945, subtext)
    return [latex]


def main() -> None:
    args = parse_args()
    setup_root(batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    positions, values = load_points(Path(args.input))
    if args.x_min is not None and args.x_max is not None:
        runs = sorted(
            (run for run, x in positions.items()
             if args.x_min <= x <= args.x_max and run in values),
            key=lambda run: positions[run],
        )
    else:
        runs = sorted((run for run in positions if run in values),
                      key=lambda run: positions[run])
    if not runs:
        raise RuntimeError("No runs with valid TOA fits found in {}".format(args.input))
    print("Selected {} run(s), X = {:.1f} .. {:.1f} mm:".format(
        len(runs), positions[runs[0]], positions[runs[-1]]))

    columns = build_columns(args.map_run)
    if len(columns) != 4:
        print("Warning: found {} 3mm tower-X columns, not 4 - "
              "check --map-run.".format(len(columns)))
    for label, branches in columns:
        print("  column '{}': {} channels".format(label, len(branches)))

    toa_min, toa_max = shared_toa_range(values, runs, columns, args)
    print("Shared TOA axis: {:.3f} .. {:.3f} ns (fixed across all pages)".format(
        toa_min, toa_max))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / "toa_column_hist_3mm.pdf"
    root_path = outdir / "toa_column_hist_3mm.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")

    canvas = ROOT.TCanvas("toa_column_hist_3mm", "", 900, 700)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.11)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    canvas.Print(str(pdf_path) + "[")

    page_summary: List[Dict[str, Any]] = []
    for run in runs:
        value_map = values.get(run, {})
        keep: List[Any] = []
        hists: List[ROOT.TH1D] = []
        for index, (label, branches) in enumerate(columns):
            hist = ROOT.TH1D(
                "h_toa_col{}_run{}".format(index, run), "",
                args.toa_bins, toa_min, toa_max,
            )
            for branch in branches:
                if branch in value_map:
                    hist.Fill(value_map[branch])
            color = COLUMN_COLORS[index % len(COLUMN_COLORS)]
            hist.SetLineColor(color)
            hist.SetLineWidth(3)
            hist.SetFillStyle(0)
            hist.SetStats(0)
            hist.GetXaxis().SetTitle("Mean fitted TOA [ns]")
            hist.GetYaxis().SetTitle("Channels")
            hist.GetXaxis().SetTitleSize(0.045)
            hist.GetYaxis().SetTitleSize(0.045)
            hist.GetXaxis().SetLabelSize(0.038)
            hist.GetYaxis().SetLabelSize(0.038)
            hist.GetYaxis().SetTitleOffset(1.30)
            hists.append(hist)

        ymax = max((h.GetMaximum() for h in hists), default=1.0)
        canvas.Clear()
        legend = ROOT.TLegend(0.68, 0.65, 0.93, 0.87)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.030)
        for index, hist in enumerate(hists):
            hist.SetMaximum(1.25 * max(ymax, 1.0))
            hist.Draw("HIST" if index == 0 else "HIST SAME")
            legend.AddEntry(hist, columns[index][0], "l")
        legend.Draw()
        keep.extend(draw_pad_label(
            "3mm channels by column",
            "run {}   X = {:.1f} mm".format(run, positions[run]),
        ))
        canvas.Update()
        canvas.Print(str(pdf_path))
        root_file.cd()
        canvas.Write("toa_column_hist_3mm_run{}".format(run))
        canvas._keepalive = (hists, legend, keep)

        page_summary.append({
            "run": run,
            "position_mm": positions[run],
            "columns": {
                columns[index][0]: {
                    "n_channels_filled": int(hists[index].GetEntries()),
                }
                for index in range(len(columns))
            },
        })
        print("  page {}/{}: run {} at X = {:.1f} mm".format(
            len(page_summary), len(runs), run, positions[run]))

    canvas.Print(str(pdf_path) + "]")
    root_file.Close()

    summary = {
        "input": args.input,
        "map_run": args.map_run,
        "position_window_mm": (
            None if args.x_min is None else [args.x_min, args.x_max]
        ),
        "columns": [label for label, _ in columns],
        "toa_range_ns": [toa_min, toa_max],
        "n_pages": len(page_summary),
        "pages": page_summary,
    }
    with (outdir / "toa_column_hist_3mm_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("")
    print("Wrote {} ({} pages, one per scan run)".format(pdf_path, len(runs)))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(outdir / "toa_column_hist_3mm_summary.json"))


if __name__ == "__main__":
    main()
