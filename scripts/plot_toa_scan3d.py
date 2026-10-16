#!/usr/bin/env python3
"""Interactive 3D "salami slice" view of the TOA position maps.

Same per-run, per-channel fitted TOA as make_toa_position_maps.py, but
instead of one flat page per scan run, every run is stacked along a length
axis at its table position X, colour still TOA. A slider walks slice-by-slice
along that axis (Plotly.js, self-contained HTML, no server needed).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

from utils.scan3d import build_channel_geometry, write_scan3d_html

DEFAULT_INPUT = (
    "results/timing_velocity_smoke/runs_1897_1943_reference_mpv/timing_summary.json"
)
DEFAULT_OUTDIR = "results/timing_velocity_smoke/TimingHeatmaps"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactive 3D (x, y, scan position) TOA view."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--output-name", default=None,
        help="Output HTML filename. Default: toa_scan3d[_relative][_ps].html.",
    )
    parser.add_argument(
        "--map-run", type=int, default=1897,
        help="Run whose channel map supplies the physical tower positions.",
    )
    parser.add_argument(
        "--x-min", type=float, default=None,
        help="Lowest table position to include, in mm. Default: all runs.",
    )
    parser.add_argument(
        "--x-max", type=float, default=None,
        help="Highest table position to include, in mm. Default: all runs.",
    )
    parser.add_argument(
        "--units", choices=("ns", "ps"), default="ns",
        help="Unit for the colour axis. Default: ns.",
    )
    parser.add_argument(
        "--relative", action="store_true",
        help=(
            "Plot each channel's mean TOA minus its average over the "
            "selected runs, instead of the raw mean."
        ),
    )
    return parser.parse_args()


def load_points(path: Path) -> Tuple[Dict[int, float], Dict[int, Dict[str, float]]]:
    """Return (position by run, mean TOA [ns] by run and channel). Mirrors
    make_toa_position_maps.py's load_points."""
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


def main() -> None:
    args = parse_args()
    positions, values = load_points(Path(args.input))

    if args.x_min is not None and args.x_max is not None:
        runs = [
            run for run, x in positions.items()
            if args.x_min <= x <= args.x_max and run in values
        ]
    else:
        runs = [run for run in positions if run in values]
    if not runs:
        raise RuntimeError("No runs with valid TOA fits found in {}".format(args.input))

    if args.relative:
        totals: Dict[str, List[float]] = {}
        for run in runs:
            for branch, value in values.get(run, {}).items():
                totals.setdefault(branch, []).append(value)
        reference = {branch: sum(vals) / len(vals) for branch, vals in totals.items()}
        values = {
            run: {
                branch: value - reference[branch]
                for branch, value in values.get(run, {}).items()
                if branch in reference
            }
            for run in runs
        }

    unit_scale = 1000.0 if args.units == "ps" else 1.0
    if unit_scale != 1.0:
        values = {
            run: {branch: value * unit_scale for branch, value in per_run.items()}
            for run, per_run in values.items()
        }

    geometry = build_channel_geometry(map_run=args.map_run)

    values_by_category: Dict[str, Dict[int, Dict[str, float]]] = {}
    for category, channels in geometry.items():
        branch_set = {c["branch"] for c in channels}
        values_by_category[category] = {
            run: {b: v for b, v in values.get(run, {}).items() if b in branch_set}
            for run in runs
        }

    value_label = (
        "mean TOA - channel average [{}]".format(args.units) if args.relative
        else "mean fitted TOA [{}]".format(args.units)
    )
    value_fmt = "{:.0f}" if args.units == "ps" else "{:.3f}"

    output_name = args.output_name
    if output_name is None:
        output_name = "toa_scan3d_relative" if args.relative else "toa_scan3d"
        if args.units != "ns":
            output_name += "_" + args.units
        output_name += ".html"
    output_path = Path(args.outdir) / output_name
    write_scan3d_html(
        output_path=output_path,
        page_title="TOA vs. face position vs. scan position",
        subtitle="runs {}-{} &middot; {} slices &middot; map from run {}".format(
            min(runs), max(runs), len(runs), args.map_run,
        ),
        geometry=geometry,
        runs=runs,
        positions=positions,
        values_by_category=values_by_category,
        value_label=value_label,
        value_fmt=value_fmt,
        z_axis_title="table scan position X [mm]",
        colorscale="Turbo",
    )
    print("Wrote {} ({} slices)".format(output_path, len(runs)))


if __name__ == "__main__":
    main()
