#!/usr/bin/env python3
"""Interactive 3D "salami slice" view of the ADC position maps.

Same per-run, per-channel mean ADC as make_adc_position_maps.py, but instead
of one flat page per scan position, every run is stacked along a length axis
at its table position X, colour still ADC. A slider walks slice-by-slice
along that axis (Plotly.js, self-contained HTML, no server needed).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple

from utils.scan3d import build_channel_geometry, write_scan3d_html

DEFAULT_INPUT = "results/attenuation/runs_1897_1942_no_pedestal/attenuation_summary.json"
DEFAULT_OUTDIR = (
    "results/attenuation/runs_1897_1942_no_pedestal/"
    "exponential_fit_per_channel/integralmaps"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Interactive 3D (x, y, scan position) ADC view."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--observable", choices=("auto", "adc", "pe"), default="auto",
        help=(
            "Quantity stored in the input points. 'auto' selects mean_pe when "
            "present, otherwise mean_adc (default: auto)."
        ),
    )
    parser.add_argument(
        "--output-name", default=None,
        help="Output HTML filename. Default: adc_scan3d[_relative].html.",
    )
    parser.add_argument(
        "--map-run", type=int, default=1897,
        help="Run whose channel map supplies the physical tower positions.",
    )
    parser.add_argument(
        "--x-min", type=float, default=None,
        help="Lowest table position X to include, in mm. Default: all runs.",
    )
    parser.add_argument(
        "--x-max", type=float, default=None,
        help="Highest table position X to include, in mm. Default: all runs.",
    )
    parser.add_argument(
        "--relative", action="store_true",
        help=(
            "Plot each channel's mean ADC divided by its average over the "
            "selected runs, instead of the raw mean."
        ),
    )
    return parser.parse_args()


def load_points(
    path: Path, observable: str = "auto",
) -> Tuple[Dict[int, float], Dict[int, Dict[str, float]], str]:
    """Return position/value maps and the resolved ADC-or-PE observable."""
    with path.open() as handle:
        source = json.load(handle)
    first_point = next(
        (
            point
            for entries in source.get("points", {}).values()
            for point in entries[:1]
        ),
        {},
    )
    if observable == "auto":
        observable = "pe" if "mean_pe" in first_point else "adc"
    value_key = "mean_pe" if observable == "pe" else "mean_adc"
    positions: Dict[int, float] = {}
    values: Dict[int, Dict[str, float]] = {}
    for branch, entries in source["points"].items():
        for entry in entries:
            run = int(entry["run"])
            mean = float(entry[value_key])
            positions.setdefault(run, float(entry["position_mm"]))
            if math.isfinite(mean):
                values.setdefault(run, {})[branch] = mean
    return positions, values, observable


def main() -> None:
    args = parse_args()
    positions, values, args.observable = load_points(
        Path(args.input), args.observable
    )
    is_pe = args.observable == "pe"

    if args.x_min is not None and args.x_max is not None:
        runs = [run for run, x in positions.items() if args.x_min <= x <= args.x_max]
    else:
        runs = list(positions)
    if not runs:
        raise RuntimeError("No runs selected from {}".format(args.input))

    if args.relative:
        totals: Dict[str, List[float]] = {}
        for run in runs:
            for branch, value in values.get(run, {}).items():
                totals.setdefault(branch, []).append(value)
        reference = {
            branch: sum(vals) / len(vals)
            for branch, vals in totals.items()
            if vals and abs(sum(vals) / len(vals)) > 1.0e-9
        }
        values = {
            run: {
                branch: value / reference[branch]
                for branch, value in values.get(run, {}).items()
                if branch in reference
            }
            for run in runs
        }

    geometry = build_channel_geometry(map_run=args.map_run)

    # Split the flat per-run branch->value map into per-category maps so each
    # family gets its own colour range, mirroring z_range_for_category.
    values_by_category: Dict[str, Dict[int, Dict[str, float]]] = {}
    for category, channels in geometry.items():
        branch_set = {c["branch"] for c in channels}
        values_by_category[category] = {
            run: {b: v for b, v in values.get(run, {}).items() if b in branch_set}
            for run in runs
        }

    value_label = (
        "mean {} / channel average".format("P.E." if is_pe else "ADC")
        if args.relative else (
            "mean photoelectrons [P.E.]" if is_pe else "mean ADC [counts]"
        )
    )
    value_fmt = "{:.3f}" if args.relative else "{:.1f}"

    output_name = args.output_name or (
        "{}_scan3d{}.html".format(
            "pe" if is_pe else "adc", "_relative" if args.relative else ""
        )
    )
    output_path = Path(args.outdir) / output_name
    write_scan3d_html(
        output_path=output_path,
        page_title="{} vs. face position vs. scan position".format(
            "P.E." if is_pe else "ADC"
        ),
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
