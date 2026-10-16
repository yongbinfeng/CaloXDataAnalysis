#!/usr/bin/env python3
"""Overlay the 3mm and 6mm velocity fits for each material family.

Reads a timing_summary.json produced with --reference-mpv-logic
--slope-family material_size (see plot_velocityfits.py / Step 1 of the
timing calibration). For each family (CerQuartz, CerPlastic) it draws, on
one axes:

  * the per-run mean TOA, averaged (error-weighted) over that family+size's
    channels, as data points with error bars;
  * the shared-family fit line y = intercept + slope * position, one per
    size, using the size's shared slope and its error-weighted mean
    intercept across channels (the same combination Step 2's
    apply_size_offset_calibration.py uses);

and records the 6mm-minus-3mm intercept difference in the legend.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FAMILIES = ("CerQuartz", "CerPlastic")
SIZES = ("3 mm", "6 mm")
COLORS = {"3 mm": "tab:blue", "6 mm": "tab:red"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Overlay 3mm/6mm velocity fits per family, with the "
                     "intercept difference recorded in the legend."
    )
    parser.add_argument(
        "--input",
        default="results/CalibrationTiming/StepOneFits/timing_summary.json",
    )
    parser.add_argument(
        "--outdir",
        default="results/CalibrationTiming/StepOneFits",
    )
    parser.add_argument(
        "--output-name", default="velocity_fits_by_size_overlay.png",
    )
    return parser.parse_args()


def weighted_mean(values: List[Tuple[float, float]]) -> Tuple[float, float]:
    """(value, err) pairs -> error-weighted mean and its combined error."""
    weight_sum = 0.0
    numerator = 0.0
    for value, err in values:
        if err is None or err <= 0.0:
            continue
        weight = 1.0 / (err * err)
        numerator += value * weight
        weight_sum += weight
    if weight_sum <= 0.0:
        return None, None
    return numerator / weight_sum, weight_sum ** -0.5


def group_channels(summary: Dict) -> Dict[str, Dict[str, List[str]]]:
    channels = summary["channels"]
    grouped: Dict[str, Dict[str, List[str]]] = {
        family: {size: [] for size in SIZES} for family in FAMILIES
    }
    for branch, info in channels.items():
        material = info.get("material")
        size = info.get("size")
        if material in FAMILIES and size in SIZES:
            grouped[material][size].append(branch)
    return grouped


def group_intercept(summary: Dict, branches: List[str]) -> Tuple[float, float]:
    velocity_fits = summary["velocity_fits"]
    pairs = []
    for branch in branches:
        fit = velocity_fits.get(branch)
        if not fit or not fit.get("fit_ok"):
            continue
        intercept = fit.get("intercept_ns")
        err = fit.get("intercept_err_ns")
        if intercept is not None and err is not None:
            pairs.append((float(intercept), float(err)))
    return weighted_mean(pairs)


def per_run_points(
    summary: Dict, branches: List[str]
) -> Tuple[List[float], List[float], List[float]]:
    """Error-weighted mean TOA per run, averaged over `branches`."""
    positions = {int(run): float(x) for run, x in summary["run_positions_mm"].items()}
    timing_fits = summary["timing_fits"]
    by_run: Dict[int, List[Tuple[float, float]]] = {}
    for branch in branches:
        per_run = timing_fits.get(branch, {})
        for run_str, entry in per_run.items():
            if not entry.get("fit_ok"):
                continue
            mean = entry.get("mean_ns")
            err = entry.get("point_error_ns")
            if mean is None or err is None or err <= 0.0:
                continue
            run = int(run_str)
            if run not in positions:
                continue
            by_run.setdefault(run, []).append((float(mean), float(err)))

    xs, ys, yerrs = [], [], []
    for run in sorted(by_run, key=lambda r: positions[r]):
        mean, err = weighted_mean(by_run[run])
        if mean is None:
            continue
        xs.append(positions[run])
        ys.append(mean)
        yerrs.append(err)
    return xs, ys, yerrs


def main() -> None:
    args = parse_args()
    with open(args.input) as handle:
        summary = json.load(handle)

    shared_slope_fits = summary["shared_slope_fits"]
    grouped_channels = group_channels(summary)

    fig, axes = plt.subplots(1, len(FAMILIES), figsize=(7 * len(FAMILIES), 6))
    if len(FAMILIES) == 1:
        axes = [axes]

    for ax, family in zip(axes, FAMILIES):
        intercepts = {}
        for size in SIZES:
            branches = grouped_channels[family][size]
            xs, ys, yerrs = per_run_points(summary, branches)
            color = COLORS[size]
            ax.errorbar(
                xs, ys, yerr=yerrs, fmt="o", ms=4, color=color, alpha=0.6,
                label=None,
            )

            slope_fit = shared_slope_fits.get(
                "{} {}".format(family, size.replace(" ", ""))
            )
            intercept, intercept_err = group_intercept(summary, branches)
            intercepts[size] = (intercept, intercept_err)
            if slope_fit is None or intercept is None or not xs:
                continue
            slope = slope_fit["slope_ns_per_mm"]
            speed = slope_fit["speed_cm_per_ns"]
            x_line = [min(xs), max(xs)]
            y_line = [intercept + slope * x for x in x_line]
            ax.plot(
                x_line, y_line, "-", color=color, linewidth=2,
                label="{}: v = {:.2f} cm/ns, intercept = {:.3f} ± {:.3f} ns".format(
                    size, speed, intercept, intercept_err
                ),
            )

        int3, err3 = intercepts.get("3 mm", (None, None))
        int6, err6 = intercepts.get("6 mm", (None, None))
        if int3 is not None and int6 is not None:
            diff = int6 - int3
            diff_err = (err3 ** 2 + err6 ** 2) ** 0.5
            ax.plot(
                [], [], " ",
                label="Δ intercept (6mm - 3mm) = {:+.3f} ± {:.3f} ns".format(
                    diff, diff_err
                ),
            )

        ax.set_title(family)
        ax.set_xlabel("Table position [mm]")
        ax.set_ylabel("Mean fitted TOA [ns]")
        ax.legend(fontsize=8.5, frameon=True, loc="best")
        ax.grid(alpha=0.3)

    fig.tight_layout()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / args.output_name
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print("Wrote {}".format(output_path))

    if output_path.suffix != ".pdf":
        pdf_path = output_path.with_suffix(".pdf")
        fig.savefig(pdf_path, bbox_inches="tight")
        print("Wrote {}".format(pdf_path))


if __name__ == "__main__":
    main()
