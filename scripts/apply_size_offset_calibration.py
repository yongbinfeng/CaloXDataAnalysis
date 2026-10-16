#!/usr/bin/env python3
"""Timing calibration Step 2: align 6 mm channels to the 3 mm y-intercept.

Reads a timing_summary.json produced with --reference-mpv-logic
--slope-family material_size (Step 1: per-(material, size) velocity fits,
one shared slope and one intercept per channel). For each material family
(CerQuartz, CerPlastic) it combines the per-channel y-intercepts within each
fiber size into one number using an error-weighted mean, then takes

    offset[family] = intercept_6mm[family] - intercept_3mm[family]

and subtracts that offset from every 6 mm channel's per-run mean TOA
(timing_fits), so the 6 mm group's timing centers land on the 3 mm group's
scale. 3 mm channels are left untouched (3 mm is the calibration anchor).

The result is a timing_summary.json copy with the same structure as the
input, suitable as --input to make_toa_position_maps.py / plot_toa_scan3d.py
to build "post-calibration" heatmaps.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Dict, List, Tuple

FAMILIES = ("CerQuartz", "CerPlastic")
SIZES = ("3 mm", "6 mm")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Shift 6 mm channel timing to match the 3 mm y-intercept, per "
            "material family, using an existing by-size velocity-fit summary."
        )
    )
    parser.add_argument(
        "--input",
        default="results/CalibrationTiming/StepOneFits/timing_summary.json",
        help="timing_summary.json from a --slope-family material_size run.",
    )
    parser.add_argument(
        "--outdir",
        default="results/CalibrationTiming/StepTwoPostCalibration",
    )
    parser.add_argument(
        "--output-name", default="timing_summary_corrected.json",
        help="Filename for the corrected timing_summary.json copy.",
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
        raise RuntimeError("No channels with a positive intercept error to combine.")
    mean = numerator / weight_sum
    mean_err = weight_sum ** -0.5
    return mean, mean_err


def group_intercepts(summary: Dict) -> Dict[str, Dict[str, List[Tuple[float, float]]]]:
    channels = summary["channels"]
    velocity_fits = summary["velocity_fits"]
    grouped: Dict[str, Dict[str, List[Tuple[float, float]]]] = {
        family: {size: [] for size in SIZES} for family in FAMILIES
    }
    for branch, info in channels.items():
        material = info.get("material")
        size = info.get("size")
        if material not in FAMILIES or size not in SIZES:
            continue
        fit = velocity_fits.get(branch)
        if not fit or not fit.get("fit_ok"):
            continue
        intercept = fit.get("intercept_ns")
        err = fit.get("intercept_err_ns")
        if intercept is None or err is None:
            continue
        grouped[material][size].append((float(intercept), float(err)))
    return grouped


def compute_offsets(summary: Dict) -> Dict[str, Dict]:
    grouped = group_intercepts(summary)
    offsets: Dict[str, Dict] = {}
    for family in FAMILIES:
        n3, n6 = len(grouped[family]["3 mm"]), len(grouped[family]["6 mm"])
        if n3 == 0 or n6 == 0:
            raise RuntimeError(
                "Family {} has {} 3mm and {} 6mm channels with valid "
                "intercepts; need at least one of each.".format(family, n3, n6)
            )
        mean3, err3 = weighted_mean(grouped[family]["3 mm"])
        mean6, err6 = weighted_mean(grouped[family]["6 mm"])
        offset = mean6 - mean3
        offset_err = (err3 ** 2 + err6 ** 2) ** 0.5
        offsets[family] = {
            "intercept_3mm_ns": mean3,
            "intercept_3mm_err_ns": err3,
            "n_channels_3mm": n3,
            "intercept_6mm_ns": mean6,
            "intercept_6mm_err_ns": err6,
            "n_channels_6mm": n6,
            "offset_6mm_minus_3mm_ns": offset,
            "offset_err_ns": offset_err,
        }
    return offsets


def apply_offsets(summary: Dict, offsets: Dict[str, Dict]) -> Dict:
    corrected = copy.deepcopy(summary)
    channels = summary["channels"]
    timing_fits = corrected["timing_fits"]
    n_shifted = 0
    for branch, info in channels.items():
        if info.get("size") != "6 mm":
            continue
        family = info.get("material")
        if family not in offsets:
            continue
        offset = offsets[family]["offset_6mm_minus_3mm_ns"]
        per_run = timing_fits.get(branch)
        if not per_run:
            continue
        for entry in per_run.values():
            mean = entry.get("mean_ns")
            if mean is not None:
                entry["mean_ns"] = mean - offset
        n_shifted += 1
    corrected["size_offset_calibration"] = {
        "anchor": "3 mm",
        "shifted": "6 mm",
        "offsets": offsets,
        "n_channels_shifted": n_shifted,
    }
    return corrected


def main() -> None:
    args = parse_args()
    with open(args.input) as handle:
        summary = json.load(handle)

    offsets = compute_offsets(summary)
    print("Per-family 3mm/6mm y-intercept offset (6mm shifted to match 3mm):")
    for family, info in offsets.items():
        print(
            "  {:<11} 3mm = {:7.4f} +/- {:.4f} ns ({} ch)   "
            "6mm = {:7.4f} +/- {:.4f} ns ({} ch)   "
            "offset = {:+.4f} +/- {:.4f} ns".format(
                family,
                info["intercept_3mm_ns"], info["intercept_3mm_err_ns"], info["n_channels_3mm"],
                info["intercept_6mm_ns"], info["intercept_6mm_err_ns"], info["n_channels_6mm"],
                info["offset_6mm_minus_3mm_ns"], info["offset_err_ns"],
            )
        )

    corrected = apply_offsets(summary, offsets)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / args.output_name
    with output_path.open("w") as handle:
        json.dump(corrected, handle, indent=2, sort_keys=True)
    print("Shifted {} 6mm channel(s).".format(corrected["size_offset_calibration"]["n_channels_shifted"]))
    print("Wrote {}".format(output_path))

    offsets_path = outdir / "size_offsets.json"
    with offsets_path.open("w") as handle:
        json.dump(offsets, handle, indent=2, sort_keys=True)
    print("Wrote {}".format(offsets_path))


if __name__ == "__main__":
    main()
