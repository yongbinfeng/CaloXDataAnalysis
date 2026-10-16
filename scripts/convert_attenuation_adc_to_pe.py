#!/usr/bin/env python3
"""Convert an existing ADC attenuation summary to photoelectrons.

This is the lightweight Step 6.5 path: it does not re-read events and does
not derive a new calibration for every attenuation run.  Instead, it applies
one existing Step 3 linear calibration, channel by channel, to the run means
and uncertainties in an attenuation summary.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer


DEFAULT_ADC_SUMMARY = Path(
    "results/attenuation/runs_1897_1942_no_pedestal/attenuation_summary.json"
)
DEFAULT_CALIBRATION = Path(
    "results/PhotoElectrons/Step3_ADCtoPEConversion/pe_calibration_summary.json"
)
DEFAULT_OUTDIR = Path(
    "results/PhotoElectrons/Step6p5AttenuationCorrections"
)


@dataclass(frozen=True)
class PEAttenuationPoint:
    run: int
    position_mm: float
    mean_pe: float
    mean_error_pe: float
    stddev_pe: float
    n_events: int
    mean_adc: float
    mean_error_adc: float
    stddev_adc: float
    slope_pe_per_adc: float
    intercept_pe: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Apply a shared Step 3 ADC-to-PE calibration to an existing "
            "per-channel ADC attenuation summary."
        )
    )
    parser.add_argument(
        "--adc-summary", type=Path, default=DEFAULT_ADC_SUMMARY,
        help="Input attenuation_summary.json containing ADC run means.",
    )
    parser.add_argument(
        "--calibration-summary", type=Path, default=DEFAULT_CALIBRATION,
        help="Step 3 pe_calibration_summary.json used for every run.",
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--channels", nargs="+", default=None,
        help="Optional exact DRS branch names to include.",
    )
    parser.add_argument(
        "--min-events", type=int, default=2,
        help="Drop points with fewer events (default: 2).",
    )
    parser.add_argument("--plots-per-row", type=int, default=2)
    args = parser.parse_args()
    if args.min_events < 1:
        parser.error("--min-events must be at least 1")
    if args.plots_per_row < 1:
        parser.error("--plots-per-row must be at least 1")
    return args


def load_json(path: Path) -> Any:
    if not path.is_file():
        raise SystemExit("Input file does not exist: {}".format(path))
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")


def collect_pe_points(
    adc_summary: Dict[str, Any], calibration_summary: Dict[str, Any],
    requested_channels: set[str] | None, min_events: int,
) -> Tuple[Dict[str, List[PEAttenuationPoint]], Dict[str, str], List[str]]:
    adc_channels = adc_summary.get("points")
    calibrations = calibration_summary.get("channels")
    if not isinstance(adc_channels, dict):
        raise SystemExit("ADC summary has no object-valued 'points' field")
    if not isinstance(calibrations, dict):
        raise SystemExit("Calibration summary has no object-valued 'channels' field")

    points: Dict[str, List[PEAttenuationPoint]] = {}
    labels: Dict[str, str] = {}
    warnings: List[str] = []

    for branch in sorted(adc_channels):
        if requested_channels is not None and branch not in requested_channels:
            continue
        calibration = calibrations.get(branch)
        if not isinstance(calibration, dict) or calibration.get("skipped"):
            warnings.append("{}: missing or skipped calibration".format(branch))
            continue
        try:
            slope = float(calibration["slope_pe_per_adc"])
            intercept = float(calibration["intercept_pe"])
        except (KeyError, TypeError, ValueError):
            warnings.append("{}: invalid slope or intercept".format(branch))
            continue
        if not math.isfinite(slope) or not math.isfinite(intercept):
            warnings.append("{}: non-finite slope or intercept".format(branch))
            continue

        converted: List[PEAttenuationPoint] = []
        for raw_point in adc_channels[branch]:
            try:
                run = int(raw_point["run"])
                position = float(raw_point["position_mm"])
                mean_adc = float(raw_point["mean_adc"])
                mean_error_adc = float(raw_point["mean_error_adc"])
                stddev_adc = float(raw_point["stddev_adc"])
                n_events = int(raw_point["n_events"])
            except (KeyError, TypeError, ValueError):
                warnings.append("{}: malformed ADC point skipped".format(branch))
                continue
            values = (position, mean_adc, mean_error_adc, stddev_adc)
            if n_events < min_events or not all(math.isfinite(x) for x in values):
                continue
            converted.append(
                PEAttenuationPoint(
                    run=run,
                    position_mm=position,
                    mean_pe=slope * mean_adc + intercept,
                    mean_error_pe=abs(slope) * mean_error_adc,
                    stddev_pe=abs(slope) * stddev_adc,
                    n_events=n_events,
                    mean_adc=mean_adc,
                    mean_error_adc=mean_error_adc,
                    stddev_adc=stddev_adc,
                    slope_pe_per_adc=slope,
                    intercept_pe=intercept,
                )
            )
        if converted:
            points[branch] = converted
            labels[branch] = str(calibration.get("short_label", branch))

    if requested_channels is not None:
        missing = sorted(requested_channels - set(adc_channels))
        warnings.extend("{}: not present in ADC summary".format(x) for x in missing)
    return points, labels, warnings


def draw_header(canvas: Any, title: str, subtitle: str) -> Any:
    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.034)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.14, 0.944, "CaloX Data")
    latex.SetTextFont(42)
    latex.SetTextSize(0.023)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.94, 0.944, subtitle)
    latex.SetTextAlign(11)
    latex.SetTextSize(0.024)
    latex.DrawLatex(0.14, 0.895, title)
    return latex


def make_canvas(
    branch: str, label: str, values: List[PEAttenuationPoint],
    position_label: str, selection: str,
) -> Tuple[Any, str]:
    ordered = sorted(values, key=lambda point: (point.position_mm, point.run))
    graph = ROOT.TGraphErrors(len(ordered))
    graph.SetName("gr_pe_attenuation_{}".format(safe_name(branch)))
    graph.SetTitle("")
    for index, point in enumerate(ordered):
        graph.SetPoint(index, point.position_mm, point.mean_pe)
        graph.SetPointError(index, 0.0, point.mean_error_pe)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.05)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)

    name = "pe_attenuation_{}".format(safe_name(branch))
    canvas = ROOT.TCanvas(name, "", 1100, 800)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    graph.Draw("AP")
    graph.GetXaxis().SetTitle("{} [mm]".format(position_label))
    graph.GetYaxis().SetTitle("Mean photoelectrons")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.GetXaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetTitleOffset(1.35)

    legend = ROOT.TLegend(0.15, 0.72, 0.57, 0.81)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.023)
    legend.AddEntry(graph, "Converted run mean (#sigma/#sqrt{N} errors)", "pe")
    legend.Draw()
    header = draw_header(canvas, "PE attenuation scan: {}".format(label), selection)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, legend, header]
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    adc_summary = load_json(args.adc_summary)
    calibration_summary = load_json(args.calibration_summary)
    requested = set(args.channels) if args.channels else None
    points, labels, warnings = collect_pe_points(
        adc_summary, calibration_summary, requested, args.min_events
    )
    if not points:
        raise SystemExit("No ADC attenuation points could be converted")
    for warning in warnings:
        print("[WARN] {}".format(warning))

    args.outdir.mkdir(parents=True, exist_ok=True)
    root_path = args.outdir / "pe_attenuation_scans.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    if root_file.IsZombie():
        raise RuntimeError("Could not create {}".format(root_path))

    selection = str(adc_summary.get("selection", "Unspecified selection"))
    position_label = str(adc_summary.get("position_key", "Position"))
    keys: List[str] = []
    canvas_jsons: Dict[str, str] = {}
    for branch in sorted(points):
        canvas, canvas_json = make_canvas(
            branch, labels[branch], points[branch], position_label, selection
        )
        root_file.cd()
        canvas.Write(canvas.GetName())
        keys.append(str(canvas.GetName()))
        canvas_jsons[str(canvas.GetName())] = canvas_json
    root_file.Close()

    source_runs = [int(run) for run in adc_summary.get("runs", [])]
    converted_runs = sorted({point.run for values in points.values() for point in values})
    summary = {
        "source_adc_summary": str(args.adc_summary),
        "source_calibration_summary": str(args.calibration_summary),
        "source_calibration_run": calibration_summary.get("run"),
        "conversion": "mean_pe = slope_pe_per_adc * mean_adc + intercept_pe",
        "uncertainty": (
            "ADC statistical errors scaled by abs(slope_pe_per_adc); Step 3 "
            "calibration coefficients are treated as fixed"
        ),
        "observable": "mean photoelectrons from converted ADC attenuation means",
        "position_key": adc_summary.get("position_key"),
        "selection": selection,
        "source_runs": source_runs,
        "runs": converted_runs,
        "channels": len(points),
        "warnings": warnings,
        "points": {
            branch: [
                asdict(point)
                for point in sorted(values, key=lambda point: point.run)
            ]
            for branch, values in points.items()
        },
    }
    json_path = args.outdir / "pe_attenuation_summary.json"
    with json_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")

    html_path = args.outdir / "pe_attenuation_scans.html"
    generate_jsroot_html(
        keys, canvas_jsons, plots_per_row=args.plots_per_row,
        output_html=str(html_path),
        title="CaloX TB2026 per-channel PE attenuation scans (Step 6.5)",
        intro_text=(
            "**Observable:** existing ADC attenuation run means converted with "
            "PE = slope x ADC + intercept.\n"
            "**ADC source:** {}.\n"
            "**Calibration source:** {} (one shared per-channel Step 3 conversion "
            "for every run).\n"
            "**Error bars:** ADC standard error scaled by |slope|; calibration "
            "coefficients are treated as fixed.\n"
            "**Position:** {} [mm].\n"
            "**Selection:** {}."
        ).format(
            args.adc_summary, args.calibration_summary, position_label, selection
        ),
    )

    print("Converted {} channels and {} distinct runs.".format(
        len(points), len(converted_runs)
    ))
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
