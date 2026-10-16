#!/usr/bin/env python3
"""Overlay existing Step 2 PE fits for runs at several scan positions.

The script does not refit anything. It reads Step 1's stored ADC-integral
histograms plus Step 2's ``pe_fits_summary.json`` and ``pe_fit_curves.json``
files. For each channel, the data spectra are drawn as thin dotted outlines
and the corresponding total multi-peak fits as thick solid curves.
"""

from __future__ import annotations

import argparse
import bisect
import json
from array import array
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import ROOT

from plot_adc_position_scan_overlay import (
    DEFAULT_INPUT_DIR,
    ROOT_COLORS,
    content_range,
    draw_header,
    load_json,
    resolve_target_runs,
    safe_name,
)
from utils.html_generator import generate_jsroot_html


DEFAULT_STEP2_DIR = Path(
    "results/PhotoElectrons/Step6_AttenuationCorrections/Step2_PEFits"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Overlay stored Step 2 PE-fit curves for runs nearest the requested "
            "table positions."
        )
    )
    parser.add_argument("--step1-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--step2-dir", type=Path, default=DEFAULT_STEP2_DIR)
    parser.add_argument(
        "--outdir", type=Path, default=None,
        help="Output directory (default: --step2-dir).",
    )
    parser.add_argument(
        "--target-positions", nargs="+", type=float,
        default=[-200.0, 87.0, 213.0, 386.0, 809.0],
        help="Target table X positions in mm (default: -200 87 213 386 809).",
    )
    parser.add_argument("--runlist", type=Path, default=Path("data/Runlist.json"))
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument(
        "--raw-yields", action="store_true",
        help="Draw raw event counts instead of unit-area-normalized spectra/fits.",
    )
    parser.add_argument(
        "--fits-only", action="store_true",
        help="Hide the dotted Step 1 data outlines and draw only fitted curves.",
    )
    parser.add_argument("--zoom-keep", type=float, default=0.995)
    parser.add_argument("--plots-per-row", type=int, default=3)
    args = parser.parse_args()
    if not args.target_positions:
        parser.error("--target-positions needs at least one value")
    if not 0.0 < args.zoom_keep <= 1.0:
        parser.error("--zoom-keep must be in (0, 1]")
    if args.plots_per_row < 1:
        parser.error("--plots-per-row must be at least 1")
    return args


def _deduplicate_curve(points: Sequence[Sequence[float]]) -> Tuple[List[float], List[float]]:
    values = {float(point[0]): float(point[1]) for point in points if len(point) >= 2}
    xs = sorted(values)
    return xs, [values[x] for x in xs]


def _interpolate(xs: Sequence[float], ys: Sequence[float], x: float) -> float:
    if not xs or x < xs[0] or x > xs[-1]:
        return 0.0
    index = bisect.bisect_left(xs, x)
    if index == 0:
        return ys[0]
    if index == len(xs):
        return ys[-1]
    if xs[index] == x:
        return ys[index]
    x0, x1 = xs[index - 1], xs[index]
    y0, y1 = ys[index - 1], ys[index]
    if x1 == x0:
        return y1
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def make_total_fit_graph(
    branch: str, run: int, component_curves: Dict[str, Any], scale: float,
) -> Optional[Any]:
    """Sum the stored RooFit component curves on their combined x grid."""
    series = []
    x_values = set()
    for peak in sorted(component_curves, key=lambda value: int(value)):
        xs, ys = _deduplicate_curve(component_curves[peak])
        if len(xs) < 2:
            continue
        series.append((xs, ys))
        x_values.update(xs)
    if not series or not x_values:
        return None

    xs = sorted(x_values)
    ys = [
        scale * sum(_interpolate(comp_x, comp_y, x) for comp_x, comp_y in series)
        for x in xs
    ]
    graph = ROOT.TGraph(len(xs), array("d", xs), array("d", ys))
    graph.SetName("total_fit_run{}_{}".format(run, safe_name(branch)))
    return graph


def make_overlay_canvas(
    branch: str, short_label: str, selection: str,
    run_positions: Sequence[Tuple[float, int, float]],
    source_hists: Dict[int, Any], channel_fits: Dict[int, Dict[str, Any]],
    channel_curves: Dict[int, Dict[str, Any]], normalise: bool,
    fits_only: bool, zoom_keep: float,
    colors: Optional[Dict[int, int]] = None,
) -> Optional[Tuple[Any, str]]:
    entries = []
    for target, run, actual in run_positions:
        source = source_hists.get(run)
        fit = channel_fits.get(run)
        curves = channel_curves.get(run)
        if source is None or fit is None or curves is None or source.GetEntries() <= 0:
            continue
        hist = source.Clone("fit_overlay_data_run{}_{}".format(run, safe_name(branch)))
        hist.SetDirectory(0)
        integral = hist.Integral()
        scale = 1.0 / integral if normalise and integral > 0.0 else 1.0
        if normalise and integral > 0.0:
            hist.Scale(scale)
        graph = make_total_fit_graph(branch, run, curves, scale)
        if graph is not None:
            entries.append((target, run, actual, hist, graph, fit))
    if not entries:
        return None

    ranges = [content_range(entry[3], zoom_keep) for entry in entries]
    ranges = [item for item in ranges if item is not None]
    if ranges:
        xmin = min(item[0] for item in ranges)
        xmax = max(item[1] for item in ranges)
    else:
        xmin = min(entry[3].GetXaxis().GetXmin() for entry in entries)
        xmax = max(entry[3].GetXaxis().GetXmax() for entry in entries)

    ymax = max(entry[3].GetMaximum() for entry in entries)
    for _, _, _, _, graph, _ in entries:
        ymax = max(ymax, max(float(graph.GetY()[i]) for i in range(graph.GetN())))
    ymin = 1.0e-5 if normalise else 0.5

    canvas_name = "pe_fit_position_overlay_{}".format(safe_name(branch))
    canvas = ROOT.TCanvas(canvas_name, "", 1050, 780)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    canvas.SetLogy(1)

    frame = canvas.DrawFrame(xmin, ymin, xmax, max(10.0 * ymin, 40.0 * ymax))
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("ADC integral [ADC]")
    frame.GetYaxis().SetTitle("Normalized events" if normalise else "Events")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.35)

    run_legend = ROOT.TLegend(0.45, 0.64, 0.95, 0.885)
    run_legend.SetBorderSize(0)
    run_legend.SetFillStyle(0)
    run_legend.SetTextFont(42)
    run_legend.SetTextSize(0.021)
    keepalive: List[Any] = [frame, run_legend]

    for index, (target, run, actual, hist, graph, fit) in enumerate(entries):
        color = (colors[run] if colors and run in colors
                 else ROOT_COLORS[index % len(ROOT_COLORS)])
        if not fits_only:
            hist.SetLineColor(color)
            hist.SetLineWidth(1)
            hist.SetLineStyle(3)
            hist.SetFillStyle(0)
            hist.SetStats(0)
            hist.SetTitle("")
            hist.Draw("HIST SAME")
            keepalive.append(hist)

        graph.SetLineColor(color)
        graph.SetLineWidth(3)
        graph.SetLineStyle(1)
        graph.Draw("L SAME")
        run_legend.AddEntry(
            graph,
            "Run {}: X={:.1f} mm, gain={:.1f} ADC/p.e.".format(
                run, actual, float(fit["dpe"])
            ),
            "l",
        )
        keepalive.append(graph)

    run_legend.Draw()
    if not fits_only:
        style_legend = ROOT.TLegend(0.69, 0.535, 0.94, 0.62)
        style_legend.SetBorderSize(0)
        style_legend.SetFillStyle(0)
        style_legend.SetTextFont(42)
        style_legend.SetTextSize(0.021)
        data_style = ROOT.TGraph()
        data_style.SetLineColor(ROOT.kGray + 2)
        data_style.SetLineWidth(1)
        data_style.SetLineStyle(3)
        fit_style = ROOT.TGraph()
        fit_style.SetLineColor(ROOT.kGray + 2)
        fit_style.SetLineWidth(3)
        style_legend.AddEntry(data_style, "Step 1 data", "l")
        style_legend.AddEntry(fit_style, "Step 2 total fit", "l")
        style_legend.Draw()
        keepalive.extend([style_legend, data_style, fit_style])

    header = draw_header(canvas, short_label, selection)
    keepalive.append(header)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def load_step2_runs(
    step2_dir: Path, runs: Sequence[int],
) -> Tuple[Dict[int, Any], Dict[int, Any]]:
    summaries: Dict[int, Any] = {}
    curves: Dict[int, Any] = {}
    missing = []
    for run in runs:
        run_dir = step2_dir / "run{}".format(run)
        summary_path = run_dir / "pe_fits_summary.json"
        curves_path = run_dir / "pe_fit_curves.json"
        if not summary_path.is_file() or not curves_path.is_file():
            missing.append(run)
            continue
        summaries[run] = load_json(summary_path)
        curves[run] = load_json(curves_path)
    if missing:
        raise SystemExit(
            "Missing completed Step 2 outputs for run(s): {}. Wait for the Step 2 "
            "fit job to finish, then rerun this overlay command.".format(
                ", ".join(str(run) for run in missing)
            )
        )
    return summaries, curves


def main() -> None:
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    step1_summary_path = args.step1_dir / "pe_integrals_summary.json"
    step1_root_path = args.step1_dir / "pe_integrals.root"
    step1_summary = load_json(step1_summary_path)
    step1_runs = step1_summary.get("channels", {})
    available_runs = [
        int(name[3:]) for name in step1_runs
        if name.startswith("run") and name[3:].isdigit()
    ]
    resolved = resolve_target_runs(
        args.target_positions, args.runlist, args.position_key, available_runs
    )
    selected_runs = [run for _, run, _ in resolved]
    fit_summaries, fit_curves = load_step2_runs(args.step2_dir, selected_runs)

    source = ROOT.TFile.Open(str(step1_root_path), "READ")
    if not source or source.IsZombie():
        raise SystemExit("Cannot open input ROOT file: {}".format(step1_root_path))

    branch_sets = []
    for run in selected_runs:
        fit_channels = fit_summaries[run].get("channels", {})
        good = {
            branch for branch, fit in fit_channels.items()
            if not fit.get("skipped") and branch in fit_curves[run]
        }
        branch_sets.append(good & set(step1_runs.get("run{}".format(run), {})))
    branches = sorted(set.intersection(*branch_sets)) if branch_sets else []
    if args.channels:
        requested = set(args.channels)
        missing = sorted(requested - set(branches))
        if missing:
            print("[WARN] Requested channel(s) unavailable: {}".format(", ".join(missing)))
        branches = [branch for branch in branches if branch in requested]
    if not branches:
        source.Close()
        raise SystemExit("No successfully fitted channels are common to all selected runs.")

    normalise = not args.raw_yields
    selection = str(step1_summary.get("selection", "Step 1 selection"))
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    skipped: List[str] = []
    print("[CHANNELS] Building {} fit overlay canvas(es).".format(len(branches)))

    for index, branch in enumerate(branches, 1):
        source_hists = {}
        channel_fits = {}
        channel_curves = {}
        short_label = branch
        for run in selected_runs:
            metadata = step1_runs["run{}".format(run)][branch]
            short_label = str(metadata.get("short_label", short_label))
            hist_name = metadata.get("hist_name")
            hist = source.Get("run{}/{}".format(run, hist_name)) if hist_name else None
            if hist:
                source_hists[run] = hist
            channel_fits[run] = fit_summaries[run]["channels"][branch]
            channel_curves[run] = fit_curves[run][branch]

        built = make_overlay_canvas(
            branch, short_label, selection, resolved, source_hists,
            channel_fits, channel_curves, normalise, args.fits_only,
            args.zoom_keep,
        )
        if built is None:
            skipped.append(branch)
            continue
        canvas, canvas_json = built
        keys.append(canvas.GetName())
        jsons[canvas.GetName()] = canvas_json
        if index % 25 == 0 or index == len(branches):
            print("[BUILD] {}/{} channels".format(index, len(branches)))
    source.Close()

    outdir = args.outdir or args.step2_dir
    outdir.mkdir(parents=True, exist_ok=True)
    html_path = outdir / "pe_fit_position_overlay.html"
    run_text = ", ".join(
        "{:g} mm -> run {} ({:g} mm)".format(target, run, actual)
        for target, run, actual in resolved
    )
    intro = "\n".join([
        "**Runs (target X -> resolved run, actual X):** " + run_text + ".",
        "**Source:** Existing Step 1 histograms and Step 2 fit summaries/curves; "
        "no integrals or fits were recomputed.",
        "**Model:** Total multi-peak RooCrystalBall fit, reconstructed by summing "
        "the stored per-photoelectron component curves.",
        "**Selection:** {}.".format(selection),
        "Curves are {}. {}".format(
            "shown with raw event counts" if args.raw_yields else "unit-area normalized",
            "Only total fits are shown." if args.fits_only else (
                "Thin dotted lines are Step 1 data; thick solid lines are Step 2 fits."
            ),
        ),
    ])
    generate_jsroot_html(
        keys, jsons, plots_per_row=args.plots_per_row,
        output_html=str(html_path),
        title="CaloX PE-calibration: PE-fit position overlays",
        intro_text=intro,
    )

    output_summary = {
        "step1_root": str(step1_root_path),
        "step2_dir": str(args.step2_dir),
        "target_positions_mm": args.target_positions,
        "position_key": args.position_key,
        "resolved_runs": [
            {"target_position_mm": target, "run": run, "actual_position_mm": actual}
            for target, run, actual in resolved
        ],
        "selection": selection,
        "normalized": normalise,
        "fits_only": args.fits_only,
        "channels": branches,
        "skipped_channels": skipped,
    }
    summary_path = outdir / "pe_fit_position_overlay_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(output_summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print("[DONE] Interactive HTML: {}".format(html_path))
    print("[DONE] JSON summary: {}".format(summary_path))


if __name__ == "__main__":
    main()
