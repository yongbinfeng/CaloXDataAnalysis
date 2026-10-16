#!/usr/bin/env python3
"""Exponential fits to per-channel attenuation scans.

Two fit models are supported via --fit-mode:

  shared       ADC_ch(x) = A_ch * exp(k_family * x)
               One slope per channel family, fitted simultaneously across all
               channels of that family; each channel keeps its own normalization.

  independent  ADC_ch(x) = A_ch * exp(k_ch * x)
               Every channel gets its own slope *and* its own normalization.
               The shared-family curve is still overlaid (dashed) for reference.

Both are weighted linear least squares in log(ADC), with sigma_log = sigma_ADC / ADC.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from array import array
from typing import Any, Dict, List, Tuple

import ROOT

from channels.channel_map import build_drs_boards
from plotting.calox_plot_helper import BoardPlotHelper
from core.plot_manager import PlotManager
from utils.visualization import visualizeDRSBoards
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.utils import get_channel_var

from plot_velocityfits import draw_header, safe_name


@dataclass
class ChannelFit:
    family: str
    fit_mode: str
    n_fit_points: int
    amplitude_at_zero_adc: float
    amplitude_error_adc: float
    slope_per_mm: float
    slope_error_per_mm: float
    attenuation_length_mm: float
    attenuation_length_error_mm: float
    chi2: float
    ndf: int
    # Shared-family reference values, kept in both modes so the two models can be
    # compared channel by channel without re-running the shared fit.
    shared_slope_per_mm: float
    shared_slope_error_per_mm: float
    shared_amplitude_at_zero_adc: float
    shared_attenuation_length_mm: float


DEFAULT_INPUT = (
    "results/attenuation/runs_1897_1942_no_pedestal/attenuation_summary.json"
)
DEFAULT_OUTDIR = {
    "shared": "results/attenuation/runs_1897_1942_no_pedestal/exponential_fit",
    "independent": (
        "results/attenuation/runs_1897_1942_no_pedestal/"
        "exponential_fit_per_channel"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit exponentials to per-channel attenuation summary points."
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument(
        "--observable", choices=("auto", "adc", "pe"), default="auto",
        help=(
            "Quantity stored in the input points. 'auto' selects mean_pe when "
            "present, otherwise mean_adc (default: auto)."
        ),
    )
    parser.add_argument(
        "--fit-mode",
        choices=("shared", "independent"),
        default="shared",
        help=(
            "shared: one slope per channel family, fitted simultaneously. "
            "independent: every channel gets its own slope. Default: shared."
        ),
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default depends on --fit-mode.",
    )
    parser.add_argument(
        "--family-key",
        choices=("material", "material_size"),
        default="material",
        help=(
            "Grouping for the shared-slope fit. 'material' pools 3 mm and 6 mm "
            "together (original behaviour). 'material_size' fits each "
            "(material, fiber size) separately. Note there are no Sci 6 mm "
            "channels, so 'material_size' yields five groups, not six."
        ),
    )
    parser.add_argument(
        "--map-run", type=int, default=1897,
        help=(
            "Run whose channel map supplies the physical (x, y) tower positions "
            "and the DRS display range. Must be a tb2026 run (>= 1896) for the "
            "PHASE_2 layout, which extends to -18 cm in X. Default: 1897."
        ),
    )
    parser.add_argument(
        "--no-lambda-map", action="store_true",
        help="Skip the 2D lambda heat maps.",
    )
    parser.add_argument(
        "--map-lambda-unit", choices=("mm", "cm", "m"), default="mm",
        help="Unit for the attenuation length on the heat-map z axis. Default: mm.",
    )
    parser.add_argument(
        "--map-log-min", type=float, default=3.0,
        help=(
            "Bottom of the heat-map z scale, as log10(lambda / mm). Rescaled "
            "automatically if --map-lambda-unit is not mm."
        ),
    )
    parser.add_argument(
        "--map-log-max", type=float, default=6.0,
        help="Top of the heat-map z scale, as log10(lambda / mm).",
    )
    parser.add_argument(
        "--lambda-bins", type=int, default=40,
        help="Bins in the per-channel lambda histograms (log-spaced).",
    )
    parser.add_argument(
        "--lambda-log-min", type=float, default=3.0,
        help="log10(lambda/mm) at the low edge of the histogram axis.",
    )
    parser.add_argument(
        "--lambda-log-max", type=float, default=7.0,
        help="log10(lambda/mm) at the high edge of the histogram axis.",
    )
    parser.add_argument("--fit-x-max", type=float, default=785.0)
    parser.add_argument("--exclude-x-min", type=float, default=-175.0)
    parser.add_argument("--exclude-x-max", type=float, default=-85.0)
    parser.add_argument(
        "--min-fit-points",
        type=int,
        default=3,
        help=(
            "Minimum surviving scan points for an independent per-channel fit. "
            "Two points give a slope with zero degrees of freedom, so the "
            "default asks for three."
        ),
    )
    args = parser.parse_args()
    if args.outdir is None:
        args.outdir = DEFAULT_OUTDIR[args.fit_mode]
        if args.family_key == "material_size":
            args.outdir += "_by_size"
    return args


def configure_observable(args: argparse.Namespace, source: Dict[str, Any]) -> None:
    """Attach input-field and display metadata for ADC or PE summaries."""
    if args.observable == "auto":
        first_point = next(
            (
                point
                for entries in source.get("points", {}).values()
                for point in entries[:1]
            ),
            {},
        )
        args.observable = "pe" if "mean_pe" in first_point else "adc"
    if args.observable == "pe":
        args.value_key = "mean_pe"
        args.error_key = "mean_error_pe"
        args.observable_symbol = "PE"
        args.observable_axis = "Mean photoelectrons"
        args.amplitude_unit = "P.E."
    else:
        args.value_key = "mean_adc"
        args.error_key = "mean_error_adc"
        args.observable_symbol = "ADC"
        args.observable_axis = "Mean ADC integral [ADC]"
        args.amplitude_unit = "ADC"


def channel_fit_dict(result: ChannelFit, args: argparse.Namespace) -> Dict[str, Any]:
    """Use PE-specific amplitude field names in PE summaries."""
    output = asdict(result)
    if args.observable == "pe":
        for old, new in (
            ("amplitude_at_zero_adc", "amplitude_at_zero_pe"),
            ("amplitude_error_adc", "amplitude_error_pe"),
            ("shared_amplitude_at_zero_adc", "shared_amplitude_at_zero_pe"),
        ):
            output[new] = output.pop(old)
    return output


def channel_metadata(run: int = 1897) -> Dict[str, Dict[str, str]]:
    metadata: Dict[str, Dict[str, str]] = {}
    for _, board in build_drs_boards(run_number=run).items():
        for channel in board:
            if channel.is_reference:
                continue
            branch = channel.get_channel_name(blsub=False)
            family = str(get_channel_var(channel))
            size = "6 mm" if channel.is6mm else "3 mm"
            metadata[branch] = {
                "family": family,
                "material": family,
                "size": size,
                "label": "B{} G{} Ch{} | {}, {}".format(
                    board.board_no,
                    channel.group_no,
                    channel.channel_no,
                    family,
                    size,
                ),
            }
    return metadata


def is_fit_point(point: Dict[str, Any], args: argparse.Namespace) -> bool:
    x = float(point["position_mm"])
    y = float(point[args.value_key])
    error = float(point[args.error_key])
    return (
        x <= args.fit_x_max
        and not (args.exclude_x_min <= x <= args.exclude_x_max)
        and math.isfinite(y)
        and math.isfinite(error)
        and y > 0.0
        and error > 0.0
    )


def simultaneous_family_fit(
    family_channels: List[str],
    points: Dict[str, List[Dict[str, Any]]],
    args: argparse.Namespace,
) -> Tuple[float, float, Dict[str, Tuple[float, float]]]:
    """Fit log(y)=alpha_channel+k*x with one shared family slope."""
    prepared: Dict[str, List[Tuple[float, float, float]]] = {}
    sxx = 0.0
    sxy = 0.0
    for branch in family_channels:
        channel_points: List[Tuple[float, float, float]] = []
        for point in points.get(branch, []):
            if not is_fit_point(point, args):
                continue
            x = float(point["position_mm"])
            y = float(point[args.value_key])
            sigma_log = float(point[args.error_key]) / y
            channel_points.append((x, math.log(y), sigma_log))
        if len(channel_points) < 2:
            continue
        prepared[branch] = channel_points
        weights = [1.0 / sigma**2 for _, _, sigma in channel_points]
        sum_w = sum(weights)
        xbar = sum(w * x for w, (x, _, _) in zip(weights, channel_points)) / sum_w
        ybar = sum(w * y for w, (_, y, _) in zip(weights, channel_points)) / sum_w
        for weight, (x, y, _) in zip(weights, channel_points):
            sxx += weight * (x - xbar) ** 2
            sxy += weight * (x - xbar) * (y - ybar)
    if sxx <= 0.0:
        raise RuntimeError("Family fit has no usable position leverage")

    slope = sxy / sxx
    slope_error = 1.0 / math.sqrt(sxx)
    intercepts: Dict[str, Tuple[float, float]] = {}
    for branch, channel_points in prepared.items():
        weights = [1.0 / sigma**2 for _, _, sigma in channel_points]
        sum_w = sum(weights)
        xbar = sum(w * x for w, (x, _, _) in zip(weights, channel_points)) / sum_w
        ybar = sum(w * y for w, (_, y, _) in zip(weights, channel_points)) / sum_w
        alpha = ybar - slope * xbar
        alpha_error = math.sqrt(1.0 / sum_w + xbar**2 / sxx)
        intercepts[branch] = (alpha, alpha_error)
    return slope, slope_error, intercepts


def independent_channel_fit(
    channel_points: List[Dict[str, Any]],
    args: argparse.Namespace,
) -> Tuple[float, float, float, float, float, int]:
    """Fit log(y)=alpha+k*x for a single channel, with its own slope.

    Returns (slope, slope_error, alpha, alpha_error, chi2, ndf).
    Same weighted least-squares algebra as the shared fit, but with both
    parameters free for this channel alone.
    """
    prepared: List[Tuple[float, float, float]] = []
    for point in channel_points:
        if not is_fit_point(point, args):
            continue
        x = float(point["position_mm"])
        y = float(point[args.value_key])
        prepared.append((x, math.log(y), float(point[args.error_key]) / y))
    if len(prepared) < max(2, args.min_fit_points):
        raise ValueError("too few points")

    weights = [1.0 / sigma**2 for _, _, sigma in prepared]
    sum_w = sum(weights)
    sum_wx = sum(w * x for w, (x, _, _) in zip(weights, prepared))
    sum_wy = sum(w * y for w, (_, y, _) in zip(weights, prepared))
    sum_wxx = sum(w * x * x for w, (x, _, _) in zip(weights, prepared))
    sum_wxy = sum(w * x * y for w, (x, y, _) in zip(weights, prepared))
    denom = sum_w * sum_wxx - sum_wx**2
    if not math.isfinite(denom) or denom <= 0.0:
        raise ValueError("degenerate position leverage")

    slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
    alpha = (sum_wy - slope * sum_wx) / sum_w
    slope_error = math.sqrt(sum_w / denom)
    alpha_error = math.sqrt(sum_wxx / denom)
    if not math.isfinite(slope) or slope == 0.0:
        raise ValueError("non-finite slope")

    chi2 = sum(
        w * (y - alpha - slope * x) ** 2
        for w, (x, y, _) in zip(weights, prepared)
    )
    return slope, slope_error, alpha, alpha_error, chi2, max(0, len(prepared) - 2)


def make_canvas(
    branch: str,
    label: str,
    family: str,
    all_points: List[Dict[str, Any]],
    fit_result: ChannelFit,
    args: argparse.Namespace,
) -> Tuple[Any, str]:
    ordered = sorted(all_points, key=lambda point: (point["position_mm"], point["run"]))
    included = [point for point in ordered if is_fit_point(point, args)]
    excluded = [point for point in ordered if not is_fit_point(point, args)]

    graph_fit = ROOT.TGraphErrors(len(included))
    graph_fit.SetName("gr_fit_points_{}".format(safe_name(branch)))
    graph_fit.SetTitle("")
    for index, point in enumerate(included):
        graph_fit.SetPoint(index, point["position_mm"], point[args.value_key])
        graph_fit.SetPointError(index, 0.0, point[args.error_key])
    graph_fit.SetMarkerStyle(20)
    graph_fit.SetMarkerSize(1.0)
    graph_fit.SetMarkerColor(ROOT.kBlue + 1)
    graph_fit.SetLineColor(ROOT.kBlue + 1)

    graph_excluded = ROOT.TGraphErrors(len(excluded))
    graph_excluded.SetName("gr_excluded_points_{}".format(safe_name(branch)))
    for index, point in enumerate(excluded):
        graph_excluded.SetPoint(index, point["position_mm"], point[args.value_key])
        graph_excluded.SetPointError(index, 0.0, point[args.error_key])
    graph_excluded.SetMarkerStyle(24)
    graph_excluded.SetMarkerSize(1.0)
    graph_excluded.SetMarkerColor(ROOT.kGray + 2)
    graph_excluded.SetLineColor(ROOT.kGray + 2)

    xmin = min(float(point["position_mm"]) for point in ordered)
    xmax = max(float(point["position_mm"]) for point in ordered)
    fit_curve = ROOT.TF1(
        "exp_fit_{}".format(safe_name(branch)),
        "[0]*exp([1]*x)",
        xmin,
        args.fit_x_max,
    )
    fit_curve.SetParameters(
        fit_result.amplitude_at_zero_adc,
        fit_result.slope_per_mm,
    )
    fit_curve.SetLineColor(ROOT.kRed + 1)
    fit_curve.SetLineWidth(3)

    # In independent mode also draw the shared-family curve, so the per-channel
    # slope can be judged against the family average on the same axes.
    shared_curve = None
    if fit_result.fit_mode == "independent":
        shared_curve = ROOT.TF1(
            "exp_shared_{}".format(safe_name(branch)),
            "[0]*exp([1]*x)",
            xmin,
            args.fit_x_max,
        )
        shared_curve.SetParameters(
            fit_result.shared_amplitude_at_zero_adc,
            fit_result.shared_slope_per_mm,
        )
        shared_curve.SetLineColor(ROOT.kAzure + 2)
        shared_curve.SetLineStyle(2)
        shared_curve.SetLineWidth(2)

    canvas_name = "attenuation_expfit_{}".format(safe_name(branch))
    canvas = ROOT.TCanvas(canvas_name, "", 1100, 800)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    graph_fit.Draw("AP")
    graph_fit.GetXaxis().SetLimits(xmin - 20.0, xmax + 20.0)
    all_y = [float(point[args.value_key]) for point in ordered]
    graph_fit.SetMinimum(max(0.0, 0.80 * min(all_y)))
    graph_fit.SetMaximum(1.18 * max(all_y))
    graph_fit.GetXaxis().SetTitle("Table position X [mm]")
    graph_fit.GetYaxis().SetTitle(args.observable_axis)
    graph_fit.GetXaxis().SetTitleSize(0.045)
    graph_fit.GetYaxis().SetTitleSize(0.045)
    graph_fit.GetXaxis().SetLabelSize(0.04)
    graph_fit.GetYaxis().SetLabelSize(0.04)
    graph_fit.GetYaxis().SetTitleOffset(1.35)
    if excluded:
        graph_excluded.Draw("P SAME")
    if shared_curve is not None:
        shared_curve.Draw("L SAME")
    fit_curve.Draw("L SAME")

    info = ROOT.TPaveText(0.15, 0.62, 0.60, 0.82, "NDC")
    info.SetBorderSize(0)
    info.SetFillStyle(0)
    info.SetTextFont(42)
    info.SetTextSize(0.023)
    info.SetTextAlign(12)
    if fit_result.fit_mode == "independent":
        info.AddText("Independent per-channel fit ({} family)".format(family))
        info.AddText("{}(x) = A_{{ch}} exp(k_{{ch}} x)".format(
            args.observable_symbol
        ))
        slope_label = "k_{ch}"
    else:
        info.AddText("Simultaneous {} family fit".format(family))
        info.AddText("{}(x) = A_{{ch}} exp(k_{{family}} x)".format(
            args.observable_symbol
        ))
        slope_label = "k_{family}"
    info.AddText("A_{{ch}} = {:.2f} #pm {:.2f} {}".format(
        fit_result.amplitude_at_zero_adc, fit_result.amplitude_error_adc,
        args.amplitude_unit,
    ))
    info.AddText("{} = ({:.6g} #pm {:.2g}) mm^{{-1}}".format(
        slope_label, fit_result.slope_per_mm, fit_result.slope_error_per_mm
    ))
    info.AddText("#lambda = 1/|k| = {:.1f} #pm {:.1f} mm".format(
        fit_result.attenuation_length_mm,
        fit_result.attenuation_length_error_mm,
    ))
    if fit_result.ndf > 0:
        info.AddText("#chi^{{2}}/ndf = {:.2f} / {}".format(
            fit_result.chi2, fit_result.ndf
        ))
    info.Draw()

    legend = ROOT.TLegend(0.63, 0.64, 0.94, 0.82)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.021)
    legend.AddEntry(graph_fit, "Points used in fit", "pe")
    if excluded:
        legend.AddEntry(graph_excluded, "Excluded points", "pe")
    if fit_result.fit_mode == "independent":
        legend.AddEntry(fit_curve, "Per-channel exponential", "l")
        legend.AddEntry(
            shared_curve,
            "Shared {} slope".format(family),
            "l",
        )
    else:
        legend.AddEntry(fit_curve, "Shared-slope exponential", "l")
    legend.Draw()

    draw_header(canvas, "Attenuation exponential fit: {}".format(label), "MCP-clean events")
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph_fit, graph_excluded, fit_curve, shared_curve, info, legend]
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def make_lambda_hist_canvas(
    family: str,
    results: List[ChannelFit],
    args: argparse.Namespace,
) -> Tuple[Any, str, Dict[str, Any]]:
    """1D distribution of per-channel attenuation lengths for one family.

    lambda = 1/|k| spans four decades across channels, so the axis is
    log-spaced. The sign of k is not in lambda, and it matters physically -
    k < 0 means light *grows* with X rather than attenuating - so those
    channels are drawn as a separate overlay instead of being silently mixed in.
    """
    edges = array(
        "d",
        [
            10.0 ** (args.lambda_log_min
                     + index * (args.lambda_log_max - args.lambda_log_min)
                     / args.lambda_bins)
            for index in range(args.lambda_bins + 1)
        ],
    )
    stem = safe_name(family)
    hist_all = ROOT.TH1D("h_lambda_{}".format(stem), "", args.lambda_bins, edges)
    hist_neg = ROOT.TH1D("h_lambda_neg_{}".format(stem), "", args.lambda_bins, edges)
    for result in results:
        hist_all.Fill(result.attenuation_length_mm)
        if result.slope_per_mm < 0.0:
            hist_neg.Fill(result.attenuation_length_mm)

    lengths = sorted(result.attenuation_length_mm for result in results)
    n = len(lengths)
    median = lengths[n // 2]
    q16, q84 = lengths[int(0.16 * n)], lengths[min(n - 1, int(0.84 * n))]
    n_neg = sum(1 for result in results if result.slope_per_mm < 0.0)
    shared_lambda = results[0].shared_attenuation_length_mm

    hist_all.SetLineColor(ROOT.kBlue + 1)
    hist_all.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)
    hist_all.SetLineWidth(2)
    hist_all.SetStats(0)
    hist_all.GetXaxis().SetTitle("Per-channel #lambda = 1/|k| [mm]")
    hist_all.GetYaxis().SetTitle("Channels")
    hist_all.GetXaxis().SetTitleSize(0.045)
    hist_all.GetYaxis().SetTitleSize(0.045)
    hist_all.GetXaxis().SetLabelSize(0.04)
    hist_all.GetYaxis().SetLabelSize(0.04)
    hist_all.GetYaxis().SetTitleOffset(1.30)
    hist_all.SetMinimum(0.0)
    hist_all.SetMaximum(max(1.0, 1.45 * hist_all.GetMaximum()))
    hist_neg.SetLineColor(ROOT.kRed + 1)
    hist_neg.SetFillColorAlpha(ROOT.kRed + 1, 0.30)
    hist_neg.SetFillStyle(3354)
    hist_neg.SetLineWidth(2)
    hist_neg.SetStats(0)

    canvas = ROOT.TCanvas("lambda_hist_{}".format(stem), "", 1000, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.13)
    canvas.SetTicks(1, 1)
    canvas.SetLogx(1)
    hist_all.Draw("HIST")
    if n_neg:
        hist_neg.Draw("HIST SAME")

    keep: List[Any] = [hist_all, hist_neg]
    ymax = hist_all.GetMaximum()
    marker_lines: List[Tuple[Any, str]] = []
    for value, colour, style, label in (
        (median, ROOT.kBlack, 2, "Median of channels"),
        (shared_lambda, ROOT.kGreen + 3, 7, "Simultaneous-fit #lambda"),
    ):
        if not (edges[0] < value < edges[-1]):
            continue
        line = ROOT.TLine(value, 0.0, value, 0.80 * ymax)
        line.SetLineColor(colour)
        line.SetLineStyle(style)
        line.SetLineWidth(3)
        line.Draw()
        keep.append(line)
        marker_lines.append((line, label))

    stats = ROOT.TPaveText(0.55, 0.58, 0.94, 0.82, "NDC")
    stats.SetBorderSize(0)
    stats.SetFillStyle(0)
    stats.SetTextFont(42)
    stats.SetTextSize(0.024)
    stats.SetTextAlign(12)
    stats.AddText("Channels: {}".format(n))
    stats.AddText("Median #lambda = {:.0f} mm".format(median))
    stats.AddText("16-84%: {:.0f} - {:.0f} mm".format(q16, q84))
    stats.AddText("Simultaneous fit #lambda = {:.0f} mm".format(shared_lambda))
    stats.AddText("k < 0 (light grows with X): {}".format(n_neg))
    if hist_all.GetBinContent(0) or hist_all.GetBinContent(args.lambda_bins + 1):
        stats.AddText("outside axis: {:.0f}".format(
            hist_all.GetBinContent(0) + hist_all.GetBinContent(args.lambda_bins + 1)
        ))
    stats.Draw()
    keep.append(stats)

    n_entries = 1 + (1 if n_neg else 0) + len(marker_lines)
    legend = ROOT.TLegend(0.55, 0.555 - 0.042 * n_entries, 0.94, 0.555)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.023)
    legend.AddEntry(hist_all, "All channels", "lf")
    if n_neg:
        legend.AddEntry(hist_neg, "k < 0 subset", "lf")
    for line, label in marker_lines:
        legend.AddEntry(line, label, "l")
    legend.Draw()
    keep.append(legend)

    draw_header(
        canvas,
        "Per-channel attenuation length: {}".format(family),
        "MCP-clean events",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keep
    summary = {
        "n_channels": n,
        "median_lambda_mm": median,
        "q16_lambda_mm": q16,
        "q84_lambda_mm": q84,
        "min_lambda_mm": lengths[0],
        "max_lambda_mm": lengths[-1],
        "n_slope_negative": n_neg,
        "simultaneous_fit_lambda_mm": shared_lambda,
    }
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data()), summary


def make_lambda_maps(
    channel_results: Dict[str, ChannelFit],
    args: argparse.Namespace,
    outdir: Path,
) -> Dict[str, Any]:
    """2D maps of per-channel lambda over the physical calorimeter face.

    Positions come from the tb2026 channel map (i_tower_x/i_tower_y scaled to
    cm by utils.visualization), and BoardPlotHelper picks the run-dependent
    display range - PHASE_2 runs (>= 1896) extend to -18 cm in X.

    z is log10(lambda/mm) because lambda spans four decades; a linear scale
    would show one saturated cell and nothing else. Cerenkov and scintillating
    fibers share tower positions, so they are drawn as separate maps rather
    than overwriting each other - that is what the repo's Cer/Sci pair does.
    Both fiber sizes appear on each map at their true positions.
    """
    drsboards = build_drs_boards(run_number=args.map_run)
    # lambda is fitted in mm; convert before taking the log so the z axis and
    # its label agree.
    per_mm = {"mm": 1.0, "cm": 0.1, "m": 0.001}[args.map_lambda_unit]
    log_shift = math.log10(per_mm)
    zmin = args.map_log_min + log_shift
    zmax = args.map_log_max + log_shift
    zlabel = "log_{{10}}( attenuation length #lambda / 1 {} )".format(
        args.map_lambda_unit
    )
    value_map = {
        branch: math.log10(result.attenuation_length_mm * per_mm)
        for branch, result in channel_results.items()
        if result.attenuation_length_mm > 0.0
    }
    n_negative = sum(
        1 for branch, result in channel_results.items()
        if branch in value_map and result.slope_per_mm < 0.0
    )

    pm = PlotManager(
        rootdir="results/root/Run{}".format(args.map_run),
        plotdir=str(outdir / "attenuation_length_maps"),
        htmldir=str(outdir / "attenuation_length_maps_html"),
        use_jsroot=False,
    )
    pm.set_output_dir("")
    helper = BoardPlotHelper(pm, run_number=args.map_run)

    # quartzOnly: 0 = all Cerenkov, 1 = quartz only, 2 = plastic only.
    produced: List[str] = []
    for tag, quartz_only in (("Quartz", 1), ("Plastic", 2)):
        cer_hists, sci_hists = visualizeDRSBoards(
            drsboards, valuemaps=value_map,
            suffix="lambda_{}".format(tag), quartzOnly=quartz_only,
        )
        helper.plot_board_map(
            cer_hists, "AttenuationLength_log10_Cer{}".format(tag),
            extra_text="Cer{}".format(tag),
            zmin=zmin, zmax=zmax, nTextDigits=2, zlabel=zlabel,
        )
        produced.append("AttenuationLength_log10_Cer{}".format(tag))
    cer_hists, sci_hists = visualizeDRSBoards(
        drsboards, valuemaps=value_map, suffix="lambda_Sci",
    )
    helper.plot_board_map(
        sci_hists, "AttenuationLength_log10_Sci", extra_text="Sci",
        zmin=zmin, zmax=zmax, nTextDigits=2, zlabel=zlabel,
    )
    produced.append("AttenuationLength_log10_Sci")

    print("[MAP] Wrote {} attenuation-length heat map(s) to {} (map run {}, "
          "z = log10(lambda / 1 {}) in [{:.2f}, {:.2f}]).".format(
              len(produced), outdir / "attenuation_length_maps", args.map_run,
              args.map_lambda_unit, zmin, zmax))
    if n_negative:
        print("[MAP] NOTE: {} of {} mapped channels have k < 0; lambda = 1/|k| "
              "hides that sign on these maps.".format(n_negative, len(value_map)))
    return {
        "map_run": args.map_run,
        "z_variable": "log10(attenuation length / 1 {})".format(args.map_lambda_unit),
        "lambda_unit": args.map_lambda_unit,
        "z_range": [zmin, zmax],
        "n_channels_mapped": len(value_map),
        "n_slope_negative": n_negative,
        "plots": produced,
        "plot_dir": str(outdir / "attenuation_length_maps"),
    }


def main() -> None:
    args = parse_args()
    setup_root(batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    with Path(args.input).open() as handle:
        source = json.load(handle)
    configure_observable(args, source)
    points: Dict[str, List[Dict[str, Any]]] = source["points"]
    metadata = channel_metadata(run=args.map_run)

    # 3 mm and 6 mm fibers of one material are different light guides, so
    # --family-key material_size gives each its own shared slope.
    for branch, meta in metadata.items():
        meta["family"] = (
            "{} {}".format(meta["material"], meta["size"].replace(" ", ""))
            if args.family_key == "material_size" else meta["material"]
        )

    families: Dict[str, List[str]] = {}
    for branch in points:
        if branch in metadata:
            families.setdefault(metadata[branch]["family"], []).append(branch)
    print("Family grouping: {} -> {} simultaneous fit(s): {}".format(
        args.family_key, len(families),
        ", ".join("{} ({} ch)".format(k, len(v)) for k, v in sorted(families.items())),
    ))

    family_results: Dict[str, Dict[str, Any]] = {}
    channel_results: Dict[str, ChannelFit] = {}
    skipped: Dict[str, str] = {}
    for family, branches in sorted(families.items()):
        # The shared-family fit is always run: in independent mode it supplies the
        # reference curve and the family-average numbers written to the summary.
        slope, slope_error, intercepts = simultaneous_family_fit(
            branches, points, args
        )
        attenuation_length = 1.0 / abs(slope)
        attenuation_length_error = slope_error / (slope * slope)
        family_results[family] = {
            "channels_in_simultaneous_fit": len(intercepts),
            "shared_slope_per_mm": slope,
            "shared_slope_error_per_mm": slope_error,
            "attenuation_length_mm": attenuation_length,
            "attenuation_length_error_mm": attenuation_length_error,
        }
        for branch in branches:
            if branch not in intercepts:
                n_usable = sum(
                    1 for point in points.get(branch, [])
                    if is_fit_point(point, args)
                )
                skipped[branch] = (
                    "too few positive finite fit points for the logarithmic fit "
                    "({} usable)".format(n_usable)
                )
        for branch, (alpha, alpha_error) in intercepts.items():
            shared_amplitude = math.exp(alpha)
            fit_points = [point for point in points[branch] if is_fit_point(point, args)]
            shared_chi2 = sum(
                (
                    (math.log(float(point[args.value_key])) - alpha
                     - slope * float(point["position_mm"]))
                    / (float(point[args.error_key]) / float(point[args.value_key]))
                ) ** 2
                for point in fit_points
            )

            if args.fit_mode == "independent":
                try:
                    (ch_slope, ch_slope_error, ch_alpha, ch_alpha_error,
                     ch_chi2, ch_ndf) = independent_channel_fit(points[branch], args)
                except ValueError as error:
                    skipped[branch] = str(error)
                    continue
                amplitude = math.exp(ch_alpha)
                fit_slope, fit_slope_error = ch_slope, ch_slope_error
                amplitude_error = amplitude * ch_alpha_error
                chi2, ndf = ch_chi2, ch_ndf
                length = 1.0 / abs(ch_slope)
                length_error = ch_slope_error / (ch_slope * ch_slope)
            else:
                amplitude = shared_amplitude
                amplitude_error = amplitude * alpha_error
                fit_slope, fit_slope_error = slope, slope_error
                chi2, ndf = shared_chi2, max(0, len(fit_points) - 1)
                length = attenuation_length
                length_error = attenuation_length_error

            channel_results[branch] = ChannelFit(
                family=family,
                fit_mode=args.fit_mode,
                n_fit_points=len(fit_points),
                amplitude_at_zero_adc=amplitude,
                amplitude_error_adc=amplitude_error,
                slope_per_mm=fit_slope,
                slope_error_per_mm=fit_slope_error,
                attenuation_length_mm=length,
                attenuation_length_error_mm=length_error,
                chi2=chi2,
                ndf=ndf,
                shared_slope_per_mm=slope,
                shared_slope_error_per_mm=slope_error,
                shared_amplitude_at_zero_adc=shared_amplitude,
                shared_attenuation_length_mm=attenuation_length,
            )

    if args.fit_mode == "independent":
        # Spread of the per-channel slopes is the headline number this mode adds:
        # it says whether one family slope was ever a fair description.
        for family in family_results:
            slopes = [
                result.slope_per_mm
                for result in channel_results.values()
                if result.family == family
            ]
            lengths = [
                result.attenuation_length_mm
                for result in channel_results.values()
                if result.family == family
            ]
            if not slopes:
                continue
            mean_slope = sum(slopes) / len(slopes)
            rms_slope = math.sqrt(
                sum((value - mean_slope) ** 2 for value in slopes) / len(slopes)
            )
            family_results[family].update({
                "channels_with_independent_fit": len(slopes),
                "independent_slope_mean_per_mm": mean_slope,
                "independent_slope_rms_per_mm": rms_slope,
                "independent_attenuation_length_median_mm": sorted(lengths)[
                    len(lengths) // 2
                ],
                "independent_attenuation_length_min_mm": min(lengths),
                "independent_attenuation_length_max_mm": max(lengths),
            })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    root_path = outdir / "attenuation_exponential_fits.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    for branch in sorted(channel_results):
        canvas, canvas_json = make_canvas(
            branch,
            metadata[branch]["label"],
            metadata[branch]["family"],
            points[branch],
            channel_results[branch],
            args,
        )
        name = canvas.GetName()
        root_file.cd()
        canvas.Write(name)
        keys.append(name)
        jsons[name] = canvas_json
    root_file.Close()

    # Per-family lambda distributions. Only meaningful in independent mode: in
    # shared mode every channel of a family carries the same lambda by
    # construction, so the histogram would be a delta function.
    lambda_summary: Dict[str, Dict[str, Any]] = {}
    lambda_keys: List[str] = []
    map_summary: Dict[str, Any] = {}
    if args.fit_mode == "independent" and channel_results:
        by_family: Dict[str, List[ChannelFit]] = {}
        for result in channel_results.values():
            by_family.setdefault(result.family, []).append(result)
        lambda_root = ROOT.TFile(
            str(outdir / "attenuation_lambda_hists.root"), "RECREATE"
        )
        for family, results in sorted(by_family.items()):
            canvas, canvas_json, stats = make_lambda_hist_canvas(family, results, args)
            name = canvas.GetName()
            lambda_root.cd()
            canvas.Write(name)
            for obj in canvas._keepalive:
                if isinstance(obj, ROOT.TH1):
                    obj.Write()
            jsons[name] = canvas_json
            lambda_keys.append(name)
            lambda_summary[family] = stats
            print(
                "[LAMBDA] {:<16} n={:3d}  median={:8.0f} mm  16-84%={:8.0f}-{:<9.0f} "
                "simultaneous={:8.0f}  k<0: {}".format(
                    family, stats["n_channels"], stats["median_lambda_mm"],
                    stats["q16_lambda_mm"], stats["q84_lambda_mm"],
                    stats["simultaneous_fit_lambda_mm"], stats["n_slope_negative"],
                )
            )
        lambda_root.Close()
        # Summary histograms first, then the per-channel fits.
        keys = lambda_keys + keys

        if not args.no_lambda_map:
            map_summary = make_lambda_maps(channel_results, args, outdir)

    if args.fit_mode == "independent":
        title = "CaloX per-channel {} attenuation exponential fits".format(
            args.observable_symbol
        )
        model_text = (
            "**Fit model:** {symbol}_ch(x) = A_ch exp(k_ch x) - every channel gets "
            "its own slope *and* its own normalization, fitted in log({symbol}) "
            "with sigma_log = sigma_{symbol} / {symbol}.\n"
            "* The dashed blue curve is the shared-family fit, drawn for comparison.\n"
            "* The first {} plot(s) are the per-family distributions of the fitted "
            "lambda = 1/|k|, on a log axis. Black dashed = median, green = the "
            "simultaneous-fit lambda, red hatched = the k < 0 subset (light grows "
            "with X rather than attenuating).\n"
            "**Minimum points per channel:** {}.\n"
        ).format(
            len(lambda_keys), args.min_fit_points,
            symbol=args.observable_symbol,
        )
        model_string = "{}_ch(x) = A_ch * exp(k_ch * x)".format(
            args.observable_symbol
        )
    else:
        title = "CaloX simultaneous-family {} attenuation exponential fits".format(
            args.observable_symbol
        )
        model_text = (
            "**Fit model:** {symbol}_ch(x) = A_ch exp(k_family x), fitted "
            "simultaneously in log({symbol}), with one shared slope per channel "
            "family and an independent normalization per channel.\n"
        ).format(symbol=args.observable_symbol)
        model_string = "{}_ch(x) = A_ch * exp(k_family * x)".format(
            args.observable_symbol
        )

    html_path = outdir / "attenuation_exponential_fits.html"
    generate_jsroot_html(
        keys,
        jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title=title,
        intro_text=(
            model_text
            + "**Fit points:** X <= {:.0f} mm, excluding {:.0f} <= X <= {:.0f} mm.\n"
            "* Filled blue markers are used; open gray markers are excluded."
        ).format(args.fit_x_max, args.exclude_x_min, args.exclude_x_max),
    )
    output = {
        "input": args.input,
        "observable": args.observable,
        "value_key": args.value_key,
        "error_key": args.error_key,
        "fit_mode": args.fit_mode,
        "family_key": args.family_key,
        "fit_model": model_string,
        "fit_x_max_mm": args.fit_x_max,
        "excluded_x_range_mm": [args.exclude_x_min, args.exclude_x_max],
        "families": family_results,
        "channels": {
            branch: channel_fit_dict(result, args)
            for branch, result in channel_results.items()
        },
        "skipped_channels": skipped,
        "lambda_distributions": lambda_summary,
        "lambda_maps": map_summary,
    }
    json_path = outdir / "attenuation_exponential_fit_summary.json"
    with json_path.open("w") as handle:
        json.dump(output, handle, indent=2, sort_keys=True)
    print("Fit mode: {}".format(args.fit_mode))
    print("Fitted {} channels; skipped {}.".format(len(channel_results), len(skipped)))
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
