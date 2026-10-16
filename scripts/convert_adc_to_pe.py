#!/usr/bin/env python3
"""Step 3 of the ADC -> photoelectron calibration study.

For each channel, uses Step 2's fitted Crystal Ball components to find the
ADC-integral value where each pair of consecutive p.e. peaks crosses (i.e.
where an event is equally likely to belong to the lower or the higher p.e.
hypothesis). Those crossings become the bin edges of an ADC -> PE lookup:
everything between crossing i-1 and crossing i is assigned i photoelectrons.

Produces, per channel:
  * a plot of the Step 1 histogram with the Step 2 fit components overlaid
    and a vertical line at each crossing, labelled by p.e. count;
  * a calibration graph - PE (y) vs. the fitted ADC-integral peak position
    (x) - with a linear fit, i.e. the ADC -> PE conversion itself.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from utils.html_generator import generate_html
from utils.root_setup import setup_root
from utils.timing import auto_timer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Segment each channel's ADC-integral spectrum into photoelectron "
            "bins from Step 2's fitted peaks, and fit the ADC -> PE "
            "calibration (Step 3)."
        )
    )
    parser.add_argument(
        "--step1-dir", default="results/PhotoElectrons/Step1_IntegralCalculation"
    )
    parser.add_argument(
        "--step2-dir", default="results/PhotoElectrons/Step2_PEFits"
    )
    parser.add_argument(
        "--runs", type=int, nargs="+", default=None,
        help="Multiple runs to convert in one process. Reads each run's Step 2 "
        "output from --step2-dir/run<N>/ and writes to --outdir/run<N>/. "
        "Default: single run, reading --step2-dir flat (legacy behaviour).",
    )
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument(
        "--max-pe-for-calibration", type=int, default=None,
        help="Only use peaks 0..N (inclusive) in the linear calibration fit. "
        "Default: all fitted peaks.",
    )
    parser.add_argument(
        "--outdir", default="results/PhotoElectrons/Step3_ADCtoPEConversion"
    )
    return parser.parse_args()


def find_crossing(
    curve_i: List[List[float]], curve_j: List[List[float]], lo: float, hi: float
) -> float:
    """First x in [lo, hi] where curve_i and curve_j (as sampled point lists)
    swap which one is larger - the decision boundary between the two p.e.
    hypotheses. Falls back to the midpoint if the curves do not cross (too
    little overlap, or one component's fitted fraction collapsed to zero).
    """
    xi = [p[0] for p in curve_i]
    yi = [p[1] for p in curve_i]
    xj = [p[0] for p in curve_j]
    yj = [p[1] for p in curve_j]

    lo_use = max(lo, min(xi), min(xj))
    hi_use = min(hi, max(xi), max(xj))
    if hi_use <= lo_use:
        return 0.5 * (lo + hi)

    def interp(xs: List[float], ys: List[float], x: float) -> float:
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        lo_idx = 0
        hi_idx = len(xs) - 1
        while hi_idx - lo_idx > 1:
            mid = (lo_idx + hi_idx) // 2
            if xs[mid] <= x:
                lo_idx = mid
            else:
                hi_idx = mid
        x0, x1 = xs[lo_idx], xs[hi_idx]
        y0, y1 = ys[lo_idx], ys[hi_idx]
        if x1 == x0:
            return y0
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    n_steps = 2000
    step = (hi_use - lo_use) / n_steps
    prev_x = lo_use
    prev_diff = interp(xi, yi, prev_x) - interp(xj, yj, prev_x)
    for k in range(1, n_steps + 1):
        cur_x = lo_use + k * step
        cur_diff = interp(xi, yi, cur_x) - interp(xj, yj, cur_x)
        if prev_diff == 0.0:
            return prev_x
        if (prev_diff > 0) != (cur_diff > 0):
            frac = prev_diff / (prev_diff - cur_diff)
            return prev_x + frac * (cur_x - prev_x)
        prev_x, prev_diff = cur_x, cur_diff
    return 0.5 * (lo + hi)


def make_segmented_canvas(
    branch: str, run: int, hist: Any, fit: Dict[str, Any],
    curves: Dict[str, Any], edges: List[float],
) -> Any:
    canvas = ROOT.TCanvas("pe_segments_{}".format(branch), "", 1000, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetLogy(1)

    plot_lo, plot_hi = fit["plot_range"]
    hist.SetTitle("")
    hist.SetStats(0)
    hist.GetXaxis().SetTitle("ADC integral [ADC]")
    hist.GetYaxis().SetTitle("Events")
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetRangeUser(plot_lo, plot_hi)
    hist.SetLineColor(ROOT.kBlack)
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(0.6)
    hist.SetMinimum(0.5)
    hist.SetMaximum(max(10.0, 40.0 * hist.GetMaximum()))
    hist.Draw("HIST P")

    keepalive = [hist]
    colors = [ROOT.kGreen + 2] + [(i % 9) + 2 for i in range(1, 20)]
    for i_str, points in sorted(curves.items(), key=lambda kv: int(kv[0])):
        i = int(i_str)
        graph = ROOT.TGraph(len(points))
        for k, (x, y) in enumerate(points):
            graph.SetPoint(k, x, y)
        graph.SetLineColor(colors[i] if i < len(colors) else ROOT.kGray + 1)
        graph.SetLineStyle(ROOT.kDotted if i == 0 else ROOT.kDashed)
        graph.SetLineWidth(2)
        graph.Draw("L SAME")
        keepalive.append(graph)

    ymin, ymax = hist.GetMinimum(), hist.GetMaximum()
    lines = []
    for edge in edges:
        line = ROOT.TLine(edge, ymin, edge, ymax)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(ROOT.kSolid)
        line.SetLineWidth(1)
        line.Draw()
        lines.append(line)
    keepalive += lines

    labels = []
    bin_edges_full = [plot_lo] + list(edges) + [plot_hi]
    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextSize(0.028)
    latex.SetTextAlign(21)
    npe_max = fit["npe_max"]
    for i in range(min(npe_max, len(bin_edges_full) - 1)):
        lo_edge, hi_edge = bin_edges_full[i], bin_edges_full[i + 1]
        if lo_edge >= plot_hi or hi_edge <= plot_lo:
            continue
        center = 0.5 * (max(lo_edge, plot_lo) + min(hi_edge, plot_hi))
        latex.DrawLatex(center, ymax * 0.6, "{} p.e.".format(i))
    keepalive.append(latex)

    from plot_velocityfits import draw_header
    draw_header(
        canvas, "Run {} ADC integral -> PE segmentation: {}".format(
            run, fit.get("short_label", branch)),
        "Crystal Ball fit boundaries",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas


def compute_data_peaks(
    hist: Any, plot_lo: float, plot_hi: float, edges: List[float],
    npe_max: int, max_pe: int,
) -> Dict[int, Tuple[float, float]]:
    """The actual data-driven ADC-integral peak position (and its error)
    within each p.e. segment, i.e. where the real histogram content between
    two consecutive segmentation edges peaks - as opposed to `peak["mu"]`,
    which is `mu_ped + i*dpe` by construction of the fit model and so is
    exactly linear in i regardless of what the data actually look like.

    The peak bin is the maximum of a 3-bin moving average inside the segment
    (so a single fluctuating bin cannot win); the position is the
    content-weighted centroid of that bin and its two neighbours, and the
    error is that centroid's standard error. This is the peak of the
    distribution, not its mean, so the long tail of a segment does not pull it.
    """
    bin_edges_full = [plot_lo] + list(edges) + [plot_hi]
    axis = hist.GetXaxis()
    nbins = hist.GetNbinsX()
    peaks: Dict[int, Tuple[float, float]] = {}
    for i in range(min(npe_max, max_pe + 1, len(bin_edges_full) - 1)):
        first = max(1, axis.FindBin(bin_edges_full[i]))
        last = min(nbins, axis.FindBin(bin_edges_full[i + 1]))
        if last < first or hist.Integral(first, last) <= 0:
            continue
        best_bin, best_val = first, -1.0
        for b in range(first, last + 1):
            lo_b, hi_b = max(1, b - 1), min(nbins, b + 1)
            smooth = sum(hist.GetBinContent(k) for k in range(lo_b, hi_b + 1)) / (hi_b - lo_b + 1)
            if smooth > best_val:
                best_bin, best_val = b, smooth
        lo_b, hi_b = max(1, best_bin - 1), min(nbins, best_bin + 1)
        weights = [(axis.GetBinCenter(k), hist.GetBinContent(k)) for k in range(lo_b, hi_b + 1)]
        total = sum(w for _, w in weights)
        if total <= 0:
            continue
        centroid = sum(x * w for x, w in weights) / total
        variance = sum(w * (x - centroid) ** 2 for x, w in weights) / total
        error = max(math.sqrt(variance / total), 0.25 * axis.GetBinWidth(best_bin))
        peaks[i] = (centroid, error)
    return peaks


def make_calibration_canvas(
    branch: str, run: int, fit: Dict[str, Any], max_pe: int,
    data_peaks: Dict[int, Tuple[float, float]],
) -> Tuple[Any, Dict[str, float]]:
    peaks = [p for p in fit["peaks"] if p["i"] <= max_pe]
    graph = ROOT.TGraphErrors(len(peaks))
    for k, peak in enumerate(peaks):
        graph.SetPoint(k, peak["mu"], peak["i"])
        # A y-error of exactly zero makes TGraphErrors::Fit's weights blow up
        # (1/0^2) and silently return a degenerate zero-parameter "fit"; the
        # PE label really is exact, so a tiny placeholder just avoids that
        # edge case while the real uncertainty (on the ADC peak position)
        # still dominates via the x-error effective-variance weighting.
        graph.SetPointError(k, peak["mu_err"], 1.0e-3)

    canvas = ROOT.TCanvas("pe_calibration_{}".format(branch), "", 850, 650)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)

    graph.SetTitle("")
    graph.GetXaxis().SetTitle("ADC integral [ADC]")
    graph.GetYaxis().SetTitle("Photoelectrons")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.Draw("AP")

    # The fit peak positions (mu_i = mu_ped + i*dpe) are exactly linear in i
    # by construction of the Step 2 model - they cannot show curvature even
    # if the real data does. Overlay the actual data peak position within each
    # segmentation bin (Step 3's own crossings) as an independent check.
    data_graph = None
    if data_peaks:
        data_graph = ROOT.TGraphErrors(len(data_peaks))
        for k, (pe, (peak_val, peak_err)) in enumerate(sorted(data_peaks.items())):
            data_graph.SetPoint(k, peak_val, pe)
            data_graph.SetPointError(k, peak_err, 1.0e-3)
        data_graph.SetMarkerStyle(24)
        data_graph.SetMarkerColor(ROOT.kBlack)
        data_graph.SetLineColor(ROOT.kBlack)
        data_graph.Draw("P SAME")

    mus = [peak["mu"] for peak in peaks]
    pad = 0.1 * max(1.0, max(mus) - min(mus))
    fit_func = ROOT.TF1(
        "calib_fit_{}".format(branch), "pol1", min(mus) - pad, max(mus) + pad
    )
    graph.Fit(fit_func, "Q0S")
    fit_func.SetLineColor(ROOT.kRed + 1)
    fit_func.SetLineWidth(2)
    fit_func.Draw("SAME")

    intercept, slope = fit_func.GetParameter(0), fit_func.GetParameter(1)
    intercept_err, slope_err = fit_func.GetParError(0), fit_func.GetParError(1)
    chi2, ndf = fit_func.GetChisquare(), fit_func.GetNDF()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(0.17, 0.85, "PE = {:.5f}(#pm{:.5f}) #times ADC + {:.3f}(#pm{:.3f})".format(
        slope, slope_err, intercept, intercept_err))
    if ndf > 0:
        latex.DrawLatexNDC(0.17, 0.80, "#chi^{{2}}/ndf = {:.2f}/{}".format(chi2, ndf))

    legend = ROOT.TLegend(0.60, 0.15, 0.94, 0.30)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.026)
    legend.AddEntry(graph, "Fit peak position (Step 2)", "p")
    if data_graph is not None:
        legend.AddEntry(data_graph, "Data peak in segment (Step 3)", "p")
    legend.Draw()

    from plot_velocityfits import draw_header
    draw_header(
        canvas, "Run {} ADC -> PE calibration: {}".format(
            run, fit.get("short_label", branch)),
        "Linear fit to Step 2 peak positions",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, fit_func, latex, legend] + ([data_graph] if data_graph else [])

    return canvas, {
        "slope_pe_per_adc": slope, "slope_err": slope_err,
        "intercept_pe": intercept, "intercept_err": intercept_err,
        "chi2": chi2, "ndf": ndf,
        "n_points": len(peaks),
    }


def make_calibration_divadc_canvas(
    branch: str, run: int, fit: Dict[str, Any], max_pe: int,
    slope: float, intercept: float,
    data_peaks: Dict[int, Tuple[float, float]],
) -> Any:
    """The calibration line divided by ADC: (slope*ADC + intercept)/ADC =
    slope + intercept/ADC. A perfectly proportional (zero-pedestal-offset)
    calibration would be flat at `slope`; the intercept/ADC term is how far
    the fit is from that, and it is largest at low ADC.
    """
    peaks = [p for p in fit["peaks"] if p["i"] <= max_pe and p["mu"] > 0]

    graph = ROOT.TGraph(len(peaks))
    for k, peak in enumerate(peaks):
        graph.SetPoint(k, peak["mu"], peak["i"] / peak["mu"])

    canvas = ROOT.TCanvas("pe_calibration_divadc_{}".format(branch), "", 850, 650)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)

    mus = [peak["mu"] for peak in peaks]
    pad = 0.1 * max(1.0, max(mus) - min(mus))
    xlo = min(mus) - pad
    if xlo <= 0.0:
        # Guard against the 1/x singularity: never let the divided curve's
        # range reach x = 0.
        xlo = max(1.0, 0.5 * min(mus))
    xhi = max(mus) + pad

    graph.SetTitle("")
    graph.GetXaxis().SetTitle("ADC integral [ADC]")
    graph.GetYaxis().SetTitle("Fitted PE / ADC integral")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.GetXaxis().SetLimits(xlo, xhi)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.Draw("AP")

    # Same actual-data check as make_calibration_canvas: the fit points
    # (peak["mu"]) are exactly linear in i by construction, so divide the
    # real per-segment ADC peak position instead to see how the data itself compares.
    data_graph = None
    data_items = [(pe, peak) for pe, (peak, _err) in data_peaks.items()
                  if pe <= max_pe and peak > 0]
    if data_items:
        data_graph = ROOT.TGraph(len(data_items))
        for k, (pe, peak) in enumerate(sorted(data_items)):
            data_graph.SetPoint(k, peak, pe / peak)
        data_graph.SetMarkerStyle(24)
        data_graph.SetMarkerColor(ROOT.kBlack)
        data_graph.SetLineColor(ROOT.kBlack)
        data_graph.Draw("P SAME")

    ratio_func = ROOT.TF1(
        "calib_fit_divadc_{}".format(branch), "[0] + [1]/x", xlo, xhi
    )
    ratio_func.SetParameters(slope, intercept)
    ratio_func.SetLineColor(ROOT.kRed + 1)
    ratio_func.SetLineWidth(2)
    ratio_func.Draw("SAME")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.03)
    latex.DrawLatexNDC(
        0.17, 0.85,
        "(PE / ADC) = {:.3g} + {:.3g} / ADC".format(slope, intercept),
    )

    legend = ROOT.TLegend(0.60, 0.15, 0.94, 0.30)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.026)
    legend.AddEntry(graph, "Fit peak position (Step 2)", "p")
    if data_graph is not None:
        legend.AddEntry(data_graph, "Data peak in segment (Step 3)", "p")
    legend.Draw()

    from plot_velocityfits import draw_header
    draw_header(
        canvas, "Run {} ADC -> PE calibration / ADC: {}".format(
            run, fit.get("short_label", branch)),
        "Fitted line divided by ADC",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, ratio_func, latex, legend] + ([data_graph] if data_graph else [])

    return canvas


def convert_run(
    run: int, step2_dir: Path, step1_summary: Dict[str, Any], root_file: Any,
    outdir: Path, args: argparse.Namespace,
) -> None:
    with (step2_dir / "pe_fits_summary.json").open() as handle:
        step2_summary = json.load(handle)
    with (step2_dir / "pe_fit_curves.json").open() as handle:
        all_curves = json.load(handle)

    run_key = "run{}".format(run)
    channels_meta = step1_summary["channels"][run_key]
    fits = step2_summary["channels"]

    requested = set(args.channels) if args.channels else None

    outdir.mkdir(parents=True, exist_ok=True)
    seg_dir = outdir / "segmented"
    calib_dir = outdir / "calibration"
    calib_divadc_dir = outdir / "calibration_divADC"
    seg_dir.mkdir(parents=True, exist_ok=True)
    calib_dir.mkdir(parents=True, exist_ok=True)
    calib_divadc_dir.mkdir(parents=True, exist_ok=True)

    bin_summary: Dict[str, Any] = {}
    calib_summary: Dict[str, Any] = {}
    seg_png: List[str] = []
    calib_png: List[str] = []
    calib_divadc_png: List[str] = []

    branches = sorted(b for b, f in fits.items() if not f.get("skipped"))
    if requested is not None:
        branches = [b for b in branches if b in requested]

    print("[STEP3] Segmenting {} channel(s) for run {}...".format(len(branches), run))
    for branch in branches:
        fit = fits[branch]
        curves = all_curves.get(branch)
        if not curves:
            continue
        npe_max = fit["npe_max"]
        max_pe = args.max_pe_for_calibration
        if max_pe is None:
            max_pe = npe_max - 1

        edges = []
        for i in range(npe_max - 1):
            ci, cj = curves.get(str(i)), curves.get(str(i + 1))
            if ci is None or cj is None:
                continue
            lo = fit["peaks"][i]["mu"]
            hi = fit["peaks"][i + 1]["mu"]
            edges.append(find_crossing(ci, cj, lo, hi))

        hist_name = channels_meta[branch]["hist_name"]
        hist = root_file.Get("{}/{}".format(run_key, hist_name))
        if not hist:
            continue
        hist = hist.Clone()
        hist.SetDirectory(0)

        seg_canvas = make_segmented_canvas(branch, run, hist, fit, curves, edges)
        seg_path = seg_dir / "{}.png".format(branch)
        seg_canvas.SaveAs(str(seg_path))
        seg_png.append(str(seg_path.relative_to(outdir)))

        plot_lo, plot_hi = fit["plot_range"]
        data_peaks = compute_data_peaks(hist, plot_lo, plot_hi, edges, npe_max, max_pe)
        calib_canvas, calib_result = make_calibration_canvas(
            branch, run, fit, max_pe, data_peaks
        )
        calib_path = calib_dir / "{}.png".format(branch)
        calib_canvas.SaveAs(str(calib_path))
        calib_png.append(str(calib_path.relative_to(outdir)))

        calib_divadc_canvas = make_calibration_divadc_canvas(
            branch, run, fit, max_pe,
            calib_result["slope_pe_per_adc"], calib_result["intercept_pe"],
            data_peaks,
        )
        calib_divadc_path = calib_divadc_dir / "{}.png".format(branch)
        calib_divadc_canvas.SaveAs(str(calib_divadc_path))
        calib_divadc_png.append(str(calib_divadc_path.relative_to(outdir)))

        bin_summary[branch] = {
            "npe_max": npe_max,
            "bin_edges_adc": edges,
            "bin_labels": ["{} p.e.".format(i) for i in range(npe_max)],
            "short_label": fit.get("short_label"),
        }
        calib_result["short_label"] = fit.get("short_label")
        calib_summary[branch] = calib_result

    bins_path = outdir / "pe_bin_boundaries.json"
    with bins_path.open("w") as handle:
        json.dump({"run": run, "channels": bin_summary}, handle, indent=2, sort_keys=True)

    calib_path_json = outdir / "pe_calibration_summary.json"
    with calib_path_json.open("w") as handle:
        json.dump({"run": run, "channels": calib_summary}, handle, indent=2, sort_keys=True)

    generate_html(
        seg_png, str(outdir), plots_per_row=3,
        output_html=str(outdir / "pe_segmentation.html"),
        title="CaloX PE-calibration Step 3: ADC segmentation",
        intro_text=(
            "**Run:** {}.\n".format(run)
            + "Vertical lines mark where consecutive Crystal Ball p.e. "
            "components (from Step 2) cross - the ADC-integral value above "
            "which an event is more consistent with the higher p.e. count "
            "than the lower one. Everything between two consecutive lines "
            "is assigned that many photoelectrons."
        ),
    )
    generate_html(
        calib_png, str(outdir), plots_per_row=3,
        output_html=str(outdir / "pe_calibration.html"),
        title="CaloX PE-calibration Step 3: ADC -> PE calibration fit",
        intro_text=(
            "**Run:** {}.\n".format(run)
            + "Linear fit of photoelectron count (y) vs. the Step 2 fitted "
            "ADC-integral peak position (x, blue), per channel. The fitted "
            "slope is the PE-per-ADC conversion factor; its inverse is the "
            "gain (ADC per p.e.). Open black points are the actual ADC-"
            "integral peak position (not the mean) of the events in each Step 3 segmentation bin - "
            "since the blue points are `mu_ped + i*dpe` by construction of "
            "the Step 2 model, they are exactly linear regardless of the "
            "data, so the black points are the real check of linearity."
        ),
    )
    generate_html(
        calib_divadc_png, str(outdir), plots_per_row=3,
        output_html=str(outdir / "pe_calibration_divADC.html"),
        title="CaloX PE-calibration Step 3: calibration fit / ADC",
        intro_text=(
            "**Run:** {}.\n".format(run)
            + "The linear ADC -> PE calibration fit (`pe_calibration.html`), "
            "divided by ADC: (slope x ADC + intercept) / ADC = slope + "
            "intercept/ADC. A perfectly proportional calibration (zero "
            "pedestal offset) would be flat at the slope; deviation from "
            "flat, especially at low ADC, shows the effect of the nonzero "
            "intercept. Open black points are the same actual per-segment "
            "ADC-integral peak positions as in `pe_calibration.html`, divided by ADC."
        ),
    )
    print("Wrote {}".format(bins_path))
    print("Wrote {}".format(calib_path_json))
    print("Wrote {}".format(outdir / "pe_segmentation.html"))
    print("Wrote {}".format(outdir / "pe_calibration.html"))
    print("Wrote {}".format(outdir / "pe_calibration_divADC.html"))


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    step1_dir = Path(args.step1_dir)
    step2_dir = Path(args.step2_dir)
    with (step1_dir / "pe_integrals_summary.json").open() as handle:
        step1_summary = json.load(handle)

    root_file = ROOT.TFile.Open(str(step1_dir / "pe_integrals.root"))
    base_outdir = Path(args.outdir)

    if args.runs:
        for run in args.runs:
            run_key = "run{}".format(run)
            run_step2_dir = step2_dir / run_key
            if not (run_step2_dir / "pe_fits_summary.json").exists():
                print("[SKIP RUN] {} - no Step 2 output at {}.".format(run, run_step2_dir))
                continue
            print("\n[RUN {}]".format(run))
            convert_run(
                run, run_step2_dir, step1_summary, root_file,
                base_outdir / run_key, args,
            )
    else:
        with (step2_dir / "pe_fits_summary.json").open() as handle:
            run = json.load(handle)["run"]
        convert_run(run, step2_dir, step1_summary, root_file, base_outdir, args)

    root_file.Close()


if __name__ == "__main__":
    main()
