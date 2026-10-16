#!/usr/bin/env python3
"""Step 4 (PE vs TOA) / Step 5 (PE vs rise time) of the ADC -> PE calibration
study.

For each channel, uses Step 3's discrete PE bin edges
(pe_bin_boundaries.json) to assign every event's ADC integral
({ch}_energy_adaptive, recomputed the same way as Step 1) to a photoelectron
bin, then plots the mean of a per-event timing quantity in each bin:

  --quantity toa       {ch}_Time_cfd_mcp [ns], the repo's standard
                       MCP/reference-corrected time of arrival.
  --quantity risetime  {ch}_risetime_adaptive [ns], the new CFD 20% -> 90%
                       leading-edge rise time (variables.drs.add_rise_time).

Same "PE on Y, mean of X per bin" style as Step 3's calibration plot - just a
different X quantity - so the two are directly comparable.

Caveat: the 0 p.e. (pedestal) bin has no real pulse to time, so its TOA/rise
time reflects wherever compute_cfd_integral's threshold-free CFD happened to
land on noise, not a physical leading edge. Treat that point with caution.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.html_generator import generate_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from variables.drs import add_cfd_to_baseline_energy, add_rise_time

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    draw_header,
    make_manager_args,
    safe_name,
)


QUANTITY_INFO = {
    "toa": {
        "label": "Time of arrival [ns]",
        "step_name": "Step4TimevsPE",
        "title": "Time vs PE",
    },
    "risetime": {
        "label": "Rise time, CFD 20% -> 90% [ns]",
        "step_name": "Step5PEvsRiseTime",
        "title": "Rise time vs PE",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Mean per-event timing quantity in each Step 3 PE bin, per channel."
    )
    parser.add_argument("--quantity", choices=sorted(QUANTITY_INFO), required=True)
    parser.add_argument(
        "--step3-dir", default="results/PhotoElectrons/Step3_ADCtoPEConversion"
    )
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument("--particle", default=None)
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--n-threads", type=int, default=10)

    parser.add_argument(
        "--integral", choices=["adaptive", "repo"], default="adaptive",
        help="ADC integral binned into PE: must match the one used in Step 1 "
        "(adaptive = {ch}_energy_adaptive, repo = {ch}_energy).",
    )

    # Must match Step 1's compute_cfd_to_baseline_integral parameters so the
    # ADC integral recomputed here lines up with Step 3's bin edges.
    parser.add_argument("--anchor-offset", type=int, default=-330)
    parser.add_argument("--search-half", type=int, default=60)
    parser.add_argument("--cfd-fraction", type=float, default=0.20)
    parser.add_argument("--stop-fraction", type=float, default=0.05)
    parser.add_argument("--stop-abs", type=float, default=4.0)
    parser.add_argument("--stop-consecutive", type=int, default=3)
    parser.add_argument("--max-width", type=int, default=250)
    parser.add_argument("--baseline-pre-lo", type=int, default=80)
    parser.add_argument("--baseline-pre-hi", type=int, default=20)
    parser.add_argument("--min-amplitude", type=float, default=0.0)
    parser.add_argument("--smooth-half", type=int, default=3)
    parser.add_argument("--min-width", type=int, default=40)

    # Rise-time-specific thresholds.
    parser.add_argument("--rise-frac-lo", type=float, default=0.20)
    parser.add_argument("--rise-frac-hi", type=float, default=0.90)

    parser.add_argument(
        "--min-entries-per-bin", type=float, default=5.0,
        help="Drop PE bins with fewer events than this from the plot/fit.",
    )
    parser.add_argument("--outdir", default=None,
                        help="Default: results/PhotoElectrons/<Step4TimevsPE|Step5PEvsRiseTime>")
    return parser.parse_args()


def pe_bin_expr(column: str, edges: List[float]) -> str:
    """Nested-ternary RDF expression assigning `column` to bin i where
    edges[i-1] <= column < edges[i] (bin 0 below edges[0], last bin above the
    final edge) - the same crossings Step 3 drew as vertical lines.
    """
    npe_max = len(edges) + 1
    expr = str(npe_max - 1)
    for i in reversed(range(len(edges))):
        expr = "(({}) < ({:.10g}) ? {} : {})".format(column, edges[i], i, expr)
    return expr


def make_pe_vs_quantity_canvas(
    branch: str, run: int, label: str, quantity_label: str,
    points: List[Tuple[int, float, float]],
) -> Tuple[Any, Dict[str, float]]:
    """points: list of (pe, mean_value, mean_error)."""
    graph = ROOT.TGraphErrors(len(points))
    for k, (pe, mean_val, mean_err) in enumerate(points):
        graph.SetPoint(k, mean_val, pe)
        graph.SetPointError(k, mean_err, 1.0e-3)

    canvas = ROOT.TCanvas("pe_vs_q_{}".format(branch), "", 850, 650)
    canvas.SetLeftMargin(0.14)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)

    graph.SetTitle("")
    graph.GetXaxis().SetTitle(quantity_label)
    graph.GetYaxis().SetTitle("Photoelectrons")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.Draw("AP")

    vals = [mean_val for _, mean_val, _ in points]
    result: Dict[str, float] = {"n_points": len(points)}
    if len(points) >= 2 and max(vals) > min(vals):
        pad = 0.1 * (max(vals) - min(vals))
        fit_func = ROOT.TF1(
            "pe_vs_q_fit_{}".format(branch), "pol1", min(vals) - pad, max(vals) + pad
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
        latex.SetTextSize(0.028)
        latex.DrawLatexNDC(0.17, 0.85, "PE = {:.4g}(#pm{:.2g}) #times X + {:.3g}(#pm{:.2g})".format(
            slope, slope_err, intercept, intercept_err))
        if ndf > 0:
            latex.DrawLatexNDC(0.17, 0.80, "#chi^{{2}}/ndf = {:.2f}/{}".format(chi2, ndf))
        canvas._fit_keepalive = [fit_func, latex]
        result.update({
            "slope_pe_per_unit": slope, "slope_err": slope_err,
            "intercept_pe": intercept, "intercept_err": intercept_err,
            "chi2": chi2, "ndf": ndf,
        })

    draw_header(
        canvas, "Run {} {}: {}".format(run, quantity_label, label),
        "Mean per PE bin (Step 3 edges)",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph]

    return canvas, result


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    info = QUANTITY_INFO[args.quantity]
    outdir = Path(args.outdir) if args.outdir else Path("results/PhotoElectrons") / info["step_name"]
    outdir.mkdir(parents=True, exist_ok=True)

    step3_dir = Path(args.step3_dir)
    with (step3_dir / "pe_bin_boundaries.json").open() as handle:
        bins_data = json.load(handle)
    run = bins_data["run"]
    bin_meta = bins_data["channels"]

    requested = set(args.channels) if args.channels else None
    selection = "MCP-clean events" if not args.no_mcp_clean else "All events"

    manager_args = argparse.Namespace(
        run=run, first_event=args.first_event, last_event=args.last_event,
        json_file=args.json_file,
    )
    manager = CaloXAnalysisManager(manager_args).prepare(do_fers=False)
    manager.apply_mcp_selection(flag_only=True)
    if not args.no_mcp_clean:
        manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
    if args.particle:
        manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

    live = active_timing_channels(manager, requested)
    branches = sorted(b for b in live if b in bin_meta)

    rdf = add_cfd_to_baseline_energy(
        manager.rdf, branches,
        anchor_offset=args.anchor_offset, search_half=args.search_half,
        cfd_fraction=args.cfd_fraction, stop_fraction=args.stop_fraction,
        stop_abs=args.stop_abs, stop_consecutive=args.stop_consecutive,
        max_width=args.max_width, baseline_pre_lo=args.baseline_pre_lo,
        baseline_pre_hi=args.baseline_pre_hi, min_amplitude=args.min_amplitude,
        smooth_half=args.smooth_half, min_width=args.min_width,
    )

    if args.quantity == "risetime":
        rdf = add_rise_time(
            rdf, branches, anchor_offset=args.anchor_offset,
            search_half=args.search_half, frac_lo=args.rise_frac_lo,
            frac_hi=args.rise_frac_hi, baseline_pre_lo=args.baseline_pre_lo,
            baseline_pre_hi=args.baseline_pre_hi, min_amplitude=args.min_amplitude,
        )

    columns = {str(c) for c in rdf.GetColumnNames()}

    # Book everything lazily, then run in one graph.
    actions: Dict[str, Dict[int, Dict[str, Any]]] = {}
    used_branches: List[str] = []
    for branch in branches:
        adc_col = "{}_energy{}".format(
            branch, "" if args.integral == "repo" else "_adaptive")
        if args.quantity == "toa":
            value_col = "{}_Time_cfd_mcp".format(branch)
            valid_expr = "std::isfinite({})".format(value_col)
        else:
            value_col = "{}_risetime_adaptive".format(branch)
            valid_col = "{}_risetime_valid_adaptive".format(branch)
            valid_expr = "{} && std::isfinite({})".format(valid_col, value_col)
        if adc_col not in columns or value_col not in columns:
            continue

        edges = bin_meta[branch]["bin_edges_adc"]
        npe_max = bin_meta[branch]["npe_max"]
        bin_col = "{}_pe_bin_step4".format(branch)
        branch_rdf = rdf.Define(bin_col, pe_bin_expr(adc_col, edges))

        actions[branch] = {}
        for pe in range(npe_max):
            filtered = branch_rdf.Filter(
                "{} == {} && {}".format(bin_col, pe, valid_expr)
            )
            actions[branch][pe] = {
                "mean": filtered.Mean(value_col),
                "stddev": filtered.StdDev(value_col),
                "count": filtered.Count(),
            }
        used_branches.append(branch)

    flat = [
        action
        for per_pe in actions.values()
        for entry in per_pe.values()
        for action in entry.values()
    ]
    print("[STEP4/5] Running {} RDF actions over {} channel(s)...".format(
        len(flat), len(used_branches)))
    ROOT.RDF.RunGraphs(flat)

    summary: Dict[str, Any] = {}
    png_files: List[str] = []
    for branch in used_branches:
        points = []
        for pe, entry in sorted(actions[branch].items()):
            count = float(entry["count"].GetValue())
            if count < args.min_entries_per_bin:
                continue
            mean_val = float(entry["mean"].GetValue())
            stddev = float(entry["stddev"].GetValue())
            mean_err = stddev / (count ** 0.5) if count > 0 else 0.0
            points.append((pe, mean_val, mean_err))
        if len(points) < 2:
            continue

        canvas, fit_result = make_pe_vs_quantity_canvas(
            branch, run, bin_meta[branch].get("short_label", branch),
            info["label"], points,
        )
        png_path = outdir / "{}.png".format(branch)
        canvas.SaveAs(str(png_path))
        png_files.append(str(png_path.relative_to(outdir)))

        fit_result["short_label"] = bin_meta[branch].get("short_label")
        fit_result["points"] = [
            {"pe": pe, "mean": mean_val, "mean_err": mean_err}
            for pe, mean_val, mean_err in points
        ]
        summary[branch] = fit_result

    json_path = outdir / "{}_summary.json".format(args.quantity)
    with json_path.open("w") as handle:
        json.dump({"run": run, "quantity": args.quantity, "channels": summary},
                  handle, indent=2, sort_keys=True)

    html_path = outdir / "{}.html".format(info["step_name"].lower())
    generate_html(
        png_files, str(outdir), plots_per_row=3, output_html=str(html_path),
        title="CaloX PE-calibration {}".format(info["title"]),
        intro_text=(
            "**Run:** {}.\n".format(run)
            + "**Y axis:** photoelectron count, from Step 3's bin edges "
            "(pe_bin_boundaries.json) applied to `{{ch}}_energy{}`.\n".format(
                "" if args.integral == "repo" else "_adaptive")
            + "**X axis:** mean `{}` of the events in each PE bin, error bars "
            "are the standard error of that mean.\n".format(info["label"])
            + "Bins with fewer than {:.0f} events are dropped.\n".format(
                args.min_entries_per_bin)
            + "**Caveat:** the 0 p.e. bin has no real pulse to time - its "
            "point reflects wherever the threshold-free CFD landed on noise, "
            "not a physical leading edge.\n"
            + "**Selection:** {}.".format(selection)
        ),
    )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
