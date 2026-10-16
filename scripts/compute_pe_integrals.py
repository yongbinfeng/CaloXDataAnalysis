#!/usr/bin/env python3
"""Step 1 of the ADC -> photoelectron calibration study: per-channel
CFD-to-baseline-crossing ADC integral histograms for one or more runs.

This uses exactly one integral definition - compute_cfd_to_baseline_integral,
via variables.drs.add_cfd_to_baseline_energy ({ch}_energy_adaptive) - because
it is the one definition in the repo that both preserves the pedestal (0 p.e.)
peak (no amplitude threshold) and follows the pulse shape on both 3 mm and
6 mm fibers alike. Step 2 resolves individual photoelectron peaks from the
histograms this script writes.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ROOT

from configs.plot_config import get_drs_energy_range
from core.analysis_manager import CaloXAnalysisManager
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from variables.drs import add_cfd_to_baseline_energy

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    draw_header,
    make_manager_args,
    safe_name,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Per-channel CFD-to-baseline-crossing ADC-integral histograms "
            "(Step 1 of the ADC -> PE calibration)."
        )
    )
    parser.add_argument("--runs", nargs="+", type=int, default=[1897])
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument("--particle", default=None)
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--n-threads", type=int, default=10)

    # compute_cfd_to_baseline_integral parameters - same defaults as
    # variables.drs.add_cfd_to_baseline_energy / the adc_integral_fixedwindow_test
    # study, so Step 1's histograms match the "CFD -> baseline crossing" curve
    # from that study exactly.
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

    parser.add_argument(
        "--integral", choices=["adaptive", "repo"], default="adaptive",
        help="adaptive: {ch}_energy_adaptive (CFD -> baseline crossing). "
        "repo: {ch}_energy, the repo's standard compute_cfd_integral fixed "
        "window (returns exactly 0 below its amplitude threshold).",
    )
    parser.add_argument("--rebin-factor", type=int, default=2)
    parser.add_argument(
        "--range-sigma", type=float, default=8.0,
        help="Half-width of the auto-derived histogram range, in RMS units.",
    )
    parser.add_argument(
        "--no-auto-range", action="store_true",
        help="Use the repo's fixed (nbins, xmin, xmax) instead of a per-channel "
        "mean/RMS-derived range.",
    )
    parser.add_argument(
        "--zoom-keep", type=float, default=0.995,
        help="Fraction of entries the displayed (zoomed) x-range must cover.",
    )
    parser.add_argument(
        "--outdir", default="results/PhotoElectrons/Step1_IntegralCalculation"
    )
    return parser.parse_args()


def model_range(is6mm: bool, args: argparse.Namespace) -> Tuple[int, float, float]:
    nbins, xmin, xmax = get_drs_energy_range(is6mm)
    return max(1, nbins) * max(1, args.rebin_factor), xmin, xmax


def auto_range(
    mean: float, rms: float, fallback: Tuple[int, float, float], args: argparse.Namespace
) -> Tuple[int, float, float]:
    if not (math.isfinite(mean) and math.isfinite(rms) and rms > 0.0):
        return fallback
    lo = min(mean - args.range_sigma * rms, 0.0)
    hi = max(mean + args.range_sigma * rms, 0.0)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return fallback
    return fallback[0], lo, hi


def content_range(hist: Any, keep: float) -> Optional[Tuple[float, float]]:
    axis = hist.GetXaxis()
    nbins = hist.GetNbinsX()
    total = hist.Integral(0, nbins + 1)
    if total <= 0.0:
        return None
    tail = max(0.0, min(0.45, 0.5 * (1.0 - keep)))
    want_lo, want_hi = tail * total, (1.0 - tail) * total
    running = hist.GetBinContent(0)
    lo, hi = axis.GetXmin(), axis.GetXmax()
    got_lo = False
    for index in range(1, nbins + 1):
        running += hist.GetBinContent(index)
        if not got_lo and running >= want_lo:
            lo = axis.GetBinLowEdge(index)
            got_lo = True
        if running >= want_hi:
            hi = axis.GetBinUpEdge(index)
            break
    if hi <= lo:
        return None
    pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


def make_canvas(
    run: int, channel: ChannelInfo, hist: Any, selection_label: str, zoom_keep: float,
    stat_label: str = "CFD #rightarrow baseline crossing",
) -> Tuple[Any, str]:
    name = "pe_integral_run{}_{}".format(run, safe_name(channel.branch))
    canvas = ROOT.TCanvas(name, "", 1000, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    canvas.SetLogy(1)

    hist.SetTitle("")
    hist.SetStats(0)
    hist.GetXaxis().SetTitle("ADC integral [ADC]")
    hist.GetYaxis().SetTitle("Events")
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetLabelSize(0.04)
    hist.GetYaxis().SetTitleOffset(1.35)
    hist.SetLineColor(ROOT.kMagenta + 2)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(ROOT.kMagenta - 9, 0.25)

    hist.SetMinimum(0.5)
    hist.SetMaximum(max(10.0, 40.0 * hist.GetMaximum()))
    zoom = content_range(hist, keep=zoom_keep)
    if zoom is not None:
        axis = hist.GetXaxis()
        lo = max(zoom[0], axis.GetXmin())
        hi = min(zoom[1], axis.GetXmax())
        if hi > lo:
            hist.GetXaxis().SetRangeUser(lo, hi)
    hist.Draw("HIST")

    stats = ROOT.TPaveText(0.55, 0.72, 0.95, 0.885, "NDC")
    stats.SetBorderSize(0)
    stats.SetFillStyle(0)
    stats.SetTextFont(42)
    stats.SetTextSize(0.022)
    stats.SetTextAlign(12)
    stats.AddText(stat_label)
    stats.AddText("  Entries {}  Mean {:.1f}  RMS {:.1f}".format(
        int(hist.GetEntries()), hist.GetMean(), hist.GetRMS()
    ))
    stats.Draw()

    draw_header(
        canvas,
        "Run {} ADC integral: {}".format(run, channel.short_label),
        selection_label,
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [hist, stats]
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    requested = set(args.channels) if args.channels else None
    selection = "MCP-clean events" if not args.no_mcp_clean else "All events"

    root_path = outdir / "pe_integrals.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    summary: Dict[str, Dict[str, Dict[str, float]]] = {}

    for run in args.runs:
        print("\n[RUN {}]".format(run))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        live = active_timing_channels(manager, requested)
        rdf = manager.rdf if args.integral == "repo" else add_cfd_to_baseline_energy(
            manager.rdf, sorted(live),
            anchor_offset=args.anchor_offset,
            search_half=args.search_half,
            cfd_fraction=args.cfd_fraction,
            stop_fraction=args.stop_fraction,
            stop_abs=args.stop_abs,
            stop_consecutive=args.stop_consecutive,
            max_width=args.max_width,
            baseline_pre_lo=args.baseline_pre_lo,
            baseline_pre_hi=args.baseline_pre_hi,
            min_amplitude=args.min_amplitude,
            smooth_half=args.smooth_half,
            min_width=args.min_width,
        )
        columns = {str(column) for column in rdf.GetColumnNames()}
        energy_columns: Dict[str, str] = {}
        for branch in sorted(live):
            name = "{}_energy{}".format(
                branch, "" if args.integral == "repo" else "_adaptive")
            if name in columns:
                energy_columns[branch] = name
        missing = sorted(set(live) - set(energy_columns))
        if missing:
            print("[WARN] No {}-integral column for {} channel(s); "
                  "skipped.".format(args.integral, len(missing)))

        # Pre-pass: mean/RMS per channel, to size the histogram range from the
        # data rather than from a fixed repo range tuned for a different
        # integral definition.
        model_ranges: Dict[str, Tuple[int, float, float]] = {}
        if not args.no_auto_range:
            stat_actions = {
                branch: (rdf.Mean(column), rdf.StdDev(column))
                for branch, column in energy_columns.items()
            }
            flat = [action for pair in stat_actions.values() for action in pair]
            print("[RANGE] Pre-pass over {} statistics actions...".format(len(flat)))
            ROOT.RDF.RunGraphs(flat)
            for branch, channel in sorted(live.items()):
                if branch not in stat_actions:
                    continue
                mean_action, rms_action = stat_actions[branch]
                model_ranges[branch] = auto_range(
                    float(mean_action.GetValue()), float(rms_action.GetValue()),
                    model_range(channel.size == "6 mm", args), args,
                )

        results: Dict[str, Any] = {}
        for branch, channel in sorted(live.items()):
            if branch not in energy_columns:
                continue
            column = energy_columns[branch]
            is6mm = channel.size == "6 mm"
            range_tuple = model_ranges.get(branch, model_range(is6mm, args))
            model = ROOT.RDF.TH1DModel(
                "hist_pe_integral_run{}_{}".format(run, safe_name(branch)),
                "", *range_tuple,
            )
            results[branch] = rdf.Filter(
                "std::isfinite({})".format(column)
            ).Histo1D(model, column)

        ROOT.RDF.RunGraphs(list(results.values()))

        run_dir = root_file.mkdir("run{}".format(run))
        run_summary: Dict[str, Dict[str, float]] = {}
        for branch, result in sorted(results.items()):
            hist = result.GetValue().Clone(
                "hist_pe_integral_run{}_{}".format(run, safe_name(branch))
            )
            hist.SetDirectory(0)
            run_dir.cd()
            hist.Write()

            axis = hist.GetXaxis()
            channel = live[branch]
            run_summary[branch] = {
                "axis_xmin_adc": float(axis.GetXmin()),
                "axis_xmax_adc": float(axis.GetXmax()),
                "nbins": int(axis.GetNbins()),
                "entries": float(hist.GetEntries()),
                "mean_adc": float(hist.GetMean()),
                "rms_adc": float(hist.GetRMS()),
                "material": channel.material,
                "size": channel.size,
                "short_label": channel.short_label,
                "hist_name": hist.GetName(),
            }

            canvas, canvas_json = make_canvas(
                run, live[branch], hist, selection, args.zoom_keep,
                "Repo integral (fixed window)" if args.integral == "repo"
                else "CFD #rightarrow baseline crossing",
            )
            name = canvas.GetName()
            root_file.cd()
            canvas.Write(name)
            keys.append(name)
            jsons[name] = canvas_json
        summary["run{}".format(run)] = run_summary

    root_file.Close()
    html_path = outdir / "pe_integrals.html"
    if args.integral == "repo":
        observable = (
            "**Observable:** `{ch}_energy` from `compute_cfd_integral` - the "
            "repo's standard integral: fixed window around the 20% CFD "
            "crossing, exactly 0 when the peak is below its amplitude "
            "threshold. Feeds Step 2 (photoelectron-peak fits).\n"
        )
    else:
        observable = (
            "**Observable:** `{ch}_energy_adaptive` from "
            "`compute_cfd_to_baseline_integral` (CFD start, "
            "baseline-crossing end) - the one integral definition that keeps "
            "the pedestal (0 p.e.) peak and follows the pulse on both 3 mm "
            "and 6 mm fibers. Feeds Step 2 (photoelectron-peak fits).\n"
        )
    intro_text = (
        "**Runs:** " + ", ".join(str(run) for run in args.runs) + ".\n"
        + observable +
        "**Selection:** " + selection + "."
    )
    generate_jsroot_html(
        keys, jsons, plots_per_row=3, output_html=str(html_path),
        title="CaloX PE-calibration Step 1: ADC integrals",
        intro_text=intro_text,
    )
    json_path = outdir / "pe_integrals_summary.json"
    with json_path.open("w") as handle:
        json.dump(
            {
                "runs": args.runs,
                "selection": selection,
                "integral_definition": (
                    "compute_cfd_integral" if args.integral == "repo"
                    else "compute_cfd_to_baseline_integral"
                ),
                "channels": summary,
            },
            handle, indent=2, sort_keys=True,
        )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
