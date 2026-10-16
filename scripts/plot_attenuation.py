#!/usr/bin/env python3
"""Per-channel DRS ADC-integral attenuation scans for TB2026."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List

import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    draw_header,
    load_run_positions,
    make_manager_args,
    safe_name,
)


@dataclass(frozen=True)
class AttenuationPoint:
    run: int
    position_mm: float
    mean_adc: float
    mean_error_adc: float
    stddev_adc: float
    n_events: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot mean per-channel DRS ADC pulse integral versus table position."
    )
    parser.add_argument("--runs", nargs="+", type=int, required=True)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--position-label", default=None)
    parser.add_argument("--positions-json", default=None)
    parser.add_argument("--positions", nargs="+", default=None, metavar="RUN:MM")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--channels", nargs="+", default=None, metavar="DRS_BRANCH")
    parser.add_argument("--particle", default=None)
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument(
        "--outdir", default="results/attenuation/runs_1897_1943"
    )
    args = parser.parse_args()
    if len(set(args.runs)) != len(args.runs):
        parser.error("--runs contains a duplicate run number")
    if args.position_label is None:
        args.position_label = args.position_key
    # load_run_positions is shared with the velocity driver.
    args.no_velocity_fits = False
    return args


def collect_run_points(
    args: argparse.Namespace,
    positions: Dict[int, float],
) -> tuple[Dict[str, ChannelInfo], Dict[str, List[AttenuationPoint]]]:
    requested = set(args.channels) if args.channels else None
    channel_info: Dict[str, ChannelInfo] = {}
    points: Dict[str, List[AttenuationPoint]] = {}

    for run in args.runs:
        if run not in positions:
            print("[WARN] Run {} has no position; skipped.".format(run))
            continue
        print("\n[RUN {}] position={:g} mm".format(run, positions[run]))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        live = active_timing_channels(manager, requested)
        for branch, meta in live.items():
            channel_info.setdefault(branch, meta)

        actions: Dict[str, tuple[Any, Any, Any]] = {}
        all_actions: List[Any] = []
        for branch in sorted(live):
            column = "{}_energy".format(branch)
            rdf = manager.rdf.Filter("std::isfinite({})".format(column))
            mean = rdf.Mean(column)
            stddev = rdf.StdDev(column)
            count = rdf.Count()
            actions[branch] = (mean, stddev, count)
            all_actions.extend((mean, stddev, count))
        if all_actions:
            ROOT.RDF.RunGraphs(all_actions)

        for branch, (mean_result, stddev_result, count_result) in actions.items():
            count = int(count_result.GetValue())
            if count < 2:
                print("[WARN] Run {} {} has only {} events.".format(run, branch, count))
                continue
            mean = float(mean_result.GetValue())
            stddev = float(stddev_result.GetValue())
            if not math.isfinite(mean) or not math.isfinite(stddev):
                continue
            point = AttenuationPoint(
                run=run,
                position_mm=positions[run],
                mean_adc=mean,
                mean_error_adc=stddev / math.sqrt(count),
                stddev_adc=stddev,
                n_events=count,
            )
            points.setdefault(branch, []).append(point)
            print(
                "[ADC] {} mean={:.3f} +/- {:.3f}, RMS={:.3f}, N={}".format(
                    live[branch].short_label, mean, point.mean_error_adc, stddev, count
                )
            )
    return channel_info, points


def make_canvas(
    branch: str,
    channel: ChannelInfo,
    points: List[AttenuationPoint],
    args: argparse.Namespace,
    selection_label: str,
) -> tuple[Any, str]:
    ordered = sorted(points, key=lambda point: (point.position_mm, point.run))
    graph = ROOT.TGraphErrors(len(ordered))
    graph.SetName("gr_attenuation_{}".format(safe_name(branch)))
    graph.SetTitle("")
    for index, point in enumerate(ordered):
        graph.SetPoint(index, point.position_mm, point.mean_adc)
        graph.SetPointError(index, 0.0, point.mean_error_adc)

    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.05)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)

    name = "attenuation_{}".format(safe_name(branch))
    canvas = ROOT.TCanvas(name, "", 1100, 800)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    graph.Draw("AP")
    graph.GetXaxis().SetTitle("{} [mm]".format(args.position_label))
    graph.GetYaxis().SetTitle("Mean ADC integral [ADC]")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.GetXaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetTitleOffset(1.35)

    legend = ROOT.TLegend(0.15, 0.72, 0.53, 0.81)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.023)
    legend.AddEntry(graph, "Run mean (#sigma/#sqrt{N} errors)", "pe")
    legend.Draw()
    draw_header(canvas, "Attenuation scan: {}".format(channel.short_label), selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, legend]
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    positions = load_run_positions(args)
    channel_info, points = collect_run_points(args, positions)
    if not points:
        raise RuntimeError("No attenuation points were produced")

    selection = "MCP-clean events" if not args.no_mcp_clean else "All events"
    root_path = outdir / "attenuation_scans.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    for branch in sorted(points):
        canvas, canvas_json = make_canvas(
            branch, channel_info[branch], points[branch], args, selection
        )
        name = "attenuation_{}".format(safe_name(branch))
        root_file.cd()
        canvas.Write(name)
        keys.append(name)
        jsons[name] = canvas_json
    root_file.Close()

    html_path = outdir / "attenuation_scans.html"
    generate_jsroot_html(
        keys,
        jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title="CaloX TB2026 per-channel ADC attenuation scans",
        intro_text=(
            "**Observable:** mean baseline-subtracted CFD pulse integral per run/channel.\n"
            "**Error bars:** standard deviation / sqrt(N).\n"
            "**Position:** {} [mm].\n"
            "**Selection:** {}."
        ).format(args.position_label, selection),
    )

    summary = {
        "runs": args.runs,
        "position_key": args.position_key,
        "selection": selection,
        "observable": "mean baseline-subtracted CFD pulse integral",
        "points": {
            branch: [asdict(point) for point in sorted(values, key=lambda p: p.run)]
            for branch, values in points.items()
        },
    }
    with open(outdir / "attenuation_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("\nWrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(outdir / "attenuation_summary.json"))


if __name__ == "__main__":
    main()
