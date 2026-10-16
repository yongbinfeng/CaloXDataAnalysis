#!/usr/bin/env python3
"""Recover reference-MPV velocity plots from a completed timing ROOT file."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List

import ROOT

from channels.channel_map import build_drs_boards
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.utils import get_channel_var, get_hist_mpv

from plot_velocityfits import (
    ChannelInfo,
    VelocityPoint,
    compute_shared_family_slopes,
    make_velocity_canvas,
    safe_name,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--runs", nargs="+", type=int, required=True)
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fit-min-events", type=int, default=25)
    parser.add_argument("--position-label", default="Table position X")
    return parser.parse_args()


def channel_metadata(run: int) -> Dict[str, ChannelInfo]:
    output: Dict[str, ChannelInfo] = {}
    for _, board in build_drs_boards(run_number=run).items():
        for channel in board:
            if channel.is_reference:
                continue
            material = str(get_channel_var(channel))
            if material == "Sci":
                continue
            branch = channel.get_channel_name(blsub=False)
            size = "6 mm" if channel.is6mm else "3 mm"
            output[branch] = ChannelInfo(
                branch=branch,
                label="{} | {}, {}".format(branch, material, size),
                short_label="B{} G{} Ch{} | {}, {}".format(
                    board.board_no, channel.group_no, channel.channel_no,
                    material, size,
                ),
                material=material,
                size=size,
            )
    return output


def main() -> None:
    args = parse_args()
    setup_root(batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    with Path(args.runlist).open() as handle:
        runlist = json.load(handle)
    positions = {
        run: float(runlist[str(run)]["Table position X"])
        for run in args.runs
    }
    channels = channel_metadata(args.runs[0])
    source = ROOT.TFile.Open(args.input_root, "READ")
    if not source or source.IsZombie():
        raise RuntimeError("Cannot open {}".format(args.input_root))

    points_by_branch: Dict[str, List[VelocityPoint]] = {}
    for branch in sorted(channels):
        points: List[VelocityPoint] = []
        for run in args.runs:
            hist = source.Get(
                "run{}/toa/{}".format(run, safe_name(branch))
            )
            if not hist:
                # The failed producer changed directory after its first channel,
                # so the remaining histograms landed in timing_fits/runN.
                hist = source.Get(
                    "timing_fits/run{}/{}".format(run, safe_name(branch))
                )
            if not hist:
                continue
            n_entries = int(round(hist.Integral(1, hist.GetNbinsX())))
            if n_entries < args.fit_min_events:
                continue
            mpv, error = get_hist_mpv(hist, window_ts=1.2)
            if error <= 0.0 or error >= 999.0:
                continue
            points.append(VelocityPoint(
                run=run,
                position_mm=positions[run],
                mean_ns=float(mpv),
                point_error_ns=float(error),
                model="CoreMPV",
            ))
        points_by_branch[branch] = points

    shared = compute_shared_family_slopes(points_by_branch, channels)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_root_path = outdir / "timing_velocity_fits_recovered.root"
    output_root = ROOT.TFile(str(output_root_path), "RECREATE")
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    summaries = {}
    for branch in sorted(points_by_branch):
        canvas, canvas_json, summary = make_velocity_canvas(
            branch,
            channels[branch],
            points_by_branch[branch],
            args,
            "MCP-clean events",
            shared_fit=shared.get(channels[branch].material),
        )
        summaries[branch] = {
            **asdict(summary),
            "points": [asdict(point) for point in points_by_branch[branch]],
        }
        if canvas is None or canvas_json is None:
            continue
        name = canvas.GetName()
        output_root.cd()
        canvas.Write(name)
        keys.append(name)
        jsons[name] = canvas_json
    output_root.Close()
    source.Close()

    html_path = outdir / "timing_velocity_fits.html"
    generate_jsroot_html(
        keys,
        jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title="CaloX reference-MPV Cherenkov velocity fits",
        intro_text=(
            "**Runs:** {}.\n"
            "* Robust +/-1.2 ns core MPV, equivalent to +/-6 TS.\n"
            "* One simultaneous weighted slope per Cherenkov material family; "
            "independent channel intercepts.\n"
            "* Sci channels excluded."
        ).format(", ".join(str(run) for run in args.runs)),
    )
    summary_path = outdir / "timing_summary.json"
    with summary_path.open("w") as handle:
        json.dump({
            "recovered_from": args.input_root,
            "runs": args.runs,
            "positions_mm": positions,
            "shared_family_fits": {
                family: {
                    "slope_ns_per_mm": values[0],
                    "slope_error_ns_per_mm": values[1],
                    "speed_cm_per_ns": 0.1 / abs(values[0]),
                    "speed_error_cm_per_ns": (
                        0.1 * abs(values[1]) / (values[0] ** 2)
                    ),
                }
                for family, values in shared.items()
            },
            "velocity_fits": summaries,
        }, handle, indent=2, sort_keys=True)
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(output_root_path))
    print("Wrote {}".format(summary_path))


if __name__ == "__main__":
    main()
