#!/usr/bin/env python3
"""Rebuild pair_timing_resolution.html from an already-produced results directory.

Re-reads the saved ROOT histograms and the JSON fit summary written by
scripts/plot_channel_pair_timing_resolution.py and regenerates the full JSROOT
gallery (per-pair fits, the resolution-vs-position overlay, and the combined
resolution-distribution overlay) without re-running the RDataFrame analysis.

Useful after changing/adding a plotting function in that script: just rerun
this against the existing --outdir instead of resubmitting the multi-hour job.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import ROOT

from scripts.plot_channel_pair_timing_resolution import (
    ChannelPair,
    PairFitSummary,
    make_pair_canvas,
    make_resolution_distribution_canvas,
    make_scan_overlay_canvas,
    safe_name,
)
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir", required=True,
        help="Directory containing pair_timing_resolution.root and _summary.json.",
    )
    return parser.parse_args()


def record_to_pair(record: Dict[str, Any]) -> ChannelPair:
    return ChannelPair(
        material=record["material"],
        bridge_no=record["bridge_no"], board_no=record["board_no"], group_no=record["group_no"],
        i_tower_y=record["i_tower_y"],
        branch_low=record["branch_low"], branch_high=record["branch_high"],
        ix_low=record["ix_low"], ix_high=record["ix_high"],
    )


def record_to_summary(record: Dict[str, Any]) -> PairFitSummary:
    return PairFitSummary(
        n_entries=record["n_entries"],
        mean_ns=record["mean_ns"], mean_err_ns=record["mean_err_ns"],
        sigma_ns=record["sigma_ns"], sigma_err_ns=record["sigma_err_ns"],
        fit_ok=record["fit_ok"],
    )


def main() -> None:
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    summary_path = outdir / "pair_timing_resolution_summary.json"
    root_path = outdir / "pair_timing_resolution.root"
    with open(summary_path) as handle:
        payload = json.load(handle)
    records: List[Dict[str, Any]] = payload["records"]
    materials: List[str] = payload["materials"]
    position_label = payload.get("position_key", "position")

    root_file = ROOT.TFile.Open(str(root_path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot open {}".format(root_path))

    canvas_keys: List[str] = []
    canvas_jsons: Dict[str, str] = {}

    records_by_run: Dict[int, List[Dict[str, Any]]] = {}
    for record in records:
        records_by_run.setdefault(record["run"], []).append(record)

    for run in sorted(records_by_run):
        for record in records_by_run[run]:
            hist = root_file.Get("run{}/{}/{}".format(run, record["material"], safe_name(record["pair_id"])))
            if not hist:
                print("[WARN] missing histogram for run {} {}".format(run, record["pair_id"]))
                continue
            hist = hist.Clone(safe_name("regen_{}_{}".format(run, record["pair_id"])))
            hist.SetDirectory(0)
            pair = record_to_pair(record)
            summary = record_to_summary(record)
            position_mm = record.get("position_mm")
            subtitle = "run {}, x={} mm".format(run, position_mm) if position_mm is not None else "run {}".format(run)
            canvas, canvas_json = make_pair_canvas(hist, pair, run, summary, subtitle)
            canvas_keys.append(canvas.GetName())
            canvas_jsons[canvas.GetName()] = canvas_json
        canvas_keys.append("NEWLINE")

    scan_canvas = make_scan_overlay_canvas(records, materials, position_label)
    if scan_canvas is not None:
        canvas, canvas_json = scan_canvas
        canvas_keys.insert(0, "NEWLINE")
        canvas_keys.insert(0, canvas.GetName())
        canvas_jsons[canvas.GetName()] = canvas_json

    distribution_canvas = make_resolution_distribution_canvas(records, materials)
    if distribution_canvas is not None:
        canvas, canvas_json = distribution_canvas
        canvas_keys.insert(0, "NEWLINE")
        canvas_keys.insert(0, canvas.GetName())
        canvas_jsons[canvas.GetName()] = canvas_json

    n_ok = sum(1 for r in records if r["fit_ok"])
    intro_lines = [
        "**Method**: for each pair of horizontally-adjacent, same-material 3 mm "
        "channels within a single DRS readout group, dt = TS_cfd_ref(high-x) - "
        "TS_cfd_ref(low-x) is histogrammed per run and fit with a Gaussian. The "
        "group's reference-channel jitter cancels exactly in this difference. "
        "Quoted resolution = sigma(dt) / sqrt(2), the single-channel-equivalent "
        "time resolution.",
        "* Electron ID: Cherenkov-counter beam-PID cut ({})".format(payload.get("particle") or "none"),
        "* Runs: {}-{} ({} runs), position key: {}".format(
            min(payload["runs"]), max(payload["runs"]), len(payload["runs"]), position_label),
        "* Materials: {}".format(", ".join(materials)),
        "* {} of {} (run, pair) histograms have a valid Gaussian fit".format(n_ok, len(records)),
        "* Not yet split by photoelectron count; each combination is one horizontal channel pair.",
    ]
    generate_jsroot_html(
        canvas_keys, canvas_jsons, plots_per_row=3,
        output_html=str(outdir / "pair_timing_resolution.html"),
        title="CaloX TB2026 adjacent-channel-pair timing resolution",
        intro_text="\n".join(intro_lines),
    )
    print("[WROTE] {}".format(outdir / "pair_timing_resolution.html"))


if __name__ == "__main__":
    main()
