#!/usr/bin/env python3
"""Step 6 (AttenuationCorrections) of the ADC -> photoelectron study.

Combines a per-run Step 1 (mean/RMS ADC integral) summary with a per-run
Step 3 (ADC -> PE calibration) summary into mean photoelectrons per channel
per run, then plots mean PE vs. beam position - the same "attenuation scan"
as scripts/plot_attenuation.py, but on PE instead of raw ADC.

No event-level pass is needed here: since each run's ADC -> PE calibration is
linear (PE = slope*ADC + intercept, from Step 3), the mean and standard error
of PE follow directly from the mean/RMS/count of ADC already in Step 1's
summary - mean(PE) = slope*mean(ADC) + intercept, err(mean PE) =
|slope|*rms(ADC)/sqrt(N) - without re-reading any waveform.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List

import ROOT

from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer

from plot_velocityfits import draw_header, load_run_positions, safe_name


@dataclass(frozen=True)
class PEAttenuationPoint:
    run: int
    position_mm: float
    mean_pe: float
    mean_pe_error: float
    n_events: int
    slope_pe_per_adc: float
    intercept_pe: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Mean photoelectrons per channel vs. beam position (Step 6)."
    )
    parser.add_argument("--runs", nargs="+", type=int, required=True)
    parser.add_argument(
        "--step1-dir",
        default="results/PhotoElectrons/Step6_AttenuationCorrections/Step1_IntegralCalculation",
    )
    parser.add_argument(
        "--step3-dir",
        default="results/PhotoElectrons/Step6_AttenuationCorrections/Step3_ADCtoPEConversion",
        help="Directory containing one run<N>/pe_calibration_summary.json per run.",
    )
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--position-label", default=None)
    parser.add_argument("--positions-json", default=None)
    parser.add_argument("--positions", nargs="+", default=None, metavar="RUN:MM")
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument(
        "--min-events", type=float, default=50.0,
        help="Drop a run/channel point if Step 1 saw fewer events than this.",
    )
    parser.add_argument(
        "--outdir", default="results/PhotoElectrons/Step6_AttenuationCorrections"
    )
    args = parser.parse_args()
    if args.position_label is None:
        args.position_label = args.position_key
    args.no_velocity_fits = False
    return args


def collect_pe_points(
    args: argparse.Namespace, positions: Dict[int, float]
) -> Dict[str, List[PEAttenuationPoint]]:
    step1_dir = Path(args.step1_dir)
    step3_dir = Path(args.step3_dir)
    with (step1_dir / "pe_integrals_summary.json").open() as handle:
        step1_summary = json.load(handle)
    step1_channels = step1_summary["channels"]

    requested = set(args.channels) if args.channels else None
    points: Dict[str, List[PEAttenuationPoint]] = {}
    labels: Dict[str, str] = {}

    for run in args.runs:
        if run not in positions:
            continue
        run_key = "run{}".format(run)
        step1_run = step1_channels.get(run_key)
        calib_path = step3_dir / run_key / "pe_calibration_summary.json"
        if step1_run is None or not calib_path.exists():
            print("[SKIP RUN] {} - missing Step 1 or Step 3 output.".format(run))
            continue
        with calib_path.open() as handle:
            calib = json.load(handle)["channels"]

        for branch, entry in step1_run.items():
            if requested is not None and branch not in requested:
                continue
            calib_entry = calib.get(branch)
            if calib_entry is None or calib_entry.get("skipped"):
                continue
            n_events = entry.get("entries", 0.0)
            if n_events < args.min_events:
                continue
            slope = calib_entry["slope_pe_per_adc"]
            intercept = calib_entry["intercept_pe"]
            mean_adc = entry["mean_adc"]
            rms_adc = entry["rms_adc"]
            mean_pe = slope * mean_adc + intercept
            mean_pe_error = abs(slope) * rms_adc / math.sqrt(n_events) if n_events > 0 else 0.0

            points.setdefault(branch, []).append(PEAttenuationPoint(
                run=run, position_mm=positions[run], mean_pe=mean_pe,
                mean_pe_error=mean_pe_error, n_events=int(n_events),
                slope_pe_per_adc=slope, intercept_pe=intercept,
            ))
            labels[branch] = entry.get("short_label", branch)

    return points, labels


def make_canvas(branch: str, label: str, pts: List[PEAttenuationPoint],
                position_label: str) -> Any:
    pts = sorted(pts, key=lambda p: (p.position_mm, p.run))
    graph = ROOT.TGraphErrors(len(pts))
    for k, p in enumerate(pts):
        graph.SetPoint(k, p.position_mm, p.mean_pe)
        graph.SetPointError(k, 0.0, p.mean_pe_error)

    canvas = ROOT.TCanvas("pe_attenuation_{}".format(safe_name(branch)), "", 1000, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    graph.SetTitle("")
    graph.GetXaxis().SetTitle("{} [mm]".format(position_label))
    graph.GetYaxis().SetTitle("Mean photoelectrons")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.Draw("AP")

    draw_header(
        canvas, "PE attenuation scan: {}".format(label),
        "Runs {}-{}".format(pts[0].run, pts[-1].run) if pts else "",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph]
    return canvas


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    positions = load_run_positions(args)
    points, labels = collect_pe_points(args, positions)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    root_path = outdir / "mean_pe_vs_position.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")

    keys: List[str] = []
    jsons: Dict[str, str] = {}
    summary: Dict[str, Any] = {}

    for branch in sorted(points):
        pts = points[branch]
        if len(pts) < 2:
            continue
        canvas = make_canvas(branch, labels[branch], pts, args.position_label)
        name = canvas.GetName()
        root_file.cd()
        canvas.Write(name)
        keys.append(name)
        jsons[name] = str(ROOT.TBufferJSON.ToJSON(canvas).Data())
        summary[branch] = {
            "short_label": labels[branch],
            "points": [asdict(p) for p in sorted(pts, key=lambda p: p.position_mm)],
        }

    root_file.Close()

    html_path = outdir / "mean_pe_vs_position.html"
    generate_jsroot_html(
        keys, jsons, plots_per_row=2, output_html=str(html_path),
        title="CaloX PE-calibration Step 6: mean PE vs beam position",
        intro_text=(
            "**Runs:** {}-{}.\n".format(min(args.runs), max(args.runs))
            + "**Observable:** mean photoelectrons per channel per run, from "
            "PE = slope x mean(ADC) + intercept using that run's own "
            "independently-fit Step 1-3 ADC -> PE calibration "
            "(results/PhotoElectrons/Step6_AttenuationCorrections/"
            "Step{{1,2,3}}_*). Error bars are the standard error of the mean, "
            "propagated from the ADC-integral RMS through the calibration "
            "slope.\n"
            + "Compare with the ADC-based attenuation scan at "
            "results/attenuation/runs_1897_1942_no_pedestal.\n"
            + "Runs/channels with fewer than {:.0f} events, or where a "
            "channel's calibration fit was skipped, are dropped.".format(
                args.min_events)
        ),
    )
    json_path = outdir / "mean_pe_vs_position_summary.json"
    with json_path.open("w") as handle:
        json.dump(
            {
                "runs": args.runs, "position_key": args.position_key,
                "channels": summary,
            },
            handle, indent=2, sort_keys=True,
        )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
