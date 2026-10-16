#!/usr/bin/env python3
"""Step 2 of the ADC -> photoelectron calibration study: fit each channel's
Step 1 ADC-integral histogram with a multi-peak RooCrystalBall model
(utils.fitter.fit_pe_spectrum) so that the position of each bump - pedestal,
1 p.e., 2 p.e., ... - is the *fitted* peak position rather than a raw bin
maximum.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import ROOT

from utils.fitter import fit_pe_spectrum
from utils.html_generator import generate_html
from utils.root_setup import setup_root
from utils.timing import auto_timer


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Multi-peak RooCrystalBall fits of Step 1's ADC-integral "
            "histograms (Step 2 of the ADC -> PE calibration)."
        )
    )
    parser.add_argument(
        "--step1-dir", default="results/PhotoElectrons/Step1_IntegralCalculation",
        help="Directory containing pe_integrals.root and pe_integrals_summary.json.",
    )
    parser.add_argument("--run", type=int, default=None,
                        help="Single run to fit, writing flat into --outdir. Default: the "
                        "only/first run in the Step 1 summary.")
    parser.add_argument(
        "--runs", type=int, nargs="+", default=None,
        help="Multiple runs to fit in one process (avoids re-paying ROOT startup per "
        "run). Each run's output goes into --outdir/run<N>/ instead of flat into "
        "--outdir. Overrides --run.",
    )
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument("--npe-max", type=int, default=8)
    parser.add_argument(
        "--min-entries", type=float, default=200.0,
        help="Skip channels with fewer histogram entries than this - too few "
        "events to resolve individual p.e. peaks.",
    )
    parser.add_argument(
        "--outdir", default="results/PhotoElectrons/Step2_PEFits"
    )
    return parser.parse_args()


def fit_run(
    run: int, root_file: Any, channels_meta: Dict[str, Any], outdir: Path,
    args: argparse.Namespace,
) -> None:
    run_key = "run{}".format(run)
    requested = set(args.channels) if args.channels else None

    outdir.mkdir(parents=True, exist_ok=True)
    fits_dir = outdir / "fits"
    fits_dir.mkdir(parents=True, exist_ok=True)

    summary: Dict[str, Any] = {}
    curves: Dict[str, Any] = {}
    png_files: List[str] = []

    branches = sorted(channels_meta)
    if requested is not None:
        branches = [b for b in branches if b in requested]

    print("[STEP2] Fitting {} channel(s) for run {}...".format(len(branches), run))
    for idx, branch in enumerate(branches):
        meta = channels_meta[branch]
        hist_name = meta.get("hist_name", "hist_pe_integral_run{}_{}".format(run, branch))
        hist = root_file.Get("{}/{}".format(run_key, hist_name))
        if not hist:
            print("[SKIP] {} - histogram {} not found.".format(branch, hist_name))
            continue
        hist = hist.Clone()
        hist.SetDirectory(0)

        if hist.GetEntries() < args.min_entries:
            print("[SKIP] {} - only {} entries (< {}).".format(
                branch, int(hist.GetEntries()), args.min_entries))
            summary[branch] = {
                "skipped": True,
                "reason": "too few entries",
                "entries": float(hist.GetEntries()),
                "material": meta.get("material"),
                "size": meta.get("size"),
                "short_label": meta.get("short_label"),
            }
            continue

        outname = branch
        title = "Run {} {}".format(run, meta.get("short_label", branch))
        print("[{}/{}] Fitting {}...".format(idx + 1, len(branches), branch))
        try:
            result = fit_pe_spectrum(
                hist, str(fits_dir), outname,
                npe_max=args.npe_max, run_number=run, title=title,
            )
        except Exception as exc:  # noqa: BLE001 - keep going on a bad channel
            print("[FAIL] {} - {}".format(branch, exc))
            summary[branch] = {
                "skipped": True,
                "reason": "fit failed: {}".format(exc),
                "entries": float(hist.GetEntries()),
                "material": meta.get("material"),
                "size": meta.get("size"),
                "short_label": meta.get("short_label"),
            }
            continue

        channel_curves = result.pop("curves")
        result["entries"] = float(hist.GetEntries())
        result["material"] = meta.get("material")
        result["size"] = meta.get("size")
        result["short_label"] = meta.get("short_label")
        result["png"] = str(Path(result["png"]).relative_to(outdir))
        summary[branch] = result
        curves[branch] = channel_curves
        png_files.append(result["png"])

    json_path = outdir / "pe_fits_summary.json"
    with json_path.open("w") as handle:
        json.dump(
            {"run": run, "npe_max": args.npe_max, "channels": summary},
            handle, indent=2, sort_keys=True,
        )
    curves_path = outdir / "pe_fit_curves.json"
    with curves_path.open("w") as handle:
        json.dump(curves, handle)

    html_path = outdir / "pe_fits.html"
    generate_html(
        png_files, str(outdir), plots_per_row=3, output_html=str(html_path),
        title="CaloX PE-calibration Step 2: photoelectron-peak fits",
        intro_text=(
            "**Run:** {}.\n".format(run)
            + "**Model:** multi-peak RooCrystalBall (pedestal + up to {} p.e. "
            "peaks), shared tail parameters, sigma scaling as sqrt(i+1) for "
            "photon-counting statistics - `utils.fitter.fit_pe_spectrum`. "
            "Initial pedestal/gain guesses are auto-detected per channel via "
            "TSpectrum peak-finding.\n".format(args.npe_max)
            + "Channels with fewer than {:.0f} histogram entries are "
            "skipped.".format(args.min_entries)
        ),
    )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(json_path))
    print("Wrote {}".format(curves_path))


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    step1_dir = Path(args.step1_dir)
    with (step1_dir / "pe_integrals_summary.json").open() as handle:
        step1_summary = json.load(handle)
    available_runs = [int(r.replace("run", "")) for r in step1_summary["channels"]]

    root_path = step1_dir / "pe_integrals.root"
    root_file = ROOT.TFile.Open(str(root_path))
    base_outdir = Path(args.outdir)

    if args.runs:
        for run in args.runs:
            run_key = "run{}".format(run)
            if run_key not in step1_summary["channels"]:
                print("[SKIP RUN] {} - not in Step 1 summary.".format(run))
                continue
            print("\n[RUN {}]".format(run))
            fit_run(
                run, root_file, step1_summary["channels"][run_key],
                base_outdir / run_key, args,
            )
    else:
        run = args.run
        if run is None:
            if len(available_runs) != 1:
                raise SystemExit(
                    "Step 1 summary has multiple runs {}; pass --run or --runs.".format(
                        available_runs)
                )
            run = available_runs[0]
        run_key = "run{}".format(run)
        fit_run(run, root_file, step1_summary["channels"][run_key], base_outdir, args)

    root_file.Close()


if __name__ == "__main__":
    main()
