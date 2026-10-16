#!/usr/bin/env python3
"""Build per-channel position-scan overlays from existing Step 1 histograms.

No event loop or ADC-integral calculation is performed. The script reads the
histograms already stored by ``compute_pe_integrals.py`` in ``pe_integrals.root``
and writes a new interactive HTML page with one unfilled overlay per channel.

Each requested target position is resolved to the closest run present in both
the Step 1 summary and the run list. By default, ``-200 87 213 386 809`` mm
selects runs 1907, 1933, 1934, 1936, and 1940.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT

from utils.html_generator import generate_jsroot_html


DEFAULT_INPUT_DIR = Path(
    "results/PhotoElectrons/Step6_AttenuationCorrections/"
    "Step1_IntegralCalculation"
)

ROOT_COLORS = [
    ROOT.kBlue + 1,
    ROOT.kRed + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kCyan + 2,
    ROOT.kViolet + 1,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Overlay existing Step 1 ADC-integral histograms for runs nearest "
            "the requested table positions."
        )
    )
    parser.add_argument(
        "--input-dir", type=Path, default=DEFAULT_INPUT_DIR,
        help=(
            "Directory containing pe_integrals.root and "
            "pe_integrals_summary.json."
        ),
    )
    parser.add_argument(
        "--input-root", type=Path, default=None,
        help="Override the input ROOT file (default: <input-dir>/pe_integrals.root).",
    )
    parser.add_argument(
        "--input-summary", type=Path, default=None,
        help=(
            "Override the input summary "
            "(default: <input-dir>/pe_integrals_summary.json)."
        ),
    )
    parser.add_argument(
        "--outdir", type=Path, default=None,
        help="Output directory (default: the Step 1 input directory).",
    )
    parser.add_argument(
        "--target-positions", nargs="+", type=float,
        default=[-200.0, 87.0, 213.0, 386.0, 809.0],
        help="Target table X positions in mm (default: -200 87 213 386 809).",
    )
    parser.add_argument("--runlist", type=Path, default=Path("data/Runlist.json"))
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument(
        "--channels", nargs="+", default=None,
        help="Optional exact channel branch names to include.",
    )
    parser.add_argument(
        "--raw-yields", action="store_true",
        help="Draw raw event counts instead of unit-area-normalized spectra.",
    )
    parser.add_argument(
        "--zoom-keep", type=float, default=0.995,
        help="Central fraction of each spectrum used to choose the shared x range.",
    )
    parser.add_argument("--plots-per-row", type=int, default=3)
    args = parser.parse_args()

    if not args.target_positions:
        parser.error("--target-positions needs at least one value")
    if not 0.0 < args.zoom_keep <= 1.0:
        parser.error("--zoom-keep must be in (0, 1]")
    if args.plots_per_row < 1:
        parser.error("--plots-per-row must be at least 1")
    return args


def parse_length_mm(value: Any) -> float:
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    if match is None:
        raise ValueError("cannot parse a numeric position from {!r}".format(value))
    return float(match.group(0))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")


def load_json(path: Path) -> Any:
    if not path.is_file():
        raise SystemExit("Input file does not exist: {}".format(path))
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def resolve_target_runs(
    targets: Sequence[float], runlist_path: Path, position_key: str,
    available_runs: Iterable[int],
) -> List[Tuple[float, int, float]]:
    """Return ``(target, run, actual_position)`` in the requested order."""
    runlist = load_json(runlist_path)
    candidates: List[Tuple[int, float]] = []
    for run in sorted(set(available_runs)):
        record = runlist.get(str(run))
        if not isinstance(record, dict) or position_key not in record:
            continue
        try:
            candidates.append((run, parse_length_mm(record[position_key])))
        except ValueError:
            continue

    if not candidates:
        raise SystemExit(
            "No Step 1 runs have a numeric {!r} in {}.".format(
                position_key, runlist_path
            )
        )

    resolved: List[Tuple[float, int, float]] = []
    used_runs = set()
    for target in targets:
        remaining = [entry for entry in candidates if entry[0] not in used_runs]
        if not remaining:
            raise SystemExit("Not enough distinct runs to match all target positions.")
        run, actual = min(
            remaining, key=lambda entry: (abs(entry[1] - target), entry[0])
        )
        used_runs.add(run)
        resolved.append((target, run, actual))
        print(
            "[POSITION] target X={:g} mm -> run {} "
            "(X={:g} mm, difference={:+.1f} mm)".format(
                target, run, actual, actual - target
            )
        )
    return resolved


def content_range(hist: Any, keep: float) -> Optional[Tuple[float, float]]:
    """Central content range, including a small visual margin."""
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
    margin = 0.06 * (hi - lo)
    return lo - margin, hi + margin


def draw_header(canvas: Any, title: str, subtitle: str) -> Any:
    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.040)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.14, 0.944, "CaloX Data")
    latex.SetTextFont(42)
    latex.SetTextSize(0.026)
    latex.DrawLatex(0.14, 0.895, title)
    latex.SetTextSize(0.024)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.94, 0.944, subtitle)
    return latex


def make_overlay_canvas(
    branch: str, short_label: str, selection: str,
    run_positions: Sequence[Tuple[float, int, float]],
    source_hists: Dict[int, Any], normalise: bool, zoom_keep: float,
    colors: Optional[Dict[int, int]] = None,
) -> Optional[Tuple[Any, str]]:
    plotted = []
    for target, run, actual in run_positions:
        source = source_hists.get(run)
        if source is None or source.GetEntries() <= 0:
            continue
        hist = source.Clone("plot_run{}_{}".format(run, safe_name(branch)))
        hist.SetDirectory(0)
        if normalise and hist.Integral() > 0.0:
            hist.Scale(1.0 / hist.Integral())
        plotted.append((target, run, actual, hist))
    if not plotted:
        return None

    ranges = [content_range(hist, zoom_keep) for _, _, _, hist in plotted]
    ranges = [item for item in ranges if item is not None]
    if ranges:
        xmin = min(item[0] for item in ranges)
        xmax = max(item[1] for item in ranges)
    else:
        xmin = min(hist.GetXaxis().GetXmin() for _, _, _, hist in plotted)
        xmax = max(hist.GetXaxis().GetXmax() for _, _, _, hist in plotted)

    ymax = max(hist.GetMaximum() for _, _, _, hist in plotted)
    ymin = 1.0e-5 if normalise else 0.5
    canvas_name = "pe_integral_position_overlay_{}".format(safe_name(branch))
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

    legend = ROOT.TLegend(0.47, 0.64, 0.95, 0.885)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.023)

    keepalive: List[Any] = [frame, legend]
    for index, (target, run, actual, hist) in enumerate(plotted):
        hist.SetLineColor(
            colors[run] if colors and run in colors
            else ROOT_COLORS[index % len(ROOT_COLORS)]
        )
        hist.SetLineWidth(2)
        hist.SetFillStyle(0)
        hist.SetStats(0)
        hist.SetTitle("")
        hist.Draw("HIST SAME")
        legend.AddEntry(
            hist,
            "Run {}: X={:.1f} mm (target {:g}), N={}".format(
                run, actual, target, int(hist.GetEntries())
            ),
            "l",
        )
        keepalive.append(hist)

    legend.Draw()
    header = draw_header(canvas, short_label, selection)
    keepalive.append(header)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def main() -> None:
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    input_root = args.input_root or args.input_dir / "pe_integrals.root"
    input_summary = args.input_summary or args.input_dir / "pe_integrals_summary.json"
    outdir = args.outdir or args.input_dir
    outdir.mkdir(parents=True, exist_ok=True)

    summary = load_json(input_summary)
    channel_runs = summary.get("channels", {})
    available_runs = [
        int(name[3:]) for name in channel_runs
        if name.startswith("run") and name[3:].isdigit()
    ]
    resolved = resolve_target_runs(
        args.target_positions, args.runlist, args.position_key, available_runs
    )

    empty = [entry for entry in resolved if not channel_runs.get("run{}".format(entry[1]))]
    for target, run, _ in empty:
        print("[WARN] run {} (target X={:g} mm) has no channels in Step 1; "
              "dropped from the overlay.".format(run, target))
    resolved = [entry for entry in resolved if entry not in empty]

    source = ROOT.TFile.Open(str(input_root), "READ")
    if not source or source.IsZombie():
        raise SystemExit("Cannot open input ROOT file: {}".format(input_root))

    selected_runs = [run for _, run, _ in resolved]
    branch_sets = [set(channel_runs.get("run{}".format(run), {})) for run in selected_runs]
    branches = sorted(set.intersection(*branch_sets)) if branch_sets else []
    if args.channels:
        requested = set(args.channels)
        missing = sorted(requested - set(branches))
        if missing:
            print("[WARN] Requested channel(s) unavailable: {}".format(", ".join(missing)))
        branches = [branch for branch in branches if branch in requested]
    if not branches:
        source.Close()
        raise SystemExit("No channels are common to all selected runs.")

    selection = str(summary.get("selection", "Step 1 selection"))
    normalise = not args.raw_yields
    keys: List[str] = []
    jsons: Dict[str, str] = {}
    skipped: List[str] = []

    print("[INPUT] {}".format(input_root))
    print("[CHANNELS] Building {} overlay canvas(es).".format(len(branches)))
    for index, branch in enumerate(branches, 1):
        source_hists: Dict[int, Any] = {}
        short_label = branch
        for run in selected_runs:
            metadata = channel_runs["run{}".format(run)][branch]
            short_label = str(metadata.get("short_label", short_label))
            hist_name = metadata.get("hist_name")
            hist = source.Get("run{}/{}".format(run, hist_name)) if hist_name else None
            if hist:
                source_hists[run] = hist
        built = make_overlay_canvas(
            branch, short_label, selection, resolved, source_hists,
            normalise, args.zoom_keep,
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

    html_path = outdir / "pe_integrals_position_overlay.html"
    run_text = ", ".join(
        "{:g} mm -> run {} ({:g} mm)".format(target, run, actual)
        for target, run, actual in resolved
    )
    if summary.get("integral_definition") == "compute_cfd_integral":
        observable = (
            "**Observable:** `{ch}_energy` from `compute_cfd_integral` "
            "(the repo's standard fixed-window integral)."
        )
    else:
        observable = (
            "**Observable:** `{ch}_energy_adaptive` from "
            "`compute_cfd_to_baseline_integral`."
        )
    intro = "\n".join([
        "**Runs (target X -> resolved run, actual X):** " + run_text + ".",
        "**Source:** Existing Step 1 histograms from `{}`; no integrals were recomputed.".format(
            input_root
        ),
        observable,
        "**Selection:** {}.".format(selection),
        "The spectra are {} and drawn as unfilled outlines.".format(
            "shown with raw event counts" if args.raw_yields else "unit-area normalized"
        ),
    ])
    generate_jsroot_html(
        keys, jsons, plots_per_row=args.plots_per_row,
        output_html=str(html_path),
        title="CaloX PE-calibration: ADC integral position overlays",
        intro_text=intro,
    )

    output_summary = {
        "source_root": str(input_root),
        "source_summary": str(input_summary),
        "target_positions_mm": args.target_positions,
        "position_key": args.position_key,
        "resolved_runs": [
            {"target_position_mm": target, "run": run, "actual_position_mm": actual}
            for target, run, actual in resolved
        ],
        "selection": selection,
        "normalized": normalise,
        "channels": branches,
        "skipped_channels": skipped,
    }
    summary_path = outdir / "pe_integrals_position_overlay_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(output_summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    print("[DONE] Interactive HTML: {}".format(html_path))
    print("[DONE] JSON summary: {}".format(summary_path))


if __name__ == "__main__":
    main()
