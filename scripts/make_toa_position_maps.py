#!/usr/bin/env python3
"""Per-run 2D maps of the mean fitted TOA over the calorimeter face.

One page per scan run: x and y are the physical channel positions from the
tb2026 channel map, and z is that run's mean fitted TOA (the same per-run,
per-channel timing center used as the input point to the velocity fits in
plot_velocityfits.py / timing_summary.json's "timing_fits"). Pages are ordered
by table position, so scrolling the PDF is a slice-by-slice walk along the
scan axis, mirroring make_adc_position_maps.py but with TOA instead of ADC.

Cerenkov and scintillating fibers share tower positions, so a page carries one
pad per category (CerQuartz | CerPlastic | Sci) rather than overwriting them.
The z scale of each category is held fixed across all pages - otherwise every
page auto-scales to itself and the position-to-position trend, which is the
whole point, becomes invisible.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from channels.channel_map import build_drs_boards
from utils.root_setup import setup_root
from utils.visualization import (
    visualizeDRSBoards, get_drs_display_x, _X_SCALE, _Y_SCALE,
)

CATEGORIES = (
    ("CerQuartz", "cer", 1),
    ("CerPlastic", "cer", 2),
    ("Sci", "sci", 0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-run mean-TOA heat maps over the calorimeter face."
    )
    parser.add_argument(
        "--input",
        default=(
            "results/timing_velocity_smoke/runs_1897_1943_reference_mpv/"
            "timing_summary.json"
        ),
        help="timing_summary.json holding per-run, per-channel fitted TOA.",
    )
    parser.add_argument(
        "--outdir",
        default="results/timing_velocity_smoke/TimingHeatmaps",
    )
    parser.add_argument(
        "--x-min", type=float, default=None,
        help="Lowest table position to include, in mm. Default: no cut (all runs).",
    )
    parser.add_argument(
        "--x-max", type=float, default=None,
        help="Highest table position to include, in mm. Default: no cut (all runs).",
    )
    parser.add_argument(
        "--map-run", type=int, default=1897,
        help=(
            "Run whose channel map supplies the physical tower positions. Must "
            "be a tb2026 run (>= 1896) for the PHASE_2 layout. Default: 1897."
        ),
    )
    parser.add_argument(
        "--z-percentile", type=float, default=2.0,
        help=(
            "Percentile trimmed from each end when fixing a category's z range "
            "across all pages. Default: 2."
        ),
    )
    parser.add_argument(
        "--z-min", type=float, default=None,
        help="Override the common z minimum for every category.",
    )
    parser.add_argument(
        "--z-max", type=float, default=None,
        help="Override the common z maximum for every category.",
    )
    parser.add_argument(
        "--combined-pdf", action="store_true",
        help=(
            "Also write the 3-pads-per-page PDF with all families side by side. "
            "By default only the per-family PDFs are written."
        ),
    )
    parser.add_argument(
        "--no-split-pdfs", action="store_true",
        help="Skip the per-family PDFs (use with --combined-pdf).",
    )
    parser.add_argument(
        "--units", choices=("ns", "ps"), default="ns",
        help=(
            "Unit for the z axis. Picoseconds are the useful scale once "
            "--relative has removed the per-channel cable offsets: the "
            "position-dependent walk is ~53 ps per 10 mm step, against an "
            "absolute TOA spread of ~8.5 ns between channels."
        ),
    )
    parser.add_argument(
        "--relative", action="store_true",
        help=(
            "Plot each channel's mean TOA minus its average over the selected "
            "runs, instead of the raw mean. Useful when the absolute scale "
            "(dominated by cable-length offsets between channels) hides the "
            "position-to-position trend."
        ),
    )
    parser.add_argument(
        "--channel-size", choices=("all", "3mm", "6mm"), default="all",
        help=(
            "Restrict to one channel size instead of drawing both on top of "
            "each other. 3mm and 6mm channels occupy the same tower slots, so "
            "'all' (the default) is what the other scripts draw. Selecting "
            "'3mm' also zooms the display to the small central block the 3mm "
            "channels actually occupy (they only exist in a few central "
            "towers), instead of the full calorimeter face where they'd be a "
            "few pixels."
        ),
    )
    parser.add_argument(
        "--text", action="store_true",
        help=(
            "Print the numeric value in each cell. Off by default: with three "
            "pads on a page the 3 mm cells are too narrow for four digits and "
            "the labels overlap into noise."
        ),
    )
    parser.add_argument(
        "--palette", choices=("kBird", "kRainBow"), default="kBird",
        help="ROOT color palette for the z axis. Default: kBird.",
    )
    return parser.parse_args()


def load_points(path: Path) -> Tuple[Dict[int, float], Dict[int, Dict[str, float]]]:
    """Return (position by run, mean TOA [ns] by run and channel).

    Reads timing_summary.json's "run_positions_mm" (run -> mm) and
    "timing_fits" (branch -> run -> {mean_ns, fit_ok, ...}), the same
    per-run timing centers fed into the velocity fits.
    """
    with path.open() as handle:
        source = json.load(handle)
    positions = {int(run): float(x) for run, x in source["run_positions_mm"].items()}
    values: Dict[int, Dict[str, float]] = {}
    for branch, per_run in source["timing_fits"].items():
        for run_str, entry in per_run.items():
            if not entry.get("fit_ok"):
                continue
            mean = entry.get("mean_ns")
            if mean is None or not math.isfinite(mean):
                continue
            run = int(run_str)
            if run not in positions:
                continue
            values.setdefault(run, {})[branch] = float(mean)
    return positions, values


def z_range_for_category(
    values: Dict[int, Dict[str, float]],
    runs: List[int],
    members: set,
    args: argparse.Namespace,
) -> Tuple[float, float]:
    """One z range per category, held fixed across every page."""
    if args.z_min is not None and args.z_max is not None:
        return args.z_min, args.z_max
    pool = sorted(
        value
        for run in runs
        for branch, value in values.get(run, {}).items()
        if branch in members
    )
    if not pool:
        return 0.0, 1.0
    frac = max(0.0, min(0.45, args.z_percentile / 100.0))
    lo = pool[int(frac * (len(pool) - 1))]
    hi = pool[int((1.0 - frac) * (len(pool) - 1))]
    if hi <= lo:
        lo, hi = pool[0], max(pool[-1], pool[0] + 1.0)
    return (args.z_min if args.z_min is not None else lo,
            args.z_max if args.z_max is not None else hi)


def style_hist(hist: Any, xmin: float, xmax: float, zmin: float, zmax: float,
               ztitle: str, ymin: float = None, ymax: float = None) -> None:
    hist.SetStats(0)
    hist.SetTitle("")
    hist.GetXaxis().SetRangeUser(xmin, xmax)
    if ymin is not None and ymax is not None:
        hist.GetYaxis().SetRangeUser(ymin, ymax)
    hist.GetXaxis().SetTitle("X [cm]")
    hist.GetYaxis().SetTitle("Y [cm]")
    hist.GetZaxis().SetTitle(ztitle)
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetZaxis().SetTitleSize(0.040)
    hist.GetZaxis().SetTitleOffset(1.60)
    hist.GetXaxis().SetLabelSize(0.038)
    hist.GetYaxis().SetLabelSize(0.038)
    hist.GetZaxis().SetLabelSize(0.032)
    hist.SetMinimum(zmin)
    hist.SetMaximum(zmax)
    hist.SetMarkerSize(0.62)


def blank_empty_bins(hist: Any, zmin: float, zmax: float) -> None:
    """Push never-filled bins below the z minimum so ROOT leaves them white.

    COL/COLZ only skips a bin when its content falls under the pad minimum.
    With a positive minimum an unfilled bin (content 0) is already under it, so
    this never mattered - but in relative mode the minimum is negative, 0 sits
    inside the range, and every empty bin gets painted as a valid measurement.
    The map then reads as a full calorimeter face at "zero difference" instead
    of a handful of instrumented towers.
    """
    # The floor must stay below the pad minimum so COL skips the bin, but it
    # must also stay positive whenever the minimum is, or a log z axis cannot
    # represent it and paints the bin at the bottom of the palette instead of
    # leaving it blank - which floods the pad and, where a second histogram is
    # overlaid, hides the cells underneath it.
    floor = (0.5 * zmin) if zmin > 0.0 else (zmin - max(1.0, abs(zmax - zmin)))
    # Include the under/overflow rims: the display range is wider than the
    # histogram, so an unblanked x-underflow column paints as a stripe down the
    # left edge of the pad.
    for bx in range(0, hist.GetNbinsX() + 2):
        for by in range(0, hist.GetNbinsY() + 2):
            if hist.GetBinContent(bx, by) == 0.0:
                hist.SetBinContent(bx, by, floor)


def draw_category(hists: List[Any], with_text: bool) -> None:
    """Draw a category's histograms, putting the palette on a filled one.

    visualizeDRSBoards always returns a (6 mm, 3 mm) pair, and one of them can
    be empty - there are no Sci 6 mm channels in tb2026. Drawing COLZ on the
    empty member gives a pad with no colour bar at all.
    """
    suffix = " TEXT" if with_text else ""
    order = sorted(range(len(hists)), key=lambda i: -hists[i].GetEntries())
    first = True
    for index in order:
        hist = hists[index]
        if hist.GetEntries() <= 0 and not first:
            continue
        hist.Draw(("COLZ" + suffix) if first else ("COL" + suffix + " SAME"))
        first = False


def draw_pad_label(text: str, subtext: str) -> List[Any]:
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.045)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.12, 0.945, text)
    latex.SetTextFont(42)
    latex.SetTextSize(0.036)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.88, 0.945, subtext)
    return [latex]


def main() -> None:
    args = parse_args()
    setup_root(batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(getattr(ROOT, args.palette))
    ROOT.gStyle.SetPaintTextFormat("5.0f" if args.units == "ps" else "5.2f")

    positions, values = load_points(Path(args.input))
    if args.x_min is not None and args.x_max is not None:
        runs = sorted(
            (run for run, x in positions.items()
             if args.x_min <= x <= args.x_max and run in values),
            key=lambda run: positions[run],
        )
    else:
        runs = sorted((run for run in positions if run in values),
                      key=lambda run: positions[run])
    if not runs:
        raise RuntimeError("No runs with valid TOA fits found in {}".format(args.input))
    print("Selected {} run(s), X = {:.1f} .. {:.1f} mm:".format(
        len(runs), positions[runs[0]], positions[runs[-1]]))
    for run in runs:
        print("   run {}  X = {:.1f} mm  ({} channels)".format(
            run, positions[run], len(values.get(run, {}))))

    if args.relative:
        # Per-channel average over the selected runs, used as the baseline.
        totals: Dict[str, List[float]] = {}
        for run in runs:
            for branch, value in values.get(run, {}).items():
                totals.setdefault(branch, []).append(value)
        reference = {branch: sum(vals) / len(vals) for branch, vals in totals.items()}
        values = {
            run: {
                branch: value - reference[branch]
                for branch, value in per_run.items()
                if branch in reference
            }
            for run, per_run in values.items()
        }
        print("Relative mode: z = mean TOA minus that channel's average over the "
              "{} selected runs.".format(len(runs)))

    unit_scale = 1000.0 if args.units == "ps" else 1.0
    if unit_scale != 1.0:
        values = {
            run: {branch: value * unit_scale for branch, value in per_run.items()}
            for run, per_run in values.items()
        }
        print("Units: ps (values scaled by {:.0f}).".format(unit_scale))

    drsboards = build_drs_boards(run_number=args.map_run)
    xmin, xmax, _ = get_drs_display_x(args.map_run)

    # Which branches belong to which category, so each gets its own fixed scale.
    # When zooming to one channel size, restrict membership (and therefore the
    # z range and per-page channel counts) to that size too.
    members: Dict[str, set] = {name: set() for name, _, _ in CATEGORIES}
    zoom_ix: List[float] = []
    zoom_iy: List[float] = []
    for _, board in drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue
            if args.channel_size == "3mm" and channel.is6mm:
                continue
            if args.channel_size == "6mm" and not channel.is6mm:
                continue
            branch = channel.get_channel_name(blsub=False)
            if not channel.isCer:
                members["Sci"].add(branch)
            elif channel.isQuartz:
                members["CerQuartz"].add(branch)
            else:
                members["CerPlastic"].add(branch)
            zoom_ix.append(channel.i_tower_x)
            zoom_iy.append(channel.i_tower_y)

    ymin_display: float = None
    ymax_display: float = None
    if args.channel_size == "3mm" and zoom_ix:
        # 3mm channels sit in a handful of central towers; zoom the display
        # to that block (plus a small margin) instead of the full face where
        # they'd be a sliver of pixels. Cell half-widths: 6mm/3mm share the X
        # binning, but 3mm rows are 4x finer in Y (see visualizeDRSBoards).
        margin = 0.3
        xmin = min(zoom_ix) * _X_SCALE - _X_SCALE / 2 - margin
        xmax = max(zoom_ix) * _X_SCALE + _X_SCALE / 2 + margin
        ymin_display = min(zoom_iy) * _Y_SCALE - (_Y_SCALE / 4) / 2 - margin
        ymax_display = max(zoom_iy) * _Y_SCALE + (_Y_SCALE / 4) / 2 + margin
        print("3mm zoom: X = {:.2f} .. {:.2f} cm, Y = {:.2f} .. {:.2f} cm".format(
            xmin, xmax, ymin_display, ymax_display))

    zranges = {
        name: z_range_for_category(values, runs, members[name], args)
        for name, _, _ in CATEGORIES
    }
    for name, (lo, hi) in sorted(zranges.items()):
        print("z range {:<11} {:9.2f} .. {:9.2f} {} (fixed across all pages)".format(
            name, lo, hi, args.units))

    ztitle = (
        "Mean TOA - channel average [{}]".format(args.units)
        if args.relative else "Mean fitted TOA [{}]".format(args.units)
    )
    stem = "toa_position_maps_relative" if args.relative else "toa_position_maps"
    if args.channel_size != "all":
        stem += "_" + args.channel_size
    if args.units != "ns":
        stem += "_" + args.units

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    root_path = outdir / "{}.root".format(stem)
    root_file = ROOT.TFile(str(root_path), "RECREATE")

    do_split = not args.no_split_pdfs
    do_combined = args.combined_pdf
    if not do_split and not do_combined:
        raise RuntimeError("Nothing to write: --no-split-pdfs without --combined-pdf")

    def style_pad(pad: Any) -> None:
        pad.SetLeftMargin(0.135)
        pad.SetRightMargin(0.215)
        pad.SetTopMargin(0.11)
        pad.SetBottomMargin(0.13)
        pad.SetTicks(1, 1)

    combined_pdf = outdir / "{}_all_families.pdf".format(stem)
    canvas = None
    if do_combined:
        canvas = ROOT.TCanvas("toa_maps", "", 1950, 700)
        canvas.Divide(3, 1)
        canvas.Print(str(combined_pdf) + "[")

    # One canvas and one PDF per family, each with a page per scan run.
    family_canvas: Dict[str, Any] = {}
    family_pdf: Dict[str, Path] = {}
    if do_split:
        for name, _, _ in CATEGORIES:
            single = ROOT.TCanvas("toa_map_{}".format(name), "", 900, 780)
            style_pad(single)
            path = outdir / "{}_{}.pdf".format(stem, name)
            single.Print(str(path) + "[")
            family_canvas[name] = single
            family_pdf[name] = path

    page_summary: List[Dict[str, Any]] = []
    for run in runs:
        value_map = values.get(run, {})
        keep: List[Any] = []
        entry = {"run": run, "position_mm": positions[run], "categories": {}}
        subtitle = "run {}   X = {:.1f} mm".format(run, positions[run])

        for index, (name, kind, quartz_only) in enumerate(CATEGORIES, start=1):
            cer_hists, sci_hists = visualizeDRSBoards(
                drsboards, valuemaps=value_map,
                suffix="toa_run{}_{}".format(run, name),
                quartzOnly=quartz_only,
            )
            hists = cer_hists if kind == "cer" else sci_hists
            if args.channel_size == "3mm":
                hists = [hists[1]]
            elif args.channel_size == "6mm":
                hists = [hists[0]]
            zmin, zmax = zranges[name]
            for hist in hists:
                style_hist(hist, xmin, xmax, zmin, zmax, ztitle,
                          ymin=ymin_display, ymax=ymax_display)
                blank_empty_bins(hist, zmin, zmax)
            keep.extend(hists)

            if do_split:
                single = family_canvas[name]
                single.cd()
                single.Clear()
                style_pad(single)
                draw_category(hists, args.text)
                keep.extend(draw_pad_label(name, subtitle))
                single.Update()
                single.Print(str(family_pdf[name]))
                root_file.cd()
                single.Write("toa_map_{}_run{}".format(name, run))

            if do_combined:
                pad = canvas.cd(index)
                style_pad(pad)
                draw_category(hists, args.text)
                keep.extend(draw_pad_label(name, subtitle))

            filled = sum(1 for b in value_map if b in members[name])
            in_cat = [value_map[b] for b in value_map if b in members[name]]
            entry["categories"][name] = {
                "n_channels": filled,
                "mean_toa_ns": sum(in_cat) / len(in_cat) if in_cat else None,
                "z_range": [zmin, zmax],
            }

        if do_combined:
            canvas.cd()
            canvas.Update()
            canvas.Print(str(combined_pdf))
            root_file.cd()
            canvas.Write("toa_maps_run{}".format(run))
            canvas._keepalive = keep
        page_summary.append(entry)
        print("  page {}/{}: run {} at X = {:.1f} mm".format(
            len(page_summary), len(runs), run, positions[run]))

    if do_combined:
        canvas.Print(str(combined_pdf) + "]")
    for name, single in family_canvas.items():
        single.Print(str(family_pdf[name]) + "]")
    root_file.Close()

    summary = {
        "input": args.input,
        "map_run": args.map_run,
        "position_window_mm": (
            None if args.x_min is None else [args.x_min, args.x_max]
        ),
        "z_variable": ztitle,
        "relative": args.relative,
        "channel_size": args.channel_size,
        "units": args.units,
        "z_ranges": {name: list(rng) for name, rng in zranges.items()},
        "n_pages": len(page_summary),
        "per_family_pdfs": {name: str(path) for name, path in family_pdf.items()},
        "pages": page_summary,
    }
    with (outdir / "{}_summary.json".format(stem)).open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("")
    for name in sorted(family_pdf):
        print("Wrote {} ({} pages, one per scan run)".format(
            family_pdf[name], len(runs)))
    if do_combined:
        print("Wrote {} ({} pages, all families side by side)".format(
            combined_pdf, len(runs)))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(outdir / "{}_summary.json".format(stem)))


if __name__ == "__main__":
    main()
