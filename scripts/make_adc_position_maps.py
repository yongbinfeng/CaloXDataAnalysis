#!/usr/bin/env python3
"""Per-run 2D maps of the mean ADC integral over the calorimeter face.

One page per scan position: x and y are the physical channel positions from the
tb2026 channel map, and z is the mean ADC integral for that run. Pages are
ordered by table position, so scrolling the PDF is a slice-by-slice walk along
the scan axis.

Default window is the band the attenuation exponential fit throws away
(-166 <= X <= -85 mm, runs 1910-1918): those points are drawn as open grey
markers in the fit plots and are otherwise never looked at.

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
from typing import Any, Dict, List, Optional, Tuple

import ROOT

from channels.channel_map import build_drs_boards
from utils.root_setup import setup_root
from utils.visualization import visualizeDRSBoards, get_drs_display_x

CATEGORIES = (
    ("CerQuartz", "cer", 1),
    ("CerPlastic", "cer", 2),
    ("Sci", "sci", 0),
)

# Quartz and plastic Cerenkov fibers never share a tower - checked against the
# tb2026 map: 96 quartz and 60 plastic cells, zero collisions - so they can be
# drawn on one pad without overwriting each other. visualizeDRSBoards already
# supports this: quartzOnly=0 returns every Cerenkov channel.
CATEGORIES_MERGED_CER = (
    ("Cer", "cer", 0),
    ("Sci", "sci", 0),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-run mean-ADC heat maps over the calorimeter face."
    )
    parser.add_argument(
        "--input",
        default=(
            "results/attenuation/runs_1897_1942_no_pedestal/"
            "attenuation_summary.json"
        ),
        help="attenuation_summary.json holding per-run, per-channel mean ADC.",
    )
    parser.add_argument(
        "--observable", choices=("auto", "adc", "pe"), default="auto",
        help=(
            "Quantity stored in the input points. 'auto' selects mean_pe when "
            "present, otherwise mean_adc (default: auto)."
        ),
    )
    parser.add_argument(
        "--outdir",
        default="results/attenuation/runs_1897_1942_no_pedestal/adc_position_maps",
    )
    parser.add_argument(
        "--x-min", type=float, default=-166.0,
        help="Lowest table position X to include, in mm. Default: -166.",
    )
    parser.add_argument(
        "--x-max", type=float, default=-85.0,
        help="Highest table position X to include, in mm. Default: -85.",
    )
    parser.add_argument(
        "--all-positions", action="store_true",
        help="Ignore --x-min/--x-max and make a page for every run in the file.",
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
        "--relative", action="store_true",
        help=(
            "Plot each channel's mean ADC divided by its average over the "
            "selected runs, instead of the raw mean. Across this window the "
            "run-to-run change is only a few percent while the absolute scale "
            "spans a factor ~25, so on the raw map every page looks identical; "
            "the ratio is what makes the position dependence visible."
        ),
    )
    parser.add_argument(
        "--merge-cer", action="store_true",
        help=(
            "Draw quartz and plastic Cerenkov channels on one pad instead of "
            "two. Their tower positions never collide, so nothing is "
            "overwritten and the whole Cerenkov face is visible at once."
        ),
    )
    parser.add_argument(
        "--logz", action="store_true",
        help=(
            "Log z scale. The absolute ADC ranges span a factor of 26 "
            "(CerQuartz, CerPlastic) to 69 (Sci), so on a linear scale the few "
            "brightest towers use most of the colour range and everything else "
            "is compressed into the bottom fifth. Log z spreads the colour over "
            "the whole range, which is what makes small tower-to-tower and "
            "position-to-position differences visible. Written with a _logz "
            "suffix so the linear maps are kept."
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
    return parser.parse_args()


def load_points(
    path: Path, observable: str = "auto",
) -> Tuple[Dict[int, float], Dict[int, Dict[str, float]], str]:
    """Return position/value maps and the resolved ADC-or-PE observable."""
    with path.open() as handle:
        source = json.load(handle)
    first_point = next(
        (
            point
            for entries in source.get("points", {}).values()
            for point in entries[:1]
        ),
        {},
    )
    if observable == "auto":
        observable = "pe" if "mean_pe" in first_point else "adc"
    value_key = "mean_pe" if observable == "pe" else "mean_adc"
    positions: Dict[int, float] = {}
    values: Dict[int, Dict[str, float]] = {}
    for branch, entries in source["points"].items():
        for entry in entries:
            run = int(entry["run"])
            mean = float(entry[value_key])
            positions.setdefault(run, float(entry["position_mm"]))
            if math.isfinite(mean):
                values.setdefault(run, {})[branch] = mean
    return positions, values, observable


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
               ztitle: str) -> None:
    hist.SetStats(0)
    hist.SetTitle("")
    hist.GetXaxis().SetRangeUser(xmin, xmax)
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
    ROOT.gStyle.SetPalette(ROOT.kBird)

    categories = CATEGORIES_MERGED_CER if args.merge_cer else CATEGORIES

    positions, values, args.observable = load_points(
        Path(args.input), args.observable
    )
    is_pe = args.observable == "pe"
    stem_prefix = "pe" if is_pe else "adc"
    unit_label = "P.E." if is_pe else "ADC"
    ROOT.gStyle.SetPaintTextFormat(
        "4.2f" if args.relative or is_pe else "4.0f"
    )
    if args.all_positions:
        runs = sorted(positions, key=lambda run: positions[run])
    else:
        runs = sorted(
            (run for run, x in positions.items() if args.x_min <= x <= args.x_max),
            key=lambda run: positions[run],
        )
    if not runs:
        raise RuntimeError(
            "No runs with {:.1f} <= X <= {:.1f} mm in {}".format(
                args.x_min, args.x_max, args.input
            )
        )
    print("Selected {} run(s) between X = {:.1f} and {:.1f} mm:".format(
        len(runs), positions[runs[0]], positions[runs[-1]]))
    for run in runs:
        print("   run {}  X = {:.1f} mm  ({} channels)".format(
            run, positions[run], len(values.get(run, {}))))

    if args.relative:
        # Per-channel average over the selected runs, used as the denominator.
        totals: Dict[str, List[float]] = {}
        for run in runs:
            for branch, value in values.get(run, {}).items():
                totals.setdefault(branch, []).append(value)
        reference = {
            branch: sum(vals) / len(vals)
            for branch, vals in totals.items()
            if vals and abs(sum(vals) / len(vals)) > 1.0e-9
        }
        values = {
            run: {
                branch: value / reference[branch]
                for branch, value in per_run.items()
                if branch in reference
            }
            for run, per_run in values.items()
        }
        print("Relative mode: z = mean {} / that channel's average over the "
              "{} selected runs.".format(unit_label, len(runs)))

    drsboards = build_drs_boards(run_number=args.map_run)
    xmin, xmax, _ = get_drs_display_x(args.map_run)

    # Which branches belong to which category, so each gets its own fixed scale.
    members: Dict[str, set] = {name: set() for name, _, _ in categories}
    for _, board in drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue
            branch = channel.get_channel_name(blsub=False)
            if not channel.isCer:
                members["Sci"].add(branch)
            elif args.merge_cer:
                members["Cer"].add(branch)
            elif channel.isQuartz:
                members["CerQuartz"].add(branch)
            else:
                members["CerPlastic"].add(branch)

    zranges = {
        name: z_range_for_category(values, runs, members[name], args)
        for name, _, _ in categories
    }
    if args.logz:
        # A log axis cannot start at or below zero; fall back to three decades
        # under the top of the range when the percentile floor is non-positive.
        for name, (lo, hi) in list(zranges.items()):
            if lo <= 0.0:
                pool = [
                    v for run in runs
                    for b, v in values.get(run, {}).items()
                    if b in members[name] and v > 0.0
                ]
                lo = min(pool) if pool else max(hi, 1.0) / 1000.0
            zranges[name] = (max(lo, hi / 1.0e4), hi)
    for name, (lo, hi) in sorted(zranges.items()):
        print("z range {:<11} {:9.1f} .. {:9.1f} {} (fixed across all pages)".format(
            name, lo, hi, unit_label))

    ztitle = (
        "Mean {} / channel average".format(unit_label)
        if args.relative else (
            "Mean photoelectrons [P.E.]" if is_pe
            else "Mean ADC integral [ADC counts]"
        )
    )
    stem = "{}_position_maps_relative".format(stem_prefix) if args.relative else (
        "{}_position_maps".format(stem_prefix)
    )
    if args.merge_cer:
        stem += "_mergedcer"
    if args.logz:
        stem += "_logz"

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    root_path = outdir / "{}.root".format(stem)
    root_file = ROOT.TFile(str(root_path), "RECREATE")

    do_split = not args.no_split_pdfs
    do_combined = args.combined_pdf
    if not do_split and not do_combined:
        raise RuntimeError("Nothing to write: --no-split-pdfs without --combined-pdf")

    def style_pad(pad: Any) -> None:
        pad.SetLogz(1 if args.logz else 0)
        pad.SetLeftMargin(0.135)
        pad.SetRightMargin(0.215)
        pad.SetTopMargin(0.11)
        pad.SetBottomMargin(0.13)
        pad.SetTicks(1, 1)

    combined_pdf = outdir / "{}_all_families.pdf".format(stem)
    canvas = None
    if do_combined:
        canvas = ROOT.TCanvas("{}_maps".format(stem_prefix), "", 1950, 700)
        canvas.Divide(len(categories), 1)
        canvas.Print(str(combined_pdf) + "[")

    # One canvas and one PDF per family, each with a page per scan position.
    family_canvas: Dict[str, Any] = {}
    family_pdf: Dict[str, Path] = {}
    if do_split:
        for name, _, _ in categories:
            single = ROOT.TCanvas(
                "{}_map_{}".format(stem_prefix, name), "", 900, 780
            )
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

        for index, (name, kind, quartz_only) in enumerate(categories, start=1):
            cer_hists, sci_hists = visualizeDRSBoards(
                drsboards, valuemaps=value_map,
                suffix="{}_run{}_{}".format(stem_prefix, run, name),
                quartzOnly=quartz_only,
            )
            hists = cer_hists if kind == "cer" else sci_hists
            zmin, zmax = zranges[name]
            for hist in hists:
                style_hist(hist, xmin, xmax, zmin, zmax, ztitle)
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
                single.Write("{}_map_{}_run{}".format(stem_prefix, name, run))

            if do_combined:
                pad = canvas.cd(index)
                style_pad(pad)
                draw_category(hists, args.text)
                keep.extend(draw_pad_label(name, subtitle))

            filled = sum(1 for b in value_map if b in members[name])
            in_cat = [value_map[b] for b in value_map if b in members[name]]
            entry["categories"][name] = {
                "n_channels": filled,
                "mean_{}".format(args.observable): (
                    sum(in_cat) / len(in_cat) if in_cat else None
                ),
                "z_range": [zmin, zmax],
            }

        if do_combined:
            canvas.cd()
            canvas.Update()
            canvas.Print(str(combined_pdf))
            root_file.cd()
            canvas.Write("{}_maps_run{}".format(stem_prefix, run))
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
        "observable": args.observable,
        "map_run": args.map_run,
        "position_window_mm": (
            None if args.all_positions else [args.x_min, args.x_max]
        ),
        "z_variable": ztitle,
        "relative": args.relative,
        "logz": args.logz,
        "z_ranges_{}".format(args.observable): {
            name: list(rng) for name, rng in zranges.items()
        },
        "n_pages": len(page_summary),
        "per_family_pdfs": {name: str(path) for name, path in family_pdf.items()},
        "pages": page_summary,
    }
    with (outdir / "{}_summary.json".format(stem)).open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print("")
    for name in sorted(family_pdf):
        print("Wrote {} ({} pages, one per scan position)".format(
            family_pdf[name], len(runs)))
    if do_combined:
        print("Wrote {} ({} pages, all families side by side)".format(
            combined_pdf, len(runs)))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(outdir / "{}_summary.json".format(stem)))


if __name__ == "__main__":
    main()
