#!/usr/bin/env python3
"""Per-event 2D TOA heatmaps over the calorimeter face - same rendering as
make_toa_position_maps.py (one ROOT COLZ page per family, z scale fixed
across all pages, extreme values clipped into the ends of the colour scale
rather than dropped), but paged by event instead of by scan run.

z is a genuine, non-negative time of arrival, on a scale comparable to
make_toa_position_maps.py's maps (a handful to a few tens of ns), not raw
_Time_cfd_mcp: for this run's channel population that quantity runs both
negative and out past +170 ns, because its group_mpv calibration (built for
the attenuation-scan runs) leaves ~100+ ns of per-channel cable-length offset
uncorrected here - which would swamp any real shower-timing structure. Each
channel's own average over this run's other clean hits is subtracted first
(pass --absolute to skip that and see the raw quantity instead), then the
whole page set is shifted so the earliest kept hit sits at 0.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import ROOT

from channels.channel_map import build_drs_boards
from core.analysis_manager import CaloXAnalysisManager
from utils.root_setup import setup_root
from utils.scan3d import compute_toa_references
from utils.visualization import visualizeDRSBoards, get_drs_display_x

from plot_velocityfits import active_timing_channels, make_manager_args

CATEGORIES = (
    ("CerQuartz", "cer", 1),
    ("CerPlastic", "cer", 2),
    ("Sci", "sci", 0),
)

DEFAULT_OUTDIR = "results/ShoweShape/EventTOAHeatmaps"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-event mean-TOA heat maps over the calorimeter face."
    )
    parser.add_argument(
        "--run", type=int, default=1831,
        help="tb2026 run (>= 1828 era). Default: 1831 (120 GeV electron).",
    )
    parser.add_argument("--particle", default=None)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=20000)
    parser.add_argument(
        "--events", type=int, nargs="+", default=None,
        help="Explicit event_n list. Default: the --n-events highest-ADC events.",
    )
    parser.add_argument(
        "--n-events", type=int, default=20,
        help="Number of highest-total-ADC events to page through. Default: 20.",
    )
    parser.add_argument(
        "--min-adc", type=float, default=15.0,
        help="Drop a channel from an event's page if its ADC is below this "
             "(channels with no real pulse have a meaningless CFD crossing).",
    )
    parser.add_argument(
        "--absolute", action="store_true",
        help=(
            "Use the raw fitted _Time_cfd_mcp instead of subtracting each "
            "channel's own run average first. This run's channels carry "
            "~100+ ns of uncorrected cable-offset spread, so the raw value "
            "swamps real shower timing - only use this to inspect that "
            "offset itself, not to compare timing across the face."
        ),
    )
    parser.add_argument(
        "--min-ref-hits", type=int, default=5,
        help="Without --absolute: minimum same-run, above-threshold hits to "
             "trust a channel's reference TOA. Ignored with --absolute.",
    )
    parser.add_argument(
        "--toa-window-ns", type=float, default=10.0,
        help=(
            "Without --absolute: drop a channel from an event's page if "
            "|TOA - its run average| exceeds this many ns. On a typical "
            "event roughly a fifth of channels have a badly-fit CFD crossing "
            "tens of ns off, even above the ADC cut; that is far more than "
            "the z-range percentile trim absorbs, so without this cut those "
            "bad fits dominate the axis range on every page. Ignored with "
            "--absolute."
        ),
    )
    parser.add_argument(
        "--map-run", type=int, default=None,
        help="Run whose channel map supplies tower positions. Default: --run.",
    )
    parser.add_argument(
        "--z-percentile", type=float, default=2.0,
        help="Percentile trimmed from each end when fixing a category's z "
             "range across all pages. Default: 2.",
    )
    parser.add_argument("--z-min", type=float, default=None)
    parser.add_argument("--z-max", type=float, default=None)
    parser.add_argument(
        "--units", choices=("ns", "ps"), default="ns",
        help="Unit for the z axis. Default: ns.",
    )
    parser.add_argument(
        "--combined-pdf", action="store_true",
        help="Also write the 3-pads-per-page PDF with all families side by side.",
    )
    parser.add_argument(
        "--no-split-pdfs", action="store_true",
        help="Skip the per-family PDFs (use with --combined-pdf).",
    )
    parser.add_argument(
        "--text", action="store_true",
        help="Print the numeric value in each cell.",
    )
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--n-threads", type=int, default=4)
    return parser.parse_args()


def z_range_for_category(
    values: Dict[int, Dict[str, float]],
    events: List[int],
    members: set,
    args: argparse.Namespace,
) -> Tuple[float, float]:
    """One z range per category, held fixed across every page."""
    if args.z_min is not None and args.z_max is not None:
        return args.z_min, args.z_max
    pool = sorted(
        value
        for event in events
        for branch, value in values.get(event, {}).items()
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
    """Push never-filled bins below the z minimum so ROOT leaves them white."""
    floor = (0.5 * zmin) if zmin > 0.0 else (zmin - max(1.0, abs(zmax - zmin)))
    for bx in range(0, hist.GetNbinsX() + 2):
        for by in range(0, hist.GetNbinsY() + 2):
            if hist.GetBinContent(bx, by) == 0.0:
                hist.SetBinContent(bx, by, floor)


def draw_category(hists: List[Any], with_text: bool) -> None:
    """Draw a category's histograms, putting the palette on a filled one."""
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


def select_events(data: Dict[str, Any], energy_columns: List[str],
                   args: argparse.Namespace) -> List[int]:
    event_numbers = np.asarray(data["event_n"], dtype=np.int64)
    if args.events is not None:
        requested = set(args.events)
        return [int(n) for n in event_numbers if int(n) in requested]
    energy_matrix = np.column_stack([
        np.nan_to_num(np.asarray(data[col], dtype=float), nan=0.0)
        for col in energy_columns
    ])
    totals = np.clip(energy_matrix, 0.0, None).sum(axis=1)
    n_pick = min(args.n_events, totals.size)
    order = np.argsort(totals)[::-1][:n_pick]
    return [int(event_numbers[i]) for i in order]


def main() -> None:
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetPaintTextFormat("5.0f" if args.units == "ps" else "5.2f")

    manager = CaloXAnalysisManager(make_manager_args(args, args.run)).prepare(do_fers=False)
    manager.apply_mcp_selection(flag_only=True)
    manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

    channels = active_timing_channels(manager, requested=None)
    branches = sorted(channels)
    if not branches:
        raise RuntimeError("No active calorimeter DRS channels found for run {}".format(args.run))

    energy_columns = ["{}_energy".format(b) for b in branches]
    toa_columns = ["{}_Time_cfd_mcp".format(b) for b in branches]
    data = manager.rdf.AsNumpy(["event_n"] + energy_columns + toa_columns)
    event_numbers = np.asarray(data["event_n"], dtype=np.int64)
    if event_numbers.size == 0:
        raise RuntimeError("No selected events found for run {}".format(args.run))

    events = select_events(data, energy_columns, args)
    if not events:
        raise RuntimeError("No events selected for run {}".format(args.run))

    use_relative = not args.absolute
    references: Dict[str, float] = {}
    if use_relative:
        # _Time_cfd_mcp's absolute scale is dominated by per-channel
        # cable-length offsets, not real shower timing - the single-run
        # analogue of make_toa_position_maps.py's --relative, since there is
        # no repeated scan of this run to average over instead.
        references = compute_toa_references(data, branches, args.min_adc, args.min_ref_hits)
        print("{}/{} channels have a stable same-run TOA reference (>= {} hits "
              "above {:.1f} ADC).".format(len(references), len(branches),
                                          args.min_ref_hits, args.min_adc))

    index_of = {int(n): i for i, n in enumerate(event_numbers)}
    values: Dict[int, Dict[str, float]] = {}
    for event in events:
        index = index_of[event]
        per_event: Dict[str, float] = {}
        for branch in branches:
            if use_relative and branch not in references:
                continue
            e = float(data["{}_energy".format(branch)][index])
            t = float(data["{}_Time_cfd_mcp".format(branch)][index])
            if not (math.isfinite(e) and math.isfinite(t)) or e < args.min_adc:
                continue
            if use_relative:
                rel_t = t - references[branch]
                if abs(rel_t) > args.toa_window_ns:
                    continue
                per_event[branch] = rel_t
            else:
                per_event[branch] = t
        values[event] = per_event
        print("  event {}: {} channels kept".format(event, len(per_event)))
    events = [event for event in events if values.get(event)]
    if not events:
        raise RuntimeError("No event survived the cuts")

    # Whichever mode, shift the whole page set so the earliest kept hit sits
    # at 0 - relative-mode deltas are naturally signed, and this run's raw
    # _Time_cfd_mcp isn't anchored near zero the way the attenuation-scan
    # runs are (TimingHeatmaps' "7.6 .. 16 ns" is specific to their
    # group_mpv calibration). A single constant shift, not a per-channel
    # one, so it never touches the relative spacing the page is meant to show.
    global_min = min(v for per_event in values.values() for v in per_event.values())
    values = {
        event: {branch: v - global_min for branch, v in per_event.items()}
        for event, per_event in values.items()
    }
    print("Zero-referenced to the earliest kept hit ({:.2f} ns) so TOA reads "
          "as a non-negative time of arrival.".format(global_min))

    unit_scale = 1000.0 if args.units == "ps" else 1.0
    if unit_scale != 1.0:
        values = {
            event: {branch: value * unit_scale for branch, value in per_event.items()}
            for event, per_event in values.items()
        }
        print("Units: ps (values scaled by {:.0f}).".format(unit_scale))

    map_run = args.map_run or args.run
    drsboards = build_drs_boards(run_number=map_run)
    xmin, xmax, _ = get_drs_display_x(map_run)

    members: Dict[str, set] = {name: set() for name, _, _ in CATEGORIES}
    for _, board in drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue
            branch = channel.get_channel_name(blsub=False)
            if not channel.isCer:
                members["Sci"].add(branch)
            elif channel.isQuartz:
                members["CerQuartz"].add(branch)
            else:
                members["CerPlastic"].add(branch)

    zranges = {
        name: z_range_for_category(values, events, members[name], args)
        for name, _, _ in CATEGORIES
    }
    for name, (lo, hi) in sorted(zranges.items()):
        print("z range {:<11} {:9.2f} .. {:9.2f} {} (fixed across all pages)".format(
            name, lo, hi, args.units))

    ztitle = "Time of arrival [{}]".format(args.units)
    stem = "event_toa_heatmaps_absolute" if args.absolute else "event_toa_heatmaps"

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
    for event in events:
        value_map = values.get(event, {})
        keep: List[Any] = []
        entry = {"run": args.run, "event": event, "categories": {}}
        subtitle = "run {}   event {}".format(args.run, event)

        for index, (name, kind, quartz_only) in enumerate(CATEGORIES, start=1):
            cer_hists, sci_hists = visualizeDRSBoards(
                drsboards, valuemaps=value_map,
                suffix="toa_evt{}_{}".format(event, name),
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
                single.Write("toa_map_{}_evt{}".format(name, event))

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
            canvas.Write("toa_maps_evt{}".format(event))
            canvas._keepalive = keep
        page_summary.append(entry)
        print("  page {}/{}: event {}".format(len(page_summary), len(events), event))

    if do_combined:
        canvas.Print(str(combined_pdf) + "]")
    for name, single in family_canvas.items():
        single.Print(str(family_pdf[name]) + "]")
    root_file.Close()

    import json
    summary = {
        "run": args.run,
        "map_run": map_run,
        "events": events,
        "z_variable": ztitle,
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
        print("Wrote {} ({} pages, one per event)".format(family_pdf[name], len(events)))
    if do_combined:
        print("Wrote {} ({} pages, all families side by side)".format(
            combined_pdf, len(events)))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(outdir / "{}_summary.json".format(stem)))


if __name__ == "__main__":
    main()
