#!/usr/bin/env python3
"""Per-channel ADC-integral distributions for selected TB2026 runs.

Two integral definitions are histogrammed on the same axis so they can be
compared directly:

  window  {ch}_energy  - the repo definition. compute_cfd_integral() finds the
          peak, walks back to the 20% CFD crossing, and sums the
          baseline-subtracted waveform over a fixed window of
          [crossing - window_pre, crossing + window_post) = 50 time slices by
          default. Returns exactly 0 when the peak is below min_amplitude.

  full    SumRange({ch}_blsub, 0, 1024) - the whole record after the same
          baseline subtraction, with no time-window selection at all.

The difference between the two is everything the fixed window throws away:
late scintillation tails, reflections, undershoot, and baseline residual
integrated over the remaining ~974 slices.
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
from variables.drs import add_fixed_window_energy, add_cfd_to_baseline_energy

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    draw_header,
    make_manager_args,
    safe_name,
)


# Rendered verbatim into the page intro: generate_jsroot_html passes through any
# line that begins with "<", so the whole table has to be a single line.
_TABLE_ROWS = [
    ("What", "Old &mdash; <code>{ch}_energy</code>",
     "Fixed &mdash; <code>{ch}_energy_fixed</code>",
     "Adaptive &mdash; <code>{ch}_energy_adaptive</code>"),
    ("C++ routine", "<code>compute_cfd_integral</code>",
     "<code>compute_fixed_window_integral</code>",
     "<code>compute_cfd_to_baseline_integral</code>"),
    ("Window start",
     "per event, channel's own 20% CFD crossing",
     "group Channel8 reference crossing + fixed offset",
     "per event, 20% CFD crossing on the locally-corrected waveform"),
    ("Window end",
     "start + 50 slices, fixed",
     "start + 39 slices, fixed",
     "<b>where the pulse returns to baseline</b> &mdash; walks forward from the "
     "peak until a 7-slice moving average stays &le; max(4&nbsp;ADC, 5% of peak) "
     "for 3 slices"),
    ("Typical width, 6&nbsp;mm", "50 (fixed)", "39 (fixed)",
     "<b>47&ndash;75, median 51</b> &mdash; follows the pulse"),
    ("Peak search",
     "whole record &mdash; on 16&ndash;24&nbsp;ADC pulses this latches onto "
     "noise, landing in the true pulse window in only <b>5.6%</b> of events",
     "not used; window position comes from the reference",
     "restricted to reference &plusmn;60 slices, so it cannot run away"),
    ("Baseline",
     "global median of slices [0,&nbsp;200), applied to the whole record",
     "local mean of [anchor&minus;80, anchor&minus;20)",
     "local mean of [peak&minus;80, peak&minus;20)"),
    ("Baseline drift",
     "uncorrected &mdash; Brg1 (6&nbsp;mm) boards drift &minus;2 to "
     "&minus;3.8&nbsp;ADC/sample, biasing a 50-slice sum by &minus;100 to &minus;190&nbsp;ADC",
     "cancels locally", "cancels locally"),
    ("Amplitude threshold",
     "peak &lt; 5&nbsp;ADC returns exactly 0",
     "none &mdash; empty events form a real pedestal peak",
     "none, but a width floor of 40 slices stops the window collapsing on "
     "low-amplitude channels"),
    ("Behaviour on weak channels",
     "integrates noise wherever the CFD happened to land",
     "pedestal centred on zero, &sigma; roughly halved vs old",
     "<b>biased</b> &mdash; an adaptive width shrinks with the signal, so small "
     "pulses are under-counted and the response is non-linear"),
    ("Best used for",
     "timing &mdash; it is the only one providing <code>_TS_cfd</code>, "
     "<code>_TS_peak</code>, needed by every velocity fit",
     "pedestal-preserving spectra, where zero light must be representable",
     "charge on channels with real signal; captures the full 6&nbsp;mm tail "
     "both fixed windows clip"),
]
def definition_table_html(include_fixed: bool, include_adaptive: bool) -> str:
    """Comparison table, restricted to the definitions actually plotted.

    Listing a definition that has no curve on the page invites the reader to
    look for it, so the columns follow the overlays.
    """
    keep = [0] + ([2] if include_fixed else []) + ([3] if include_adaptive else [])
    rows = [[row[i] for i in keep] for row in _TABLE_ROWS]
    # Column 1 is the old definition and is always the reference point.
    rows = [[row[0], orig[1]] + row[1:] for row, orig in zip(rows, _TABLE_ROWS)]
    return (
        '<table style="border-collapse:collapse;margin:10px 0;font-size:13px;">'
        + "".join(
            '<tr>' + "".join(
                '<{tag} style="border:1px solid #ccd;padding:4px 9px;'
                'text-align:left;vertical-align:top;{bg}">{cell}</{tag}>'.format(
                    tag="th" if i == 0 else "td",
                    bg="background:#eef2f7;" if i == 0 else "",
                    cell=cell,
                )
                for cell in row
            ) + '</tr>'
            for i, row in enumerate(rows)
        )
        + '</table>'
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Make per-channel CFD ADC-integral histograms for selected runs."
    )
    parser.add_argument("--runs", nargs="+", type=int, required=True)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--channels", nargs="+", default=None)
    parser.add_argument("--particle", default=None)
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument(
        "--no-full-waveform",
        action="store_true",
        help="Only histogram the repo (windowed CFD) integral.",
    )
    parser.add_argument(
        "--fixed-window", action="store_true",
        help=(
            "Overlay a third definition: the reference-anchored fixed-window "
            "integral with a local pre-pulse baseline "
            "(variables.drs.add_fixed_window_energy). Built for the 6 mm "
            "channels, whose low amplitude makes the per-event CFD anchor "
            "unreliable and whose boards drift enough to bias the global "
            "median(0,200) baseline."
        ),
    )
    parser.add_argument(
        "--fixed-anchor-offset", type=int, default=-336,
        help=(
            "Slices from the group Channel8 reference crossing to the start of "
            "the fixed window. Measured on run 1917; re-check for other runs. "
            "Default: -336."
        ),
    )
    parser.add_argument("--fixed-window-pre", type=int, default=4)
    parser.add_argument("--fixed-window-post", type=int, default=35)
    parser.add_argument("--fixed-baseline-lo", type=int, default=80)
    parser.add_argument("--fixed-baseline-hi", type=int, default=20)
    parser.add_argument(
        "--adaptive-window", action="store_true",
        help=(
            "Overlay a fourth definition: CFD start, baseline-crossing end "
            "(variables.drs.add_cfd_to_baseline_energy). The window length "
            "follows the pulse instead of being fixed, which is what the 6 mm "
            "fibers need - their pulse is ~42-45 slices against ~13 for 3 mm."
        ),
    )
    parser.add_argument("--adaptive-stop-abs", type=float, default=4.0,
                        help="Absolute floor of the baseline-crossing stop level [ADC].")
    parser.add_argument("--adaptive-stop-fraction", type=float, default=0.05)
    parser.add_argument("--adaptive-stop-consecutive", type=int, default=3)
    parser.add_argument("--adaptive-smooth-half", type=int, default=3,
                        help="Half-width of the moving average used for the stop test.")
    parser.add_argument("--adaptive-min-width", type=int, default=40,
                        help="Floor on the adaptive window width [slices].")
    parser.add_argument(
        "--write-adaptive-hist", action="store_true",
        help="Also write the adaptive histograms into the output .root. Off by default.",
    )
    parser.add_argument(
        "--write-fixed-hist", action="store_true",
        help=(
            "Also write the fixed-window histograms into the output .root. Off "
            "by default: the definition is still under test, so it shows up in "
            "the plots without landing in the ROOT file."
        ),
    )
    parser.add_argument(
        "--waveform-slices",
        type=int,
        default=1024,
        help="Record length summed for the full-waveform integral. Default: 1024.",
    )
    parser.add_argument(
        "--range-low-pad",
        type=float,
        default=0.6,
        help=(
            "Extend the low edge of the shared histogram axis by this fraction of "
            "the repo energy range, to catch undershoot in the full-waveform "
            "integral. Default: 0.6."
        ),
    )
    parser.add_argument(
        "--range-high-scale",
        type=float,
        default=1.6,
        help=(
            "Multiply the high edge of the shared histogram axis by this factor. "
            "Default: 1.6."
        ),
    )
    parser.add_argument("--rebin-factor", type=int, default=2)
    parser.add_argument(
        "--logy", action="store_true",
        help=(
            "Log y axis. These distributions are a tall pedestal peak plus a "
            "long beam tail spanning a decade or more in x, so no linear range "
            "shows both: on a linear y the pedestal dominates and the tail is "
            "invisible, and zooming to the pedestal clips the tail. Log y shows "
            "both at once."
        ),
    )
    parser.add_argument(
        "--zoom-keep", type=float, default=0.995,
        help=(
            "Fraction of each histogram's entries the displayed x-range must "
            "cover. Default 0.995. Set to 1 to show the full stored range."
        ),
    )
    parser.add_argument(
        "--no-auto-range",
        action="store_true",
        help=(
            "Use the widened repo energy range instead of deriving a per-channel "
            "axis from a mean/RMS pre-pass. The fixed range is far too narrow for "
            "the full-record integral on 6 mm channels, so auto-range is default."
        ),
    )
    parser.add_argument(
        "--range-ignore-full", action="store_true",
        help=(
            "Leave the full-record integral out of the auto-range calculation. "
            "Its RMS is ~4x the windowed one, so when all three definitions "
            "share an axis it stretches the range and squashes the two "
            "comparable distributions into a spike. The full-record histogram "
            "still draws; it just spills into the overflow."
        ),
    )
    parser.add_argument(
        "--range-sigma",
        type=float,
        default=8.0,
        help=(
            "Half-width of the *stored* histogram range in RMS units. Kept "
            "generous (8) so the tails are not pushed into the overflow bin, "
            "which showed up as a spike at the last bin and clipped content at "
            "both edges. The displayed range is narrowed separately by "
            "--zoom-keep."
        ),
    )
    parser.add_argument(
        "--outdir", default="results/attenuation/adc_distributions_representative_runs"
    )
    return parser.parse_args()


def shared_energy_model_range(
    is6mm: bool, args: argparse.Namespace
) -> Tuple[int, float, float]:
    """Common (nbins, xmin, xmax) for both integral definitions.

    Both observables must share binning, otherwise the overlay compares
    different bin widths. The repo range is the anchor; it is widened downward
    (undershoot) and upward (late tails) since the full-record sum is not
    bounded by the CFD window.
    """
    nbins, xmin, xmax = get_drs_energy_range(is6mm)
    span = xmax - xmin
    return (
        max(1, nbins) * max(1, args.rebin_factor),
        xmin - args.range_low_pad * span,
        xmax * args.range_high_scale,
    )


def auto_model_range(
    stats: List[Tuple[float, float]],
    fallback: Tuple[int, float, float],
    args: argparse.Namespace,
) -> Tuple[int, float, float]:
    """Shared axis covering mean +/- range_sigma * RMS of every observable.

    `stats` holds one (mean, rms) pair per observable. The full-record integral
    picks up baseline noise from ~1000 slices, so its width is set by the data,
    not by the repo's CFD-window range; deriving the axis is the only way both
    definitions land on screen for 3 mm and 6 mm channels alike.
    """
    usable = [
        (mean, rms)
        for mean, rms in stats
        if math.isfinite(mean) and math.isfinite(rms) and rms > 0.0
    ]
    if not usable:
        return fallback
    lo = min(mean - args.range_sigma * rms for mean, rms in usable)
    hi = max(mean + args.range_sigma * rms for mean, rms in usable)
    # The windowed integral piles up at exactly 0 for sub-threshold peaks; keep
    # that bin visible.
    lo, hi = min(lo, 0.0), max(hi, 0.0)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return fallback
    return fallback[0], lo, hi


def content_range(
    hists: List[Any], keep: float = 0.995
) -> Optional[Tuple[float, float]]:
    """Display range covering `keep` of each histogram's entries, unioned.

    Taking the first and last non-empty bin instead (the previous behaviour)
    lets a handful of far-out events stretch the view across the whole axis,
    while a range set from mean +/- 4 RMS clipped real content into the
    overflow bin. A quantile keeps the bulk on screen without either failure.
    """
    lo, hi = None, None
    tail = max(0.0, min(0.45, 0.5 * (1.0 - keep)))
    for hist in hists:
        axis = hist.GetXaxis()
        nbins = hist.GetNbinsX()
        total = hist.Integral(0, nbins + 1)
        if total <= 0.0:
            continue
        want_lo, want_hi = tail * total, (1.0 - tail) * total
        running = hist.GetBinContent(0)
        low_edge, high_edge = axis.GetXmin(), axis.GetXmax()
        got_lo = False
        for index in range(1, nbins + 1):
            running += hist.GetBinContent(index)
            if not got_lo and running >= want_lo:
                low_edge = axis.GetBinLowEdge(index)
                got_lo = True
            if running >= want_hi:
                high_edge = axis.GetBinUpEdge(index)
                break
        lo = low_edge if lo is None else min(lo, low_edge)
        hi = high_edge if hi is None else max(hi, high_edge)
    if lo is None or hi is None or hi <= lo:
        return None
    pad = 0.06 * (hi - lo)
    return lo - pad, hi + pad


def make_canvas(
    run: int,
    channel: ChannelInfo,
    hist_window: Any,
    hist_full: Optional[Any],
    hist_fixed: Optional[Any],
    hist_adaptive: Optional[Any],
    selection_label: str,
    zoom_hists: Optional[List[Any]] = None,
    zoom_keep: float = 0.995,
    logy: bool = False,
) -> tuple[Any, str]:
    name = "adc_integral_run{}_{}".format(run, safe_name(channel.branch))
    canvas = ROOT.TCanvas(name, "", 1000, 750)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    drawn = (
        [hist_window]
        + ([hist_full] if hist_full is not None else [])
        + ([hist_fixed] if hist_fixed is not None else [])
        + ([hist_adaptive] if hist_adaptive is not None else [])
    )
    for hist in drawn:
        hist.SetTitle("")
        hist.SetStats(0)
        hist.GetXaxis().SetTitle("ADC integral [ADC]")
        hist.GetYaxis().SetTitle("Events")
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetTitleOffset(1.35)

    hist_window.SetLineColor(ROOT.kBlue + 1)
    hist_window.SetFillColorAlpha(ROOT.kAzure - 9, 0.35)
    hist_window.SetLineWidth(2)
    if hist_full is not None:
        hist_full.SetLineColor(ROOT.kRed + 1)
        hist_full.SetFillStyle(0)
        hist_full.SetLineWidth(2)
    if hist_fixed is not None:
        hist_fixed.SetLineColor(ROOT.kBlack)
        hist_fixed.SetFillStyle(0)
        hist_fixed.SetLineWidth(2)
        hist_fixed.SetLineStyle(2)
    if hist_adaptive is not None:
        hist_adaptive.SetLineColor(ROOT.kMagenta + 2)
        hist_adaptive.SetFillStyle(0)
        hist_adaptive.SetLineWidth(2)

    ymax = max(hist.GetMaximum() for hist in drawn)
    # Extra headroom so the stats box and legend sit above the distribution
    # rather than on top of it.
    if logy:
        canvas.SetLogy(1)
        hist_window.SetMinimum(0.5)
        hist_window.SetMaximum(max(10.0, 60.0 * ymax))
    else:
        hist_window.SetMinimum(0.0)
        hist_window.SetMaximum(max(1.0, 1.95 * ymax))
    zoom = content_range(zoom_hists or drawn, keep=zoom_keep)
    if zoom is not None:
        hist_window.GetXaxis().SetRangeUser(*zoom)

    hist_window.Draw("HIST")
    if hist_full is not None:
        hist_full.Draw("HIST SAME")
    if hist_fixed is not None:
        hist_fixed.Draw("HIST SAME")
    if hist_adaptive is not None:
        hist_adaptive.Draw("HIST SAME")

    stats = ROOT.TPaveText(0.55, 0.555, 0.95, 0.885, "NDC")
    stats.SetBorderSize(0)
    stats.SetFillStyle(0)
    stats.SetTextFont(42)
    stats.SetTextSize(0.021)
    stats.SetTextAlign(12)
    stats.AddText("Windowed CFD integral (repo)")
    stats.AddText("  Entries {}  Mean {:.1f}  RMS {:.1f}".format(
        int(hist_window.GetEntries()), hist_window.GetMean(), hist_window.GetRMS()
    ))
    if hist_full is not None:
        stats.AddText("Full-waveform integral")
        stats.AddText("  Entries {}  Mean {:.1f}  RMS {:.1f}".format(
            int(hist_full.GetEntries()), hist_full.GetMean(), hist_full.GetRMS()
        ))
        window_mean = hist_window.GetMean()
        if abs(window_mean) > 1.0e-9:
            stats.AddText("  Mean_{{full}} / Mean_{{window}} = {:.3f}".format(
                hist_full.GetMean() / window_mean
            ))
    if hist_fixed is not None:
        stats.AddText("Fixed window, local baseline")
        stats.AddText("  Entries {}  Mean {:.1f}  RMS {:.1f}".format(
            int(hist_fixed.GetEntries()), hist_fixed.GetMean(), hist_fixed.GetRMS()
        ))
    if hist_adaptive is not None:
        stats.AddText("CFD #rightarrow baseline crossing")
        stats.AddText("  Entries {}  Mean {:.1f}  RMS {:.1f}".format(
            int(hist_adaptive.GetEntries()), hist_adaptive.GetMean(),
            hist_adaptive.GetRMS()
        ))
    stats.Draw()

    legend = ROOT.TLegend(0.55, 0.395, 0.95, 0.545)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.020)
    legend.AddEntry(hist_window, "Repo: CFD window [-5, +45] TS", "lf")
    if hist_full is not None:
        legend.AddEntry(hist_full, "Whole record, baseline-subtracted", "l")
    if hist_fixed is not None:
        legend.AddEntry(hist_fixed, "Fixed window @ ref, local baseline", "l")
    if hist_adaptive is not None:
        legend.AddEntry(hist_adaptive, "CFD #rightarrow baseline crossing", "l")
    legend.Draw()

    draw_header(
        canvas,
        "Run {} ADC integral: {}".format(run, channel.short_label),
        selection_label,
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = drawn + [stats, legend]
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
    do_full = not args.no_full_waveform
    root_path = outdir / "adc_integral_distributions.root"
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
        columns = {str(column) for column in manager.rdf.GetColumnNames()}
        rdf = manager.rdf

        # Full-record integral: same baseline-subtracted waveform the CFD sees,
        # summed over every slice, with no time-window selection.
        full_columns: Dict[str, str] = {}
        if do_full:
            for branch in sorted(live):
                blsub = "{}_blsub".format(branch)
                if blsub not in columns:
                    continue
                full_column = "{}_energy_fullwf".format(branch)
                if full_column not in columns:
                    rdf = rdf.Define(
                        full_column,
                        "SumRange({}, 0, {})".format(blsub, args.waveform_slices),
                    )
                full_columns[branch] = full_column
            missing = sorted(set(live) - set(full_columns))
            if missing:
                print("[WARN] No _blsub column for {} channel(s); "
                      "full-waveform integral skipped there.".format(len(missing)))

        # Reference-anchored fixed window with a local pre-pulse baseline.
        fixed_columns: Dict[str, str] = {}
        if args.fixed_window:
            before = {str(c) for c in rdf.GetColumnNames()}
            rdf = add_fixed_window_energy(
                rdf, sorted(live),
                anchor_offset=args.fixed_anchor_offset,
                window_pre=args.fixed_window_pre,
                window_post=args.fixed_window_post,
                baseline_pre_lo=args.fixed_baseline_lo,
                baseline_pre_hi=args.fixed_baseline_hi,
            )
            after = {str(c) for c in rdf.GetColumnNames()}
            for branch in sorted(live):
                name = "{}_energy_fixed".format(branch)
                if name in after:
                    fixed_columns[branch] = name
            print("[FIXED] Defined {} fixed-window column(s); anchor = ref_TS {:+d}, "
                  "window [{:+d}, {:+d}), baseline [{:+d}, {:+d}).".format(
                      len(fixed_columns), args.fixed_anchor_offset,
                      -args.fixed_window_pre, args.fixed_window_post,
                      -args.fixed_baseline_lo, -args.fixed_baseline_hi))
            skipped = sorted(set(live) - set(fixed_columns))
            if skipped:
                print("[WARN] No reference-TS column for {} channel(s); "
                      "fixed-window integral skipped there.".format(len(skipped)))

        # CFD start, baseline-crossing end: window length follows the pulse.
        adaptive_columns: Dict[str, str] = {}
        if args.adaptive_window:
            rdf = add_cfd_to_baseline_energy(
                rdf, sorted(live),
                anchor_offset=args.fixed_anchor_offset + 6,
                stop_fraction=args.adaptive_stop_fraction,
                stop_abs=args.adaptive_stop_abs,
                stop_consecutive=args.adaptive_stop_consecutive,
                smooth_half=args.adaptive_smooth_half,
                min_width=args.adaptive_min_width,
            )
            after = {str(c) for c in rdf.GetColumnNames()}
            for branch in sorted(live):
                name = "{}_energy_adaptive".format(branch)
                if name in after:
                    adaptive_columns[branch] = name
            print("[ADAPT] Defined {} adaptive column(s); peak search @ ref {:+d} "
                  "+/-60, stop <= max({:.1f}, {:.2f}*peak) for {} slices "
                  "(smoothed +/-{}), min width {}.".format(
                      len(adaptive_columns), args.fixed_anchor_offset + 6,
                      args.adaptive_stop_abs, args.adaptive_stop_fraction,
                      args.adaptive_stop_consecutive, args.adaptive_smooth_half,
                      args.adaptive_min_width))

        # Pre-pass: mean/RMS of both observables, so the shared axis can be sized
        # from the data rather than from the CFD-window-tuned repo range.
        model_ranges: Dict[str, Tuple[int, float, float]] = {}
        if not args.no_auto_range:
            stat_actions: Dict[str, List[Any]] = {}
            for branch in sorted(live):
                actions = [
                    rdf.Mean("{}_energy".format(branch)),
                    rdf.StdDev("{}_energy".format(branch)),
                ]
                if branch in full_columns and not args.range_ignore_full:
                    actions += [
                        rdf.Mean(full_columns[branch]),
                        rdf.StdDev(full_columns[branch]),
                    ]
                if branch in fixed_columns:
                    actions += [
                        rdf.Mean(fixed_columns[branch]),
                        rdf.StdDev(fixed_columns[branch]),
                    ]
                if branch in adaptive_columns:
                    actions += [
                        rdf.Mean(adaptive_columns[branch]),
                        rdf.StdDev(adaptive_columns[branch]),
                    ]
                stat_actions[branch] = actions
            flat = [action for actions in stat_actions.values() for action in actions]
            print("[RANGE] Pre-pass over {} statistics actions...".format(len(flat)))
            ROOT.RDF.RunGraphs(flat)
            for branch, channel in sorted(live.items()):
                values = [float(action.GetValue()) for action in stat_actions[branch]]
                pairs = list(zip(values[0::2], values[1::2]))
                model_ranges[branch] = auto_model_range(
                    pairs,
                    shared_energy_model_range(channel.size == "6 mm", args),
                    args,
                )

        results: Dict[str, Any] = {}
        results_full: Dict[str, Any] = {}
        results_fixed: Dict[str, Any] = {}
        results_adaptive: Dict[str, Any] = {}
        for branch, channel in sorted(live.items()):
            column = "{}_energy".format(branch)
            is6mm = channel.size == "6 mm"
            model_range = model_ranges.get(
                branch, shared_energy_model_range(is6mm, args)
            )
            model = ROOT.RDF.TH1DModel(
                "hist_adc_integral_run{}_{}".format(run, safe_name(branch)),
                "",
                *model_range,
            )
            results[branch] = rdf.Filter(
                "std::isfinite({})".format(column)
            ).Histo1D(model, column)

            if branch in full_columns:
                full_column = full_columns[branch]
                model_full = ROOT.RDF.TH1DModel(
                    "hist_adc_integral_fullwf_run{}_{}".format(
                        run, safe_name(branch)
                    ),
                    "",
                    *model_range,
                )
                results_full[branch] = rdf.Filter(
                    "std::isfinite({})".format(full_column)
                ).Histo1D(model_full, full_column)

            if branch in fixed_columns:
                fixed_column = fixed_columns[branch]
                model_fixed = ROOT.RDF.TH1DModel(
                    "hist_adc_integral_fixed_run{}_{}".format(
                        run, safe_name(branch)
                    ),
                    "",
                    *model_range,
                )
                results_fixed[branch] = rdf.Filter(
                    "std::isfinite({})".format(fixed_column)
                ).Histo1D(model_fixed, fixed_column)

            if branch in adaptive_columns:
                adaptive_column = adaptive_columns[branch]
                model_adaptive = ROOT.RDF.TH1DModel(
                    "hist_adc_integral_adaptive_run{}_{}".format(
                        run, safe_name(branch)
                    ),
                    "",
                    *model_range,
                )
                results_adaptive[branch] = rdf.Filter(
                    "std::isfinite({})".format(adaptive_column)
                ).Histo1D(model_adaptive, adaptive_column)

        ROOT.RDF.RunGraphs(
            list(results.values()) + list(results_full.values())
            + list(results_fixed.values()) + list(results_adaptive.values())
        )

        run_dir = root_file.mkdir("run{}".format(run))
        run_summary: Dict[str, Dict[str, float]] = {}
        for branch, result in sorted(results.items()):
            hist = result.GetValue().Clone(
                "hist_adc_integral_run{}_{}".format(run, safe_name(branch))
            )
            hist.SetDirectory(0)
            run_dir.cd()
            hist.Write()

            hist_full = None
            if branch in results_full:
                hist_full = results_full[branch].GetValue().Clone(
                    "hist_adc_integral_fullwf_run{}_{}".format(run, safe_name(branch))
                )
                hist_full.SetDirectory(0)
                run_dir.cd()
                hist_full.Write()

            hist_fixed = None
            if branch in results_fixed:
                hist_fixed = results_fixed[branch].GetValue().Clone(
                    "hist_adc_integral_fixed_run{}_{}".format(run, safe_name(branch))
                )
                hist_fixed.SetDirectory(0)
                # Under test: kept out of the ROOT file unless asked for.
                if args.write_fixed_hist:
                    run_dir.cd()
                    hist_fixed.Write()

            axis = hist.GetXaxis()
            hist_adaptive = None
            if branch in results_adaptive:
                hist_adaptive = results_adaptive[branch].GetValue().Clone(
                    "hist_adc_integral_adaptive_run{}_{}".format(run, safe_name(branch))
                )
                hist_adaptive.SetDirectory(0)
                if args.write_adaptive_hist:
                    run_dir.cd()
                    hist_adaptive.Write()

            entry = {
                "axis_xmin_adc": float(axis.GetXmin()),
                "axis_xmax_adc": float(axis.GetXmax()),
                "window_entries": float(hist.GetEntries()),
                "window_mean_adc": float(hist.GetMean()),
                "window_rms_adc": float(hist.GetRMS()),
            }
            if hist_full is not None:
                entry.update({
                    "full_entries": float(hist_full.GetEntries()),
                    "full_mean_adc": float(hist_full.GetMean()),
                    "full_rms_adc": float(hist_full.GetRMS()),
                })
                if abs(hist.GetMean()) > 1.0e-9:
                    entry["full_over_window_mean"] = (
                        float(hist_full.GetMean()) / float(hist.GetMean())
                    )
            if hist_fixed is not None:
                entry.update({
                    "fixed_entries": float(hist_fixed.GetEntries()),
                    "fixed_mean_adc": float(hist_fixed.GetMean()),
                    "fixed_rms_adc": float(hist_fixed.GetRMS()),
                })
                if abs(hist.GetMean()) > 1.0e-9:
                    entry["fixed_over_window_mean"] = (
                        float(hist_fixed.GetMean()) / float(hist.GetMean())
                    )
            if hist_adaptive is not None:
                entry.update({
                    "adaptive_entries": float(hist_adaptive.GetEntries()),
                    "adaptive_mean_adc": float(hist_adaptive.GetMean()),
                    "adaptive_rms_adc": float(hist_adaptive.GetRMS()),
                })
                if abs(hist.GetMean()) > 1.0e-9:
                    entry["adaptive_over_window_mean"] = (
                        float(hist_adaptive.GetMean()) / float(hist.GetMean())
                    )
            run_summary[branch] = entry

            # The full-record definition is ~4x wider; if it drives the zoom
            # the two comparable distributions collapse into a spike.
            zoom_hists = [hist] + ([hist_fixed] if hist_fixed is not None else [])
            if hist_adaptive is not None:
                zoom_hists.append(hist_adaptive)
            if hist_full is not None and not args.range_ignore_full:
                zoom_hists.append(hist_full)
            canvas, canvas_json = make_canvas(
                run, live[branch], hist, hist_full, hist_fixed, hist_adaptive,
                selection,
                zoom_hists=zoom_hists, zoom_keep=args.zoom_keep,
                logy=args.logy,
            )
            name = canvas.GetName()
            root_file.cd()
            canvas.Write(name)
            keys.append(name)
            jsons[name] = canvas_json
        summary["run{}".format(run)] = run_summary

    root_file.Close()
    html_path = outdir / "adc_integral_distributions.html"
    observable_text = (
        "**Observables (same binning, overlaid):**\n"
        "* Blue, filled - the repo integral, `{ch}_energy`. `compute_cfd_integral()` "
        "locates the peak, walks back to the 20% CFD crossing, and sums the "
        "baseline-subtracted waveform over a fixed 50-slice window "
        "`[crossing - 5, crossing + 45)`. It returns exactly 0 when the peak is "
        "below the 5-ADC amplitude threshold.\n"
        "* Red, open - the full-record integral, `SumRange({ch}_blsub, 0, "
        + str(args.waveform_slices)
        + ")`: the same baseline-subtracted waveform summed over every slice, with "
        "no time-window selection.\n"
        "* Baseline for both is the median of slices [0, 200), subtracted from the "
        "whole record (sign-flipped for inverted channels).\n"
    ) if do_full else (
        "**Observable:** baseline-subtracted CFD pulse integral.\n"
    )
    if args.fixed_window:
        observable_text += (
            "* Black, dashed - **under test**: fixed window anchored on the group "
            "Channel8 reference crossing at ref_TS {:+d}, spanning [{:+d}, {:+d}) "
            "slices, rebaselined on [{:+d}, {:+d}) immediately ahead of the pulse. "
            "No amplitude threshold, so empty events form a pedestal peak at zero "
            "instead of being forced to exactly 0. Built for the 6 mm channels; "
            "the anchor offset was measured on run 1917.\n"
        ).format(args.fixed_anchor_offset, -args.fixed_window_pre,
                 args.fixed_window_post, -args.fixed_baseline_lo,
                 -args.fixed_baseline_hi)
        if not args.write_fixed_hist:
            observable_text += (
                "* The fixed-window histograms are **not** written to the .root "
                "file (pass --write-fixed-hist to persist them).\n"
            )
    if args.adaptive_window:
        observable_text += (
            "* Magenta - **under test**: CFD start, baseline-crossing end. The "
            "window runs from the 20% CFD crossing to where a 7-slice moving "
            "average of the pulse falls back to max({:.0f} ADC, {:.0%} of peak) "
            "for {} consecutive slices, with a {}-slice width floor. Peak search "
            "restricted to the reference crossing +/- 60 slices.\n"
        ).format(args.adaptive_stop_abs, args.adaptive_stop_fraction,
                 args.adaptive_stop_consecutive, args.adaptive_min_width)
        if not args.write_adaptive_hist:
            observable_text += (
                "* The adaptive histograms are **not** written to the .root file "
                "(pass --write-adaptive-hist to persist them).\n"
            )
    if args.fixed_window or args.adaptive_window:
        observable_text += definition_table_html(
            args.fixed_window, args.adaptive_window) + "\n"

    generate_jsroot_html(
        keys,
        jsons,
        plots_per_row=3,
        output_html=str(html_path),
        title="CaloX ADC-integral distributions",
        intro_text=(
            "**Runs:** {}.\n".format(", ".join(str(run) for run in args.runs))
            + observable_text
            + "**Selection:** {}.".format(selection)
        ),
    )
    json_path = outdir / "adc_integral_summary.json"
    with json_path.open("w") as handle:
        json.dump(
            {
                "runs": args.runs,
                "selection": selection,
                "waveform_slices": args.waveform_slices,
                "channels": summary,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(root_path))
    print("Wrote {}".format(json_path))


if __name__ == "__main__":
    main()
