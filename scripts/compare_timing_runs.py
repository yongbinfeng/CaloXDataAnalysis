
"""Overlay MCP-corrected time-of-arrival distributions across TB2026 runs.

This is a native CaloXDataAnalysis multi-run driver.  Each run gets its own
CaloXAnalysisManager, so the active DRS channel map, DRS polarity handling,
reference-channel correction, MCP correction, and the group-MPV time
calibration are all taken from the current tb2026 code for that run.

For every active non-reference calorimeter DRS channel it writes:
  * one raw TH1D per run to a combined ROOT file;
  * one interactive JSROOT HTML overlay with all requested runs, normalized only for drawing.

It also writes a run-mapped MCP timing-difference histogram per run, displays it in the same HTML gallery, and fits a
Gaussian to each.  The reported single-MCP resolution assumes two equal,
independent MCP resolutions: sigma_single = sigma_delta_t / sqrt(2).

Examples
--------
  python scripts/compare_timing_runs.py --runs 1776 1777 1778

  python scripts/compare_timing_runs.py \
      --runs 1776 1777 \
      --json-file data/datafiles_local.json \
      --particle pion \
      --outdir results/timing_compare/runs_1776_1777

The default MCP clean-pulse selection is intentional.  Do not replace it with
apply_mcp_diff_selection: that selection would bias the MCP delta-t width.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import ROOT

from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from utils.data_loader import getRunInfo
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var
from variables.drs import MCP_REF, _TS_TO_NS


# Keep a stable run-to-colour mapping in every output plot.
ROOT_COLORS = [
    ROOT.kBlue + 1,
    ROOT.kRed + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
    ROOT.kCyan + 2,
    ROOT.kViolet + 1,
    ROOT.kTeal + 3,
    ROOT.kPink + 7,
    ROOT.kAzure + 2,
]


@dataclass(frozen=True)
class ChannelInfo:
    """Display metadata taken from the current run's CaloX channel object."""

    branch: str
    label: str


@dataclass(frozen=True)
class MCPPair:
    """Logical MCP pair selected from the live TB2026 channel map."""

    first: str
    second: str

    @property
    def label(self) -> str:
        return f"{self.second} - {self.first}"


@dataclass
class FitSummary:
    """Gaussian MCP delta-t fit result in ns."""

    n_entries: int
    mean_ns: float | None
    mean_err_ns: float | None
    sigma_ns: float | None
    sigma_err_ns: float | None
    fit_ok: bool

    @property
    def single_mcp_sigma_ns(self) -> float | None:
        return None if self.sigma_ns is None else self.sigma_ns / math.sqrt(2.0)

    @property
    def single_mcp_sigma_err_ns(self) -> float | None:
        return None if self.sigma_err_ns is None else self.sigma_err_ns / math.sqrt(2.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make per-channel MCP-corrected TOA overlays and MCP2-MCP1 "
            "timing-resolution overlays for multiple TB2026 runs."
        )
    )
    parser.add_argument(
        "--runs",
        nargs="+",
        type=int,
        required=True,
        help="Run numbers to compare.  Their files are looked up in --json-file.",
    )
    parser.add_argument(
        "--json-file",
        default="data/datafiles_local.json",
        help="Run-to-ROOT-file JSON map used by CaloXAnalysisManager.",
    )
    parser.add_argument(
        "--first-event",
        type=int,
        default=0,
        help="First event index, inclusive (default: 0).",
    )
    parser.add_argument(
        "--last-event",
        type=int,
        default=-1,
        help="Last event index, exclusive; -1 means the manager default (default: -1).",
    )
    parser.add_argument(
        "--particle",
        default=None,
        help=(
            "Optional CaloX particle selection, e.g. pion/electron/proton. "
            "Applied after the MCP clean selection."
        ),
    )
    parser.add_argument(
        "--no-mcp-clean",
        action="store_true",
        help="Disable the native MCP clean-pulse selection (not recommended for a resolution fit).",
    )
    parser.add_argument(
        "--channels",
        nargs="+",
        default=None,
        metavar="DRS_BRANCH",
        help=(
            "Optional exact DRS branch-name subset. By default all mapped, "
            "non-reference calorimeter DRS channels are used."
        ),
    )
    parser.add_argument(
        "--toa-bins",
        type=int,
        default=400,
        help="Number of bins for each TOA histogram (default: 400).",
    )
    parser.add_argument(
        "--toa-range",
        nargs=2,
        type=float,
        default=(-30.0, 30.0),
        metavar=("MIN_NS", "MAX_NS"),
        help="Shared TOA axis range in ns (default: -30 30).",
    )
    parser.add_argument(
        "--mcp-dt-bins",
        type=int,
        default=400,
        help="Number of bins for MCP2-MCP1 delta-t (default: 400).",
    )
    parser.add_argument(
        "--mcp-dt-range",
        nargs=2,
        type=float,
        default=(-5.0, 5.0),
        metavar=("MIN_NS", "MAX_NS"),
        help="Shared MCP2-MCP1 axis range in ns (default: -5 5).",
    )
    parser.add_argument(
        "--mcp-a",
        default=None,
        help=(
            "Optional logical first-MCP override. Use only with --mcp-b. "
            "Default: resolve the native TB2026 MCP pair from the run map."
        ),
    )
    parser.add_argument(
        "--mcp-b",
        default=None,
        help=(
            "Optional logical second-MCP override. Use only with --mcp-a. "
            "Default: resolve the native TB2026 MCP pair from the run map."
        ),
    )
    parser.add_argument(
        "--raw-yields",
        action="store_true",
        help="Draw raw event counts instead of unit-area-normalized overlays.",
    )
    parser.add_argument(
        "--n-threads",
        type=int,
        default=10,
        help="ROOT implicit-MT thread count (default: 10).",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help=(
            "Combined output directory. Default: results/timing_compare/"
            "runs_<run>_<run>_..."
        ),
    )
    args = parser.parse_args()

    if len(set(args.runs)) != len(args.runs):
        parser.error("--runs contains a duplicate run number")
    if (args.mcp_a is None) != (args.mcp_b is None):
        parser.error("--mcp-a and --mcp-b must be supplied together")

    if args.toa_range[0] >= args.toa_range[1]:
        parser.error("--toa-range requires MIN_NS < MAX_NS")
    if args.mcp_dt_range[0] >= args.mcp_dt_range[1]:
        parser.error("--mcp-dt-range requires MIN_NS < MAX_NS")
    if args.toa_bins <= 0 or args.mcp_dt_bins <= 0:
        parser.error("Histogram bin counts must be positive")

    if args.outdir is None:
        run_tag = "_".join(f"run{run}" for run in args.runs)
        args.outdir = os.path.join("results", "timing_compare", run_tag)

    return args


def safe_name(value: str) -> str:
    """Make an arbitrary branch name safe for a ROOT object/file name."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def make_manager_args(args: argparse.Namespace, run: int) -> argparse.Namespace:
    """Build exactly the small args interface CaloXAnalysisManager needs."""
    return argparse.Namespace(
        run=run,
        first_event=args.first_event,
        last_event=args.last_event,
        json_file=args.json_file,
    )


def active_timing_channels(manager: CaloXAnalysisManager, requested: set[str] | None) -> dict[str, ChannelInfo]:
    """Return mapped non-reference calorimeter channels with Time_cfd_mcp available.

    No DRS branch, board number, 3-mm/6-mm layout, or material family is
    hard-coded here.  The current run's build_drs_boards() map supplies all of
    that, which is why this remains correct when TB2026 channel assignments
    change from one run period to another.
    """
    columns = {str(column) for column in manager.rdf.GetColumnNames()}
    selected: dict[str, ChannelInfo] = {}

    for _, board in manager.drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue

            branch = channel.get_channel_name(blsub=False)
            if requested is not None and branch not in requested:
                continue

            toa_column = f"{branch}_Time_cfd_mcp"
            if toa_column not in columns:
                # This protects runs where a mapped channel is physically absent
                # or lacks a valid time-MPV group calibration.
                continue

            try:
                material = get_channel_var(channel)
            except Exception:
                material = "DRS"
            size = "6 mm" if getattr(channel, "is6mm", False) else "3 mm"
            selected[branch] = ChannelInfo(branch=branch, label=f"{branch} | {material}, {size}")

    return selected


def resolve_mcp_pair(
    run: int,
    args: argparse.Namespace,
    columns: set[str],
) -> MCPPair | None:
    """Choose the MCP pair from the run's live TB2026 map.

    The native selection manager uses keys 0 and 4 from
    ``get_mcp_channels(run)`` for its MCP time-difference selection.  We
    preserve that convention when both timing columns exist.  If a run has a
    smaller valid map, use the first and last *live* MCP timing columns rather
    than silently reverting to names from a previous run period.
    """
    mapped = list(get_mcp_channels(run).keys())
    live = [det for det in mapped if f"{det}_TS_cfd_ref" in columns]

    if args.mcp_a is not None:
        pair = MCPPair(args.mcp_a, args.mcp_b)
        missing = [
            f"{det}_TS_cfd_ref"
            for det in (pair.first, pair.second)
            if f"{det}_TS_cfd_ref" not in columns
        ]
        if missing:
            print(
                f"[WARN] Run {run}: requested MCP pair {pair.label} is unavailable; "
                f"missing {', '.join(missing)}"
            )
            return None
        print(f"[MCP] Run {run}: using user-requested pair {pair.label}")
        return pair

    # This is the pair used by SelectionManager.apply_mcp_diff_selection().
    if len(mapped) >= 5:
        native_pair = MCPPair(mapped[0], mapped[4])
        if (
            f"{native_pair.first}_TS_cfd_ref" in columns
            and f"{native_pair.second}_TS_cfd_ref" in columns
        ):
            print(f"[MCP] Run {run}: using native mapped pair {native_pair.label}")
            return native_pair

    # Some runs may expose only two usable MCP timing columns after the raw
    # branch availability check in variables/drs.py.
    if len(live) >= 2:
        fallback_pair = MCPPair(live[0], live[-1])
        print(
            f"[WARN] Run {run}: native MCP pair is unavailable; "
            f"using live mapped fallback {fallback_pair.label}"
        )
        return fallback_pair

    mapped_text = ", ".join(mapped) if mapped else "none"
    live_text = ", ".join(live) if live else "none"
    print(
        f"[WARN] Run {run}: no usable MCP delta-t pair. "
        f"Mapped MCP labels: {mapped_text}; live timing labels: {live_text}"
    )
    return None


def mcp_reference_is_live(run: int, columns: set[str]) -> bool:
    """Whether the timing pipeline can apply its configured MCP reference."""
    return (
        MCP_REF in get_mcp_channels(run)
        and f"{MCP_REF}_TS_cfd_ref" in columns
    )


def book_run_histograms(
    manager: CaloXAnalysisManager,
    run: int,
    args: argparse.Namespace,
    channel_info: dict[str, ChannelInfo],
) -> tuple[
    dict[str, ROOT.RDF.RResultPtr],
    ROOT.RDF.RResultPtr | None,
    MCPPair | None,
]:
    """Book all lazy actions for one run, without evaluating the dataframe yet."""
    toa_results: dict[str, ROOT.RDF.RResultPtr] = {}

    for branch in sorted(channel_info):
        column = f"{branch}_Time_cfd_mcp"
        model = ROOT.RDF.TH1DModel(
            f"toa_run{run}_{safe_name(branch)}",
            f"Run {run};Time of Arrival [ns];Events",
            args.toa_bins,
            args.toa_range[0],
            args.toa_range[1],
        )
        # Reject NaN values before filling.  The formula is valid in cling and
        # does not impose a physics timing window.
        toa_results[branch] = manager.rdf.Filter(
            f"std::isfinite({column})"
        ).Histo1D(model, column)

    columns = {str(column) for column in manager.rdf.GetColumnNames()}
    pair = resolve_mcp_pair(run, args, columns)
    if pair is None:
        return toa_results, None, None

    mcp_dt_col = "__mcp_delta_t_ns"
    mcp_rdf = manager.rdf.Define(
        mcp_dt_col,
        f"({pair.second}_TS_cfd_ref - {pair.first}_TS_cfd_ref) * {_TS_TO_NS}",
    )
    model = ROOT.RDF.TH1DModel(
        f"mcp_delta_t_run{run}",
        f"Run {run};{pair.label} [ns];Events",
        args.mcp_dt_bins,
        args.mcp_dt_range[0],
        args.mcp_dt_range[1],
    )
    mcp_dt_result = mcp_rdf.Filter(
        f"std::isfinite({mcp_dt_col})"
    ).Histo1D(model, mcp_dt_col)
    return toa_results, mcp_dt_result, pair

def clone_hist(result: ROOT.RDF.RResultPtr, name: str) -> ROOT.TH1:
    """Detach an RDataFrame histogram from its owning dataframe/tree."""
    source = result.GetValue()
    hist = source.Clone(name)
    hist.SetDirectory(0)
    return hist


def fit_gaussian(hist: ROOT.TH1, stem: str) -> FitSummary:
    """Fit the populated central peak of an MCP delta-t histogram."""
    n_entries = int(hist.GetEntries())
    if n_entries < 25 or hist.Integral() <= 0:
        return FitSummary(n_entries, None, None, None, None, False)

    axis = hist.GetXaxis()
    peak = axis.GetBinCenter(hist.GetMaximumBin())
    rms = hist.GetRMS()
    if not math.isfinite(rms) or rms <= 0.0:
        return FitSummary(n_entries, None, None, None, None, False)

    fit_half_width = min(2.5 * rms, 0.45 * (axis.GetXmax() - axis.GetXmin()))
    fit_lo = max(axis.GetXmin(), peak - fit_half_width)
    fit_hi = min(axis.GetXmax(), peak + fit_half_width)
    if fit_hi <= fit_lo:
        return FitSummary(n_entries, None, None, None, None, False)

    fit = ROOT.TF1(f"fit_{safe_name(stem)}", "gaus", fit_lo, fit_hi)
    fit.SetParameters(hist.GetMaximum(), peak, max(rms, axis.GetBinWidth(1)))
    fit.SetLineWidth(2)

    status = int(hist.Fit(fit, "Q0RS"))
    if status != 0 or fit.GetParameter(2) <= 0.0:
        return FitSummary(n_entries, None, None, None, None, False)

    return FitSummary(
        n_entries=n_entries,
        mean_ns=float(fit.GetParameter(1)),
        mean_err_ns=float(fit.GetParError(1)),
        sigma_ns=float(abs(fit.GetParameter(2))),
        sigma_err_ns=float(fit.GetParError(2)),
        fit_ok=True,
    )


def draw_label(canvas: ROOT.TCanvas, title: str, subtitle: str) -> None:
    """Draw a compact CaloX-style plot header without a Preliminary tag."""
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.14, 0.94, "CaloX Data")
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    latex.DrawLatex(0.14, 0.895, title)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.94, 0.94, subtitle)
    canvas.Update()


def histogram_for_plot(hist: ROOT.TH1, normalise: bool) -> ROOT.TH1:
    out = hist.Clone(f"plot_{hist.GetName()}")
    out.SetDirectory(0)
    if normalise and out.Integral() > 0.0:
        out.Scale(1.0 / out.Integral())
    return out


def make_overlay_canvas(
    histograms: Iterable[tuple[int, ROOT.TH1]],
    canvas_name: str,
    title: str,
    x_title: str,
    y_title: str,
    legend_lines: dict[int, str],
    normalise: bool,
    fit_summaries: dict[int, FitSummary] | None = None,
    selection_label: str = "MCP-clean events",
) -> tuple[ROOT.TCanvas, str] | None:
    """Build one interactive JSROOT canvas; do not write PNG or PDF files."""
    entries = [
        (run, histogram_for_plot(hist, normalise))
        for run, hist in histograms
        if hist.GetEntries() > 0
    ]
    if not entries:
        return None

    canvas = ROOT.TCanvas(canvas_name, "", 1200, 850)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    ymax = max(hist.GetMaximum() for _, hist in entries)
    first = True
    legend = ROOT.TLegend(0.52, 0.62, 0.94, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.027)

    draw_hists: list[ROOT.TH1] = []
    fits: list[ROOT.TF1] = []
    for index, (run, hist) in enumerate(entries):
        color = ROOT_COLORS[index % len(ROOT_COLORS)]
        hist.SetLineColor(color)
        hist.SetLineWidth(3)
        hist.SetFillStyle(0)
        hist.SetStats(0)
        hist.GetXaxis().SetTitle(x_title)
        hist.GetYaxis().SetTitle(y_title)
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetLabelSize(0.04)
        hist.GetYaxis().SetTitleOffset(1.25)
        hist.SetMaximum(1.28 * ymax)
        hist.Draw("HIST" if first else "HIST SAME")
        first = False
        draw_hists.append(hist)
        legend.AddEntry(hist, legend_lines[run], "l")

        # For Delta-t only, draw the fit curve in the same run colour.
        if fit_summaries is not None:
            fit_summary = fit_summaries.get(run)
            if (
                fit_summary
                and fit_summary.fit_ok
                and fit_summary.mean_ns is not None
                and fit_summary.sigma_ns is not None
            ):
                func = ROOT.TF1(
                    f"plot_fit_{safe_name(hist.GetName())}",
                    "[0]*exp(-0.5*((x-[1])/[2])^2)",
                    hist.GetXaxis().GetXmin(),
                    hist.GetXaxis().GetXmax(),
                )
                func.SetParameters(hist.GetMaximum(), fit_summary.mean_ns, fit_summary.sigma_ns)
                func.SetLineColor(color)
                func.SetLineStyle(2)
                func.SetLineWidth(2)
                func.Draw("SAME")
                fits.append(func)

    draw_label(canvas, title, selection_label)
    legend.Draw()
    canvas.RedrawAxis()
    canvas.Update()

    # Serialize while all draw objects are alive.  The repository's
    # generate_jsroot_html() embeds this JSON into one interactive HTML gallery.
    canvas_json = str(ROOT.TBufferJSON.ToJSON(canvas).Data())
    canvas._keepalive = (draw_hists, fits, legend)  # prevent PyROOT GC before Write()
    return canvas, canvas_json


def write_jsroot_gallery(
    outdir: Path,
    canvas_keys: list[str],
    canvas_jsons: dict[str, str],
    args: argparse.Namespace,
    selection_label: str,
) -> Path:
    """Use the tb2026 JSROOT HTML generator; it writes no static image files."""
    html_path = outdir / "timing_overlays.html"
    intro = "\n".join(
        [
            f"**Runs:** {', '.join(map(str, args.runs))}",
            f"**Selection:** {selection_label}" + (f" + {args.particle}" if args.particle else ""),
            f"**TOA axis:** [{args.toa_range[0]:g}, {args.toa_range[1]:g}] ns",
            f"**MCP Delta-t axis:** [{args.mcp_dt_range[0]:g}, {args.mcp_dt_range[1]:g}] ns",
            "* Each overlay is unit-area normalized unless `--raw-yields` was used.",
            "* MCP resolution reports use sigma(MCP) = sigma(Delta-t) / sqrt(2).",
            "* Raw TH1 histograms and the drawn canvases are also stored in `timing_overlays.root`.",
        ]
    )
    generate_jsroot_html(
        canvas_keys,
        canvas_jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title="CaloX TB2026 multi-run timing overlays",
        intro_text=intro,
    )
    return html_path


def json_fit_summary(summary: FitSummary) -> dict[str, float | int | bool | None]:
    return {
        "n_entries": summary.n_entries,
        "fit_ok": summary.fit_ok,
        "mean_ns": summary.mean_ns,
        "mean_err_ns": summary.mean_err_ns,
        "sigma_delta_t_ns": summary.sigma_ns,
        "sigma_delta_t_err_ns": summary.sigma_err_ns,
        "single_mcp_sigma_ns": summary.single_mcp_sigma_ns,
        "single_mcp_sigma_err_ns": summary.single_mcp_sigma_err_ns,
    }


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    requested = set(args.channels) if args.channels else None
    root_path = outdir / "timing_overlays.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot create {root_path}")

    print("\n" + "=" * 78)
    print("CaloX TB2026 multi-run timing overlay")
    print(f"Runs              : {args.runs}")
    print(f"JSON file         : {args.json_file}")
    print(f"MCP clean         : {not args.no_mcp_clean}")
    print(f"Particle selection: {args.particle or 'none'}")
    print(f"Output            : {outdir}")
    print("=" * 78)

    # {branch: {run: raw detached TH1}}
    toa_hists: dict[str, dict[int, ROOT.TH1]] = {}
    channel_labels: dict[str, str] = {}
    mcp_dt_hists: dict[int, ROOT.TH1] = {}
    mcp_pairs: dict[int, MCPPair] = {}
    mcp_reference_applied: dict[int, bool] = {}
    fit_summaries: dict[int, FitSummary] = {}

    for run in args.runs:
        beam_type, beam_energy = getRunInfo(run)
        print(f"\n[RUN {run}] beam={beam_type}, energy={beam_energy} GeV")
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)

        # Match the native service-DRS workflow: define all per-MCP clean
        # flags first, then optionally filter on the mcp_clean particle recipe.
        # Deliberately do NOT call apply_mcp_diff_selection(), because it cuts
        # on Delta-t and would artificially narrow the fitted resolution.
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        live_columns = {str(column) for column in manager.rdf.GetColumnNames()}
        mcp_reference_applied[run] = mcp_reference_is_live(run, live_columns)
        if not mcp_reference_applied[run]:
            print(
                f"[WARN] Run {run}: configured MCP reference {MCP_REF} is not live. "
                "*_Time_cfd_mcp is present but has zero MCP correction for this run."
            )

        info = active_timing_channels(manager, requested)
        if not info:
            print(f"[WARN] Run {run}: no mapped channels expose *_Time_cfd_mcp. Skipping.")
            continue

        for branch, meta in info.items():
            channel_labels.setdefault(branch, meta.label)

        results, mcp_result, mcp_pair = book_run_histograms(manager, run, args, info)
        all_actions = list(results.values())
        if mcp_result is not None:
            all_actions.append(mcp_result)
        ROOT.RDF.RunGraphs(all_actions)

        run_dir = root_file.mkdir(f"run{run}")
        toa_dir = run_dir.mkdir("toa")
        toa_dir.cd()
        for branch, result in results.items():
            hist = clone_hist(result, safe_name(branch))
            hist.SetTitle(f"Run {run}: {channel_labels[branch]};Time of Arrival [ns];Events")
            hist.Write()
            toa_hists.setdefault(branch, {})[run] = hist

        if mcp_result is not None and mcp_pair is not None:
            run_dir.cd()
            delta_hist = clone_hist(mcp_result, "mcp_delta_t")
            delta_hist.SetTitle(f"Run {run}: {mcp_pair.label};#Delta t [ns];Events")
            delta_hist.Write()
            mcp_dt_hists[run] = delta_hist
            mcp_pairs[run] = mcp_pair
            fit_summaries[run] = fit_gaussian(delta_hist, f"run{run}_mcp_delta_t")

            summary = fit_summaries[run]
            if summary.fit_ok and summary.sigma_ns is not None:
                print(
                    f"[RUN {run}] MCP #Delta t: #sigma = {1000.0 * summary.sigma_ns:.1f} ps; "
                    f"single-MCP = {1000.0 * summary.single_mcp_sigma_ns:.1f} ps "
                    f"(N={summary.n_entries})"
                )
            else:
                print(f"[WARN] Run {run}: MCP delta-t Gaussian fit failed (N={summary.n_entries}).")

    if not toa_hists:
        raise RuntimeError("No TOA histograms were produced. Check runs, JSON paths, and MCP availability.")

    normalise = not args.raw_yields
    selection_label = "MCP-clean events" if not args.no_mcp_clean else "All events"
    channel_order = sorted(toa_hists)

    # Canvases are written to the ROOT file and serialized into a single native
    # JSROOT gallery.  No PNG/PDF files or plots/ directory are produced.
    canvas_dir = root_file.mkdir("canvases")
    canvas_keys: list[str] = []
    canvas_jsons: dict[str, str] = {}

    mcp_plot_exists = bool(mcp_dt_hists)
    if mcp_plot_exists:
        overlays = sorted(mcp_dt_hists.items())
        legend_lines: dict[int, str] = {}
        for run, hist in overlays:
            summary = fit_summaries.get(run)
            pair = mcp_pairs[run]
            pair_prefix = f"{pair.label}: "
            if summary and summary.fit_ok and summary.sigma_ns is not None:
                legend_lines[run] = (
                    f"Run {run}: {pair_prefix}#sigma_{{#Delta t}}="
                    f"{1000.0 * summary.sigma_ns:.1f} ps, "
                    f"#sigma_{{MCP}}={1000.0 * summary.single_mcp_sigma_ns:.1f} ps"
                )
            else:
                legend_lines[run] = (
                    f"Run {run}: {pair_prefix}N={int(hist.GetEntries())}, fit unavailable"
                )

        pair_labels = {pair.label for pair in mcp_pairs.values()}
        mcp_title = (
            f"{next(iter(pair_labels))} timing difference"
            if len(pair_labels) == 1
            else "MCP timing difference (native per-run pairs)"
        )
        mcp_x_title = (
            f"{next(iter(pair_labels))} [ns]"
            if len(pair_labels) == 1
            else "MCP #Delta t [ns]"
        )

        canvas_name = "mcp_delta_t_overlay"
        built = make_overlay_canvas(
            histograms=overlays,
            canvas_name=canvas_name,
            title=mcp_title,
            x_title=mcp_x_title,
            y_title="Normalized Events" if normalise else "Events",
            legend_lines=legend_lines,
            normalise=normalise,
            fit_summaries=fit_summaries,
            selection_label=selection_label,
        )
        if built is not None:
            canvas, canvas_json = built
            canvas_dir.cd()
            canvas.Write(canvas_name)
            canvas_keys.extend([canvas_name, "NEWLINE"])
            canvas_jsons[canvas_name] = canvas_json

    for branch in channel_order:
        overlays = sorted(toa_hists[branch].items())
        legend_lines = {
            run: (
                f"Run {run}: N={int(hist.GetEntries())}, "
                f"#mu={hist.GetMean():.3f} ns, RMS={1000.0 * hist.GetRMS():.0f} ps"
            )
            for run, hist in overlays
        }
        canvas_name = f"toa_{safe_name(branch)}"
        built = make_overlay_canvas(
            histograms=overlays,
            canvas_name=canvas_name,
            title=channel_labels.get(branch, branch),
            x_title="Time of Arrival [ns]",
            y_title="Normalized Events" if normalise else "Events",
            legend_lines=legend_lines,
            normalise=normalise,
            selection_label=selection_label,
        )
        if built is not None:
            canvas, canvas_json = built
            canvas_dir.cd()
            canvas.Write(canvas_name)
            canvas_keys.append(canvas_name)
            canvas_jsons[canvas_name] = canvas_json

    root_file.Close()
    html_path = write_jsroot_gallery(outdir, canvas_keys, canvas_jsons, args, selection_label)

    report = {
        "runs": args.runs,
        "json_file": args.json_file,
        "mcp_clean_selection": not args.no_mcp_clean,
        "particle_selection": args.particle,
        "toa_column": "<mapped DRS branch>_Time_cfd_mcp",
        "mcp_reference_name": MCP_REF,
        "mcp_reference_applied_by_run": {
            str(run): mcp_reference_applied.get(run, False) for run in args.runs
        },
        "mcp_delta_t_pairs": {
            str(run): [pair.first, pair.second] for run, pair in mcp_pairs.items()
        },
        "mcp_delta_t_definition_by_run": {
            str(run): (
                f"({pair.second}_TS_cfd_ref - {pair.first}_TS_cfd_ref) * {_TS_TO_NS}"
            )
            for run, pair in mcp_pairs.items()
        },
        "toa_range_ns": list(args.toa_range),
        "mcp_delta_t_range_ns": list(args.mcp_dt_range),
        "channels": {branch: channel_labels.get(branch, branch) for branch in channel_order},
        "mcp_delta_t_fits": {str(run): json_fit_summary(summary) for run, summary in fit_summaries.items()},
    }
    (outdir / "timing_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("\n[DONE]")
    print(f"  ROOT histograms : {root_path}")
    print(f"  Interactive HTML: {html_path}")
    print(f"  JSON summary    : {outdir / 'timing_summary.json'}")


if __name__ == "__main__":
    main()