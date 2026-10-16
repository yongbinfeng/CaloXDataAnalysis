#!/usr/bin/env python3
"""TB2026 timing resolution from horizontally-adjacent same-material channel pairs.

Reproduces the "timing resolution vs. beam position" method used in the
reference position-scan study: for every pair of 3 mm channels that (a) sit in
the same DRS readout group, (b) are physically adjacent along x (same tower
row), and (c) share the same scintillator family (CerQuartz-CerQuartz,
CerPlastic-CerPlastic, or Sci-Sci), the time difference

    dt = (TS_cfd_ref[higher-x channel] - TS_cfd_ref[lower-x channel]) * 0.2 ns

TS_cfd_ref is the CFD (constant-fraction) threshold-crossing time, linearly
interpolated between DRS samples (utils/functions.cc:compute_cfd_integral) --
this is the "interpolated peak position" of the reference method. (TS_peak is
the raw un-interpolated argmax sample index and quantizes dt to multiples of
one 0.2 ns DRS sample, so it is not used here.)

This is histogrammed per run and fit with a Gaussian. Because both channels
share the same DRS board+group, the group's LED reference-channel jitter
cancels exactly in the difference (see variables/drs.py:get_ref_ts_name), and
no MCP channel is live for this run range (data/Runlist.json / channel_map.py),
so no further correction is needed. The reported per-pair time resolution is
sigma(dt)/sqrt(2), i.e. the single-channel-equivalent resolution.

Electron identification uses the Cherenkov-counter beam-PID cut
(CaloXAnalysisManager.apply_beam_pid_selection), applied by default.

Outputs (under --outdir):
  pair_timing_resolution.root    raw dt histograms, Gaussian fits, and
                                  per-run/per-material resolution-spread
                                  histograms
  pair_timing_resolution.html    interactive JSROOT gallery
  pair_timing_resolution_summary.json  one record per (run, channel pair)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.data_loader import getRunInfo
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var
from variables.drs import _TS_TO_NS

MATERIAL_DISPLAY = {"CerQuartz": "CerQuartz", "CerPlastic": "CerPlastic", "Sci": "Sci"}
MATERIAL_COLOR = {
    "CerQuartz": ROOT.kOrange + 7,
    "CerPlastic": ROOT.kAzure + 2,
    "Sci": ROOT.kGreen + 2,
}


@dataclass(frozen=True)
class ChannelPair:
    material: str
    bridge_no: Any
    board_no: Any
    group_no: Any
    i_tower_y: float
    branch_low: str
    branch_high: str
    ix_low: float
    ix_high: float

    @property
    def pair_id(self) -> str:
        return "Brg{}_Board{}_Group{}_y{}_{}_x{}to{}".format(
            self.bridge_no, self.board_no, self.group_no, self.i_tower_y,
            self.material, self.ix_low, self.ix_high,
        )

    @property
    def label(self) -> str:
        return "{} pair: B{} G{} y={} (x {} vs {})".format(
            MATERIAL_DISPLAY.get(self.material, self.material),
            self.board_no, self.group_no, self.i_tower_y, self.ix_low, self.ix_high,
        )


@dataclass
class PairFitSummary:
    n_entries: int
    mean_ns: Optional[float]
    mean_err_ns: Optional[float]
    sigma_ns: Optional[float]
    sigma_err_ns: Optional[float]
    fit_ok: bool

    @property
    def resolution_ps(self) -> Optional[float]:
        if self.sigma_ns is None:
            return None
        return 1000.0 * self.sigma_ns / math.sqrt(2.0)

    @property
    def resolution_err_ps(self) -> Optional[float]:
        if self.sigma_err_ns is None:
            return None
        return 1000.0 * self.sigma_err_ns / math.sqrt(2.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Time-difference histograms and Gaussian fits for horizontally "
            "adjacent, same-material DRS channel pairs (TB2026 3 mm SiPMs)."
        )
    )
    parser.add_argument(
        "--runs", nargs="+", type=int,
        default=list(range(1897, 1941)),
        help="Run numbers (default: 1897-1940, the electron 40 GeV 90deg position scan).",
    )
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--position-label", default=None)
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument(
        "--particle", default="electron",
        help="Beam-PID particle selection to apply, or 'none' to disable (default: electron).",
    )
    parser.add_argument(
        "--no-mcp-clean", action="store_true",
        help="Disable the MCP clean-pulse event-quality selection.",
    )
    parser.add_argument(
        "--skip-sci", action="store_true",
        help="Skip Sci-Sci pairs (CerQuartz/CerPlastic pairs only).",
    )
    parser.add_argument(
        "--materials", nargs="+", default=["CerQuartz", "CerPlastic", "Sci"],
        choices=["CerQuartz", "CerPlastic", "Sci"],
        help="Material families to include.",
    )
    parser.add_argument("--dt-bins", type=int, default=400)
    parser.add_argument(
        "--dt-range", nargs=2, type=float, default=(-4.0, 4.0),
        metavar=("MIN_NS", "MAX_NS"),
    )
    parser.add_argument("--fit-min-events", type=int, default=25)
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument(
        "--outdir", default=None,
        help="Output directory (default: results/timing_resolution_pairs/runs_<first>_<last>).",
    )
    args = parser.parse_args()

    if len(set(args.runs)) != len(args.runs):
        parser.error("--runs contains a duplicate run number")
    if args.dt_range[0] >= args.dt_range[1]:
        parser.error("--dt-range requires MIN_NS < MAX_NS")
    if args.skip_sci:
        args.materials = [m for m in args.materials if m != "Sci"]
    if str(args.particle).strip().lower() in ("none", ""):
        args.particle = None
    if args.position_label is None:
        args.position_label = args.position_key
    if args.outdir is None:
        args.outdir = os.path.join(
            "results", "timing_resolution_pairs",
            "runs_{}_{}".format(min(args.runs), max(args.runs)),
        )
    return args


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def cpp_identifier(value: str) -> str:
    """Sanitize to a valid C++ variable name (RDataFrame column names)."""
    ident = re.sub(r"[^A-Za-z0-9_]+", "_", value)
    if ident and ident[0].isdigit():
        ident = "_" + ident
    return ident


def load_run_positions(args: argparse.Namespace) -> Dict[int, float]:
    positions: Dict[int, float] = {}
    if not os.path.exists(args.runlist):
        return positions
    with open(args.runlist) as handle:
        runlist = json.load(handle)
    for run in args.runs:
        info = runlist.get(str(run))
        if info is None or args.position_key not in info:
            continue
        try:
            positions[run] = float(str(info[args.position_key]).split()[0])
        except (ValueError, IndexError):
            continue
    return positions


def find_channel_pairs(
    manager: CaloXAnalysisManager, columns: set, materials: List[str],
) -> List[ChannelPair]:
    """Enumerate x-adjacent, same-material 3 mm channel pairs within each DRS group."""
    groups: Dict[Tuple[Any, Any, Any], List[Any]] = {}
    for _, board in manager.drsboards.items():
        for channel in board:
            if channel.is_reference or getattr(channel, "is6mm", False):
                continue
            branch = channel.get_channel_name(blsub=False)
            if "{}_TS_cfd_ref".format(branch) not in columns:
                continue
            key = (channel.bridge_no, board.board_no, channel.group_no)
            groups.setdefault(key, []).append(channel)

    pairs: List[ChannelPair] = []
    for (bridge_no, board_no, group_no), channels in groups.items():
        rows: Dict[float, List[Any]] = {}
        for channel in channels:
            rows.setdefault(channel.i_tower_y, []).append(channel)
        for i_tower_y, row in rows.items():
            cer = sorted((c for c in row if c.isCer), key=lambda c: c.i_tower_x)
            sci = sorted((c for c in row if not c.isCer), key=lambda c: c.i_tower_x)
            for lst in (cer, sci):
                for i in range(len(lst) - 1):
                    low, high = lst[i], lst[i + 1]
                    material_low = get_channel_var(low)
                    material_high = get_channel_var(high)
                    if material_low != material_high or material_low not in materials:
                        continue
                    pairs.append(ChannelPair(
                        material=material_low,
                        bridge_no=bridge_no, board_no=board_no, group_no=group_no,
                        i_tower_y=i_tower_y,
                        branch_low=low.get_channel_name(blsub=False),
                        branch_high=high.get_channel_name(blsub=False),
                        ix_low=low.i_tower_x, ix_high=high.i_tower_x,
                    ))
    return pairs


def book_pair_histograms(
    manager: CaloXAnalysisManager, run: int, args: argparse.Namespace,
    pairs: List[ChannelPair],
) -> Dict[str, Any]:
    results: Dict[str, Any] = {}
    rdf = manager.rdf
    for pair in pairs:
        column = "__dt_{}".format(cpp_identifier(pair.pair_id))
        rdf = rdf.Define(
            column,
            "({}_TS_cfd_ref - {}_TS_cfd_ref) * {}".format(
                pair.branch_high, pair.branch_low, _TS_TO_NS,
            ),
        )
        model = ROOT.RDF.TH1DModel(
            "dt_run{}_{}".format(run, safe_name(pair.pair_id)),
            "Run {};{} [ns];Events".format(run, pair.label),
            args.dt_bins, args.dt_range[0], args.dt_range[1],
        )
        results[pair.pair_id] = rdf.Filter(
            "std::isfinite({})".format(column)
        ).Histo1D(model, column)
    return results


def in_range_entries(hist: Any) -> int:
    return int(round(hist.Integral(1, hist.GetNbinsX())))


def estimate_peak_sigma(hist: Any) -> Tuple[float, float]:
    """Return peak position and a robust width estimate from local FWHM."""
    axis = hist.GetXaxis()
    peak_bin = hist.GetMaximumBin()
    peak_x = axis.GetBinCenter(peak_bin)
    max_y = hist.GetBinContent(peak_bin)
    bin_width = axis.GetBinWidth(1)

    if max_y <= 0.0:
        return peak_x, max(bin_width, 0.25)

    half = 0.5 * max_y
    left = peak_bin
    while left > 1 and hist.GetBinContent(left) >= half:
        left -= 1
    right = peak_bin
    while right < hist.GetNbinsX() and hist.GetBinContent(right) >= half:
        right += 1

    fwhm = axis.GetBinCenter(right) - axis.GetBinCenter(left)
    if not math.isfinite(fwhm) or fwhm <= bin_width:
        rms = hist.GetRMS()
        if math.isfinite(rms) and rms > 0.0:
            return peak_x, max(bin_width, min(rms, 3.0))
        return peak_x, max(bin_width, 0.25)
    return peak_x, max(bin_width, fwhm / 2.355)


def fit_pair_gaussian(hist: Any, stem: str, min_events: int) -> PairFitSummary:
    n_entries = in_range_entries(hist)
    if n_entries < min_events:
        return PairFitSummary(n_entries, None, None, None, None, False)

    axis = hist.GetXaxis()
    peak, sigma_guess = estimate_peak_sigma(hist)
    half_width = min(2.5 * sigma_guess, 0.45 * (axis.GetXmax() - axis.GetXmin()))
    lo = max(axis.GetXmin(), peak - half_width)
    hi = min(axis.GetXmax(), peak + half_width)
    if hi <= lo:
        return PairFitSummary(n_entries, None, None, None, None, False)

    fit = ROOT.TF1("fit_pair_{}".format(safe_name(stem)), "gaus", lo, hi)
    fit.SetParameters(hist.GetMaximum(), peak, max(sigma_guess, axis.GetBinWidth(1)))
    status = int(hist.Fit(fit, "Q0RS"))
    sigma = abs(float(fit.GetParameter(2)))
    if status != 0 or not math.isfinite(sigma) or sigma <= 0.0:
        return PairFitSummary(n_entries, None, None, None, None, False)
    return PairFitSummary(
        n_entries=n_entries,
        mean_ns=float(fit.GetParameter(1)), mean_err_ns=float(fit.GetParError(1)),
        sigma_ns=sigma, sigma_err_ns=float(fit.GetParError(2)), fit_ok=True,
    )


def draw_header(canvas: Any, title: str, subtitle: str) -> None:
    canvas.cd()
    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(62)
    latex.SetTextSize(0.034)
    latex.SetTextAlign(11)
    latex.DrawLatex(0.14, 0.944, "CaloX Data")
    latex.SetTextFont(42)
    latex.SetTextSize(0.023)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.94, 0.944, subtitle)
    latex.SetTextAlign(11)
    latex.SetTextSize(0.024)
    latex.DrawLatex(0.14, 0.895, title)
    canvas.Update()


def make_pair_canvas(
    hist: Any, pair: ChannelPair, run: int, summary: PairFitSummary, subtitle: str,
) -> Tuple[Any, str]:
    canvas_name = "canvas_run{}_{}".format(run, safe_name(pair.pair_id))
    canvas = ROOT.TCanvas(canvas_name, "", 900, 700)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    hist.SetLineColor(MATERIAL_COLOR.get(pair.material, ROOT.kBlack))
    hist.SetLineWidth(2)
    hist.SetStats(0)
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetLabelSize(0.038)
    hist.GetYaxis().SetLabelSize(0.038)
    hist.GetYaxis().SetTitleOffset(1.35)
    hist.Draw("HIST")

    keepalive: List[Any] = []
    if summary.fit_ok and summary.mean_ns is not None and summary.sigma_ns is not None:
        func = ROOT.TF1(
            "curve_{}".format(canvas_name),
            "[0]*exp(-0.5*((x-[1])/[2])^2)",
            hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax(),
        )
        func.SetParameters(hist.GetMaximum(), summary.mean_ns, summary.sigma_ns)
        func.SetLineColor(ROOT.kRed + 1)
        func.SetLineWidth(2)
        func.Draw("SAME")
        keepalive.append(func)

        info = ROOT.TPaveText(0.15, 0.63, 0.55, 0.87, "NDC")
        info.SetBorderSize(0)
        info.SetFillStyle(0)
        info.SetTextFont(42)
        info.SetTextSize(0.026)
        info.SetTextAlign(12)
        info.AddText("N = {}".format(summary.n_entries))
        info.AddText("mean = {:.4f} ns".format(summary.mean_ns))
        info.AddText("#sigma_{dt} = " + "{:.4f} ns".format(summary.sigma_ns))
        info.AddText("resolution = #sigma_{dt}/#sqrt{2} = " + "{:.1f} ps".format(summary.resolution_ps))
        info.Draw()
        keepalive.append(info)
    else:
        info = ROOT.TPaveText(0.15, 0.75, 0.55, 0.87, "NDC")
        info.SetBorderSize(0)
        info.SetFillStyle(0)
        info.SetTextFont(42)
        info.SetTextSize(0.026)
        info.AddText("N = {} (fit failed or insufficient stats)".format(summary.n_entries))
        info.Draw()
        keepalive.append(info)

    draw_header(canvas, pair.label, subtitle)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def make_resolution_meta_hist(
    name: str, title: str, resolutions_ps: List[float],
) -> Optional[Any]:
    if not resolutions_ps:
        return None
    lo = max(0.0, min(resolutions_ps) - 20.0)
    hi = max(resolutions_ps) + 20.0
    hist = ROOT.TH1D(name, title, 30, lo, hi)
    for value in resolutions_ps:
        hist.Fill(value)
    hist.SetDirectory(0)
    return hist


def make_resolution_distribution_canvas(
    records: List[Dict[str, Any]], materials: List[str], bins: int = 60,
) -> Optional[Tuple[Any, str]]:
    """One overlaid histogram per material: every fitted (run, pair) resolution as one entry."""
    values_by_material = {
        material: [r["resolution_ps"] for r in records if r["fit_ok"] and r["material"] == material]
        for material in materials
    }
    all_values = [v for values in values_by_material.values() for v in values]
    if not all_values:
        return None

    hi = max(all_values) * 1.05
    canvas = ROOT.TCanvas("canvas_resolution_distribution", "", 1000, 750)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    keepalive: List[Any] = []
    legend = ROOT.TLegend(0.55, 0.68, 0.94, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.026)

    ymax = 0.0
    hists = []
    for material in materials:
        values = values_by_material[material]
        if not values:
            continue
        hist = ROOT.TH1D(
            "resdist_{}".format(safe_name(material)),
            ";Time resolution per channel pair [ps];Channel-pair combinations",
            bins, 0.0, hi,
        )
        for value in values:
            hist.Fill(value)
        hist.SetDirectory(0)
        hist.SetLineColor(MATERIAL_COLOR.get(material, ROOT.kBlack))
        hist.SetLineWidth(3)
        hist.SetStats(0)
        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetXaxis().SetLabelSize(0.038)
        hist.GetYaxis().SetLabelSize(0.038)
        hist.GetYaxis().SetTitleOffset(1.3)
        ymax = max(ymax, hist.GetMaximum())
        hists.append((material, hist))

    for index, (material, hist) in enumerate(hists):
        hist.SetMaximum(1.15 * ymax)
        hist.Draw("HIST" if index == 0 else "HIST SAME")
        keepalive.append(hist)
        legend.AddEntry(
            hist, "{} (N={})".format(MATERIAL_DISPLAY.get(material, material), len(values_by_material[material])), "l",
        )

    legend.Draw()
    keepalive.append(legend)
    draw_header(
        canvas, "Distribution of per-pair time resolutions (all runs, all combinations)",
        "e+ 40 GeV, 90 deg",
    )
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def make_scan_overlay_canvas(
    records: List[Dict[str, Any]], materials: List[str], position_label: str,
) -> Optional[Tuple[Any, str]]:
    entries = [r for r in records if r["fit_ok"] and r["resolution_ps"] is not None]
    if not entries:
        return None

    canvas = ROOT.TCanvas("canvas_resolution_vs_position", "", 1100, 800)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    xs = [r["position_mm"] for r in entries]
    ys = [r["resolution_ps"] for r in entries]
    frame = ROOT.TH1F(
        "frame_resolution_vs_position", ";{} [mm];Time resolution per channel pair [ps]".format(position_label),
        100, min(xs) - 20.0, max(xs) + 20.0,
    )
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.15 * max(ys))
    frame.SetStats(0)
    frame.GetXaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.3)
    frame.Draw("AXIS")

    keepalive: List[Any] = [frame]
    legend = ROOT.TLegend(0.6, 0.68, 0.94, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.026)

    for material in materials:
        mat_entries = [r for r in entries if r["material"] == material]
        if not mat_entries:
            continue
        color = MATERIAL_COLOR.get(material, ROOT.kBlack)

        scatter = ROOT.TGraph(len(mat_entries))
        for i, r in enumerate(mat_entries):
            scatter.SetPoint(i, r["position_mm"], r["resolution_ps"])
        scatter.SetMarkerStyle(24)
        scatter.SetMarkerSize(0.8)
        scatter.SetMarkerColor(color)
        scatter.Draw("P SAME")
        keepalive.append(scatter)

        by_position: Dict[float, List[float]] = {}
        for r in mat_entries:
            by_position.setdefault(r["position_mm"], []).append(r["resolution_ps"])
        positions = sorted(by_position)
        trend = ROOT.TGraphErrors(len(positions))
        for i, pos in enumerate(positions):
            values = by_position[pos]
            mean = sum(values) / len(values)
            spread = (sum((v - mean) ** 2 for v in values) / len(values)) ** 0.5 if len(values) > 1 else 0.0
            trend.SetPoint(i, pos, mean)
            trend.SetPointError(i, 0.0, spread / max(1, len(values)) ** 0.5)
        trend.SetLineColor(color)
        trend.SetLineWidth(3)
        trend.SetMarkerColor(color)
        trend.SetMarkerStyle(20)
        trend.SetMarkerSize(1.1)
        trend.Draw("LP SAME")
        keepalive.append(trend)

        legend.AddEntry(trend, "{} (mean over pairs)".format(MATERIAL_DISPLAY.get(material, material)), "lp")
        legend.AddEntry(scatter, "{} (per pair)".format(MATERIAL_DISPLAY.get(material, material)), "p")

    legend.Draw()
    draw_header(canvas, "Time resolution vs. beam position (per adjacent-pair, no p.e. binning)", "e+ 40 GeV, 90 deg")
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def make_manager_args(args: argparse.Namespace, run: int) -> argparse.Namespace:
    return argparse.Namespace(
        run=run, first_event=args.first_event, last_event=args.last_event,
        json_file=args.json_file,
    )


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    positions = load_run_positions(args)

    root_path = outdir / "pair_timing_resolution.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot create {}".format(root_path))

    print("\n" + "=" * 78)
    print("CaloX TB2026 adjacent-channel-pair timing resolution")
    print("Runs              : {}".format(args.runs))
    print("Materials         : {}".format(args.materials))
    print("Particle selection: {}".format(args.particle or "none"))
    print("Output            : {}".format(outdir))
    print("=" * 78)

    records: List[Dict[str, Any]] = []
    canvas_keys: List[str] = []
    canvas_jsons: Dict[str, str] = {}

    for run in args.runs:
        beam_type, beam_energy = getRunInfo(run)
        position_mm = positions.get(run)
        print("\n[RUN {}] beam={}, energy={} GeV, position={} mm".format(
            run, beam_type, beam_energy, position_mm))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)

        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        columns = {str(column) for column in manager.rdf.GetColumnNames()}
        pairs = find_channel_pairs(manager, columns, args.materials)
        if not pairs:
            print("[WARN] Run {}: no qualifying channel pairs found; skipped.".format(run))
            continue
        print("[INFO] Run {}: booking {} channel-pair histograms.".format(run, len(pairs)))

        results = book_pair_histograms(manager, run, args, pairs)
        ROOT.RDF.RunGraphs(list(results.values()))

        run_dir = root_file.mkdir("run{}".format(run))
        resolutions_by_material: Dict[str, List[float]] = {m: [] for m in args.materials}

        for pair in pairs:
            hist = results[pair.pair_id].GetValue().Clone(safe_name(pair.pair_id))
            hist.SetDirectory(0)
            mat_dir_name = "run{}/{}".format(run, pair.material)
            mat_dir = root_file.GetDirectory(mat_dir_name) or run_dir.mkdir(pair.material)
            mat_dir.cd()
            hist.Write()

            summary = fit_pair_gaussian(hist, "{}_{}".format(run, pair.pair_id), args.fit_min_events)
            if summary.fit_ok:
                mat_dir.cd()
                gaus_fit = hist.GetFunction("fit_pair_{}".format(safe_name("{}_{}".format(run, pair.pair_id))))
                if gaus_fit is not None:
                    gaus_fit.Write("fit_{}".format(safe_name(pair.pair_id)))
                resolutions_by_material[pair.material].append(summary.resolution_ps)
                print("[FIT] Run {} {}: mean={:.4f} ns, sigma={:.4f} ns, resolution={:.1f} ps, N={}".format(
                    run, pair.label, summary.mean_ns, summary.sigma_ns, summary.resolution_ps, summary.n_entries))
            else:
                print("[WARN] Run {} {}: fit failed (N={}).".format(run, pair.label, summary.n_entries))

            records.append({
                "run": run, "position_mm": position_mm, "pair_id": pair.pair_id,
                "material": pair.material,
                "bridge_no": pair.bridge_no, "board_no": pair.board_no, "group_no": pair.group_no,
                "i_tower_y": pair.i_tower_y,
                "branch_low": pair.branch_low, "branch_high": pair.branch_high,
                "ix_low": pair.ix_low, "ix_high": pair.ix_high,
                "n_entries": summary.n_entries, "fit_ok": summary.fit_ok,
                "mean_ns": summary.mean_ns, "mean_err_ns": summary.mean_err_ns,
                "sigma_ns": summary.sigma_ns, "sigma_err_ns": summary.sigma_err_ns,
                "resolution_ps": summary.resolution_ps, "resolution_err_ps": summary.resolution_err_ps,
            })

            canvas, canvas_json = make_pair_canvas(
                hist, pair, run, summary,
                "run {}, x={} mm".format(run, position_mm) if position_mm is not None else "run {}".format(run),
            )
            canvas_keys.append(canvas.GetName())
            canvas_jsons[canvas.GetName()] = canvas_json

        canvas_keys.append("NEWLINE")

        meta_dir = root_file.GetDirectory("resolution_meta") or root_file.mkdir("resolution_meta")
        meta_run_dir = meta_dir.mkdir("run{}".format(run))
        for material in args.materials:
            values = resolutions_by_material[material]
            hist = make_resolution_meta_hist(
                "res_run{}_{}".format(run, material),
                "Run {} {} pair resolutions;Resolution [ps];Pairs".format(run, material),
                values,
            )
            if hist is not None:
                meta_run_dir.cd()
                hist.Write()

    if not records:
        raise RuntimeError("No channel-pair histograms were produced for any requested run.")

    all_meta_dir = root_file.GetDirectory("resolution_meta").mkdir("all_runs")
    for material in args.materials:
        values = [r["resolution_ps"] for r in records if r["fit_ok"] and r["material"] == material]
        hist = make_resolution_meta_hist(
            "res_allruns_{}".format(material),
            "All runs {} pair resolutions;Resolution [ps];Pairs".format(material),
            values,
        )
        if hist is not None:
            all_meta_dir.cd()
            hist.Write()

    scan_canvas = make_scan_overlay_canvas(records, args.materials, args.position_label)
    if scan_canvas is not None:
        canvas, canvas_json = scan_canvas
        canvas_keys.insert(0, "NEWLINE")
        canvas_keys.insert(0, canvas.GetName())
        canvas_jsons[canvas.GetName()] = canvas_json

    distribution_canvas = make_resolution_distribution_canvas(records, args.materials)
    if distribution_canvas is not None:
        canvas, canvas_json = distribution_canvas
        canvas_keys.insert(0, "NEWLINE")
        canvas_keys.insert(0, canvas.GetName())
        canvas_jsons[canvas.GetName()] = canvas_json

    root_file.Write()
    root_file.Close()

    summary_path = outdir / "pair_timing_resolution_summary.json"
    with open(summary_path, "w") as handle:
        json.dump({
            "runs": args.runs, "materials": args.materials,
            "position_key": args.position_key, "particle": args.particle,
            "records": records,
        }, handle, indent=2)
    print("\n[WROTE] {}".format(summary_path))

    n_ok = sum(1 for r in records if r["fit_ok"])
    intro_lines = [
        "**Method**: for each pair of horizontally-adjacent, same-material 3 mm "
        "channels within a single DRS readout group, dt = TS_peak_ref(high-x) - "
        "TS_peak_ref(low-x) is histogrammed per run and fit with a Gaussian. The "
        "group's reference-channel jitter cancels exactly in this difference. "
        "Quoted resolution = sigma(dt) / sqrt(2), the single-channel-equivalent "
        "time resolution.",
        "* Electron ID: Cherenkov-counter beam-PID cut ({})".format(args.particle or "none"),
        "* Runs: {}-{} ({} runs), position key: {}".format(
            min(args.runs), max(args.runs), len(args.runs), args.position_key),
        "* Materials: {}".format(", ".join(args.materials)),
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
    print("[WROTE] {}".format(root_path))


if __name__ == "__main__":
    main()
