#!/usr/bin/env python3
"""TB2026 multi-run timing overlays plus per-channel velocity fits.

For every active calorimeter DRS channel this script:
  1. builds one MCP-corrected TOA histogram per requested run;
  2. overlays the runs in a single interactive JSROOT HTML gallery;
  3. fits each run histogram before the overlay is drawn:
       * CerQuartz / CerPlastic: Gaussian;
       * Sci: Landau convolved with a Gaussian (Langau);
  4. uses the fitted distribution mean as the timing point and
     sigma_eff / sqrt(N) as its statistical uncertainty;
  5. fits TOA versus the run scan coordinate per channel and reports
     v = 0.1 / |dt/dz| in cm/ns when the coordinate is supplied in mm.

The scan-coordinate values are read from data/Runlist.json by default.  New
TB2026 runs are not assumed to exist in the upstream repository: this script
uses the user's local Runlist.json.  A small --positions-json file can override
that source when needed.

No PNGs or PDFs are written.  The outputs are only:
  timing_overlays.html       TOA and MCP interactive JSROOT gallery
  timing_velocity_fits.html   per-channel velocity interactive JSROOT gallery
  timing_overlays.root        raw histograms, fitted objects, and canvases
  timing_summary.json         machine-readable fit and velocity summaries
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import ROOT

from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from utils.data_loader import getRunInfo
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var
from variables.drs import MCP_REF, _TS_TO_NS


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
    ROOT.kSpring + 5,
]

MATERIAL_DISPLAY = {
    "CerQuartz": "CerQuartz",
    "CerPlastic": "CerPlastic",
    "Sci": "Sci",
}

_LANGAU_DECLARED = False


@dataclass(frozen=True)
class ChannelInfo:
    branch: str
    label: str
    short_label: str
    material: str
    size: str


@dataclass(frozen=True)
class MCPPair:
    first: str
    second: str

    @property
    def label(self) -> str:
        return f"{self.second} - {self.first}"


@dataclass
class MCPFitSummary:
    n_entries: int
    mean_ns: Optional[float]
    mean_err_ns: Optional[float]
    sigma_ns: Optional[float]
    sigma_err_ns: Optional[float]
    fit_ok: bool

    @property
    def single_mcp_sigma_ns(self) -> Optional[float]:
        if self.sigma_ns is None:
            return None
        return self.sigma_ns / math.sqrt(2.0)

    @property
    def single_mcp_sigma_err_ns(self) -> Optional[float]:
        if self.sigma_err_ns is None:
            return None
        return self.sigma_err_ns / math.sqrt(2.0)


@dataclass
class TimingFitSummary:
    model: str
    material: str
    n_entries: int
    mean_ns: Optional[float]
    point_error_ns: Optional[float]
    sigma_ns: Optional[float]
    sigma_err_ns: Optional[float]
    chi2_ndf: Optional[float]
    fit_lo_ns: Optional[float]
    fit_hi_ns: Optional[float]
    fit_ok: bool
    function: Any = None


@dataclass
class VelocityFitSummary:
    n_points: int
    intercept_ns: Optional[float]
    intercept_err_ns: Optional[float]
    slope_ns_per_mm: Optional[float]
    slope_err_ns_per_mm: Optional[float]
    speed_cm_per_ns: Optional[float]
    speed_err_cm_per_ns: Optional[float]
    chi2_ndf: Optional[float]
    fit_ok: bool


@dataclass(frozen=True)
class VelocityPoint:
    run: int
    position_mm: float
    mean_ns: float
    point_error_ns: float
    model: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make per-channel multi-run TOA overlays, timing fits, and "
            "per-channel TOA-versus-position velocity fits for TB2026."
        )
    )
    parser.add_argument(
        "--runs", nargs="+", type=int, required=True,
        help="Run numbers to compare. Files are read through --json-file.",
    )
    parser.add_argument(
        "--json-file", default="data/datafiles_local.json",
        help="Run-to-ROOT-file JSON used by CaloXAnalysisManager.",
    )
    parser.add_argument(
        "--first-event", type=int, default=0,
        help="First event index, inclusive (default: 0).",
    )
    parser.add_argument(
        "--last-event", type=int, default=-1,
        help="Last event index, exclusive; -1 means all events (default: -1).",
    )
    parser.add_argument(
        "--particle", default=None,
        help="Optional CaloX particle selection, for example electron or pion.",
    )
    parser.add_argument(
        "--no-mcp-clean", action="store_true",
        help="Disable the MCP clean-pulse selection.",
    )
    parser.add_argument(
        "--channels", nargs="+", default=None, metavar="DRS_BRANCH",
        help="Optional exact DRS branch subset. Default: all mapped calorimeter channels.",
    )
    parser.add_argument(
        "--toa-bins", type=int, default=400,
        help="TOA histogram bins (default: 400).",
    )
    parser.add_argument(
        "--toa-range", nargs=2, type=float, default=(-30.0, 30.0),
        metavar=("MIN_NS", "MAX_NS"),
        help="TOA axis range in ns (default: -30 30).",
    )
    parser.add_argument(
        "--mcp-dt-bins", type=int, default=400,
        help="MCP delta-t histogram bins (default: 400).",
    )
    parser.add_argument(
        "--mcp-dt-range", nargs=2, type=float, default=(-5.0, 5.0),
        metavar=("MIN_NS", "MAX_NS"),
        help="MCP delta-t axis range in ns (default: -5 5).",
    )
    parser.add_argument(
        "--mcp-a", default=None,
        help="Optional logical first-MCP override; must be passed with --mcp-b.",
    )
    parser.add_argument(
        "--mcp-b", default=None,
        help="Optional logical second-MCP override; must be passed with --mcp-a.",
    )
    parser.add_argument(
        "--raw-yields", action="store_true",
        help="Draw raw counts instead of unit-area-normalized TOA overlays.",
    )
    parser.add_argument(
        "--runlist", default="data/Runlist.json",
        help="Local run metadata JSON used to obtain scan positions.",
    )
    parser.add_argument(
        "--position-key", default="Table position Z",
        help=(
            "Field name in --runlist containing the scan coordinate in mm "
            "(default: 'Table position Z')."
        ),
    )
    parser.add_argument(
        "--positions-json", default=None,
        help=(
            "Optional override JSON of the form {\"1898\": -50.0, ...}; "
            "values are interpreted as mm."
        ),
    )
    parser.add_argument(
        "--positions", nargs="+", default=None, metavar="RUN:MM",
        help=(
            "Direct run-position override, for example "
            "--positions 1898:-100 1903:-50 1908:0 1913:50 1918:100. "
            "Values are interpreted as mm and override --runlist."
        ),
    )
    parser.add_argument(
        "--position-label", default=None,
        help="Axis label before '[mm]' on velocity plots. Default: --position-key.",
    )
    parser.add_argument(
        "--no-velocity-fits", action="store_true",
        help="Skip TOA-versus-position velocity fits.",
    )
    parser.add_argument(
        "--fit-min-events", type=int, default=25,
        help="Minimum in-range entries needed for a per-run timing fit (default: 25).",
    )
    parser.add_argument(
        "--fit-max-half-width", type=float, default=5.0,
        help="Largest symmetric Gaussian fit half-width in ns (default: 5).",
    )
    parser.add_argument(
        "--n-threads", type=int, default=10,
        help="ROOT implicit-MT thread count (default: 10).",
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Output directory. Default: results/timing_compare/runs_<run>_...",
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
    if args.fit_min_events < 1:
        parser.error("--fit-min-events must be positive")
    if args.fit_max_half_width <= 0:
        parser.error("--fit-max-half-width must be positive")

    if args.outdir is None:
        tag = "_".join("run{}".format(run) for run in args.runs)
        args.outdir = os.path.join("results", "timing_compare", tag)
    if args.position_label is None:
        args.position_label = args.position_key
    return args


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def display_material(material: str) -> str:
    return MATERIAL_DISPLAY.get(material, material)


def make_manager_args(args: argparse.Namespace, run: int) -> argparse.Namespace:
    return argparse.Namespace(
        run=run,
        first_event=args.first_event,
        last_event=args.last_event,
        json_file=args.json_file,
    )


def active_timing_channels(
    manager: CaloXAnalysisManager,
    requested: Optional[set],
) -> Dict[str, ChannelInfo]:
    """Return all mapped non-reference calorimeter channels with TOA columns."""
    columns = {str(column) for column in manager.rdf.GetColumnNames()}
    selected: Dict[str, ChannelInfo] = {}

    for _, board in manager.drsboards.items():
        for channel in board:
            if channel.is_reference:
                continue

            branch = channel.get_channel_name(blsub=False)
            if requested is not None and branch not in requested:
                continue
            if "{}_Time_cfd_mcp".format(branch) not in columns:
                continue

            try:
                material = str(get_channel_var(channel))
            except Exception:
                material = "DRS"
            size = "6 mm" if getattr(channel, "is6mm", False) else "3 mm"
            board_no = getattr(board, "board_no", "?")
            group_no = getattr(channel, "group_no", "?")
            channel_no = getattr(channel, "channel_no", "?")
            short_label = "B{} G{} Ch{} | {}, {}".format(
                board_no, group_no, channel_no, display_material(material), size
            )
            selected[branch] = ChannelInfo(
                branch=branch,
                label="{} | {}, {}".format(branch, display_material(material), size),
                short_label=short_label,
                material=material,
                size=size,
            )
    return selected


def resolve_mcp_pair(
    run: int,
    args: argparse.Namespace,
    columns: set,
) -> Optional[MCPPair]:
    """Use the native run-mapped MCP pair used by SelectionManager."""
    mapped = list(get_mcp_channels(run).keys())
    live = [name for name in mapped if "{}_TS_cfd_ref".format(name) in columns]

    if args.mcp_a is not None:
        pair = MCPPair(args.mcp_a, args.mcp_b)
        missing = [
            "{}_TS_cfd_ref".format(name)
            for name in (pair.first, pair.second)
            if "{}_TS_cfd_ref".format(name) not in columns
        ]
        if missing:
            print("[WARN] Run {}: requested MCP pair {} unavailable: {}".format(
                run, pair.label, ", ".join(missing)
            ))
            return None
        print("[MCP] Run {}: using requested pair {}".format(run, pair.label))
        return pair

    if len(mapped) >= 5:
        native_pair = MCPPair(mapped[0], mapped[4])
        if (
            "{}_TS_cfd_ref".format(native_pair.first) in columns
            and "{}_TS_cfd_ref".format(native_pair.second) in columns
        ):
            print("[MCP] Run {}: using native mapped pair {}".format(run, native_pair.label))
            return native_pair

    if len(live) >= 2:
        pair = MCPPair(live[0], live[-1])
        print("[WARN] Run {}: native pair unavailable; using {}".format(run, pair.label))
        return pair

    print("[WARN] Run {}: no usable MCP delta-t pair. mapped=[{}], live=[{}]".format(
        run, ", ".join(mapped) or "none", ", ".join(live) or "none"
    ))
    return None


def mcp_reference_is_live(run: int, columns: set) -> bool:
    return MCP_REF in get_mcp_channels(run) and "{}_TS_cfd_ref".format(MCP_REF) in columns


def book_run_histograms(
    manager: CaloXAnalysisManager,
    run: int,
    args: argparse.Namespace,
    channel_info: Dict[str, ChannelInfo],
) -> Tuple[Dict[str, Any], Optional[Any], Optional[MCPPair]]:
    toa_results: Dict[str, Any] = {}
    for branch in sorted(channel_info):
        column = "{}_Time_cfd_mcp".format(branch)
        model = ROOT.RDF.TH1DModel(
            "toa_run{}_{}".format(run, safe_name(branch)),
            "Run {};Time of Arrival [ns];Events".format(run),
            args.toa_bins,
            args.toa_range[0],
            args.toa_range[1],
        )
        toa_results[branch] = manager.rdf.Filter(
            "std::isfinite({})".format(column)
        ).Histo1D(model, column)

    columns = {str(column) for column in manager.rdf.GetColumnNames()}
    pair = resolve_mcp_pair(run, args, columns)
    if pair is None:
        return toa_results, None, None

    dt_column = "__mcp_delta_t_ns"
    dt_rdf = manager.rdf.Define(
        dt_column,
        "({}_TS_cfd_ref - {}_TS_cfd_ref) * {}".format(
            pair.second, pair.first, _TS_TO_NS
        ),
    )
    model = ROOT.RDF.TH1DModel(
        "mcp_delta_t_run{}".format(run),
        "Run {};{} [ns];Events".format(run, pair.label),
        args.mcp_dt_bins,
        args.mcp_dt_range[0],
        args.mcp_dt_range[1],
    )
    result = dt_rdf.Filter("std::isfinite({})".format(dt_column)).Histo1D(model, dt_column)
    return toa_results, result, pair


def clone_hist(result: Any, name: str) -> Any:
    hist = result.GetValue().Clone(name)
    hist.SetDirectory(0)
    return hist


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


def fit_gaussian_hist(hist: Any, stem: str, args: argparse.Namespace) -> TimingFitSummary:
    n_entries = in_range_entries(hist)
    if n_entries < args.fit_min_events:
        return TimingFitSummary(
            "Gaussian", "", n_entries, None, None, None, None,
            None, None, None, False
        )

    axis = hist.GetXaxis()
    peak, sigma_guess = estimate_peak_sigma(hist)
    half_width = min(args.fit_max_half_width, max(0.5, 2.5 * sigma_guess))
    lo = max(axis.GetXmin(), peak - half_width)
    hi = min(axis.GetXmax(), peak + half_width)
    if hi <= lo:
        return TimingFitSummary("Gaussian", "", n_entries, None, None, None, None, None, lo, hi, False)

    fit = ROOT.TF1("fit_gaus_{}".format(safe_name(stem)), "gaus", lo, hi)
    fit.SetParameters(hist.GetMaximum(), peak, max(sigma_guess, axis.GetBinWidth(1)))
    fit.SetParLimits(1, lo, hi)
    fit.SetParLimits(2, 0.25 * axis.GetBinWidth(1), max(axis.GetBinWidth(1), hi - lo))
    fit.SetLineWidth(2)

    status = int(hist.Fit(fit, "Q0RS"))
    sigma = abs(float(fit.GetParameter(2)))
    mean = float(fit.GetParameter(1))
    edge_guard = max(2.0 * axis.GetBinWidth(1), 0.01 * (axis.GetXmax() - axis.GetXmin()))
    valid = (
        status == 0
        and math.isfinite(mean)
        and math.isfinite(sigma)
        and sigma > 0.0
        and (axis.GetXmin() + edge_guard) < mean < (axis.GetXmax() - edge_guard)
        and (lo + 0.5 * axis.GetBinWidth(1)) < mean < (hi - 0.5 * axis.GetBinWidth(1))
    )
    if not valid:
        return TimingFitSummary("Gaussian", "", n_entries, None, None, None, None, None, lo, hi, False)

    ndf = int(fit.GetNDF())
    chi2_ndf = float(fit.GetChisquare() / ndf) if ndf > 0 else None
    return TimingFitSummary(
        model="Gaussian",
        material="",
        n_entries=n_entries,
        mean_ns=mean,
        point_error_ns=sigma / math.sqrt(float(n_entries)),
        sigma_ns=sigma,
        sigma_err_ns=float(fit.GetParError(2)),
        chi2_ndf=chi2_ndf,
        fit_lo_ns=lo,
        fit_hi_ns=hi,
        fit_ok=True,
        function=fit,
    )


def declare_langau() -> None:
    """Declare the Landau⊗Gaussian C++ function once for use as a TF1 pointer."""
    global _LANGAU_DECLARED
    if _LANGAU_DECLARED:
        return

    ROOT.gInterpreter.Declare(
        r'''
        #include "TMath.h"
        double calox_langau(double *x, double *par) {
            // par[0] = Landau width, par[1] = Landau MPV,
            // par[2] = normalization, par[3] = Gaussian sigma.
            const double mpshift = -0.22278298;
            const int np = 80;
            const double sc = 5.0;
            if (par[0] <= 0.0 || par[3] <= 0.0 || par[2] < 0.0) return 0.0;

            const double mpc = par[1] - mpshift * par[0];
            const double xlow = x[0] - sc * par[3];
            const double xupp = x[0] + sc * par[3];
            const double step = (xupp - xlow) / np;
            double sum = 0.0;

            for (int i = 1; i <= np / 2; ++i) {
                double xx = xlow + (i - 0.5) * step;
                sum += TMath::Landau(xx, mpc, par[0], true)
                       * TMath::Gaus(x[0], xx, par[3], true);
                xx = xupp - (i - 0.5) * step;
                sum += TMath::Landau(xx, mpc, par[0], true)
                       * TMath::Gaus(x[0], xx, par[3], true);
            }
            return par[2] * step * sum;
        }
        '''
    )
    if not hasattr(ROOT, "calox_langau"):
        raise RuntimeError("Cling did not expose calox_langau to PyROOT")
    _LANGAU_DECLARED = True

def function_moments(func: Any, lo: float, hi: float, nsteps: int = 1600) -> Tuple[Optional[float], Optional[float]]:
    """Finite-range mean and sigma of the fitted model, evaluated numerically."""
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return None, None

    dx = (hi - lo) / float(nsteps)
    sum_w = 0.0
    sum_x = 0.0
    sum_x2 = 0.0
    for index in range(nsteps):
        x = lo + (index + 0.5) * dx
        y = float(func.Eval(x))
        if not math.isfinite(y) or y <= 0.0:
            continue
        weight = y * dx
        sum_w += weight
        sum_x += x * weight
        sum_x2 += x * x * weight

    if sum_w <= 0.0:
        return None, None
    mean = sum_x / sum_w
    variance = max(0.0, sum_x2 / sum_w - mean * mean)
    return mean, math.sqrt(variance)


def fit_langau_hist(hist: Any, stem: str, args: argparse.Namespace) -> TimingFitSummary:
    """Fit Sci timing with Landau⊗Gaussian and use finite-range moments."""
    n_entries = in_range_entries(hist)
    if n_entries < args.fit_min_events:
        return TimingFitSummary("Landau#otimesGaussian", "", n_entries, None, None, None, None, None, None, None, False)

    declare_langau()
    axis = hist.GetXaxis()
    peak, sigma_guess = estimate_peak_sigma(hist)
    left_width = min(args.fit_max_half_width, max(0.7, 2.5 * sigma_guess))
    right_width = min(1.8 * args.fit_max_half_width, max(1.2, 6.0 * sigma_guess))
    lo = max(axis.GetXmin(), peak - left_width)
    hi = min(axis.GetXmax(), peak + right_width)
    if hi <= lo:
        return TimingFitSummary("Landau#otimesGaussian", "", n_entries, None, None, None, None, None, lo, hi, False)

    fit = ROOT.TF1(
        "fit_langau_{}".format(safe_name(stem)),
        ROOT.calox_langau,
        lo,
        hi,
        4,
    )
    fit.SetParNames("Landau width", "Landau MPV", "Area", "Gaussian sigma")
    bin_width = axis.GetBinWidth(1)
    fit.SetParameters(
        max(0.5 * sigma_guess, 0.5 * bin_width),
        peak,
        max(1.0, float(n_entries)),
        max(0.5 * sigma_guess, 0.5 * bin_width),
    )
    fit.SetParLimits(0, 0.15 * bin_width, max(bin_width, hi - lo))
    fit.SetParLimits(1, lo, hi)
    fit.SetParLimits(2, 0.0, max(100.0, 100.0 * n_entries))
    fit.SetParLimits(3, 0.15 * bin_width, max(bin_width, hi - lo))
    fit.SetNpx(400)
    fit.SetLineWidth(2)

    status = int(hist.Fit(fit, "Q0RS"))
    mean, sigma_eff = function_moments(fit, axis.GetXmin(), axis.GetXmax())
    valid = (
        status == 0
        and mean is not None
        and sigma_eff is not None
        and math.isfinite(mean)
        and math.isfinite(sigma_eff)
        and sigma_eff > 0.0
    )
    if not valid:
        return TimingFitSummary("Landau#otimesGaussian", "", n_entries, None, None, None, None, None, lo, hi, False)

    ndf = int(fit.GetNDF())
    chi2_ndf = float(fit.GetChisquare() / ndf) if ndf > 0 else None
    return TimingFitSummary(
        model="Landau#otimesGaussian",
        material="",
        n_entries=n_entries,
        mean_ns=mean,
        point_error_ns=sigma_eff / math.sqrt(float(n_entries)),
        sigma_ns=sigma_eff,
        sigma_err_ns=None,
        chi2_ndf=chi2_ndf,
        fit_lo_ns=lo,
        fit_hi_ns=hi,
        fit_ok=True,
        function=fit,
    )


def fit_timing_hist(hist: Any, channel: ChannelInfo, run: int, args: argparse.Namespace) -> TimingFitSummary:
    """Apply the requested material-dependent timing fit."""
    if channel.material == "Sci":
        summary = fit_langau_hist(hist, "{}_run{}".format(channel.branch, run), args)
    else:
        summary = fit_gaussian_hist(hist, "{}_run{}".format(channel.branch, run), args)
    summary.material = channel.material
    return summary


def fit_mcp_gaussian(hist: Any, stem: str) -> MCPFitSummary:
    n_entries = in_range_entries(hist)
    if n_entries < 25:
        return MCPFitSummary(n_entries, None, None, None, None, False)

    axis = hist.GetXaxis()
    peak, sigma_guess = estimate_peak_sigma(hist)
    half_width = min(2.5 * sigma_guess, 0.45 * (axis.GetXmax() - axis.GetXmin()))
    lo = max(axis.GetXmin(), peak - half_width)
    hi = min(axis.GetXmax(), peak + half_width)
    if hi <= lo:
        return MCPFitSummary(n_entries, None, None, None, None, False)

    fit = ROOT.TF1("fit_mcp_{}".format(safe_name(stem)), "gaus", lo, hi)
    fit.SetParameters(hist.GetMaximum(), peak, max(sigma_guess, axis.GetBinWidth(1)))
    fit.SetLineWidth(2)
    status = int(hist.Fit(fit, "Q0RS"))
    sigma = abs(float(fit.GetParameter(2)))
    if status != 0 or not math.isfinite(sigma) or sigma <= 0.0:
        return MCPFitSummary(n_entries, None, None, None, None, False)

    return MCPFitSummary(
        n_entries=n_entries,
        mean_ns=float(fit.GetParameter(1)),
        mean_err_ns=float(fit.GetParError(1)),
        sigma_ns=sigma,
        sigma_err_ns=float(fit.GetParError(2)),
        fit_ok=True,
    )


def histogram_for_plot(hist: Any, normalise: bool) -> Any:
    out = hist.Clone("plot_{}".format(hist.GetName()))
    out.SetDirectory(0)
    integral = out.Integral(1, out.GetNbinsX())
    if normalise and integral > 0.0:
        out.Scale(1.0 / integral)
    out.SetTitle("")
    return out


def make_fit_graph(
    hist: Any,
    summary: TimingFitSummary,
    normalise: bool,
    name: str,
) -> Optional[Any]:
    """Sample a ROOT fit as a TGraph so JSROOT can render custom Langau curves."""
    if not summary.fit_ok or summary.function is None:
        return None
    if summary.fit_lo_ns is None or summary.fit_hi_ns is None:
        return None

    lo = summary.fit_lo_ns
    hi = summary.fit_hi_ns
    if hi <= lo:
        return None

    integral = hist.Integral(1, hist.GetNbinsX())
    scale = 1.0 / integral if normalise and integral > 0.0 else 1.0
    graph = ROOT.TGraph(300)
    graph.SetName(name)
    for index in range(300):
        x = lo + (hi - lo) * index / 299.0
        y = max(0.0, float(summary.function.Eval(x))) * scale
        graph.SetPoint(index, x, y)
    return graph


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


def make_overlay_canvas(
    histograms: Iterable[Tuple[int, Any]],
    canvas_name: str,
    title: str,
    x_title: str,
    y_title: str,
    legend_lines: Dict[int, str],
    normalise: bool,
    selection_label: str,
    timing_fits: Optional[Dict[int, TimingFitSummary]] = None,
    mcp_fits: Optional[Dict[int, MCPFitSummary]] = None,
) -> Optional[Tuple[Any, str]]:
    entries = [
        (run, histogram_for_plot(hist, normalise))
        for run, hist in histograms
        if in_range_entries(hist) > 0
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
    nlegend = len(entries)
    legend_y2 = 0.835
    legend_y1 = max(0.42, legend_y2 - 0.047 * nlegend)
    legend = ROOT.TLegend(0.49, legend_y1, 0.94, legend_y2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.019)

    keepalive: List[Any] = [legend]
    first = True
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
        hist.SetMinimum(0.0)
        hist.SetMaximum(max(1.0e-12, 1.30 * ymax))
        hist.Draw("HIST" if first else "HIST SAME")
        first = False
        legend.AddEntry(hist, legend_lines[run], "l")
        keepalive.append(hist)

        if timing_fits is not None and run in timing_fits:
            graph = make_fit_graph(
                hist, timing_fits[run], normalise,
                "curve_{}_{}".format(canvas_name, run),
            )
            if graph is not None:
                graph.SetLineColor(color)
                graph.SetLineStyle(2)
                graph.SetLineWidth(2)
                graph.Draw("L SAME")
                keepalive.append(graph)

        if mcp_fits is not None:
            summary = mcp_fits.get(run)
            if summary is not None and summary.fit_ok and summary.mean_ns is not None and summary.sigma_ns is not None:
                func = ROOT.TF1(
                    "curve_mcp_{}_{}".format(canvas_name, run),
                    "[0]*exp(-0.5*((x-[1])/[2])^2)",
                    hist.GetXaxis().GetXmin(),
                    hist.GetXaxis().GetXmax(),
                )
                func.SetParameters(hist.GetMaximum(), summary.mean_ns, summary.sigma_ns)
                func.SetLineColor(color)
                func.SetLineStyle(2)
                func.SetLineWidth(2)
                func.Draw("SAME")
                keepalive.append(func)

    draw_header(canvas, title, selection_label)
    legend.Draw()
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keepalive
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def parse_length_mm(value: Any) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(value))
    if match is None:
        raise ValueError("Cannot parse a numeric position from {!r}".format(value))
    return float(match.group(0))


def _parse_inline_positions(tokens: Sequence[str], requested_runs: Sequence[int]) -> Dict[int, float]:
    """Parse --positions tokens written as RUN:MM or RUN=MM."""
    parsed: Dict[int, float] = {}
    for token in tokens:
        if ":" in token:
            run_text, position_text = token.split(":", 1)
        elif "=" in token:
            run_text, position_text = token.split("=", 1)
        else:
            raise ValueError(
                "Invalid --positions token {!r}; use RUN:MM, for example 1898:-50".format(token)
            )
        run = int(run_text.strip())
        parsed[run] = parse_length_mm(position_text)
    requested = set(requested_runs)
    extra = sorted(set(parsed) - requested)
    if extra:
        print("[WARN] --positions includes runs not requested: {}".format(extra))
    return {run: parsed[run] for run in requested_runs if run in parsed}


def _print_position_summary(positions: Dict[int, float], args: argparse.Namespace, source_label: str) -> None:
    if len(positions) < 2 or len(set(positions.values())) < 2:
        print(
            "[WARN] Need at least two distinct scan positions for velocity plots. "
            "Resolved positions from {}: {}".format(source_label, positions)
        )
        return
    print("[POSITION] Using {} scan positions from {}:".format(len(positions), source_label))
    for run in args.runs:
        if run in positions:
            print("  run {:>4d} -> {:g} mm".format(run, positions[run]))
        else:
            print("  run {:>4d} -> MISSING".format(run))


def load_run_positions(args: argparse.Namespace) -> Dict[int, float]:
    """Load scan positions in mm; direct --positions overrides every file source."""
    if args.no_velocity_fits:
        print("[POSITION] Velocity fits disabled by --no-velocity-fits.")
        return {}

    if args.positions:
        try:
            positions = _parse_inline_positions(args.positions, args.runs)
        except (ValueError, TypeError) as exc:
            raise SystemExit("[FATAL] {}".format(exc))
        _print_position_summary(positions, args, "--positions")
        return positions

    if args.positions_json:
        path = Path(args.positions_json)
        if not path.exists():
            print("[WARN] Position JSON {} does not exist; velocity fits disabled.".format(path))
            return {}
        with path.open("r", encoding="utf-8") as handle:
            source = json.load(handle)
        if isinstance(source, dict) and isinstance(source.get("positions"), dict):
            source = source["positions"]

        positions: Dict[int, float] = {}
        for run in args.runs:
            value = source.get(str(run), source.get(run))
            if value is None:
                print("[WARN] Position JSON: run {} is missing.".format(run))
                continue
            try:
                positions[run] = parse_length_mm(value)
            except ValueError as exc:
                print("[WARN] Position JSON: run {} skipped: {}".format(run, exc))
        _print_position_summary(positions, args, str(path))
        return positions

    path = Path(args.runlist)
    if not path.exists():
        print("[WARN] Runlist {} does not exist; velocity fits disabled.".format(path))
        return {}
    with path.open("r", encoding="utf-8") as handle:
        source = json.load(handle)

    positions: Dict[int, float] = {}
    for run in args.runs:
        record = source.get(str(run))
        if not isinstance(record, dict):
            print("[WARN] {}: run {} is missing; no velocity point.".format(path, run))
            continue
        if args.position_key not in record:
            keys = ", ".join(sorted(record.keys()))
            print(
                "[WARN] {}: run {} has no {!r}. Available keys: {}".format(
                    path, run, args.position_key, keys
                )
            )
            continue
        try:
            positions[run] = parse_length_mm(record[args.position_key])
        except ValueError as exc:
            print("[WARN] {}: run {} skipped: {}".format(path, run, exc))

    _print_position_summary(positions, args, "{} [{}]".format(path, args.position_key))
    return positions

def make_velocity_canvas(
    branch: str,
    channel: ChannelInfo,
    points: Sequence[VelocityPoint],
    args: argparse.Namespace,
    selection_label: str,
) -> Tuple[Optional[Any], Optional[str], VelocityFitSummary]:
    """Make a weighted TOA-versus-position graph and compute v in cm/ns."""
    ordered = sorted(points, key=lambda point: (point.position_mm, point.run))
    invalid = VelocityFitSummary(0, None, None, None, None, None, None, None, False)
    if len(ordered) < 2:
        return None, None, invalid
    if len({point.position_mm for point in ordered}) < 2:
        return None, None, invalid

    graph = ROOT.TGraphErrors(len(ordered))
    graph.SetName("gr_velocity_{}".format(safe_name(branch)))
    graph.SetTitle("")
    for index, point in enumerate(ordered):
        graph.SetPoint(index, point.position_mm, point.mean_ns)
        graph.SetPointError(index, 0.0, max(point.point_error_ns, 1.0e-6))

    xmin = min(point.position_mm for point in ordered)
    xmax = max(point.position_mm for point in ordered)
    xpad = max(1.0, 0.08 * (xmax - xmin))
    fit = ROOT.TF1(
        "fit_velocity_{}".format(safe_name(branch)),
        "pol1", xmin - xpad, xmax + xpad,
    )
    result = graph.Fit(fit, "Q0S")
    status = int(result)
    slope = float(fit.GetParameter(1))
    slope_err = float(fit.GetParError(1))
    intercept = float(fit.GetParameter(0))
    intercept_err = float(fit.GetParError(0))
    valid = status == 0 and math.isfinite(slope) and abs(slope) > 0.0

    ndf = int(fit.GetNDF())
    chi2_ndf = float(fit.GetChisquare() / ndf) if ndf > 0 else None
    speed = 0.1 / abs(slope) if valid else None
    speed_err = 0.1 * abs(slope_err) / (slope * slope) if valid else None
    summary = VelocityFitSummary(
        n_points=len(ordered),
        intercept_ns=intercept if valid else None,
        intercept_err_ns=intercept_err if valid else None,
        slope_ns_per_mm=slope if valid else None,
        slope_err_ns_per_mm=slope_err if valid else None,
        speed_cm_per_ns=speed,
        speed_err_cm_per_ns=speed_err,
        chi2_ndf=chi2_ndf if valid else None,
        fit_ok=valid,
    )
    if not valid:
        return None, None, summary

    canvas_name = "velocity_{}".format(safe_name(branch))
    canvas = ROOT.TCanvas(canvas_name, "", 1100, 800)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)

    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.15)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)
    graph.Draw("AP")
    graph.GetXaxis().SetTitle("{} [mm]".format(args.position_label))
    graph.GetYaxis().SetTitle("Fitted TOA mean [ns]")
    graph.GetXaxis().SetTitleSize(0.045)
    graph.GetYaxis().SetTitleSize(0.045)
    graph.GetXaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetLabelSize(0.04)
    graph.GetYaxis().SetTitleOffset(1.25)

    fit.SetLineColor(ROOT.kRed + 1)
    fit.SetLineWidth(3)
    fit.SetLineStyle(2)
    fit.Draw("SAME")

    info = ROOT.TPaveText(0.15, 0.64, 0.53, 0.83, "NDC")
    info.SetBorderSize(0)
    info.SetFillStyle(0)
    info.SetTextFont(42)
    info.SetTextSize(0.024)
    info.SetTextAlign(12)
    info.AddText("{} timing centers".format(channel.model_label if hasattr(channel, "model_label") else channel.material))
    info.AddText("dt/dz = ({:.5f} #pm {:.5f}) ns/mm".format(slope, slope_err))
    info.AddText("v = {:.3f} #pm {:.3f} cm/ns".format(speed, speed_err))
    if chi2_ndf is not None:
        info.AddText("#chi^{2}/ndf = {:.2f}".format(chi2_ndf))
    info.Draw()

    legend = ROOT.TLegend(0.15, 0.53, 0.48, 0.61)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.022)
    legend.AddEntry(graph, "Run timing centers (#sigma/#sqrt{N} errors)", "pe")
    legend.AddEntry(fit, "Weighted linear fit", "l")
    legend.Draw()

    draw_header(canvas, "Velocity fit: {}".format(channel.short_label), selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, fit, info, legend]
    return canvas, str(ROOT.TBufferJSON.ToJSON(canvas).Data()), summary


def mcp_summary_json(summary: MCPFitSummary) -> Dict[str, Any]:
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


def timing_summary_json(summary: TimingFitSummary) -> Dict[str, Any]:
    return {
        "model": summary.model,
        "material": summary.material,
        "n_entries": summary.n_entries,
        "fit_ok": summary.fit_ok,
        "mean_ns": summary.mean_ns,
        "point_error_ns": summary.point_error_ns,
        "sigma_ns": summary.sigma_ns,
        "sigma_err_ns": summary.sigma_err_ns,
        "chi2_ndf": summary.chi2_ndf,
        "fit_range_ns": [summary.fit_lo_ns, summary.fit_hi_ns],
    }


def velocity_summary_json(summary: VelocityFitSummary, points: Sequence[VelocityPoint]) -> Dict[str, Any]:
    payload = asdict(summary)
    payload["points"] = [asdict(point) for point in points]
    return payload


def write_jsroot_gallery(
    outdir: Path,
    canvas_keys: List[str],
    canvas_jsons: Dict[str, str],
    args: argparse.Namespace,
    selection_label: str,
    velocity_html_exists: bool,
) -> Path:
    """Write the TOA/MCP HTML gallery only."""
    html_path = outdir / "timing_overlays.html"
    intro_lines = [
        "**Runs:** {}".format(", ".join(str(run) for run in args.runs)),
        "**Selection:** {}{}".format(
            selection_label,
            " + {}".format(args.particle) if args.particle else "",
        ),
        "**TOA axis:** [{:g}, {:g}] ns".format(args.toa_range[0], args.toa_range[1]),
        "**MCP #Delta t axis:** [{:g}, {:g}] ns".format(args.mcp_dt_range[0], args.mcp_dt_range[1]),
        "* CerQuartz and CerPlastic run histograms use Gaussian timing fits.",
        "* Sci run histograms use Landau#otimesGaussian timing fits.",
        "* Raw histograms, fits, and canvases are also stored in timing_overlays.root.",
    ]
    if velocity_html_exists:
        intro_lines.append("* Per-channel timing-velocity fits are in `timing_velocity_fits.html`.")
    else:
        intro_lines.append("* No velocity page was made: at least two distinct run positions and valid timing fits are required.")

    generate_jsroot_html(
        canvas_keys,
        canvas_jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title="CaloX TB2026 timing overlays",
        intro_text="\n".join(intro_lines),
    )
    return html_path


def write_velocity_gallery(
    outdir: Path,
    canvas_keys: List[str],
    canvas_jsons: Dict[str, str],
    args: argparse.Namespace,
    selection_label: str,
    positions: Dict[int, float],
) -> Path:
    """Write a separate, easy-to-find HTML page containing only velocity fits."""
    html_path = outdir / "timing_velocity_fits.html"
    position_text = ", ".join(
        "run {} = {:g} mm".format(run, positions[run])
        for run in args.runs if run in positions
    )
    intro_lines = [
        "**Runs and positions:** {}".format(position_text),
        "**Selection:** {}{}".format(
            selection_label,
            " + {}".format(args.particle) if args.particle else "",
        ),
        "* CerQuartz and CerPlastic centers: Gaussian #mu, with #sigma/#sqrt{N} y-errors.",
        "* Sci centers: mean of the fitted Landau#otimesGaussian, with #sigma_eff/#sqrt{N} y-errors.",
        "* Each channel is fitted as TOA(z) = p0 + p1 z; v = 0.1/|p1| in cm/ns because z is in mm.",
    ]
    generate_jsroot_html(
        canvas_keys,
        canvas_jsons,
        plots_per_row=2,
        output_html=str(html_path),
        title="CaloX TB2026 per-channel timing velocity fits",
        intro_text="\n".join(intro_lines),
    )
    return html_path

def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    positions = load_run_positions(args)

    root_path = outdir / "timing_overlays.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot create {}".format(root_path))

    print("\n" + "=" * 78)
    print("CaloX TB2026 multi-run timing overlay + velocity fit")
    print("Runs              : {}".format(args.runs))
    print("JSON file         : {}".format(args.json_file))
    print("MCP clean         : {}".format(not args.no_mcp_clean))
    print("Particle selection: {}".format(args.particle or "none"))
    print("Position key      : {}".format(args.position_key))
    print("Output            : {}".format(outdir))
    print("=" * 78)

    requested = set(args.channels) if args.channels else None
    toa_hists: Dict[str, Dict[int, Any]] = {}
    channel_info: Dict[str, ChannelInfo] = {}
    timing_fits: Dict[str, Dict[int, TimingFitSummary]] = {}
    mcp_dt_hists: Dict[int, Any] = {}
    mcp_pairs: Dict[int, MCPPair] = {}
    mcp_fits: Dict[int, MCPFitSummary] = {}
    mcp_reference_applied: Dict[int, bool] = {}

    for run in args.runs:
        beam_type, beam_energy = getRunInfo(run)
        print("\n[RUN {}] beam={}, energy={} GeV".format(run, beam_type, beam_energy))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)

        # Define all clean flags first; use mcp_clean only as the event-quality filter.
        # Do NOT call apply_mcp_diff_selection() because it would bias delta-t width.
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        live_columns = {str(column) for column in manager.rdf.GetColumnNames()}
        mcp_reference_applied[run] = mcp_reference_is_live(run, live_columns)
        if not mcp_reference_applied[run]:
            print(
                "[WARN] Run {}: configured MCP reference {} is not live; "
                "*_Time_cfd_mcp has no MCP shift for this run.".format(run, MCP_REF)
            )

        info = active_timing_channels(manager, requested)
        if not info:
            print("[WARN] Run {}: no mapped channels expose *_Time_cfd_mcp; skipped.".format(run))
            continue
        for branch, meta in info.items():
            channel_info.setdefault(branch, meta)

        results, mcp_result, mcp_pair = book_run_histograms(manager, run, args, info)
        actions: List[Any] = list(results.values())
        if mcp_result is not None:
            actions.append(mcp_result)
        ROOT.RDF.RunGraphs(actions)

        run_dir = root_file.mkdir("run{}".format(run))
        toa_dir = run_dir.mkdir("toa")
        toa_dir.cd()
        for branch, result in results.items():
            hist = clone_hist(result, safe_name(branch))
            meta = channel_info[branch]
            hist.SetTitle("Run {}: {};Time of Arrival [ns];Events".format(run, meta.label))
            hist.Write()
            toa_hists.setdefault(branch, {})[run] = hist

            summary = fit_timing_hist(hist, meta, run, args)
            timing_fits.setdefault(branch, {})[run] = summary
            fit_dir_name = "timing_fits"
            if not root_file.GetDirectory(fit_dir_name):
                root_file.mkdir(fit_dir_name)
            fit_dir = root_file.GetDirectory(fit_dir_name)
            fit_run_dir = fit_dir.GetDirectory("run{}".format(run))
            if not fit_run_dir:
                fit_run_dir = fit_dir.mkdir("run{}".format(run))
            fit_run_dir.cd()
            if summary.function is not None:
                summary.function.Write("fit_{}".format(safe_name(branch)))

            if summary.fit_ok:
                print(
                    "[FIT] Run {} {}: {} mean={:.4f} ns, sigma={:.4f} ns, "
                    "error={:.4f} ns, N={}".format(
                        run, meta.short_label, summary.model.replace("#otimes", "x"),
                        summary.mean_ns, summary.sigma_ns, summary.point_error_ns,
                        summary.n_entries,
                    )
                )
            else:
                print("[WARN] Run {} {}: timing fit failed (N={}).".format(
                    run, meta.short_label, summary.n_entries
                ))

        if mcp_result is not None and mcp_pair is not None:
            run_dir.cd()
            hist = clone_hist(mcp_result, "mcp_delta_t")
            hist.SetTitle("Run {}: {};#Delta t [ns];Events".format(run, mcp_pair.label))
            hist.Write()
            mcp_dt_hists[run] = hist
            mcp_pairs[run] = mcp_pair
            mcp_fits[run] = fit_mcp_gaussian(hist, "run{}_mcp_delta_t".format(run))
            summary = mcp_fits[run]
            if summary.fit_ok and summary.sigma_ns is not None:
                print(
                    "[MCP] Run {}: sigma(delta-t)={:.1f} ps, single MCP={:.1f} ps (N={})".format(
                        run, 1000.0 * summary.sigma_ns,
                        1000.0 * summary.single_mcp_sigma_ns,
                        summary.n_entries,
                    )
                )
            else:
                print("[WARN] Run {}: MCP delta-t fit failed (N={}).".format(run, summary.n_entries))

    if not toa_hists:
        root_file.Close()
        raise RuntimeError("No TOA histograms were produced. Check runs, JSON paths, and selections.")

    normalise = not args.raw_yields
    selection_label = "MCP-clean events" if not args.no_mcp_clean else "All events"
    ordered_channels = sorted(toa_hists)

    canvas_dir = root_file.mkdir("canvases")
    canvas_keys: List[str] = []
    canvas_jsons: Dict[str, str] = {}
    velocity_canvas_keys: List[str] = []
    velocity_canvas_jsons: Dict[str, str] = {}

    if mcp_dt_hists:
        overlays = sorted(mcp_dt_hists.items())
        legend_lines: Dict[int, str] = {}
        for run, hist in overlays:
            summary = mcp_fits.get(run)
            pair = mcp_pairs[run]
            if summary is not None and summary.fit_ok and summary.sigma_ns is not None:
                legend_lines[run] = "Run {}: #sigma_{{#Delta t}}={:.1f} ps, #sigma_{{MCP}}={:.1f} ps".format(
                    run, 1000.0 * summary.sigma_ns,
                    1000.0 * summary.single_mcp_sigma_ns,
                )
            else:
                legend_lines[run] = "Run {}: N={}, fit unavailable".format(run, in_range_entries(hist))

        labels = {pair.label for pair in mcp_pairs.values()}
        title = "{} timing difference".format(next(iter(labels))) if len(labels) == 1 else "MCP timing difference"
        xtitle = "{} [ns]".format(next(iter(labels))) if len(labels) == 1 else "MCP #Delta t [ns]"
        built = make_overlay_canvas(
            overlays,
            "mcp_delta_t_overlay",
            title,
            xtitle,
            "Normalized events" if normalise else "Events",
            legend_lines,
            normalise,
            selection_label,
            mcp_fits=mcp_fits,
        )
        if built is not None:
            canvas, canvas_json = built
            canvas_dir.cd()
            canvas.Write("mcp_delta_t_overlay")
            canvas_keys.extend(["mcp_delta_t_overlay", "NEWLINE"])
            canvas_jsons["mcp_delta_t_overlay"] = canvas_json

    for branch in ordered_channels:
        overlays = sorted(toa_hists[branch].items())
        legend_lines = {}
        for run, hist in overlays:
            summary = timing_fits.get(branch, {}).get(run)
            if summary is not None and summary.fit_ok:
                model_tag = "G" if summary.model == "Gaussian" else "L#otimesG"
                legend_lines[run] = "Run {}: {} #mu={:.3f} ns, #sigma={:.0f} ps, N={}".format(
                    run, model_tag, summary.mean_ns, 1000.0 * summary.sigma_ns, summary.n_entries
                )
            else:
                legend_lines[run] = "Run {}: N={}, fit unavailable".format(run, in_range_entries(hist))

        canvas_name = "toa_{}".format(safe_name(branch))
        built = make_overlay_canvas(
            overlays,
            canvas_name,
            channel_info[branch].short_label,
            "Time of Arrival [ns]",
            "Normalized events" if normalise else "Events",
            legend_lines,
            normalise,
            selection_label,
            timing_fits=timing_fits.get(branch, {}),
        )
        if built is not None:
            canvas, canvas_json = built
            canvas_dir.cd()
            canvas.Write(canvas_name)
            canvas_keys.append(canvas_name)
            canvas_jsons[canvas_name] = canvas_json

    # Build velocity figures only after all per-run timing fits are available.
    velocity_summaries: Dict[str, Tuple[VelocityFitSummary, List[VelocityPoint]]] = {}
    velocity_keys: List[str] = []
    if len(positions) >= 2 and len(set(positions.values())) >= 2 and not args.no_velocity_fits:
        print("\n[VELOCITY] Building per-channel TOA-versus-position fits...")
        for branch in ordered_channels:
            points: List[VelocityPoint] = []
            for run, summary in sorted(timing_fits.get(branch, {}).items()):
                if run not in positions or not summary.fit_ok:
                    continue
                if summary.mean_ns is None or summary.point_error_ns is None:
                    continue
                points.append(VelocityPoint(
                    run=run,
                    position_mm=positions[run],
                    mean_ns=summary.mean_ns,
                    point_error_ns=summary.point_error_ns,
                    model=summary.model,
                ))

            canvas, canvas_json, velocity = make_velocity_canvas(
                branch, channel_info[branch], points, args, selection_label
            )
            velocity_summaries[branch] = (velocity, points)
            if canvas is None or canvas_json is None:
                continue
            canvas_name = "velocity_{}".format(safe_name(branch))
            canvas_dir.cd()
            canvas.Write(canvas_name)
            velocity_canvas_keys.append(canvas_name)
            velocity_keys.append(canvas_name)
            velocity_canvas_jsons[canvas_name] = canvas_json
            if velocity.fit_ok:
                print(
                    "[VELOCITY] {}: v={:.3f} +/- {:.3f} cm/ns, slope={:.5f} +/- {:.5f} ns/mm".format(
                        channel_info[branch].short_label,
                        velocity.speed_cm_per_ns,
                        velocity.speed_err_cm_per_ns,
                        velocity.slope_ns_per_mm,
                        velocity.slope_err_ns_per_mm,
                    )
                )
        print("[VELOCITY] Created {} per-channel velocity canvases.".format(len(velocity_keys)))
    else:
        print(
            "\n[WARN] No velocity canvases created: supply at least two distinct run positions "
            "via data/Runlist.json, --positions-json, or --positions."
        )

    root_file.Close()
    velocity_html_path = None
    if velocity_keys:
        velocity_html_path = write_velocity_gallery(
            outdir,
            velocity_canvas_keys,
            velocity_canvas_jsons,
            args,
            selection_label,
            positions,
        )
    html_path = write_jsroot_gallery(
        outdir, canvas_keys, canvas_jsons, args, selection_label, velocity_html_path is not None
    )

    report = {
        "runs": args.runs,
        "json_file": args.json_file,
        "mcp_clean_selection": not args.no_mcp_clean,
        "particle_selection": args.particle,
        "toa_range_ns": list(args.toa_range),
        "mcp_delta_t_range_ns": list(args.mcp_dt_range),
        "position_key": args.position_key,
        "position_label": args.position_label,
        "run_positions_mm": {str(run): positions.get(run) for run in args.runs},
        "velocity_html": str(velocity_html_path) if velocity_html_path is not None else None,
        "timing_point_definition": {
            "CerQuartz": "Gaussian mean; error = fitted Gaussian sigma / sqrt(N)",
            "CerPlastic": "Gaussian mean; error = fitted Gaussian sigma / sqrt(N)",
            "Sci": "finite-range mean of fitted Landau⊗Gaussian; error = fitted sigma_eff / sqrt(N)",
        },
        "mcp_reference_name": MCP_REF,
        "mcp_reference_applied_by_run": {
            str(run): mcp_reference_applied.get(run, False) for run in args.runs
        },
        "mcp_delta_t_pairs": {
            str(run): [pair.first, pair.second] for run, pair in mcp_pairs.items()
        },
        "mcp_delta_t_fits": {
            str(run): mcp_summary_json(summary) for run, summary in mcp_fits.items()
        },
        "channels": {
            branch: {
                "label": channel_info[branch].label,
                "short_label": channel_info[branch].short_label,
                "material": channel_info[branch].material,
                "size": channel_info[branch].size,
            }
            for branch in ordered_channels
        },
        "timing_fits": {
            branch: {
                str(run): timing_summary_json(summary)
                for run, summary in sorted(timing_fits.get(branch, {}).items())
            }
            for branch in ordered_channels
        },
        "velocity_fits": {
            branch: velocity_summary_json(summary, points)
            for branch, (summary, points) in velocity_summaries.items()
        },
    }
    (outdir / "timing_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("\n[DONE]")
    print("  ROOT histograms / fits : {}".format(root_path))
    print("  Timing HTML            : {}".format(html_path))
    if velocity_html_path is not None:
        print("  Velocity HTML          : {}".format(velocity_html_path))
    print("  JSON summary           : {}".format(outdir / "timing_summary.json"))


if __name__ == "__main__":
    main()