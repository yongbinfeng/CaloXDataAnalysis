#!/usr/bin/env python3
"""TB2026 adjacent-pair timing resolution vs. scan position, split by photoelectrons.

Same pair method as scripts/plot_channel_pair_timing_resolution.py (dt between
x-adjacent, same-material 3 mm channels in one DRS group, resolution =
sigma(dt)/sqrt(2)), but dt is filled only for events in a given photoelectron
bin, giving one curve per (material, p.e.) - the "Time Resolution - position
scan" summary plot, in p.e. instead of ADC.

PE definition (--pe-definition), from the single calibration run in --step3-dir,
applied to each channel's adaptive ADC integral {ch}_energy_adaptive (recomputed
with Step 1's parameters) and to every run in the scan:
  edges   Step 3's discrete ADC bin edges (pe_bin_boundaries.json) - exactly the
          binning of scripts/plot_pe_vs_event_quantity.py (default)
  linear  Step 3's linear calibration PE = slope*ADC + intercept
          (pe_calibration_summary.json), rounded to the nearest integer

A pair "has N p.e." per --pair-pe-mode:
  both  both channels fall in PE bin N (default; cleanest, fewest events)
  min   the lower of the two channels' PE bins is N

Individual pairs have too few events per p.e. bin for a stable fit, so for each
run the pairs of a material are pooled: every pair's dt histogram is re-centred on
that pair's own all-p.e. Gaussian mean (removing fixed cable/channel offsets,
rounded to one dt bin), the re-centred histograms are summed per (material, N),
and the pooled histogram is fit once with a Gaussian + flat background (low-p.e. events timed on noise form a
flat pedestal). The plotted resolution is sigma_core/sqrt(2).

Outputs (under --outdir):
  pair_timing_resolution_vs_pe_pooled_hists.json  pooled dt histogram counts (for replotting)
  pair_timing_resolution_vs_pe_summary.json  pooled points + fit parameters + per-pair entry counts
Plots are made from these by scripts/plot_pair_timing_resolution_vs_pe_plots.py.
  pair_timing_resolution_vs_pe.root          pooled dt histograms per (run, material, N)
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List

import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.timing import auto_timer
from variables.drs import _TS_TO_NS, add_cfd_to_baseline_energy

from plot_channel_pair_timing_resolution import (
    cpp_identifier,
    find_channel_pairs,
    PairFitSummary,
    load_run_positions,
    make_manager_args,
    safe_name,
)
from plot_pe_vs_event_quantity import pe_bin_expr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Adjacent-pair timing resolution vs. scan position per photoelectron bin."
    )
    parser.add_argument("--runs", nargs="+", type=int, default=list(range(1897, 1941)))
    parser.add_argument(
        "--step3-dir", default="results/PhotoElectrons/Step3_ADCtoPEConversion",
        help="Directory with pe_bin_boundaries.json (single calibration run).",
    )
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--position-label", default=None)
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--particle", default="electron",
                        help="Beam-PID particle selection, or 'none' (default: electron).")
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--materials", nargs="+", default=["CerQuartz", "CerPlastic", "Sci"],
                        choices=["CerQuartz", "CerPlastic", "Sci"])
    parser.add_argument("--pe-values", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    parser.add_argument("--pe-definition", choices=["edges", "linear"], default="edges")
    parser.add_argument("--pair-pe-mode", choices=["both", "min"], default="both")
    parser.add_argument("--integral", choices=["adaptive", "repo"], default="adaptive",
                        help="Must match the integral used for the Step 3 edges.")
    parser.add_argument("--dt-bins", type=int, default=1200)
    parser.add_argument("--dt-range", nargs=2, type=float, default=(-12.0, 12.0),
                        metavar=("MIN_NS", "MAX_NS"))
    parser.add_argument("--fit-min-events", type=int, default=30)
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument("--outdir", default="results/PhotoElectrons/PairTimingResolutionvsPE")

    # Must match Step 1's compute_cfd_to_baseline_integral parameters.
    parser.add_argument("--anchor-offset", type=int, default=-330)
    parser.add_argument("--search-half", type=int, default=60)
    parser.add_argument("--cfd-fraction", type=float, default=0.20)
    parser.add_argument("--stop-fraction", type=float, default=0.05)
    parser.add_argument("--stop-abs", type=float, default=4.0)
    parser.add_argument("--stop-consecutive", type=int, default=3)
    parser.add_argument("--max-width", type=int, default=250)
    parser.add_argument("--baseline-pre-lo", type=int, default=80)
    parser.add_argument("--baseline-pre-hi", type=int, default=20)
    parser.add_argument("--min-amplitude", type=float, default=0.0)
    parser.add_argument("--smooth-half", type=int, default=3)
    parser.add_argument("--min-width", type=int, default=40)

    args = parser.parse_args()
    if len(set(args.runs)) != len(args.runs):
        parser.error("--runs contains a duplicate run number")
    if str(args.particle).strip().lower() in ("none", ""):
        args.particle = None
    if args.position_label is None:
        args.position_label = args.position_key
    return args


def linear_pe_expr(column: str, slope: float, intercept: float) -> str:
    """Nearest-integer PE from the Step 3 linear calibration (-1 if ADC is NaN)."""
    pe = "({:.10g}) * ({}) + ({:.10g})".format(slope, column, intercept)
    return "(std::isfinite({0}) ? (int)std::floor({1} + 0.5) : -1)".format(column, pe)


def book_run(rdf, pairs, bin_meta, calib, args, run):
    """Book one dt histogram per (pair, N). Returns {(pair_id, N): RResultPtr}."""
    columns = {str(c) for c in rdf.GetColumnNames()}
    bin_cols: Dict[str, str] = {}
    for pair in pairs:
        for branch in (pair.branch_low, pair.branch_high):
            if branch in bin_cols:
                continue
            adc_col = "{}_energy{}".format(branch, "" if args.integral == "repo" else "_adaptive")
            col = "__pebin_{}".format(cpp_identifier(branch))
            if args.pe_definition == "edges":
                expr = pe_bin_expr(adc_col, bin_meta[branch]["bin_edges_adc"])
            else:
                expr = linear_pe_expr(
                    adc_col, calib[branch]["slope_pe_per_adc"], calib[branch]["intercept_pe"])
            rdf = rdf.Define(col, expr)
            bin_cols[branch] = col

    results = {}
    for pair in pairs:
        pid = cpp_identifier(pair.pair_id)
        dt_col = "__dt_{}".format(pid)
        pe_col = "__pe_{}".format(pid)
        low, high = bin_cols[pair.branch_low], bin_cols[pair.branch_high]
        rdf = rdf.Define(dt_col, "({}_TS_cfd_ref - {}_TS_cfd_ref) * {}".format(
            pair.branch_high, pair.branch_low, _TS_TO_NS))
        if args.pair_pe_mode == "both":
            rdf = rdf.Define(pe_col, "({0} == {1} ? {0} : -1)".format(low, high))
        else:
            rdf = rdf.Define(pe_col, "std::min({}, {})".format(low, high))
        results[(pair.pair_id, "all")] = rdf.Filter(
            "std::isfinite({})".format(dt_col)
        ).Histo1D(ROOT.RDF.TH1DModel(
            "dt_run{}_{}_all".format(run, safe_name(pair.pair_id)),
            "Run {} {} all p.e.;dt [ns];Events".format(run, pair.label),
            args.dt_bins, args.dt_range[0], args.dt_range[1]), dt_col)
        for n in args.pe_values:
            model = ROOT.RDF.TH1DModel(
                "dt_run{}_{}_pe{}".format(run, safe_name(pair.pair_id), n),
                "Run {} {} {} p.e.;dt [ns];Events".format(run, pair.label, n),
                args.dt_bins, args.dt_range[0], args.dt_range[1])
            results[(pair.pair_id, n)] = rdf.Filter(
                "{} == {} && std::isfinite({})".format(pe_col, n, dt_col)
            ).Histo1D(model, dt_col)
    return results


def fit_core_gaussian(hist, stem: str, min_events: int) -> PairFitSummary:
    """Gaussian + flat-background fit of the core of a pooled dt histogram.

    Low-p.e. pulses are sometimes timed on noise, which puts a flat pedestal of
    mistimed events under the real dt peak across the whole dt range. A plain
    Gaussian (or a MAD-seeded window) is dragged wide by it, so the model is
    gaus + constant. The seed is the peak position and its half-maximum width
    on a ~0.08 ns rebinning; the histogram is then rebinned so a bin is ~sigma/4
    wide, and a binned-likelihood fit is made in +-6 sigma, refined once. The
    reported sigma is that of the Gaussian core only.
    """
    nbins = hist.GetNbinsX()
    n_entries = int(round(hist.Integral(1, nbins)))
    if n_entries < min_events:
        return PairFitSummary(n_entries, None, None, None, None, False)

    axis = hist.GetXaxis()
    bin_width = axis.GetBinWidth(1)

    def rebin_to(width):
        factor = max(1, int(round(width / bin_width)))
        while nbins % factor:
            factor -= 1
        out = hist.Clone("rebin_{}_{}".format(safe_name(stem), factor))
        out.SetDirectory(0)
        if factor > 1:
            out.Rebin(factor)
        return out, factor

    # Seed: peak of a lightly rebinned histogram, width from its half maximum.
    coarse, _ = rebin_to(0.08)
    cn = coarse.GetNbinsX()
    values = [coarse.GetBinContent(i) for i in range(1, cn + 1)]
    smooth = [sum(values[max(0, i - 1):i + 2]) / len(values[max(0, i - 1):i + 2])
              for i in range(cn)]
    peak_i = max(range(cn), key=lambda i: smooth[i])
    background = sorted(values)[cn // 2]
    half = background + 0.5 * (smooth[peak_i] - background)
    left, right = peak_i, peak_i
    while left > 0 and smooth[left] >= half:
        left -= 1
    while right < cn - 1 and smooth[right] >= half:
        right += 1
    cbw = coarse.GetXaxis().GetBinWidth(1)
    mean = coarse.GetXaxis().GetBinCenter(peak_i + 1)
    sigma = max(cbw, (right - left) * cbw / 2.355)

    mean_seed_range = axis.GetXmax() - axis.GetXmin()
    fit, fit_status, rebinned, factor = None, -1, None, 1
    for _ in range(2):
        rebinned, factor = rebin_to(sigma / 4.0)
        lo = max(axis.GetXmin(), mean - 6.0 * sigma)
        hi = min(axis.GetXmax(), mean + 6.0 * sigma)
        fit = ROOT.TF1("fit_core_{}".format(safe_name(stem)), "gaus(0)+pol0(3)", lo, hi)
        bg_per_bin = background * factor / max(1, round(0.08 / bin_width))
        fit.SetParameters(max(rebinned.GetMaximum() - bg_per_bin, 1.0), mean, sigma, bg_per_bin)
        fit.SetParLimits(1, mean - 3.0 * sigma, mean + 3.0 * sigma)
        fit.SetParLimits(2, bin_width, mean_seed_range)
        fit.SetParLimits(3, 0.0, max(rebinned.GetMaximum(), 1.0))
        fit_status = int(rebinned.Fit(fit, "Q0RSL"))
        if fit_status != 0:
            break
        mean, sigma = float(fit.GetParameter(1)), abs(float(fit.GetParameter(2)))
        background = fit.GetParameter(3) * round(0.08 / bin_width) / factor
    if fit_status != 0 or not math.isfinite(sigma) or sigma <= 0.0:
        return PairFitSummary(n_entries, None, None, None, None, False)
    summary = PairFitSummary(
        n_entries=n_entries, mean_ns=mean, mean_err_ns=float(fit.GetParError(1)),
        sigma_ns=sigma, sigma_err_ns=float(fit.GetParError(2)), fit_ok=True)
    # Heights are per rebinned bin, for redrawing the fit.
    summary.amplitude = float(fit.GetParameter(0))
    summary.background = float(fit.GetParameter(3))
    summary.rebin = factor
    summary.fit_lo, summary.fit_hi = lo, hi
    return summary


def shifted_add(pooled, hist, shift_bins: int) -> None:
    """pooled += hist with content moved from bin i to bin i - shift_bins."""
    nbins = hist.GetNbinsX()
    for i in range(1, nbins + 1):
        j = i - shift_bins
        if 1 <= j <= nbins:
            content = hist.GetBinContent(i)
            if content:
                pooled.SetBinContent(j, pooled.GetBinContent(j) + content)
                pooled.SetBinError(j, math.hypot(pooled.GetBinError(j), hist.GetBinError(i)))


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    positions = load_run_positions(args)

    with (Path(args.step3_dir) / "pe_bin_boundaries.json").open() as handle:
        bins_data = json.load(handle)
    bin_meta = bins_data["channels"]
    with (Path(args.step3_dir) / "pe_calibration_summary.json").open() as handle:
        calib = json.load(handle)["channels"]
    print("PE calibration run: {}; runs: {}".format(bins_data["run"], args.runs))

    root_file = ROOT.TFile(str(outdir / "pair_timing_resolution_vs_pe.root"), "RECREATE")
    records: List[Dict[str, Any]] = []
    hist_counts: Dict[str, List[int]] = {}
    combined: List[Dict[str, Any]] = []
    beam_desc = None

    for run in args.runs:
        beam_type, beam_energy = getRunInfo(run)
        beam_desc = beam_desc or "{} {} GeV".format(beam_type, beam_energy)
        position_mm = positions.get(run)
        print("\n[RUN {}] beam={}, energy={} GeV, position={} mm".format(
            run, beam_type, beam_energy, position_mm))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        columns = {str(c) for c in manager.rdf.GetColumnNames()}
        pairs = [p for p in find_channel_pairs(manager, columns, args.materials)
                 if all(b in bin_meta and b in calib and not calib[b].get("skipped")
                        for b in (p.branch_low, p.branch_high))]
        if not pairs:
            print("[WARN] Run {}: no qualifying channel pairs; skipped.".format(run))
            continue

        branches = sorted({b for p in pairs for b in (p.branch_low, p.branch_high)})
        rdf = manager.rdf
        if args.integral == "adaptive":
            rdf = add_cfd_to_baseline_energy(
                rdf, branches,
                anchor_offset=args.anchor_offset, search_half=args.search_half,
                cfd_fraction=args.cfd_fraction, stop_fraction=args.stop_fraction,
                stop_abs=args.stop_abs, stop_consecutive=args.stop_consecutive,
                max_width=args.max_width, baseline_pre_lo=args.baseline_pre_lo,
                baseline_pre_hi=args.baseline_pre_hi, min_amplitude=args.min_amplitude,
                smooth_half=args.smooth_half, min_width=args.min_width)

        results = book_run(rdf, pairs, bin_meta, calib, args, run)
        print("[INFO] Run {}: {} pairs, running {} histograms.".format(
            run, len(pairs), len(results)))
        ROOT.RDF.RunGraphs(list(results.values()))

        run_dir = root_file.mkdir("run{}".format(run))
        pooled: Dict[Any, Any] = {}
        n_pairs: Dict[str, int] = {}
        for pair in pairs:
            all_hist = results[(pair.pair_id, "all")].GetValue()
            offset = fit_core_gaussian(
                all_hist, "{}_{}_all".format(run, pair.pair_id), args.fit_min_events)
            if not offset.fit_ok:
                print("[WARN] Run {} {}: no offset fit (N={}); pair skipped.".format(
                    run, pair.label, offset.n_entries))
                continue
            shift = int(round(offset.mean_ns / all_hist.GetXaxis().GetBinWidth(1)))
            n_pairs[pair.material] = n_pairs.get(pair.material, 0) + 1
            for n in args.pe_values:
                hist = results[(pair.pair_id, n)].GetValue()
                key = (pair.material, n)
                if key not in pooled:
                    pooled[key] = hist.Clone("pooled_run{}_{}_pe{}".format(run, pair.material, n))
                    pooled[key].Reset()
                    pooled[key].SetDirectory(0)
                    pooled[key].SetTitle("Run {} {} pooled {} p.e.;dt [ns];Events".format(
                        run, pair.material, n))
                shifted_add(pooled[key], hist, shift)
                records.append({
                    "run": run, "pair_id": pair.pair_id, "material": pair.material,
                    "pe": n, "n_entries": int(hist.GetEntries()),
                    "offset_ns": offset.mean_ns,
                })

        for (material, n), hist in sorted(pooled.items()):
            summary = fit_core_gaussian(
                hist, "{}_{}_pe{}".format(run, material, n), args.fit_min_events)
            run_dir.cd()
            hist.Write()
            hist_counts["{}|{}|{}".format(run, material, n)] = [
                int(round(hist.GetBinContent(i))) for i in range(1, hist.GetNbinsX() + 1)]
            if summary.fit_ok:
                combined.append({
                    "material": material, "pe": n, "run": run,
                    "position_mm": position_mm, "n_pairs": n_pairs[material],
                    "n_entries": summary.n_entries, "mean_ns": summary.mean_ns,
                    "sigma_ns": summary.sigma_ns, "amplitude": summary.amplitude,
                    "background": summary.background,
                    "rebin": summary.rebin, "fit_lo": summary.fit_lo, "fit_hi": summary.fit_hi,
                    "resolution_ps": summary.resolution_ps,
                    "resolution_err_ps": summary.resolution_err_ps,
                })
                print("[FIT] {} {} p.e.: N={} res={:.1f} +- {:.1f} ps".format(
                    material, n, summary.n_entries, summary.resolution_ps,
                    summary.resolution_err_ps))
            else:
                print("[WARN] {} {} p.e.: fit failed/too few events (N={}).".format(
                    material, n, summary.n_entries))

    root_file.Close()
    if not records:
        raise RuntimeError("No pair histograms were produced for any requested run.")

    with (outdir / "pair_timing_resolution_vs_pe_summary.json").open("w") as handle:
        json.dump({
            "runs": args.runs, "pe_calibration_run": bins_data["run"],
            "materials": args.materials, "pe_values": args.pe_values,
            "pair_pe_mode": args.pair_pe_mode, "pe_definition": args.pe_definition, "integral": args.integral,
            "position_key": args.position_key, "position_label": args.position_label,
            "particle": args.particle, "beam": beam_desc,
            "combined": combined, "records": records,
        }, handle, indent=2)

    with (outdir / "pair_timing_resolution_vs_pe_pooled_hists.json").open("w") as handle:
        json.dump({"dt_range": list(args.dt_range), "dt_bins": args.dt_bins,
                   "counts": hist_counts}, handle)

    print("\n[WROTE] {}".format(outdir))


if __name__ == "__main__":
    main()
