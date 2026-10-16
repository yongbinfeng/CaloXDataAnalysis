#!/usr/bin/env python3
"""Compare the PhotoElectrons study done with two integral definitions.

  adaptive  results/PhotoElectrons              ({ch}_energy_adaptive)
  repo      results/PhotoElectrons_RepoIntegral ({ch}_energy)

Two pages are written to --outdir:

  position_scan_overlay_comparison.html
      Per channel, the position-scan ADC-integral overlay of both studies side
      by side (adaptive left, repo right). Both panels are drawn by
      plot_adc_position_scan_overlay.make_overlay_canvas, so styling is identical.

  pe_attenuation_scans_comparison.html
      Per channel, the Step 6.5 mean-p.e. vs position graphs of both studies
      overlaid on one canvas (adaptive blue, repo red).

No integrals or fits are recomputed; everything is read from stored ROOT files.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ROOT

from utils.html_generator import generate_jsroot_html

from plot_adc_position_scan_overlay import ROOT_COLORS, draw_header, make_overlay_canvas
from plot_pe_fit_position_overlay import make_overlay_canvas as make_fit_overlay_canvas

ADAPTIVE_ROOT = Path("results/PhotoElectrons")
REPO_ROOT = Path("results/PhotoElectrons_RepoIntegral")

ADAPTIVE_COLOR = ROOT.kBlue + 1
REPO_COLOR = ROOT.kRed + 1


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--adaptive-dir", type=Path, default=ADAPTIVE_ROOT)
    parser.add_argument("--repo-dir", type=Path, default=REPO_ROOT)
    parser.add_argument(
        "--outdir", type=Path, default=Path("results/PhotoElectrons_Comparison")
    )
    parser.add_argument("--zoom-keep", type=float, default=0.995)
    parser.add_argument("--plots-per-row", type=int, default=2)
    parser.add_argument(
        "--no-fits", action="store_true",
        help="Position-scan page with histograms only (no Step 2 fits).",
    )
    parser.add_argument(
        "--skip-overlay", action="store_true", help="Skip the position-scan overlay page."
    )
    parser.add_argument(
        "--skip-sigma", action="store_true", help="Skip the fitted-sigma histogram page."
    )
    parser.add_argument(
        "--skip-attenuation", action="store_true", help="Skip the attenuation page."
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")


def to_json(canvas: Any) -> str:
    return str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def open_root(path: Path) -> Any:
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie():
        raise SystemExit("Cannot open ROOT file: {}".format(path))
    return handle


# ---------------------------------------------------------------------------
# Position-scan overlay, side by side
# ---------------------------------------------------------------------------

def load_resolved_runs(summary_path: Path) -> List[Tuple[float, int, float]]:
    with summary_path.open() as handle:
        summary = json.load(handle)
    runs = summary["resolved_runs"]
    if isinstance(runs, dict):  # older layout: {target: {run, actual_position_mm}}
        return [
            (float(target), int(info["run"]), float(info["actual_position_mm"]))
            for target, info in runs.items()
        ]
    return [
        (float(r["target_position_mm"]), int(r["run"]), float(r["actual_position_mm"]))
        for r in runs
    ]


def adaptive_overlay_hists(root_file: Any, runs: List[int]) -> Dict[int, Dict[str, Any]]:
    """Old overlay file layout: run<R>/hist_adc_scan_run<R>_<branch>."""
    out: Dict[int, Dict[str, Any]] = {}
    for run in runs:
        directory = root_file.Get("run{}".format(run))
        hists: Dict[str, Any] = {}
        if directory:
            prefix = "hist_adc_scan_run{}_".format(run)
            for key in directory.GetListOfKeys():
                name = key.GetName()
                if name.startswith(prefix):
                    hist = directory.Get(name)
                    hist.SetDirectory(0)
                    hists[name[len(prefix):]] = hist
        out[run] = hists
    return out


def repo_overlay_hists(
    root_file: Any, summary: Dict[str, Any], runs: List[int]
) -> Dict[int, Dict[str, Any]]:
    """Step 1 layout: run<R>/<hist_name from pe_integrals_summary.json>."""
    out: Dict[int, Dict[str, Any]] = {}
    for run in runs:
        hists: Dict[str, Any] = {}
        for branch, meta in summary["channels"].get("run{}".format(run), {}).items():
            hist = root_file.Get("run{}/{}".format(run, meta["hist_name"]))
            if hist:
                hist.SetDirectory(0)
                hists[branch] = hist
        out[run] = hists
    return out


def build_position_overlay(args: argparse.Namespace) -> None:
    ad_dir = args.adaptive_dir / "PositionScanOverlay"
    rp_dir = args.repo_dir / "PositionScanOverlay"
    ad_runs = load_resolved_runs(ad_dir / "adc_integral_position_overlay_summary.json")
    rp_runs = load_resolved_runs(rp_dir / "pe_integrals_position_overlay_summary.json")

    ad_file = open_root(ad_dir / "adc_integral_position_overlay.root")
    ad_hists = adaptive_overlay_hists(ad_file, [r for _, r, _ in ad_runs])
    ad_file.Close()

    rp_step1 = rp_dir / "Step1_IntegralCalculation"
    with (rp_step1 / "pe_integrals_summary.json").open() as handle:
        rp_summary = json.load(handle)
    rp_file = open_root(rp_step1 / "pe_integrals.root")
    rp_hists = repo_overlay_hists(rp_file, rp_summary, [r for _, r, _ in rp_runs])
    rp_file.Close()

    # Drop runs with no histograms (e.g. run 1777 has no live channels).
    ad_runs = [e for e in ad_runs if ad_hists.get(e[1])]
    rp_runs = [e for e in rp_runs if rp_hists.get(e[1])]
    print("[adaptive] runs: {}".format([e[1] for e in ad_runs]))
    print("[repo]     runs: {}".format([e[1] for e in rp_runs]))

    # One colour per run, assigned in order of target position, so a run has
    # the same colour in the adaptive and the repo panel.
    by_target = sorted({(t, r) for t, r, _ in ad_runs + rp_runs})
    run_colors = {
        run: ROOT_COLORS[i % len(ROOT_COLORS)] for i, (_, run) in enumerate(by_target)
    }
    ad_runs = sorted(ad_runs)
    rp_runs = sorted(rp_runs)

    def common(hists: Dict[int, Dict[str, Any]], entries) -> set:
        sets = [set(hists[e[1]]) for e in entries]
        return set.intersection(*sets) if sets else set()

    branches = sorted(common(ad_hists, ad_runs) & common(rp_hists, rp_runs))
    labels = {
        branch: rp_summary["channels"]["run{}".format(rp_runs[0][1])][branch].get(
            "short_label", branch)
        for branch in branches
    }

    keys: List[str] = []
    jsons: Dict[str, str] = {}
    skipped: List[str] = []
    for index, branch in enumerate(branches, 1):
        pair = []
        for tag, entries, hists, text in (
            ("adaptive", ad_runs, ad_hists, "adaptive integral"),
            ("repo", rp_runs, rp_hists, "repo integral"),
        ):
            built = make_overlay_canvas(
                branch, "{}  [{}]".format(labels[branch], text),
                "MCP-clean events", entries,
                {run: hists[run][branch] for _, run, _ in entries},
                True, args.zoom_keep, run_colors,
            )
            if built is None:
                pair = []
                break
            canvas, _ = built
            canvas.SetName("pe_overlay_{}_{}".format(tag, safe_name(branch)))
            pair.append((canvas, to_json(canvas)))
        if not pair:
            skipped.append(branch)
            continue
        for canvas, js in pair:
            keys.append(canvas.GetName())
            jsons[canvas.GetName()] = js
        if index % 25 == 0 or index == len(branches):
            print("[OVERLAY] {}/{} channels".format(index, len(branches)))

    def run_text(entries) -> str:
        return ", ".join(
            "{:g} mm -> run {} ({:g} mm)".format(t, r, a) for t, r, a in entries
        )

    intro = "\n".join([
        "**Left:** adaptive integral `{ch}_energy_adaptive` (`results/PhotoElectrons`). "
        "**Right:** repo integral `{ch}_energy` (`results/PhotoElectrons_RepoIntegral`).",
        "**Adaptive runs:** " + run_text(ad_runs) + ".",
        "**Repo runs:** " + run_text(rp_runs) + ".",
        "Spectra are unit-area normalized, drawn as unfilled outlines. The ADC "
        "axes are not on a common scale: the two integrals have different "
        "definitions.",
    ])
    out_html = args.outdir / "position_scan_overlay_comparison.html"
    generate_jsroot_html(
        keys, jsons, plots_per_row=args.plots_per_row, output_html=str(out_html),
        title="CaloX PE study: position-scan overlay, adaptive vs repo integral",
        intro_text=intro,
    )
    print("[DONE] {} ({} channels, {} skipped)".format(out_html, len(branches) - len(skipped), len(skipped)))


def build_position_fit_overlay(args: argparse.Namespace) -> bool:
    """Position-scan overlay with the Step 2 total fit drawn on the data.

    Both studies are read the same way: <study>/PositionScanOverlay holds a
    Step1_IntegralCalculation (histograms) and Step2_PEFits (fits, one run<R>
    directory per run, written by fit_pe_peaks.py --runs). Returns False when
    either study has no fits yet.
    """
    studies = []
    for tag, study_dir, text in (
        ("adaptive", args.adaptive_dir, "adaptive integral"),
        ("repo", args.repo_dir, "repo integral"),
    ):
        base = study_dir / "PositionScanOverlay"
        step1, step2 = base / "Step1_IntegralCalculation", base / "Step2_PEFits"
        if not (step1 / "pe_integrals_summary.json").is_file() or not any(step2.glob("run*/pe_fit_curves.json")):
            print("[FITS] {}: no Step 1 / Step 2 outputs under {}".format(tag, base))
            return False
        with (step1 / "pe_integrals_summary.json").open() as handle:
            summary = json.load(handle)
        studies.append((tag, text, step1, step2, summary))

    resolved = load_resolved_runs(
        args.repo_dir / "PositionScanOverlay" / "pe_integrals_position_overlay_summary.json"
    )

    loaded = []
    for tag, text, step1, step2, summary in studies:
        entries, fits, curves = [], {}, {}
        for target, run, actual in resolved:
            run_dir = step2 / "run{}".format(run)
            if not (run_dir / "pe_fit_curves.json").is_file():
                continue
            with (run_dir / "pe_fits_summary.json").open() as handle:
                fits[run] = json.load(handle)["channels"]
            with (run_dir / "pe_fit_curves.json").open() as handle:
                curves[run] = json.load(handle)
            entries.append((target, run, actual))
        root_file = open_root(step1 / "pe_integrals.root")
        hists: Dict[int, Dict[str, Any]] = {}
        for _, run, _ in entries:
            hists[run] = {}
            for branch, meta in summary["channels"].get("run{}".format(run), {}).items():
                hist = root_file.Get("run{}/{}".format(run, meta["hist_name"]))
                if hist:
                    hist.SetDirectory(0)
                    hists[run][branch] = hist
        root_file.Close()
        print("[{}] runs with fits: {}".format(tag, [e[1] for e in entries]))
        loaded.append((tag, text, sorted(entries), hists, fits, curves, summary))

    # One colour per run (by target position), identical in both panels.
    by_target = sorted({(t, r) for _, _, entries, *_ in loaded for t, r, _ in entries})
    run_colors = {
        run: ROOT_COLORS[i % len(ROOT_COLORS)] for i, (_, run) in enumerate(by_target)
    }

    def branches_of(item) -> set:
        _, _, entries, hists, fits, curves, _ = item
        sets = []
        for _, run, _ in entries:
            good = {b for b, f in fits[run].items()
                    if not f.get("skipped") and b in curves[run]}
            sets.append(good & set(hists[run]))
        return set.intersection(*sets) if sets else set()

    fitted = [branches_of(item) for item in loaded]
    all_branches = sorted(set(loaded[0][3][loaded[0][2][0][1]]) & set(loaded[1][3][loaded[1][2][0][1]]))
    branches = all_branches  # fall back to histogram-only for channels without fits
    labels = {
        branch: loaded[1][6]["channels"]["run{}".format(loaded[1][2][0][1])][branch].get(
            "short_label", branch)
        for branch in branches
    }

    keys: List[str] = []
    jsons: Dict[str, str] = {}
    n_nofit = {"adaptive": 0, "repo": 0}
    for index, branch in enumerate(branches, 1):
        for item, good in zip(loaded, fitted):
            tag, text, entries, hists, fits, curves, _ = item
            source = {run: hists[run][branch] for _, run, _ in entries if branch in hists[run]}
            label = "{}  [{}]".format(labels[branch], text)
            built = None
            if branch in good:
                built = make_fit_overlay_canvas(
                    branch, label, "MCP-clean events", entries, source,
                    {run: fits[run][branch] for _, run, _ in entries},
                    {run: curves[run][branch] for _, run, _ in entries},
                    True, False, args.zoom_keep, run_colors,
                )
            if built is None:
                n_nofit[tag] += 1
                built = make_overlay_canvas(
                    branch, label + " (no fit)", "MCP-clean events", entries, source,
                    True, args.zoom_keep, run_colors,
                )
            if built is None:
                continue
            canvas, _ = built
            canvas.SetName("pe_fit_overlay_{}_{}".format(tag, safe_name(branch)))
            keys.append(canvas.GetName())
            jsons[canvas.GetName()] = to_json(canvas)
        if index % 25 == 0 or index == len(branches):
            print("[OVERLAY] {}/{} channels".format(index, len(branches)))

    def run_text(entries) -> str:
        return ", ".join("{:g} mm -> run {} ({:g} mm)".format(t, r, a) for t, r, a in entries)

    intro = "\n".join([
        "**Left:** adaptive integral `{ch}_energy_adaptive`. **Right:** repo "
        "integral `{ch}_energy` (`results/PhotoElectrons_RepoIntegral`). Each run "
        "has the same colour in both panels.",
        "**Adaptive runs:** " + run_text(loaded[0][2]) + ".",
        "**Repo runs:** " + run_text(loaded[1][2]) + ".",
        "Thin dotted lines are the Step 1 ADC-integral spectra; thick solid lines "
        "are the Step 2 total multi-peak RooCrystalBall fits (sum of the stored "
        "per-p.e. components). The legend gives the fitted gain in ADC/p.e. per run. "
        "Spectra are unit-area normalized; the ADC axes are not on a common scale.",
        "Channels drawn without a fit (fit skipped or failed for some run): "
        "adaptive {}, repo {}.".format(n_nofit["adaptive"], n_nofit["repo"]),
    ])
    out_html = args.outdir / "position_scan_overlay_comparison.html"
    generate_jsroot_html(
        keys, jsons, plots_per_row=args.plots_per_row, output_html=str(out_html),
        title="CaloX PE study: position-scan overlay with fits, adaptive vs repo integral",
        intro_text=intro,
    )
    print("[DONE] {}".format(out_html))
    return True


def load_overlay_fits(study_dir: Path, resolved) -> Dict[int, Dict[str, Any]]:
    """Step 2 per-channel fit summaries for the overlay runs of one study."""
    step2 = study_dir / "PositionScanOverlay" / "Step2_PEFits"
    out: Dict[int, Dict[str, Any]] = {}
    for _, run, _ in resolved:
        path = step2 / "run{}".format(run) / "pe_fits_summary.json"
        if path.is_file():
            with path.open() as handle:
                out[run] = json.load(handle)["channels"]
    return out


def percentile_range(values: List[float], lo: float = 0.01, hi: float = 0.99) -> Tuple[float, float]:
    ordered = sorted(values)
    a = ordered[int(lo * (len(ordered) - 1))]
    b = ordered[int(hi * (len(ordered) - 1))]
    pad = 0.08 * (b - a if b > a else abs(b) or 1.0)
    return a - pad, b + pad


def build_sigma_hists(args: argparse.Namespace) -> bool:
    """Distributions across channels of the fitted peak width at each scan position.

    sigma = (sigmaL + sigmaR) / 2 of the Step 2 Crystal Ball fit, either in ADC
    or divided by the fitted gain dpe (ADC/p.e.) so it reads as a width in p.e.
    One canvas per (channel family, variable, study); the studies are paired
    side by side and each run keeps the colour it has on the overlay page.
    """
    resolved = load_resolved_runs(
        args.repo_dir / "PositionScanOverlay" / "pe_integrals_position_overlay_summary.json"
    )
    studies = [
        ("adaptive", "adaptive integral", load_overlay_fits(args.adaptive_dir, resolved)),
        ("repo", "repo integral", load_overlay_fits(args.repo_dir, resolved)),
    ]
    for tag, _, fits in studies:
        if not fits:
            print("[SIGMA] {}: no Step 2 fits yet under PositionScanOverlay".format(tag))
            return False

    runs = [(t, r, a) for t, r, a in resolved if all(r in fits for _, _, fits in studies)]
    by_target = sorted((t, r) for t, r, _ in runs)
    colors = {run: ROOT_COLORS[i % len(ROOT_COLORS)] for i, (_, run) in enumerate(by_target)}
    runs = sorted(runs)

    # (family, run) -> list of (sigma_adc, sigma_pe) for fitted channels.
    def collect(fits: Dict[int, Dict[str, Any]]):
        table: Dict[str, Dict[int, List[Tuple[float, float]]]] = {}
        for _, run, _ in runs:
            for branch, fit in fits[run].items():
                if fit.get("skipped") or "sigmaL" not in fit or not fit.get("dpe"):
                    continue
                sigma = 0.5 * (float(fit["sigmaL"]) + float(fit["sigmaR"]))
                gain = float(fit["dpe"])
                family = "{} {}".format(fit.get("material", "?"), fit.get("size", "?"))
                table.setdefault(family, {}).setdefault(run, []).append((sigma, sigma / gain))
        return table

    tables = {tag: collect(fits) for tag, _, fits in studies}
    families = sorted(set().union(*[set(t) for t in tables.values()]))
    variables = (
        (0, "sigma_adc", "Fitted peak width #sigma [ADC]"),
        (1, "sigma_pe", "Fitted peak width #sigma / gain [p.e.]"),
    )

    keys: List[str] = []
    jsons: Dict[str, str] = {}
    for family in families:
        for var_index, var_name, axis_title in variables:
            # Shared range and binning across the two studies for this variable.
            pooled = [
                v[var_index] for tag in tables for runvals in tables[tag].get(family, {}).values()
                for v in runvals
            ]
            if len(pooled) < 5:
                continue
            for tag, text, _ in studies:
                xmin, xmax = (
                    percentile_range([v[var_index] for runvals in tables[tag].get(family, {}).values() for v in runvals])
                    if tables[tag].get(family) else percentile_range(pooled)
                )
                hists = []
                for _, run, actual in runs:
                    vals = [v[var_index] for v in tables[tag].get(family, {}).get(run, [])]
                    hist = ROOT.TH1D(
                        "h_{}_{}_{}_{}".format(tag, var_name, safe_name(family), run), "", 30, xmin, xmax)
                    hist.SetDirectory(0)
                    for value in vals:
                        hist.Fill(value)
                    hists.append((run, actual, hist, len(vals)))
                if not any(n for *_, n in hists):
                    continue
                name = "sigma_{}_{}_{}".format(tag, var_name, safe_name(family))
                canvas = ROOT.TCanvas(name, "", 1050, 780)
                canvas.SetLeftMargin(0.13)
                canvas.SetRightMargin(0.04)
                canvas.SetTopMargin(0.12)
                canvas.SetBottomMargin(0.12)
                canvas.SetTicks(1, 1)
                ymax = max(h.GetMaximum() for _, _, h, _ in hists)
                frame = canvas.DrawFrame(xmin, 0.0, xmax, 1.9 * max(ymax, 1.0))
                frame.GetXaxis().SetTitle(axis_title)
                frame.GetYaxis().SetTitle("Channels")
                for axis in (frame.GetXaxis(), frame.GetYaxis()):
                    axis.SetTitleSize(0.045)
                    axis.SetLabelSize(0.04)
                frame.GetYaxis().SetTitleOffset(1.35)
                legend = ROOT.TLegend(0.47, 0.64, 0.95, 0.885)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextFont(42)
                legend.SetTextSize(0.022)
                for run, actual, hist, n in hists:
                    hist.SetLineColor(colors[run])
                    hist.SetLineWidth(2)
                    hist.SetStats(0)
                    hist.Draw("HIST SAME")
                    legend.AddEntry(
                        hist,
                        "Run {} (X={:.0f} mm): N={}, mean {:.3g} #pm {:.2g}".format(
                            run, actual, n, hist.GetMean(), hist.GetRMS()),
                        "l",
                    )
                legend.Draw()
                header = draw_header(
                    canvas, "{}  [{}]".format(family, text), "MCP-clean events")
                canvas.RedrawAxis()
                canvas.Update()
                canvas._keepalive = [frame, legend, header, hists]
                keys.append(name)
                jsons[name] = to_json(canvas)

    # Order so each row is adaptive | repo for one (family, variable).
    ordered = []
    for family in families:
        for _, var_name, _ in variables:
            for tag in ("adaptive", "repo"):
                key = "sigma_{}_{}_{}".format(tag, var_name, safe_name(family))
                if key in jsons:
                    ordered.append(key)

    run_text = ", ".join("{:g} mm -> run {} ({:g} mm)".format(t, r, a) for t, r, a in runs)
    intro = "\n".join([
        "**Left:** adaptive integral. **Right:** repo integral. One row per channel "
        "family and variable; each run keeps its overlay-page colour.",
        "**Runs:** " + run_text + ".",
        "**Quantity:** sigma = (sigmaL + sigmaR)/2 of the Step 2 Crystal Ball fit of "
        "each channel's ADC-integral spectrum, one entry per successfully fitted "
        "channel. The sigma/gain rows divide by the fitted ADC/p.e. gain so the "
        "width is in photoelectron units. N is the number of fitted channels.",
        "The legend gives the mean and RMS of each distribution. Only successfully "
        "fitted channels enter, so the two studies can have different N.",
    ])
    out_html = args.outdir / "sigma_hists_comparison.html"
    generate_jsroot_html(
        ordered, jsons, plots_per_row=2, output_html=str(out_html),
        title="CaloX PE study: fitted peak width by position, adaptive vs repo integral",
        intro_text=intro,
    )
    print("[DONE] {} ({} panels)".format(out_html, len(ordered)))
    return True


# ---------------------------------------------------------------------------
# Step 6.5 attenuation scans, overlaid
# ---------------------------------------------------------------------------

def find_graph(canvas: Any) -> Optional[Any]:
    for obj in canvas.GetListOfPrimitives():
        if obj.InheritsFrom("TGraphErrors"):
            return obj
    return None


def find_title(canvas: Any) -> str:
    for obj in canvas.GetListOfPrimitives():
        if obj.InheritsFrom("TLatex") and "PE attenuation scan" in obj.GetTitle():
            return str(obj.GetTitle()).replace("PE attenuation scan: ", "")
    return ""


def load_fit_summary(study_dir: Path) -> Dict[str, Any]:
    path = (study_dir / "Step6p5AttenuationCorrections" /
            "exponential_fit_per_channel" / "attenuation_exponential_fit_summary.json")
    with path.open() as handle:
        return json.load(handle)


def fit_curve(branch: str, fits: Dict[str, Any], x_min: float, color: int, tag: str) -> Optional[Any]:
    """TF1 for PE(x) = A * exp(k x) over the fitted range, or None if unfitted."""
    fit = fits["channels"].get(branch)
    if fit is None:
        return None
    curve = ROOT.TF1(
        "fit_{}_{}".format(tag, safe_name(branch)), "[0]*exp([1]*x)",
        x_min, float(fits["fit_x_max_mm"]),
    )
    curve.SetParameters(fit["amplitude_at_zero_pe"], fit["slope_per_mm"])
    curve.SetLineColor(color)
    curve.SetLineWidth(3)
    curve.SetLineStyle(1)
    curve.SetNpx(400)
    return curve


def fit_label(fit: Optional[Dict[str, Any]], name: str) -> str:
    if fit is None:
        return "{}: fit skipped".format(name)
    return "{}: #lambda = {:.0f} #pm {:.0f} mm, #chi^{{2}}/ndf = {:.1f}/{}".format(
        name, fit["attenuation_length_mm"], fit["attenuation_length_error_mm"],
        fit["chi2"], fit["ndf"],
    )


def build_attenuation_overlay(args: argparse.Namespace) -> None:
    ad_fits = load_fit_summary(args.adaptive_dir)
    rp_fits = load_fit_summary(args.repo_dir)
    ad_file = open_root(args.adaptive_dir / "Step6p5AttenuationCorrections" / "pe_attenuation_scans.root")
    rp_file = open_root(args.repo_dir / "Step6p5AttenuationCorrections" / "pe_attenuation_scans.root")
    ad_names = {k.GetName() for k in ad_file.GetListOfKeys()}
    rp_names = {k.GetName() for k in rp_file.GetListOfKeys()}
    names = sorted(ad_names & rp_names)
    only_ad, only_rp = sorted(ad_names - rp_names), sorted(rp_names - ad_names)
    print("[ATTEN] {} common canvases; only adaptive: {}, only repo: {}".format(
        len(names), len(only_ad), len(only_rp)))

    args.outdir.mkdir(parents=True, exist_ok=True)
    out_root = ROOT.TFile(str(args.outdir / "pe_attenuation_scans_comparison.root"), "RECREATE")
    keys: List[str] = []
    jsons: Dict[str, str] = {}

    for index, name in enumerate(names, 1):
        ad_canvas, rp_canvas = ad_file.Get(name), rp_file.Get(name)
        ad_graph, rp_graph = find_graph(ad_canvas), find_graph(rp_canvas)
        if ad_graph is None or rp_graph is None:
            continue
        label = find_title(rp_canvas) or find_title(ad_canvas) or name

        graphs = []
        for source, color, marker in (
            (ad_graph, ADAPTIVE_COLOR, 20), (rp_graph, REPO_COLOR, 24),
        ):
            graph = source.Clone("{}_{}".format(source.GetName(), color))
            graph.SetTitle("")
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.05)
            graph.SetMarkerColor(color)
            graph.SetLineColor(color)
            graph.SetLineWidth(2)
            graphs.append(graph)

        canvas_name = name + "_comparison"
        canvas = ROOT.TCanvas(canvas_name, "", 1100, 800)
        canvas.SetLeftMargin(0.13)
        canvas.SetRightMargin(0.04)
        canvas.SetTopMargin(0.12)
        canvas.SetBottomMargin(0.12)
        canvas.SetTicks(1, 1)

        multi = ROOT.TMultiGraph()
        for graph in graphs:
            multi.Add(graph, "P")
        multi.Draw("A")
        multi.GetXaxis().SetTitle(ad_graph.GetXaxis().GetTitle() or "Table position X [mm]")
        multi.GetYaxis().SetTitle("Mean photoelectrons")
        for axis in (multi.GetXaxis(), multi.GetYaxis()):
            axis.SetTitleSize(0.045)
            axis.SetLabelSize(0.04)
        multi.GetYaxis().SetTitleOffset(1.35)

        # Headroom above the data so the four-line legend clears the curves.
        lows = [g.GetY()[i] - g.GetEY()[i] for g in graphs for i in range(g.GetN())]
        highs = [g.GetY()[i] + g.GetEY()[i] for g in graphs for i in range(g.GetN())]
        span = max(highs) - min(lows)
        multi.SetMinimum(min(lows) - 0.05 * span)
        multi.SetMaximum(max(highs) + 0.45 * span)

        branch = name[len("pe_attenuation_"):]
        x_min = min(g.GetX()[i] for g in graphs for i in range(g.GetN()))
        curves = [
            fit_curve(branch, ad_fits, x_min, ADAPTIVE_COLOR, "adaptive"),
            fit_curve(branch, rp_fits, x_min, REPO_COLOR, "repo"),
        ]
        for curve in curves:
            if curve is not None:
                curve.Draw("SAME")

        legend = ROOT.TLegend(0.15, 0.68, 0.80, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.026)
        legend.AddEntry(graphs[0], "Adaptive integral (#sigma/#sqrt{N} errors)", "pe")
        legend.AddEntry(graphs[1], "Repo integral (#sigma/#sqrt{N} errors)", "pe")
        for curve, fits, label_name in (
            (curves[0], ad_fits, "Adaptive fit"), (curves[1], rp_fits, "Repo fit"),
        ):
            if curve is not None:
                legend.AddEntry(curve, fit_label(fits["channels"].get(branch), label_name), "l")
            else:
                legend.AddEntry(ROOT.nullptr, fit_label(None, label_name), "")
        legend.Draw()
        header = draw_header(canvas, "PE attenuation scan: {}".format(label), "MCP-clean events")
        canvas.RedrawAxis()
        canvas.Update()

        out_root.cd()
        canvas.Write(canvas_name)
        keys.append(canvas_name)
        jsons[canvas_name] = to_json(canvas)
        canvas._keepalive = [multi, graphs, curves, legend, header]
        if index % 25 == 0 or index == len(names):
            print("[ATTEN] {}/{} channels".format(index, len(names)))

    out_root.Close()
    ad_file.Close()
    rp_file.Close()

    intro = "\n".join([
        "**Blue (filled):** adaptive integral, `results/PhotoElectrons/Step6p5AttenuationCorrections`. "
        "**Red (open):** repo integral, `results/PhotoElectrons_RepoIntegral/Step6p5AttenuationCorrections`.",
        "Each study applies its own run-1897 ADC-to-PE calibration to the same "
        "attenuation ADC means, so the points can differ in scale as well as shape.",
        "**Fits:** the stored per-channel exponential `PE(x) = A exp(kx)` of each "
        "study (`exponential_fit_per_channel`), drawn in that study's colour, "
        "with its attenuation length #lambda = 1/k. Fit range x <= {:g} mm, "
        "excluding {:g} to {:g} mm. Adaptive study: {} channels fitted, {} skipped; "
        "repo study: {} fitted, {} skipped.".format(
            ad_fits["fit_x_max_mm"], *ad_fits["excluded_x_range_mm"],
            len(ad_fits["channels"]), len(ad_fits["skipped_channels"]),
            len(rp_fits["channels"]), len(rp_fits["skipped_channels"]),
        ),
        "Channels in only one study: adaptive {}, repo {}.".format(len(only_ad), len(only_rp)),
    ])
    out_html = args.outdir / "pe_attenuation_scans_comparison.html"
    generate_jsroot_html(
        keys, jsons, plots_per_row=3, output_html=str(out_html),
        title="CaloX PE study: attenuation scans, adaptive vs repo integral",
        intro_text=intro,
    )
    print("[DONE] {}".format(out_html))


def main() -> None:
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not args.skip_overlay:
        if args.no_fits or not build_position_fit_overlay(args):
            build_position_overlay(args)
    if not args.skip_sigma:
        build_sigma_hists(args)
    if not args.skip_attenuation:
        build_attenuation_overlay(args)


if __name__ == "__main__":
    main()
