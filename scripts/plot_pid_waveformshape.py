#!/usr/bin/env python3
"""Waveform-shape PID study, integrated with the CaloXDataAnalysis framework.

Ported from the standalone WaveFormAnalysis_PIDv2.py script
(https://github.com/Akshat-Shrivastava11/Dream-Timing_analysis), adapted to
this repository's conventions:

  * Runs are looked up from data/Runlist.json instead of hardcoded file
    paths: give a beam energy and this script finds one run per particle
    type at that energy (tb2026 runs, run_number > 1700), or you can pin
    exact runs with --runs.
  * Particle-ID cuts reuse configs/selection_config.get_particle_selection
    via CaloXAnalysisManager.apply_beam_pid_selection() / SelectionManager,
    the same cuts used everywhere else in this framework, instead of
    re-deriving Cherenkov/veto thresholds by hand.
  * The comparison runs over every mapped, non-reference calorimeter DRS
    channel common to all requested runs (channels/channel_map.py), not
    just three representative "family" channels.

For each channel and particle type, per-pulse rise time (10%->90%,
backward search from the peak), fall time (90%->10%, forward search
clamped to --fall-search-ns), and pulse width (FWHM, 50%->50%) are computed
directly from the baseline-subtracted waveform (`{branch}_blsub`). Outputs:

  01_histograms_<tag>.pdf   overlaid rise/fall/width distributions per channel
  02_summary_<tag>.pdf      mean +/- std table per channel/particle
  03_waveforms_<tag>.pdf    annotated example waveforms (--waveform-examples)
  shapes_<tag>.csv          per-pulse metrics
  pid_waveformshape_<tag>.html   index linking the PDFs

Example:
  python scripts/plot_pid_waveformshape.py --energy 40 \\
      --particles electron pion muon --last-event 20000

Known caveat (found 2026-08-08): tb2026 runs with 1748 < run_number <= 1827
resolve to 0 active channels here (and in plot_velocityfits.py /
compare_timing_runs.py). channels/channel_map.py hardcodes DRS bridge_no=0
for run_number >= 1748, but the raw files in that window actually use
DRS_Brg1_*  branches, so every `{branch}_blsub` Define silently no-ops.
Runs with run_number > 1828 use the FERS-CSV-driven mapping and are fine
(verified for 1831/1832/1848). Prefer energies where all requested
particles resolve to runs > 1828 (e.g. 120 GeV: electron=1831, pion=1832,
muon=1848) until that mapping is fixed upstream.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

from utils.html_generator import generate_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from variables.drs import _TS_TO_NS

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    make_manager_args,
)
from core.analysis_manager import CaloXAnalysisManager

ALL_PARTICLES = ["electron", "pion", "muon", "proton"]
DEFAULT_PARTICLES = ["electron", "pion", "muon"]

# data/Runlist.json "beam type" strings seen across tb2025/tb2026 -> our
# canonical particle names (configs/selection_config.py keys).
BEAM_TYPE_TO_PARTICLE = {
    "e+": "electron", "positron": "electron", "positrons": "electron",
    "pi+": "pion", "pion": "pion", "pions": "pion", "pions (+)": "pion",
    "mu+": "muon", "muon": "muon", "muons": "muon",
}

PARTICLE_COLORS = {
    "electron": "#E45C3A", "pion": "#3A8FE4", "muon": "#3AC46E", "proton": "#B23AE4",
}

METRICS_INFO = {
    "rise_time_ns": ("Rise Time (10-90%) [ns]", 40),
    "fall_time_ns": ("Fall Time (90-10%) [ns]", 40),
    "width_ns":     ("Pulse Width FWHM [ns]", 40),
}


# =============================================================================
# RUN DISCOVERY
# =============================================================================

def _parse_energy_gev(value: str) -> Optional[float]:
    try:
        return float(str(value).replace("GeV", "").strip())
    except (TypeError, ValueError):
        return None


def _parse_event_count(value: str) -> int:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return 0


def discover_runs(energy_gev: float, particles: Sequence[str], tolerance_gev: float,
                  runlist_path: str, min_run: int = 1828) -> Dict[str, int]:
    """Pick the highest-statistics tb2026 run per particle at a given beam energy.

    min_run defaults to 1828: run_number 1748-1827 uses the FERS-CSV-driven
    DRS mapping's predecessor, which hardcodes the wrong bridge number for
    that window (see module docstring) and silently resolves to 0 active
    channels. run_number > 1828 is unaffected. Pass a lower --min-run once
    that mapping bug is fixed upstream, or to deliberately target that range
    with an explicit --runs override.
    """
    with open(runlist_path) as handle:
        runlist = json.load(handle)

    best: Dict[str, Tuple[int, int]] = {}  # particle -> (n_events, run_number)
    for run_str, info in runlist.items():
        try:
            run_number = int(run_str)
        except ValueError:
            continue
        if run_number <= min_run:
            continue

        beam_type = str(info.get("beam type", "")).strip().lower()
        particle = BEAM_TYPE_TO_PARTICLE.get(beam_type)
        if particle is None or particle not in particles:
            continue

        run_energy = _parse_energy_gev(info.get("beam energy", ""))
        if run_energy is None or abs(run_energy - energy_gev) > tolerance_gev:
            continue

        n_events = _parse_event_count(info.get("number of evt.", "0"))
        if n_events <= 0:
            continue

        if particle not in best or n_events > best[particle][0]:
            best[particle] = (n_events, run_number)

    return {particle: run for particle, (_, run) in best.items()}


def resolve_runs(args: argparse.Namespace) -> Dict[str, List[int]]:
    runs: Dict[str, List[int]] = {}
    if args.runs:
        for item in args.runs:
            key, _, val = item.partition("=")
            if not val:
                raise ValueError(
                    "--runs entries must look like particle=run[,run...]: {!r}".format(item))
            runs[key.strip().lower()] = [int(x) for x in val.split(",")]

    remaining = [p for p in args.particles if p not in runs]
    if remaining:
        if args.energy is None:
            raise ValueError(
                "--energy is required to auto-discover runs for: {}. "
                "Pass --runs {}=RUN instead, or --energy.".format(
                    remaining, remaining[0]))
        discovered = discover_runs(args.energy, remaining, args.tolerance, args.runlist, args.min_run)
        missing = [p for p in remaining if p not in discovered]
        if missing:
            print("[WARN] No tb2026 run within {:.1f} GeV of {:.1f} GeV for: {}".format(
                args.tolerance, args.energy, ", ".join(missing)))
        for particle, run in discovered.items():
            runs[particle] = [run]

    return {p: runs[p] for p in args.particles if p in runs}


# =============================================================================
# WAVEFORM SHAPE METRICS
# =============================================================================

def _crossing_backward(t: np.ndarray, wf: np.ndarray, level: float, i_start: int) -> float:
    """First upward crossing of `level`, searching backward from i_start (the peak).

    Scanning from the peak toward the baseline avoids latching onto noise
    crossings far from the pulse.
    """
    for i in range(i_start, 0, -1):
        if wf[i] >= level > wf[i - 1]:
            frac = (level - wf[i - 1]) / (wf[i] - wf[i - 1])
            return t[i - 1] + frac * (t[i] - t[i - 1])
    return float("nan")


def _crossing_forward(t: np.ndarray, wf: np.ndarray, level: float, i_start: int, i_stop: int) -> float:
    """First downward crossing of `level`, searching forward from i_start, clamped to i_stop."""
    for i in range(i_start, i_stop):
        if wf[i] >= level > wf[i + 1]:
            frac = (wf[i] - level) / (wf[i] - wf[i + 1])
            return t[i] + frac * (t[i + 1] - t[i])
    return float("nan")


def compute_shape_metrics(wf: np.ndarray, t_abs: np.ndarray, *, min_amplitude: float,
                          min_trough: float, fall_search_ns: float,
                          bin_width_ns: float) -> Optional[dict]:
    """Rise/fall/width metrics for one baseline-subtracted waveform.

    Falling-edge crossings are only searched within
    [t_peak, t_peak + fall_search_ns] so late shower sub-pulses (common in
    scintillator-family channels) cannot be mistaken for the main pulse's
    90%/50%/10% crossings.
    """
    peak_val = float(np.max(wf))
    trough = float(np.min(wf))
    if peak_val < min_amplitude or trough < min_trough:
        return None
    peak_idx = int(np.argmax(wf))

    lvl10, lvl50, lvl90 = 0.10 * peak_val, 0.50 * peak_val, 0.90 * peak_val

    t10_rise = _crossing_backward(t_abs, wf, lvl10, peak_idx)
    t50_rise = _crossing_backward(t_abs, wf, lvl50, peak_idx)
    t90_rise = _crossing_backward(t_abs, wf, lvl90, peak_idx)

    fall_stop = min(len(wf) - 1, peak_idx + int(math.ceil(fall_search_ns / bin_width_ns)))
    t90_fall = _crossing_forward(t_abs, wf, lvl90, peak_idx, fall_stop)
    t10_fall = _crossing_forward(t_abs, wf, lvl10, peak_idx, fall_stop)
    t50_fall = _crossing_forward(t_abs, wf, lvl50, peak_idx, fall_stop)

    rise = (t90_rise - t10_rise) if np.isfinite(t90_rise) and np.isfinite(t10_rise) else float("nan")
    fall = (t10_fall - t90_fall) if np.isfinite(t10_fall) and np.isfinite(t90_fall) else float("nan")
    width = (t50_fall - t50_rise) if np.isfinite(t50_fall) and np.isfinite(t50_rise) else float("nan")

    return {
        "rise_time_ns": rise, "fall_time_ns": fall, "width_ns": width,
        "peak_adc": peak_val, "t_peak_ns": float(t_abs[peak_idx]),
        "t_10_rise": t10_rise, "t_50_rise": t50_rise, "t_90_rise": t90_rise,
        "t_90_fall": t90_fall, "t_10_fall": t10_fall, "t_50_fall": t50_fall,
    }


def collect_channel_metrics(manager: CaloXAnalysisManager, particle: str, run: int,
                            branches: Sequence[str], args: argparse.Namespace,
                            ) -> Dict[str, List[dict]]:
    """Pull `{branch}_blsub` waveforms for the selected events and compute shape
    metrics, keeping up to --n-events valid pulses per channel."""
    columns = ["event_n"] + ["{}_blsub".format(b) for b in branches]
    data = manager.rdf.AsNumpy(columns)
    n_available = len(data["event_n"])

    results: Dict[str, List[dict]] = {}
    for branch in branches:
        wf_col = data["{}_blsub".format(branch)]
        collected: List[dict] = []
        for row in range(n_available):
            if len(collected) >= args.n_events:
                break
            wf = np.asarray(wf_col[row], dtype=float)
            t_abs = np.arange(wf.size) * _TS_TO_NS
            metrics = compute_shape_metrics(
                wf, t_abs, min_amplitude=args.min_amplitude, min_trough=args.min_trough,
                fall_search_ns=args.fall_search_ns, bin_width_ns=_TS_TO_NS)
            if metrics is None:
                continue
            metrics["event_idx"] = int(data["event_n"][row])
            metrics["run"] = run
            metrics["particle"] = particle
            metrics["branch"] = branch
            if args.waveform_examples:
                metrics["_waveform"] = wf
                metrics["_t_abs"] = t_abs
            collected.append(metrics)
        results[branch] = collected
    return results


# =============================================================================
# PLOT HELPERS
# =============================================================================

def gauss(x, amplitude, mu, sigma):
    return amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def _tick_style(ax):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)


def _auto_lim(values_list: Sequence[np.ndarray], pad_frac: float = 0.10) -> Tuple[float, float]:
    finite = [v[np.isfinite(v)] for v in values_list if len(v)]
    all_v = np.concatenate(finite) if finite else np.array([])
    if all_v.size == 0:
        return (0.0, 1.0)
    lo, hi = np.percentile(all_v, 1), np.percentile(all_v, 99)
    pad = max((hi - lo) * pad_frac, 0.05)
    return (lo - pad, hi + pad)


def particle_label(particle: str, runs: Sequence[int]) -> str:
    return "{} (run {})".format(particle.capitalize(), ",".join(str(r) for r in runs))


# =============================================================================
# PDF 1 - HISTOGRAMS (one page per channel x metric, particles overlaid)
# =============================================================================

def make_histograms_pdf(all_metrics, particle_list, channel_info, branches, outdir, tag, labels):
    pdf_path = outdir / "01_histograms_{}.pdf".format(tag)
    with PdfPages(pdf_path) as pdf:
        for branch in branches:
            info = channel_info[branch]
            for metric, (xlabel, nbins) in METRICS_INFO.items():
                fig, ax = plt.subplots(figsize=(9, 6.5))
                all_vals = [
                    np.array([m[metric] for m in all_metrics[p][branch]
                             if np.isfinite(m.get(metric, np.nan))])
                    for p in particle_list
                ]
                xlim = _auto_lim(all_vals)

                for particle, vals in zip(particle_list, all_vals):
                    if len(vals) == 0:
                        continue
                    color = PARTICLE_COLORS.get(particle, "#888888")
                    counts, edges = np.histogram(vals, bins=nbins, range=xlim)
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    peak_count = counts.max() if counts.max() > 0 else 1
                    norm_counts = counts / peak_count

                    ax.fill_between(centers, norm_counts, step="mid", color=color, alpha=0.25)
                    ax.step(centers, norm_counts, where="mid", color=color, lw=1.8, alpha=0.90)

                    try:
                        popt, _ = curve_fit(
                            gauss, centers, norm_counts,
                            p0=[1.0, float(np.nanmean(vals)), max(float(np.nanstd(vals)), 1e-3)],
                            maxfev=3000)
                        x_fit = np.linspace(*xlim, 300)
                        ax.plot(x_fit, gauss(x_fit, *popt), color=color, lw=2.2, ls="--")
                        mu_s, sig_s = "{:.3f}".format(popt[1]), "{:.3f}".format(abs(popt[2]))
                    except Exception:
                        mu_s = "{:.3f}".format(np.nanmean(vals))
                        sig_s = "{:.3f}".format(np.nanstd(vals))

                    ax.plot([], [], color=color, lw=2,
                           label="{}\n  N={}, mu={} ns, sigma={} ns".format(
                               labels[particle], len(vals), mu_s, sig_s))

                ax.set_xlabel(xlabel, fontsize=13)
                ax.set_ylabel("Normalised pulses", fontsize=12)
                ax.set_xlim(xlim)
                ax.set_ylim(0, 1.12)
                ax.legend(fontsize=9, frameon=True, loc="upper right")
                ax.set_title("{}\n{}".format(info.short_label, xlabel), fontsize=12)
                _tick_style(ax)
                fig.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

    print("  [OK] Histograms PDF -> {}".format(pdf_path))
    return pdf_path


# =============================================================================
# PDF 2 - SUMMARY TABLE
# =============================================================================

def make_summary_pdf(all_metrics, particle_list, channel_info, branches, outdir, tag, run_map, labels):
    pdf_path = outdir / "02_summary_{}.pdf".format(tag)

    rows = []
    for branch in branches:
        info = channel_info[branch]
        for particle in particle_list:
            mlist = all_metrics[particle][branch]
            row = [info.short_label if particle == particle_list[0] else "", labels[particle]]
            for metric in ("rise_time_ns", "fall_time_ns", "width_ns"):
                vals = np.array([m[metric] for m in mlist if np.isfinite(m.get(metric, np.nan))])
                row.append("{:.3f} +/- {:.3f}".format(np.mean(vals), np.std(vals))
                          if len(vals) >= 5 else "-")
            row.append(str(len(mlist)))
            rows.append(row)

    header = ["Channel", "Particle", "Rise mu+/-sigma [ns]", "Fall mu+/-sigma [ns]",
              "FWHM mu+/-sigma [ns]", "N"]
    rows_per_page = 24

    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(14, 6))
        ax = fig.add_subplot(111)
        ax.axis("off")
        lines = [
            r"$\bf{Waveform\ Shape\ PID\ Study}$", "",
            "Runs: " + "; ".join(
                "{}={}".format(p, ",".join(str(r) for r in run_map[p])) for p in particle_list),
            "Channels: {}  |  Particles: {}".format(len(branches), ", ".join(labels[p] for p in particle_list)),
            "Metrics: Rise Time (10->90%) | Fall Time (90->10%) | Pulse Width (FWHM 50->50%)",
        ]
        ax.text(0.5, 0.5, "\n".join(lines), ha="center", va="center", fontsize=12,
               transform=ax.transAxes,
               bbox=dict(fc="#f7f7f7", ec="#aaaaaa", boxstyle="round,pad=0.9"))
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        for start in range(0, len(rows), rows_per_page):
            chunk = rows[start:start + rows_per_page]
            fig = plt.figure(figsize=(16, 0.34 * len(chunk) + 1.4))
            ax = fig.add_subplot(111)
            ax.axis("off")
            tbl = ax.table(cellText=chunk, colLabels=header, loc="center", cellLoc="center")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(9)
            tbl.scale(1, 1.7)
            for (r, _c), cell in tbl.get_celld().items():
                if r == 0:
                    cell.set_facecolor("#2c3e50")
                    cell.set_text_props(color="white", fontweight="bold")
                elif r % 2 == 0:
                    cell.set_facecolor("#f0f4f8")
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print("  [OK] Summary PDF -> {}".format(pdf_path))
    return pdf_path


# =============================================================================
# PDF 3 - EXAMPLE WAVEFORMS (optional, --waveform-examples)
# =============================================================================

def make_waveforms_pdf(all_metrics, particle_list, channel_info, branches, outdir, tag,
                       n_examples, fall_search_ns, labels):
    pdf_path = outdir / "03_waveforms_{}.pdf".format(tag)
    with PdfPages(pdf_path) as pdf:
        for branch in branches:
            info = channel_info[branch]
            for particle in particle_list:
                color = PARTICLE_COLORS.get(particle, "#888888")
                mlist = [m for m in all_metrics[particle][branch] if "_waveform" in m]
                if not mlist:
                    continue
                peaks = np.array([m["peak_adc"] for m in mlist])
                med_pk = np.median(peaks)
                mlist.sort(key=lambda m: abs(m["peak_adc"] - med_pk))
                examples = mlist[:n_examples]

                for ex_num, m in enumerate(examples):
                    fig, ax = plt.subplots(figsize=(10, 6))
                    wf, t = m["_waveform"], m["_t_abs"]
                    t_peak = m["t_peak_ns"]
                    lo, hi = max(0.0, t_peak - 15.0), min(t[-1], t_peak + 25.0)
                    win = (t >= lo) & (t <= hi)
                    tw, ww = t[win], wf[win]
                    tr = tw - t_peak

                    ax.plot(tr, ww, color=color, lw=1.8, zorder=3, label="Waveform")
                    pk = m["peak_adc"]

                    t10r, t90r = m.get("t_10_rise", np.nan), m.get("t_90_rise", np.nan)
                    if np.isfinite(t10r) and np.isfinite(t90r):
                        shade = (tw >= t10r) & (tw <= t90r)
                        ax.fill_between(tw[shade] - t_peak, 0, ww[shade], color="#3A8FE4",
                                       alpha=0.30, label="Rise 10-90%", zorder=2)

                    t90f, t10f = m.get("t_90_fall", np.nan), m.get("t_10_fall", np.nan)
                    fall_limit = t_peak + fall_search_ns
                    if np.isfinite(t90f) and np.isfinite(t10f):
                        shade = (tw >= t90f) & (tw <= min(t10f, fall_limit))
                        ax.fill_between(tw[shade] - t_peak, 0, ww[shade], color="#E45C3A",
                                       alpha=0.30, label="Fall 90-10%", zorder=2)
                    ax.axvline(fall_limit - t_peak, color="#E45C3A", lw=0.9, ls=(0, (3, 4)),
                              alpha=0.6, label="Fall window (+{:.0f} ns)".format(fall_search_ns))

                    for tv, lbl in [(m.get("t_50_rise", np.nan), "50% rise"),
                                    (m.get("t_50_fall", np.nan), "50% fall")]:
                        if np.isfinite(tv):
                            ax.axvline(tv - t_peak, color="purple", lw=1.0, ls=":",
                                      alpha=0.7, label=lbl)

                    for frac, lbl in [(0.10, "10% pk"), (0.50, "50% pk"), (0.90, "90% pk")]:
                        ax.axhline(frac * pk, color="gray", lw=0.8, ls=":", alpha=0.5,
                                  label=lbl if frac == 0.10 else None)
                    ax.axvline(0, color="gray", lw=0.9, ls="--", alpha=0.55, label="t_peak")

                    ann = ("Rise (10->90%) = {:.3f} ns\nFall (90->10%) = {:.3f} ns\n"
                          "Width (FWHM)  = {:.3f} ns").format(
                        m.get("rise_time_ns", float("nan")),
                        m.get("fall_time_ns", float("nan")),
                        m.get("width_ns", float("nan")))
                    ax.text(0.97, 0.97, ann, transform=ax.transAxes, fontsize=9, va="top",
                           ha="right", bbox=dict(fc="white", ec="#cccccc",
                                                 boxstyle="round,pad=0.45", alpha=0.85))

                    ax.set_xlabel("Time - t_peak [ns]", fontsize=12)
                    ax.set_ylabel("ADC (baseline subtracted)", fontsize=12)
                    ax.set_title("{}  |  {}\nExample {}/{}  -  event {}  |  peak = {:.0f} ADC".format(
                        labels[particle], info.short_label, ex_num + 1, len(examples),
                        m["event_idx"], pk), fontsize=11)
                    ax.set_xlim(lo - t_peak, hi - t_peak)
                    ax.set_ylim(bottom=-0.06 * pk)
                    ax.legend(fontsize=8, frameon=True, loc="upper left", ncol=2)
                    _tick_style(ax)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

    print("  [OK] Waveforms PDF -> {}".format(pdf_path))
    return pdf_path


# =============================================================================
# CSV
# =============================================================================

def write_csv(all_metrics, particle_list, channel_info, outdir, tag):
    import csv
    csv_path = outdir / "shapes_{}.csv".format(tag)
    fieldnames = ["particle", "branch", "material", "size", "run", "event_idx",
                 "peak_adc", "t_peak_ns", "rise_time_ns", "fall_time_ns", "width_ns",
                 "t_10_rise", "t_50_rise", "t_90_rise", "t_90_fall", "t_10_fall", "t_50_fall"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for particle in particle_list:
            for branch, mlist in all_metrics[particle].items():
                info = channel_info[branch]
                for m in mlist:
                    row = {k: v for k, v in m.items() if not k.startswith("_")}
                    row["material"] = info.material
                    row["size"] = info.size
                    writer.writerow(row)
    print("  [OK] CSV -> {}".format(csv_path))
    return csv_path


# =============================================================================
# DRIVER
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Waveform-shape PID study: compare rise/fall/FWHM pulse "
                    "shapes across particle types at a fixed beam energy, for "
                    "every mapped calorimeter DRS channel.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--energy", type=float, default=None,
                       help="Beam energy in GeV; auto-discovers one run per "
                            "particle from --runlist (tb2026 runs only).")
    parser.add_argument("--tolerance", type=float, default=2.0,
                       help="Allowed |run energy - --energy| window in GeV.")
    parser.add_argument("--particles", nargs="+", default=DEFAULT_PARTICLES,
                       choices=ALL_PARTICLES)
    parser.add_argument("--runs", nargs="+", default=None, metavar="PARTICLE=RUN[,RUN...]",
                       help="Explicit run override, e.g. --runs electron=2016 pion=2018 "
                            "muon=2041. Bypasses --energy discovery for the given "
                            "particles; comma-separated runs are concatenated.")
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--min-run", type=int, default=1828,
                       help="Lowest run_number considered by --energy auto-discovery. "
                            "Default skips the 1748-1827 window, which resolves to 0 "
                            "active channels due to a DRS bridge-mapping bug in "
                            "channels/channel_map.py (see module docstring). Ignored "
                            "for particles pinned via --runs.")
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=20000,
                       help="Upper bound on events scanned per run so this "
                            "all-channel study stays tractable. -1 for the full run.")
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument("--channels", nargs="+", default=None,
                       help="Restrict to specific DRS branch names; default is every "
                            "mapped, non-reference calorimeter channel common to all "
                            "requested runs.")
    parser.add_argument("--thickness", choices=["3mm", "6mm", "both"], default="both")
    parser.add_argument("--require-mcp-clean", dest="require_mcp_clean",
                       action="store_true", default=True,
                       help="Also require a clean MCP timing reference (recommended).")
    parser.add_argument("--no-mcp-clean", dest="require_mcp_clean", action="store_false")
    parser.add_argument("--n-events", type=int, default=500,
                       help="Max valid pulses collected per particle per channel.")
    parser.add_argument("--min-amplitude", type=float, default=100.0,
                       help="Minimum baseline-subtracted peak ADC to accept a pulse.")
    parser.add_argument("--min-trough", type=float, default=-100.0,
                       help="Minimum allowed undershoot (ADC) before rejecting a pulse.")
    parser.add_argument("--fall-search-ns", type=float, default=2.0,
                       help="Falling-edge crossings are only searched within this many "
                            "ns after the peak, to reject shower sub-pulses.")
    parser.add_argument("--waveform-examples", action="store_true",
                       help="Also render a page of annotated example waveforms per "
                            "channel/particle (slower, much larger PDF).")
    parser.add_argument("--n-examples", type=int, default=2)
    parser.add_argument("--outdir", default=None)
    return parser.parse_args()


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)

    run_map = resolve_runs(args)
    if len(run_map) < 2:
        raise RuntimeError(
            "Need resolved runs for at least two particle types to compare; got: {}".format(run_map))

    tag = "energy{:g}GeV".format(args.energy) if args.energy is not None else "custom_runs"
    outdir = Path(args.outdir) if args.outdir else Path("results/pid_waveformshape") / tag
    outdir.mkdir(parents=True, exist_ok=True)

    particle_list = [p for p in args.particles if p in run_map]
    labels = {p: particle_label(p, run_map[p]) for p in particle_list}
    print("Resolved runs: {}".format(", ".join("{}={}".format(p, run_map[p]) for p in particle_list)))

    requested = set(args.channels) if args.channels else None
    wanted_size = {"3mm": "3 mm", "6mm": "6 mm"}.get(args.thickness)

    per_particle_channels: Dict[str, Dict[str, ChannelInfo]] = {}
    per_particle_metrics: Dict[str, Dict[str, List[dict]]] = {}

    for particle in particle_list:
        combined_metrics: Dict[str, List[dict]] = defaultdict(list)
        channel_map: Dict[str, ChannelInfo] = {}

        for run in run_map[particle]:
            manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare()
            manager.apply_mcp_selection(flag_only=True)
            if args.require_mcp_clean:
                manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
            manager.apply_beam_pid_selection(flag_only=False, particle=particle)

            channels = active_timing_channels(manager, requested=requested)
            if wanted_size is not None:
                channels = {b: c for b, c in channels.items() if c.size == wanted_size}
            channel_map.update(channels)

            branches = sorted(channels)
            print("[{}] run {}: {} channels, scanning up to {} events".format(
                particle, run, len(branches),
                args.last_event if args.last_event >= 0 else "all"))

            run_metrics = collect_channel_metrics(manager, particle, run, branches, args)
            for branch, rows in run_metrics.items():
                combined_metrics[branch].extend(rows)
            del manager

        per_particle_channels[particle] = channel_map
        per_particle_metrics[particle] = combined_metrics

    common_branches = set.intersection(*[set(d) for d in per_particle_channels.values()])
    if not common_branches:
        raise RuntimeError("No calorimeter channel is common to every requested run/particle.")

    channel_info: Dict[str, ChannelInfo] = {}
    for d in per_particle_channels.values():
        channel_info.update({b: c for b, c in d.items() if b in common_branches})
    branches = sorted(common_branches, key=lambda b: channel_info[b].short_label)

    all_metrics = {
        p: {b: per_particle_metrics[p].get(b, []) for b in branches} for p in particle_list
    }
    for particle in particle_list:
        for branch in branches:
            n = len(all_metrics[particle][branch])
            print("  {} | {}: {} valid pulses".format(particle, channel_info[branch].short_label, n))

    print("\nGenerating outputs in {}".format(outdir))
    write_csv(all_metrics, particle_list, channel_info, outdir, tag)
    pdf_hist = make_histograms_pdf(all_metrics, particle_list, channel_info, branches, outdir, tag, labels)
    pdf_summary = make_summary_pdf(all_metrics, particle_list, channel_info, branches, outdir, tag, run_map, labels)
    pdf_names = [pdf_hist.stem, pdf_summary.stem]

    if args.waveform_examples:
        pdf_wave = make_waveforms_pdf(all_metrics, particle_list, channel_info, branches, outdir,
                                      tag, args.n_examples, args.fall_search_ns, labels)
        pdf_names.append(pdf_wave.stem)

    html_path = outdir / "pid_waveformshape_{}.html".format(tag)
    intro = (
        "**Runs:** " + ", ".join("{}={}".format(p, ",".join(str(r) for r in run_map[p]))
                                 for p in particle_list) + ".\n"
        "**Channels:** {} common calorimeter channels.\n"
        "**Selection:** {}beam-PID selection (configs/selection_config.py).\n"
        "Metrics: rise time (10-90%, backward search from peak), fall time (90-10%, "
        "forward search clamped to {:.1f} ns after the peak), pulse width "
        "(FWHM, 50-50%)."
    ).format(len(branches), "MCP-clean + " if args.require_mcp_clean else "", args.fall_search_ns)

    generate_html(
        ["PDF_LINK:{}".format(name) for name in pdf_names], str(outdir),
        plots_per_row=1, output_html=str(html_path),
        title="Waveform-shape PID study ({})".format(tag), intro_text=intro,
    )
    print("Wrote {}".format(html_path))


if __name__ == "__main__":
    main()
