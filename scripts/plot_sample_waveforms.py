#!/usr/bin/env python3
"""Presentation-quality example waveforms, one particle type at a time.

Picks a handful of individual, real (non-noise) pulses per particle type and
plots each as its own clean, single-pulse figure -- meant for dropping
straight into a talk, not for quantitative shape analysis (see
plot_pid_waveformshape.py for that).

Run/particle selection:
  * Default is 80 GeV, tb2026, "nominal" table position (no rotation scans,
    no 90 deg run) -- see NOMINAL_RUNS below. Only electron and pion have
    80 GeV tb2026 data; muon/proton beams were never taken at 80 GeV in
    data/Runlist.json, so those particles are silently skipped unless
    --runs pins them explicitly.
  * --energy switches to auto-discovery (nominal-angle only) at a different
    energy; --runs pins exact run numbers per particle.

"Not noise" selection, applied before any waveform is eligible:
  1. Beam PID cut (configs/selection_config.py via apply_beam_pid_selection):
     the event is actually consistent with the requested particle type.
  2. MCP-clean timing cut: a real, in-time beam particle triggered the event.
  3. Per-pulse amplitude/trough cut on the baseline-subtracted waveform
     itself: rejects pedestal-only / sub-threshold pulses.
  4. Shape sanity cut: rise/fall/width crossings must all resolve to finite
     values (compute_shape_metrics, reused from plot_pid_waveformshape.py) --
     rejects garbled/multi-pulse traces.
  5. Single-pulse cut: outside a window around the main peak, the trace must
     stay below --max-secondary-frac of the peak -- rejects pileup (a second,
     unrelated particle/afterpulse in the same trace) so each example shows
     one clean pulse.

Within the surviving pool, examples are chosen spread across the peak-ADC
distribution (not just near the median) so the "bunch" of examples shows
real pulse-to-pulse variation rather than n near-duplicates.

Output: results/Samplewaveforms/<particle>/*.png (one pulse per file) plus
an overview grid per particle and an index.html gallery.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.html_generator import generate_html
from utils.root_setup import setup_root
from variables.drs import _TS_TO_NS

from plot_velocityfits import (
    ChannelInfo,
    active_timing_channels,
    display_material,
    make_manager_args,
)
from plot_pid_waveformshape import compute_shape_metrics
from core.analysis_manager import CaloXAnalysisManager

ALL_PARTICLES = ["electron", "pion", "muon", "proton"]
DEFAULT_PARTICLES = ["electron", "pion"]

# Hand-picked nominal (no table rotation), tb2026, 80 GeV runs -- vetted
# against data/Runlist.json on 2026-08-10: highest-statistics run per
# particle with a blank "Table Degree" field (i.e. not part of the
# angle-scan series 1968-1992, and not the 90 deg run 1997) and
# run_number > 1828 (avoids the DRS bridge-mapping bug documented in
# plot_pid_waveformshape.py for 1748 < run_number <= 1827). No muon/proton
# beam was ever taken at 80 GeV in tb2026. pion's Sci channel only has ~24
# clean pulses in the whole run (vs. hundreds for CerQuartz/CerPlastic) --
# thin, but it's the only nominal 80 GeV pion run available.
NOMINAL_RUNS = {
    80.0: {"electron": 1841, "pion": 1834},
}


def _style_presentation():
    plt.rcParams.update({
        "font.size": 13,
        "font.family": "sans-serif",
        "axes.linewidth": 1.1,
        "axes.edgecolor": "#333333",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.6,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "#cccccc",
    })


def _tick_style(ax):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=6)
    ax.tick_params(which="minor", length=3)


# =============================================================================
# RUN / CHANNEL SELECTION
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


def _parse_table_degree(value) -> Optional[float]:
    """'' -> 0.0 (nominal); '3.3' -> 3.3; '90 (88.2)' -> 90.0 (scan/rotated)."""
    s = str(value).strip()
    if not s or s == "-":
        return 0.0
    head = s.split("(")[0].strip()
    try:
        return float(head)
    except ValueError:
        return None


BEAM_TYPE_TO_PARTICLE = {
    "e+": "electron", "positron": "electron", "positrons": "electron",
    "pi+": "pion", "pion": "pion", "pions": "pion", "pions (+)": "pion",
    "mu+": "muon", "muon": "muon", "muons": "muon",
}


def discover_nominal_runs(energy_gev: float, particles: Sequence[str], tolerance_gev: float,
                          runlist_path: str, min_run: int, max_angle_deg: float) -> Dict[str, int]:
    """Highest-statistics tb2026 run per particle at energy_gev, restricted to
    a nominal (unrotated) table position -- excludes angle-scan series and the
    90 deg run."""
    with open(runlist_path) as handle:
        runlist = json.load(handle)

    best: Dict[str, tuple] = {}
    for run_str, info in runlist.items():
        try:
            run_number = int(run_str)
        except ValueError:
            continue
        if run_number <= min_run:
            continue

        particle = BEAM_TYPE_TO_PARTICLE.get(str(info.get("beam type", "")).strip().lower())
        if particle is None or particle not in particles:
            continue

        run_energy = _parse_energy_gev(info.get("beam energy", ""))
        if run_energy is None or abs(run_energy - energy_gev) > tolerance_gev:
            continue

        angle = _parse_table_degree(info.get("Table Degree", ""))
        if angle is None or abs(angle) > max_angle_deg:
            continue

        n_events = _parse_event_count(info.get("number of evt.", "0"))
        if n_events <= 0:
            continue

        if particle not in best or n_events > best[particle][0]:
            best[particle] = (n_events, run_number)

    return {particle: run for particle, (_, run) in best.items()}


def resolve_runs(args: argparse.Namespace) -> Dict[str, int]:
    runs: Dict[str, int] = {}
    if args.runs:
        for item in args.runs:
            key, _, val = item.partition("=")
            if not val:
                raise ValueError("--runs entries must look like particle=run: {!r}".format(item))
            runs[key.strip().lower()] = int(val)

    remaining = [p for p in args.particles if p not in runs]
    if remaining and args.energy is not None:
        table = NOMINAL_RUNS.get(args.energy)
        if table is not None:
            for p in remaining:
                if p in table:
                    runs[p] = table[p]
        remaining = [p for p in args.particles if p not in runs]
        if remaining:
            discovered = discover_nominal_runs(
                args.energy, remaining, args.tolerance, args.runlist,
                args.min_run, args.max_angle)
            missing = [p for p in remaining if p not in discovered]
            if missing:
                print("[WARN] No nominal 80-GeV-window tb2026 run for: {} "
                      "(no beam of that type was taken at this energy?)".format(missing))
            runs.update(discovered)

    return {p: runs[p] for p in args.particles if p in runs}


def pick_channel(channel_map: Dict[str, ChannelInfo], material: str, size: str = "3 mm") -> List[str]:
    return sorted(b for b, c in channel_map.items() if c.material == material and c.size == size)


# =============================================================================
# EVENT SELECTION
# =============================================================================

def collect_valid_pulses(manager: CaloXAnalysisManager, branch: str, args: argparse.Namespace) -> List[dict]:
    """Pull `{branch}_blsub` waveforms for PID+MCP-selected events and keep
    only those passing the amplitude/trough/shape "not noise" cuts."""
    data = manager.rdf.AsNumpy(["event_n", "{}_blsub".format(branch)])
    wf_col = data["{}_blsub".format(branch)]

    valid: List[dict] = []
    for row in range(len(data["event_n"])):
        wf = np.asarray(wf_col[row], dtype=float)
        t_abs = np.arange(wf.size) * _TS_TO_NS
        metrics = compute_shape_metrics(
            wf, t_abs, min_amplitude=args.min_amplitude, min_trough=args.min_trough,
            fall_search_ns=args.fall_search_ns, bin_width_ns=_TS_TO_NS)
        if metrics is None:
            continue
        if not all(np.isfinite(metrics[k]) for k in ("rise_time_ns", "fall_time_ns", "width_ns")):
            continue

        peak_idx = int(np.argmax(wf))
        t_peak = t_abs[peak_idx]
        outside = np.abs(t_abs - t_peak) > args.pileup_window_ns
        if outside.any() and np.max(wf[outside]) > args.max_secondary_frac * metrics["peak_adc"]:
            continue

        metrics["event_idx"] = int(data["event_n"][row])
        metrics["_waveform"] = wf
        metrics["_t_abs"] = t_abs
        valid.append(metrics)
    return valid


def pick_spread_examples(valid: List[dict], n_examples: int) -> List[dict]:
    """Spread picks across the peak-ADC distribution (roughly evenly spaced
    percentiles) so examples show real pulse-to-pulse variation instead of
    n near-duplicates clustered at the median. The top decile is excluded
    from the spread: the single largest pulses in a run are disproportionately
    likely to carry extra shower substructure (more photoelectrons, more
    secondary activity), which reads as a messier trace than is typical."""
    ordered = sorted(valid, key=lambda m: m["peak_adc"])
    if len(ordered) > n_examples:
        cap = max(n_examples, int(round(0.90 * len(ordered))))
        ordered = ordered[:cap]
    if len(ordered) <= n_examples:
        return ordered
    idx = np.linspace(0, len(ordered) - 1, n_examples).round().astype(int)
    idx = sorted(set(int(i) for i in idx))
    return [ordered[i] for i in idx]


# =============================================================================
# PLOTTING
# =============================================================================

def plot_single_waveform(ax, m: dict, particle: str, info: ChannelInfo, run: int,
                         shade: bool = True, absolute_time: bool = False):
    color = "black"
    wf, t = m["_waveform"], m["_t_abs"]
    t_peak = m["t_peak_ns"]
    lo, hi = max(0.0, t_peak - 15.0), min(t[-1], t_peak + 25.0)
    win = (t >= lo) & (t <= hi)
    shift = 0.0 if absolute_time else t_peak
    tw, ww = t[win] - shift, wf[win]

    if shade:
        ax.fill_between(tw, 0, ww, color=color, alpha=0.22, zorder=1)
    ax.plot(tw, ww, color=color, lw=2.2, zorder=3)
    ax.axhline(0, color="#999999", lw=0.9, zorder=0)

    pk = m["peak_adc"]
    ann = "Peak = {:.0f} ADC\nRise (10-90%) = {:.2f} ns\nWidth (FWHM) = {:.2f} ns".format(
        pk, m["rise_time_ns"], m["width_ns"])
    ax.text(0.97, 0.95, ann, transform=ax.transAxes, fontsize=10, va="top", ha="right",
            color="black",
            bbox=dict(fc="white", ec="#cccccc", boxstyle="round,pad=0.4", alpha=0.92))

    ax.set_xlabel("Time [ns]" if absolute_time else "Time - t$_{peak}$ [ns]",
                 fontsize=12, loc="right")
    ax.set_ylabel("A.U.", fontsize=12, color="black", loc="top", labelpad=8)
    ax.set_xlim(lo - shift, hi - shift)
    ax.set_ylim(bottom=-100, top=1.25 * pk)
    _tick_style(ax)
    return ax


def make_individual_plot(m: dict, particle: str, info: ChannelInfo, run: int, out_path: Path,
                         shade: bool = True, absolute_time: bool = False):
    fig, ax = plt.subplots(figsize=(8, 5.2))
    plot_single_waveform(ax, m, particle, info, run, shade=shade, absolute_time=absolute_time)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_overview_grid(examples_by_channel: Dict[str, List[dict]], particle: str,
                       channel_info: Dict[str, ChannelInfo], run: int, out_path: Path,
                       shade: bool = True, absolute_time: bool = False):
    branches = list(examples_by_channel.keys())
    n_cols = max((len(v) for v in examples_by_channel.values()), default=0)
    n_rows = len(branches)
    if n_rows == 0 or n_cols == 0 or n_cols > 8:
        return None

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.6 * n_cols, 3.6 * n_rows), squeeze=False)
    for r, branch in enumerate(branches):
        info = channel_info[branch]
        examples = examples_by_channel[branch]
        for c in range(n_cols):
            ax = axes[r][c]
            if c >= len(examples):
                ax.axis("off")
                continue
            m = examples[c]
            plot_single_waveform(ax, m, particle, info, run, shade=shade, absolute_time=absolute_time)
            ax.set_xlabel("")
            ax.text(0.97, 0.7, "evt {}".format(m["event_idx"]), transform=ax.transAxes,
                    fontsize=9, va="top", ha="right", color="#555555")
            if c == 0:
                ax.text(-0.32, 0.5, display_material(info.material), transform=ax.transAxes,
                       fontsize=12, fontweight="bold", va="center", ha="center", rotation=90)
        for a in axes[r]:
            pass
    xlabel = "t [ns]" if absolute_time else "t - t$_{peak}$ [ns]"
    for c in range(n_cols):
        axes[-1][c].set_xlabel(xlabel, fontsize=10, loc="right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


# =============================================================================
# DRIVER
# =============================================================================

_energy_tag_value = "80"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Save a handful of presentation-quality, individually "
                    "selected (non-noise) example waveforms per particle type.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--energy", type=float, default=80.0,
                       help="Beam energy in GeV. Runs are the vetted nominal-angle "
                            "tb2026 picks in NOMINAL_RUNS when available, else "
                            "auto-discovered (nominal angle only).")
    parser.add_argument("--tolerance", type=float, default=2.0)
    parser.add_argument("--max-angle", type=float, default=5.0,
                       help="Reject runs with |Table Degree| beyond this (excludes "
                            "the angle-scan series and the 90 deg run).")
    parser.add_argument("--particles", nargs="+", default=DEFAULT_PARTICLES, choices=ALL_PARTICLES)
    parser.add_argument("--runs", nargs="+", default=None, metavar="PARTICLE=RUN",
                       help="Pin exact runs, e.g. --runs electron=1841 pion=1834.")
    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--min-run", type=int, default=1828)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=20000)
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument("--materials", nargs="+", default=["Sci", "CerQuartz", "CerPlastic"],
                       help="Channel materials to sample (one representative "
                            "channel per material).")
    parser.add_argument("--size", default="3 mm", choices=["3 mm", "6 mm"])
    parser.add_argument("--n-examples", type=int, default=4,
                       help="Number of example pulses per particle per channel.")
    parser.add_argument("--min-amplitude", type=float, default=150.0,
                       help="Minimum baseline-subtracted peak ADC to accept a pulse as signal (not noise).")
    parser.add_argument("--min-trough", type=float, default=-100.0)
    parser.add_argument("--fall-search-ns", type=float, default=2.0)
    parser.add_argument("--pileup-window-ns", type=float, default=6.0,
                       help="Half-width around the main peak excluded from the "
                            "secondary-pulse check.")
    parser.add_argument("--max-secondary-frac", type=float, default=0.2,
                       help="Reject the pulse if anything outside the pileup "
                            "window exceeds this fraction of the peak (pileup cut).")
    parser.add_argument("--require-mcp-clean", dest="require_mcp_clean",
                       action="store_true", default=True)
    parser.add_argument("--no-mcp-clean", dest="require_mcp_clean", action="store_false")
    parser.add_argument("--shade", dest="shade", action="store_true", default=True,
                       help="Fill the area under each pulse.")
    parser.add_argument("--no-shade", dest="shade", action="store_false",
                       help="Plot the waveform as a bare line, no fill.")
    parser.add_argument("--absolute-time", action="store_true", default=False,
                       help="Plot the native waveform time axis instead of shifting so "
                            "the peak sits at t=0.")
    parser.add_argument("--outdir", default="results/Samplewaveforms")
    return parser.parse_args()


def main() -> None:
    global _energy_tag_value
    args = parse_args()
    _energy_tag_value = "{:g}".format(args.energy) if args.energy is not None else "?"
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    _style_presentation()

    run_map = resolve_runs(args)
    if not run_map:
        raise RuntimeError("No runs resolved for particles {}".format(args.particles))
    print("Resolved runs: {}".format(run_map))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gallery_files: List[str] = []
    summary_lines = []

    for particle in args.particles:
        if particle not in run_map:
            continue
        run = run_map[particle]
        pdir = outdir / particle
        pdir.mkdir(parents=True, exist_ok=True)

        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare()
        manager.apply_mcp_selection(flag_only=True)
        if args.require_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        manager.apply_beam_pid_selection(flag_only=False, particle=particle)

        channel_map = active_timing_channels(manager, requested=None)

        examples_by_channel: Dict[str, List[dict]] = {}
        channel_info: Dict[str, ChannelInfo] = {}
        for material in args.materials:
            candidates = pick_channel(channel_map, material, args.size)
            chosen_branch = None
            chosen_valid: List[dict] = []
            for branch in candidates:
                valid = collect_valid_pulses(manager, branch, args)
                if len(valid) >= min(args.n_examples, 2):
                    chosen_branch, chosen_valid = branch, valid
                    break
                if len(valid) > len(chosen_valid):
                    chosen_branch, chosen_valid = branch, valid
            if chosen_branch is None:
                print("[WARN] {}: no usable {} channel found".format(particle, material))
                continue

            examples = pick_spread_examples(chosen_valid, args.n_examples)
            examples_by_channel[chosen_branch] = examples
            channel_info[chosen_branch] = channel_map[chosen_branch]
            print("  {} | {} ({}): {} valid pulses, {} selected as examples".format(
                particle, chosen_branch, display_material(material), len(chosen_valid), len(examples)))

            for i, m in enumerate(examples, start=1):
                fname = "{}_{}_ex{}_evt{}.png".format(
                    particle, display_material(material), i, m["event_idx"])
                make_individual_plot(m, particle, channel_map[chosen_branch], run, pdir / fname,
                                     shade=args.shade, absolute_time=args.absolute_time)
                gallery_files.append("{}/{}".format(particle, fname))

        overview_path = outdir / "{}_overview.png".format(particle)
        made = make_overview_grid(examples_by_channel, particle, channel_info, run, overview_path,
                                  shade=args.shade, absolute_time=args.absolute_time)
        if made is not None:
            gallery_files.insert(0, overview_path.name)

        summary_lines.append("* **{}**: run {}, channels {}".format(
            particle.capitalize(), run,
            ", ".join(channel_map[b].short_label for b in examples_by_channel)))

        del manager

    html_path = outdir / "index.html"
    intro = (
        "**Energy:** {} GeV, tb2026, nominal table position (angle-scan runs and the "
        "90 deg run excluded).\n"
        "**Selection:** beam-PID cut + MCP-clean timing cut + per-pulse amplitude/trough "
        "cut (min {:.0f} ADC peak) + finite rise/fall/width shape check -- rejects "
        "noise/pedestal/garbled pulses.\n"
        + "\n".join(summary_lines)
    ).format(_energy_tag_value, args.min_amplitude)
    generate_html(gallery_files, str(outdir), plots_per_row=4, output_html=str(html_path),
                 title="Sample waveforms ({} GeV)".format(_energy_tag_value), intro_text=intro)
    print("Wrote {}".format(html_path))


if __name__ == "__main__":
    main()
