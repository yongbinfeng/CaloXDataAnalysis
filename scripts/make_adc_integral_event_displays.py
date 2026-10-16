#!/usr/bin/env python3
"""Proof-of-concept pulse displays with the CFD ADC integral highlighted."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.html_generator import generate_html
from utils.root_setup import setup_root
from utils.timing import auto_timer

from variables.drs import add_cfd_to_baseline_energy, get_ref_ts_name
from plot_velocityfits import active_timing_channels, make_manager_args


WINDOW_PRE = 5
WINDOW_POST = 45


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Highlight the CFD ADC integration window on individual DRS pulses."
    )
    parser.add_argument("--run", type=int, default=1917)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=1000)
    parser.add_argument("--n-events", type=int, default=3)
    parser.add_argument("--channels-per-event", type=int, default=8)
    parser.add_argument(
        "--size",
        choices=("all", "3mm", "6mm"),
        default="all",
        help=(
            "Restrict the candidate channel pool by fiber size. Event ranking "
            "uses the same pool, so 6mm picks events with real 6 mm signal. "
            "Default: all."
        ),
    )
    parser.add_argument(
        "--pin-channels",
        nargs="+",
        default=None,
        help=(
            "Channels always shown first, whatever their rank and size. Accepts "
            "a full branch name (DRS_Brg1_Board1_Group1_Channel3) or a short "
            "spec (\"B1 G1 Ch3\", \"Brg1B1G1Ch3\"). Short specs that match on "
            "more than one bridge are rejected - qualify them with Brg<n>."
        ),
    )
    parser.add_argument(
        "--compare-windows", action="store_true",
        help=(
            "Overlay the competing integration windows on each pulse instead of "
            "just the repo one: the fixed 50-slice CFD window, the "
            "reference-anchored fixed window, and the adaptive CFD-to-baseline "
            "window. Each is annotated with the charge it collects."
        ),
    )
    parser.add_argument(
        "--no-fixed-overlay", action="store_true",
        help=(
            "Drop the reference-anchored fixed window from the overlay, leaving "
            "the repo CFD window and the adaptive CFD-to-baseline window."
        ),
    )
    parser.add_argument(
        "--all-channels", action="store_true",
        help=(
            "Draw every channel in the pool rather than the strongest N. The "
            "figure is split across several PNGs per event; see "
            "--panels-per-figure."
        ),
    )
    parser.add_argument(
        "--panels-per-figure", type=int, default=12,
        help="Panels on one PNG when --all-channels is used. Default: 12.",
    )
    parser.add_argument("--anchor-offset", type=int, default=-330,
                        help="Slices from the group reference crossing to the pulse.")
    parser.add_argument("--fixed-pre", type=int, default=4)
    parser.add_argument("--fixed-post", type=int, default=35)
    parser.add_argument("--stop-abs", type=float, default=4.0,
                        help="Absolute floor of the baseline-crossing stop level [ADC].")
    parser.add_argument("--stop-fraction", type=float, default=0.05)
    parser.add_argument("--stop-consecutive", type=int, default=3)
    parser.add_argument("--min-width", type=int, default=40,
                        help="Floor on the adaptive window width [slices].")
    parser.add_argument("--smooth-half", type=int, default=3,
                        help="Half-width of the moving average used for the stop test.")
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument(
        "--outdir", default="results/attenuation/adc_integral_event_displays_run1917"
    )
    return parser.parse_args()


SHORT_SPEC = re.compile(
    r"^(?:brg(?P<bridge>\d+))?\s*b(?P<board>\d+)\s*g(?P<group>\d+)\s*ch(?P<channel>\d+)$"
)


def branch_key(branch: str) -> str:
    """Bridge/board/group/channel of a branch name, as a comparable key."""
    match = re.search(
        r"DRS_(?:Brg(\d+)_)?Board(\d+)_Group(\d+)_Channel(\d+)$", branch
    )
    if match is None:
        return branch
    bridge = match.group(1) if match.group(1) is not None else "0"
    return "brg{}b{}g{}ch{}".format(bridge, *match.group(2, 3, 4))


def resolve_channel_spec(spec: str, branches: List[str]) -> str:
    """Map a user channel spec onto exactly one branch name.

    Short specs like "B1 G1 Ch3" omit the bridge, and that label is not unique:
    Brg0 and Brg1 both carry a Board1/Group1/Channel3. Rather than silently
    picking one, an ambiguous spec is an error naming both candidates.
    """
    if spec in branches:
        return spec
    match = SHORT_SPEC.match(spec.strip().lower().replace("_", " "))
    if match is None:
        raise ValueError(
            "Cannot parse channel spec {!r}. Use a branch name or 'B1 G1 Ch3'.".format(spec)
        )
    if match.group("bridge") is not None:
        wanted = "brg{bridge}b{board}g{group}ch{channel}".format(**match.groupdict())
        matches = [b for b in branches if branch_key(b) == wanted]
    else:
        tail = "b{board}g{group}ch{channel}".format(**match.groupdict())
        matches = [b for b in branches if branch_key(b).endswith(tail)]
    if not matches:
        raise ValueError("No active channel matches {!r}".format(spec))
    if len(matches) > 1:
        raise ValueError(
            "{!r} is ambiguous across bridges: {}. Qualify it, e.g. 'Brg1{}'.".format(
                spec, ", ".join(sorted(matches)), spec.strip().replace(" ", "")
            )
        )
    return matches[0]


def panel_title(branch: str, channel) -> str:
    """Panel heading that keeps the bridge, unlike ChannelInfo.short_label."""
    match = re.search(r"DRS_(?:Brg(\d+)_)?Board", branch)
    bridge = match.group(1) if match and match.group(1) is not None else "0"
    return "Brg{} {}".format(bridge, channel.short_label)


def event_mask(values: np.ndarray, event_number: int) -> np.ndarray:
    return np.asarray(values, dtype=np.int64) == event_number


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    manager = CaloXAnalysisManager(make_manager_args(args, args.run)).prepare(do_fers=False)
    manager.apply_mcp_selection(flag_only=True)
    manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
    channels = active_timing_channels(manager, requested=None)
    branches = sorted(channels)
    if not branches:
        raise RuntimeError("No active calorimeter DRS channels found")

    pinned: List[str] = []
    for spec in (args.pin_channels or []):
        resolved = resolve_channel_spec(spec, branches)
        if resolved not in pinned:
            pinned.append(resolved)
    if pinned:
        print("Pinned channels: {}".format(
            ", ".join("{} ({})".format(b, panel_title(b, channels[b])) for b in pinned)
        ))

    size_label = {"3mm": "3 mm", "6mm": "6 mm"}.get(args.size)
    if size_label is None:
        pool = list(branches)
    else:
        pool = [b for b in branches if channels[b].size == size_label]
        if not pool:
            raise RuntimeError("No active {} channels found".format(size_label))
        print("Channel pool restricted to {}: {} channels".format(size_label, len(pool)))
    # Pinned channels ride along even when the size filter excludes them.
    pool_set = set(pool)
    ranking_pool = [b for b in pool if b not in pinned]

    # Rank clean events by summed positive channel integral over the pool, so a
    # 6 mm request picks events with 6 mm signal rather than Sci-dominated ones.
    energy_columns = ["{}_energy".format(branch) for branch in branches]
    scalar_data = manager.rdf.AsNumpy(["event_n"] + energy_columns)
    event_numbers = np.asarray(scalar_data["event_n"], dtype=np.int64)
    if event_numbers.size == 0:
        raise RuntimeError("No MCP-clean events in the requested range")
    energy_matrix = np.column_stack([
        np.nan_to_num(np.asarray(scalar_data[column], dtype=float), nan=0.0)
        for column in energy_columns
    ])
    branch_index = {branch: index for index, branch in enumerate(branches)}
    pool_columns = [branch_index[branch] for branch in pool]
    totals = np.clip(energy_matrix[:, pool_columns], 0.0, None).sum(axis=1)
    selected_indices = np.argsort(totals)[::-1][: args.n_events]
    selected_events = [int(event_numbers[index]) for index in selected_indices]
    print("Selected events: {}".format(selected_events))

    event_filter = " || ".join(
        "event_n == {}".format(event) for event in selected_events
    )
    display_rdf = manager.rdf.Filter(event_filter)
    if args.compare_windows:
        display_rdf = add_cfd_to_baseline_energy(
            display_rdf, branches,
            anchor_offset=args.anchor_offset,
            stop_fraction=args.stop_fraction, stop_abs=args.stop_abs,
            stop_consecutive=args.stop_consecutive,
            smooth_half=args.smooth_half, min_width=args.min_width,
        )
    waveform_columns: List[str] = ["event_n"]
    for branch in branches:
        waveform_columns.extend([
            "{}_blsub".format(branch),
            "{}_energy".format(branch),
            "{}_TS_cfd".format(branch),
        ])
        if args.compare_windows:
            waveform_columns.extend([
                "{}_energy_adaptive".format(branch),
                "{}_start_adaptive".format(branch),
                "{}_end_adaptive".format(branch),
                "{}_baseline_adaptive".format(branch),
                "{}_tscfd_adaptive".format(branch),
                get_ref_ts_name(branch),
            ])
    waveform_columns = list(dict.fromkeys(waveform_columns))
    pulse_data = display_rdf.AsNumpy(waveform_columns)

    png_names: List[str] = []
    metadata: Dict[str, object] = {
        "run": args.run,
        "selection": "MCP-clean",
        "integration_window": {"samples_before_anchor": WINDOW_PRE,
                               "samples_after_anchor": WINDOW_POST},
        "channel_pool": args.size,
        "pinned_channels": pinned,
        "events": {},
    }
    for event_number in selected_events:
        scalar_index = int(np.flatnonzero(event_numbers == event_number)[0])
        per_event_energy = energy_matrix[scalar_index]
        n_slots = max(args.channels_per_event, len(pinned))
        ranked = sorted(
            ranking_pool,
            key=lambda branch: per_event_energy[branch_index[branch]],
            reverse=True,
        )
        if args.all_channels:
            chosen = pinned + [b for b in ranked]
        else:
            chosen = pinned + ranked[: max(0, n_slots - len(pinned))]

        mask = event_mask(pulse_data["event_n"], event_number)
        row = int(np.flatnonzero(mask)[0])
        chunk = max(1, args.panels_per_figure) if args.all_channels else len(chosen)
        pages = [chosen[i:i + chunk] for i in range(0, len(chosen), chunk)] or [[]]
        event_records: List[Dict[str, object]] = []
        for page_no, page_channels in enumerate(pages, start=1):
            ncols = 2
            nrows = int(math.ceil(len(page_channels) / ncols))
            fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.6 * nrows),
                                     squeeze=False)
            for axis, branch in zip(axes.flat, page_channels):
                waveform = np.asarray(pulse_data["{}_blsub".format(branch)][row], dtype=float)
                stored_integral = float(pulse_data["{}_energy".format(branch)][row])
                cfd_time = float(pulse_data["{}_TS_cfd".format(branch)][row])
                anchor = int(math.floor(cfd_time)) + 1
                start = max(0, anchor - WINDOW_PRE)
                stop = min(waveform.size, anchor + WINDOW_POST)
                recomputed = float(waveform[start:stop].sum())
                samples = np.arange(waveform.size)

                extra_lo, extra_hi = start, stop
                ad_start = ad_stop = None
                fx_start = fx_stop = None
                if args.compare_windows:
                    ad_start = int(pulse_data["{}_start_adaptive".format(branch)][row])
                    ad_stop = int(pulse_data["{}_end_adaptive".format(branch)][row])
                    ref_ts = float(pulse_data[get_ref_ts_name(branch)][row])
                    fx_anchor = int(round(ref_ts)) + args.anchor_offset - 6
                    fx_start = max(0, fx_anchor - args.fixed_pre)
                    fx_stop = min(waveform.size, fx_anchor + args.fixed_post)
                    spans = [(ad_start, ad_stop)]
                    if not args.no_fixed_overlay:
                        spans.append((fx_start, fx_stop))
                    for lo_, hi_ in spans:
                        if lo_ is not None and hi_ is not None and hi_ > lo_:
                            extra_lo = min(extra_lo, lo_)
                            extra_hi = max(extra_hi, hi_)

                view_lo = max(0, extra_lo - 25)
                view_hi = min(waveform.size, extra_hi + 30)
                view = slice(view_lo, view_hi)
                axis.plot(samples[view], waveform[view], color="#225ea8", lw=1.7)
                integration_samples = samples[start:stop]
                axis.fill_between(
                    integration_samples,
                    0.0,
                    waveform[start:stop],
                    color="#fd8d3c",
                    alpha=0.40,
                    step="mid",
                    label="repo: CFD +/- fixed 50",
                )
                if args.compare_windows and ad_stop is not None and ad_stop > ad_start:
                    axis.axvspan(ad_start, ad_stop, color="#2ca25f", alpha=0.16,
                                 label="adaptive: CFD -> baseline crossing")
                    for edge in (ad_start, ad_stop):
                        axis.axvline(edge, color="#006d2c", lw=1.6)
                if (args.compare_windows and not args.no_fixed_overlay
                        and fx_stop is not None and fx_stop > fx_start):
                    axis.axvspan(fx_start, fx_stop, color="#54278f", alpha=0.10,
                                 label="fixed window @ reference")
                    for edge in (fx_start, fx_stop):
                        axis.axvline(edge, color="#54278f", lw=1.2, ls=":")
                axis.axvline(anchor, color="#cb181d", ls="--", lw=1.4, label="CFD anchor")
                axis.axhline(0.0, color="black", lw=0.7, alpha=0.6)
                title = panel_title(branch, channels[branch])
                if branch in pinned:
                    title += "  [pinned]"
                elif branch not in pool_set:
                    title += "  [outside pool]"
                axis.set_title(title, fontsize=10)
                axis.set_xlabel("Time sample")
                axis.set_ylabel("Baseline-subtracted ADC")
                axis.grid(alpha=0.18)
                note = "repo   = {:>8.1f} ADC   [{}, {})  w={}".format(
                    stored_integral, start, stop, stop - start)
                if args.compare_windows:
                    ad_e = float(pulse_data["{}_energy_adaptive".format(branch)][row])
                    fx_e = float(waveform[fx_start:fx_stop].sum()) if fx_stop > fx_start else 0.0
                    note += "\nadapt  = {:>8.1f} ADC   [{}, {})  w={}".format(
                        ad_e, ad_start, ad_stop, ad_stop - ad_start)
                    if not args.no_fixed_overlay:
                        note += "\nfixed  = {:>8.1f} ADC   [{}, {})  w={}".format(
                            fx_e, fx_start, fx_stop, fx_stop - fx_start)
                axis.text(
                    0.98, 0.95, note,
                    transform=axis.transAxes,
                    ha="right", va="top", fontsize=8, family="monospace",
                    bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.9},
                )
                record = {
                    "channel": branch,
                    "stored_integral_adc": stored_integral,
                    "recomputed_integral_adc": recomputed,
                    "cfd_time_sample": cfd_time,
                    "integration_start": start,
                    "integration_stop_exclusive": stop,
                }
                if args.compare_windows:
                    record.update({
                        "adaptive_integral_adc": float(
                            pulse_data["{}_energy_adaptive".format(branch)][row]),
                        "adaptive_start": ad_start,
                        "adaptive_stop_exclusive": ad_stop,
                        "adaptive_width": ad_stop - ad_start,
                        "fixed_ref_start": fx_start,
                        "fixed_ref_stop_exclusive": fx_stop,
                    })
                event_records.append(record)
            for axis in axes.flat[len(page_channels):]:
                axis.set_visible(False)
            handles, labels = axes.flat[0].get_legend_handles_labels()
            # Legend below the title, not on top of it.
            fig.legend(
                handles, labels, loc="upper center", ncol=4, frameon=False,
                bbox_to_anchor=(0.5, 0.972),
            )
            pool_text = ("all channels" if size_label is None
                         else "{} channels".format(size_label))
            scope = "all" if args.all_channels else "strongest"
            part = "  (part {}/{})".format(page_no, len(pages)) if len(pages) > 1 else ""
            fig.suptitle(
                "CaloX run {} event {} \u2014 {} {}{}".format(
                    args.run, event_number, scope, pool_text, part
                ),
                fontsize=15, y=0.998,
            )
            fig.tight_layout(rect=(0, 0, 1, 0.955))
            png_name = "adc_integral_run{}_event{}{}.png".format(
                args.run, event_number,
                "_part{:02d}".format(page_no) if len(pages) > 1 else "",
            )
            fig.savefig(outdir / png_name, dpi=150)
            plt.close(fig)
            png_names.append(png_name)
        metadata["events"][str(event_number)] = event_records

    with (outdir / "event_display_summary.json").open("w") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
    html_path = outdir / "adc_integral_event_displays.html"
    generate_html(
        png_names,
        str(outdir),
        plots_per_row=1,
        output_html=str(html_path),
        title="ADC-integral pulse event displays",
        intro_text=(
            "**Run:** {}. **Events:** {}.\n"
            "**Channel pool:** {}{}\n"
            "Orange shading is the exact 50-sample CFD energy-integration window: "
            "5 samples before and 45 after the 20% CFD anchor."
        ).format(
            args.run,
            ", ".join(str(event) for event in selected_events),
            "all sizes" if size_label is None else size_label,
            (". **Pinned:** " + ", ".join(
                panel_title(b, channels[b]) for b in pinned
            ) + ".") if pinned else ".",
        ),
    )
    print("Wrote {}".format(html_path))
    print("Wrote {}".format(outdir / "event_display_summary.json"))


if __name__ == "__main__":
    main()
