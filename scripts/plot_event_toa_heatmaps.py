#!/usr/bin/env python3
"""2D per-event calorimeter-face heatmaps, colour = TOA, one event per slide.

Same channel geometry and TOA-relative-to-its-own-run-average trick as
plot_event_shower_scan3d.py, but flattened to a plain 2D face map (no TOA
stacking) and built for many events at once: a slider pages from one event's
heatmap to the next, so a whole batch of events can be browsed in one page.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from core.analysis_manager import CaloXAnalysisManager
from utils.root_setup import setup_root
from utils.scan3d import build_channel_geometry, compute_toa_references, write_event_heatmap_gallery_html

from plot_velocityfits import active_timing_channels, make_manager_args

DEFAULT_OUTDIR = "results/ShoweShape"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "2D calorimeter-face heatmaps coloured by TOA, one event per "
            "slide, for a batch of events in one run."
        )
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
        help="Number of highest-total-ADC events to show. Default: 20.",
    )
    parser.add_argument(
        "--min-adc", type=float, default=15.0,
        help="Drop a channel from an event's slide if its ADC is below this.",
    )
    parser.add_argument(
        "--min-ref-hits", type=int, default=5,
        help="Minimum same-run, above-threshold hits to trust a channel's "
             "reference TOA; channels with fewer are dropped entirely.",
    )
    parser.add_argument(
        "--toa-window-ns", type=float, default=10.0,
        help="Drop a channel from a slide if |TOA - its run average| exceeds "
             "this many ns (mistimed CFD fits, not real shower structure).",
    )
    parser.add_argument(
        "--map-run", type=int, default=None,
        help="Run whose channel map supplies tower positions. Default: --run.",
    )
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    parser.add_argument("--output-name", default=None)
    parser.add_argument("--n-threads", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)

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

    if args.events is not None:
        requested = set(args.events)
        indices = [i for i, n in enumerate(event_numbers) if int(n) in requested]
        if not indices:
            raise RuntimeError("None of --events {} found for run {}".format(args.events, args.run))
    else:
        energy_matrix = np.column_stack([
            np.nan_to_num(np.asarray(data[col], dtype=float), nan=0.0)
            for col in energy_columns
        ])
        totals = np.clip(energy_matrix, 0.0, None).sum(axis=1)
        n_pick = min(args.n_events, totals.size)
        indices = list(np.argsort(totals)[::-1][:n_pick])
        print("Auto-selected top {} events by total ADC (run {}).".format(n_pick, args.run))

    # Same trick as plot_event_shower_scan3d.py: _Time_cfd_mcp's absolute
    # scale is dominated by per-channel cable-length offsets, so each
    # channel's own typical arrival time (over this run's other clean,
    # above-threshold events) is subtracted before using TOA as colour.
    references = compute_toa_references(data, branches, args.min_adc, args.min_ref_hits)
    print("{}/{} channels have a stable same-run TOA reference (>= {} hits "
          "above {:.1f} ADC).".format(len(references), len(branches),
                                      args.min_ref_hits, args.min_adc))

    geometry = build_channel_geometry(map_run=args.map_run or args.run)
    values_by_category: Dict[str, Dict[int, Dict[str, Tuple[float, float]]]] = {
        cat: {} for cat in geometry
    }
    branch_category = {
        branch: cat
        for cat, chlist in geometry.items()
        for branch in (c["branch"] for c in chlist)
    }

    events: List[int] = []
    for index in indices:
        event_number = int(event_numbers[index])
        per_event: Dict[str, Tuple[float, float]] = {}
        for branch, reference in references.items():
            e = float(data["{}_energy".format(branch)][index])
            t = float(data["{}_Time_cfd_mcp".format(branch)][index])
            if not (math.isfinite(e) and math.isfinite(t)):
                continue
            if e < args.min_adc:
                continue
            rel_t = t - reference
            if abs(rel_t) > args.toa_window_ns:
                continue
            per_event[branch] = (rel_t, e)
        if not per_event:
            print("  event {}: no channel survived the cuts, skipping".format(event_number))
            continue
        events.append(event_number)
        for branch, value in per_event.items():
            cat = branch_category.get(branch)
            if cat is None:
                continue
            values_by_category[cat].setdefault(event_number, {})[branch] = value
        print("  event {}: {} channels kept".format(event_number, len(per_event)))

    if not events:
        raise RuntimeError("No event survived the cuts")
    events.sort()

    output_name = args.output_name or "event_toa_heatmaps_run{}.html".format(args.run)
    output_path = Path(args.outdir) / output_name
    write_event_heatmap_gallery_html(
        output_path=output_path,
        page_title="Per-event TOA heatmaps",
        subtitle=(
            "run {} &middot; {} events &middot; TOA relative to each "
            "channel's own run average"
        ).format(args.run, len(events)),
        geometry=geometry,
        events=events,
        values_by_category=values_by_category,
        value_label="TOA - channel run average [ns]",
        value_fmt="{:.2f}",
        colorscale="Turbo",
    )
    print("Wrote {} ({} events)".format(output_path, len(events)))


if __name__ == "__main__":
    main()
