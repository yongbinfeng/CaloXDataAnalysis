#!/usr/bin/env python3
"""Interactive 3D "stacked heatmap" shower display for a single event.

Same rendering as plot_adc_scan3d.py / plot_toa_scan3d.py (utils/scan3d.py),
but for one event instead of a whole scan: every DRS channel becomes one
flat, ADC-coloured tile at its physical (x, y) tower position, and the tiles
are stacked along a TOA axis instead of the table-position axis - so the
length dimension is now "when" instead of "where", and colour is ADC instead
of TOA. TOA is binned (--time-bin-ns) so channels that fired together form
one flat plane, the same way each scan run formed one plane before.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List

import numpy as np

from core.analysis_manager import CaloXAnalysisManager
from utils.root_setup import setup_root
from utils.scan3d import build_channel_geometry, compute_toa_references, write_scan3d_html

from plot_velocityfits import active_timing_channels, make_manager_args

DEFAULT_OUTDIR = "results/ShoweShape"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Single-event 3D shower display: calorimeter face x, y, "
            "colour = ADC, stacked along TOA."
        )
    )
    parser.add_argument(
        "--run", type=int, default=1831,
        help="tb2026 run (>= 1828 era). Default: 1831 (120 GeV electron).",
    )
    parser.add_argument(
        "--particle", default=None,
        help=(
            "Override the beam-type PID selection (electron/pion/proton). "
            "Default: whatever the run's beam type resolves to."
        ),
    )
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=20000)
    parser.add_argument(
        "--event", type=int, default=None,
        help=(
            "Specific event_n to display. Default: the highest total-ADC "
            "event among clean, PID-selected events in the range."
        ),
    )
    parser.add_argument(
        "--time-bin-ns", type=float, default=0.2,
        help="TOA bin width, in ns, used to stack channels into planes.",
    )
    parser.add_argument(
        "--toa-window-ns", type=float, default=10.0,
        help=(
            "Drop channels whose relative TOA falls outside +/- this many ns. "
            "The real shower front lands almost all channels within a few ns "
            "of each other; a handful of badly-timed CFD fits otherwise land "
            "tens of ns away and stretch the z axis until the real cluster "
            "is squeezed into what looks like a single slice. Default: 10."
        ),
    )
    parser.add_argument(
        "--min-adc", type=float, default=15.0,
        help=(
            "Drop channels with less ADC than this for the chosen event, and "
            "use only same-run events at or above it to estimate each "
            "channel's own reference TOA (see --min-ref-hits). Channels with "
            "no real pulse have a CFD crossing set by noise. Default: 15."
        ),
    )
    parser.add_argument(
        "--min-ref-hits", type=int, default=5,
        help=(
            "Minimum same-run, above-threshold hits required to trust a "
            "channel's reference TOA; channels with fewer are dropped. "
            "Default: 5."
        ),
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

    if args.event is None:
        energy_matrix = np.column_stack([
            np.nan_to_num(np.asarray(data[col], dtype=float), nan=0.0)
            for col in energy_columns
        ])
        totals = np.clip(energy_matrix, 0.0, None).sum(axis=1)
        event_index = int(np.argmax(totals))
        event_number = int(event_numbers[event_index])
        print("Auto-selected event {} (highest total ADC: {:.0f})".format(
            event_number, totals[event_index]))
    else:
        matches = np.flatnonzero(event_numbers == args.event)
        if matches.size == 0:
            raise RuntimeError(
                "Event {} not found among selected events for run {}".format(
                    args.event, args.run)
            )
        event_index = int(matches[0])
        event_number = args.event

    # _Time_cfd_mcp's absolute scale is dominated by per-channel cable-length
    # offsets (the existing scan scripts see ~8.5 ns of that spread and call
    # it out explicitly), not real shower timing. There is no repeated scan
    # of this run to average over, so instead each channel's own "typical"
    # arrival time is estimated from every other clean event in this run
    # where that channel actually had a real pulse, and only the deviation
    # from that is kept - the single-event analogue of --relative.
    adc_all: Dict[str, float] = {}
    toa_all: Dict[str, float] = {}
    for branch, energy_col, toa_col in zip(branches, energy_columns, toa_columns):
        e = float(data[energy_col][event_index])
        t = float(data[toa_col][event_index])
        if math.isfinite(e) and math.isfinite(t):
            adc_all[branch] = e
            toa_all[branch] = t
    if not adc_all:
        raise RuntimeError("No channel had both a finite ADC and TOA for this event")

    all_adc_values = sorted(adc_all.values())
    print("ADC range this event: {:.1f} .. {:.1f} ({} channels); "
          "dropping < {:.1f}".format(
              all_adc_values[0], all_adc_values[-1], len(all_adc_values),
              args.min_adc))

    adc = {b: v for b, v in adc_all.items() if v >= args.min_adc}
    if not adc:
        raise RuntimeError(
            "No channel survived the --min-adc={} cut for this event".format(
                args.min_adc)
        )

    references = compute_toa_references(data, list(adc), args.min_adc, args.min_ref_hits)
    n_no_reference = len(adc) - len(references)
    toa = {b: toa_all[b] - ref for b, ref in references.items()}
    adc = {b: v for b, v in adc.items() if b in toa}
    if not adc:
        raise RuntimeError(
            "No channel had a stable per-channel TOA reference in this run"
        )
    if n_no_reference > 0:
        print("Dropped {} channel(s) with fewer than {} same-run reference "
              "hits (no stable TOA baseline)".format(n_no_reference, args.min_ref_hits))

    n_before_window = len(toa)
    toa = {b: t for b, t in toa.items() if abs(t) <= args.toa_window_ns}
    adc = {b: v for b, v in adc.items() if b in toa}
    n_outliers = n_before_window - len(toa)
    if not adc:
        raise RuntimeError(
            "No channel survived the +/-{} ns TOA window".format(args.toa_window_ns)
        )
    if n_outliers:
        print("Dropped {} channel(s) outside +/-{:.1f} ns as mistimed outliers "
              "(bad CFD fit, not real shower structure)".format(
                  n_outliers, args.toa_window_ns))
    print("Kept {}/{} channels above the ADC threshold; TOA (relative to each "
          "channel's own run average) now spans {:.2f} .. {:.2f} ns".format(
              len(adc), len(adc_all), min(toa.values()), max(toa.values())))

    bin_width = args.time_bin_ns
    bin_of: Dict[str, int] = {b: round(t / bin_width) for b, t in toa.items()}
    bins: List[int] = sorted(set(bin_of.values()))
    bin_positions = {b: b * bin_width for b in bins}

    geometry = build_channel_geometry(map_run=args.map_run or args.run)
    values_by_category: Dict[str, Dict[int, Dict[str, float]]] = {}
    for category, chlist in geometry.items():
        branch_set = {c["branch"] for c in chlist}
        per_bin: Dict[int, Dict[str, float]] = {b: {} for b in bins}
        for branch, value in adc.items():
            if branch not in branch_set:
                continue
            per_bin[bin_of[branch]][branch] = value
        values_by_category[category] = per_bin

    n_channels = len(adc)
    output_name = args.output_name or "event_shower_run{}_evt{}.html".format(
        args.run, event_number)
    output_path = Path(args.outdir) / output_name
    write_scan3d_html(
        output_path=output_path,
        page_title="Event shower: ADC vs. face position vs. TOA",
        subtitle=(
            "run {} &middot; event {} &middot; {} channels &middot; {} ns bins "
            "&middot; TOA relative to each channel's own run average"
        ).format(args.run, event_number, n_channels, bin_width),
        geometry=geometry,
        runs=bins,
        positions=bin_positions,
        values_by_category=values_by_category,
        value_label="ADC [counts]",
        value_fmt="{:.1f}",
        z_axis_title="TOA - channel run average [ns]",
        colorscale="Turbo",
    )
    print("Wrote {} ({} channels, {} time bins)".format(
        output_path, n_channels, len(bins)))


if __name__ == "__main__":
    main()
