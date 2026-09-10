"""Build the interactive DRS channel-overlay page from existing histograms.

This is a post-processing step: it reads a ROOT file that a histogram-booking
script has already written, pulls out the per-channel profiles, and writes a
self-contained HTML page on which the curves can be switched on and off. It
runs no event loop and needs neither the raw data files nor the channel map,
so it is quick to re-run while iterating on the plot.

Usage
-----
  python scripts/make_drs_overlay.py --run 1994 \
      --channels data/channel_maps/testingfibers.json

  # a different histogram from the same file
  python scripts/make_drs_overlay.py --run 1994 \
      --channels data/channel_maps/testingfibers.json \
      --hist-pattern "prof_{channel}_blsub_VS_ts" --output DRS/overlay_raw.html

If the ROOT file or the histograms are missing, the run has not been produced
yet -- run scripts/check_drs_mcp.py for that run first.
"""
import argparse
import json
import os
import sys

import ROOT

from utils.overlay import build_channel_overlay
from utils.plot_helper import get_run_paths

ROOT.gROOT.SetBatch(True)

_RUN_FIRST = ("  Run it first, e.g.:\n"
              "    python3 scripts/check_drs_mcp.py --run {run}"
              " --channels <channels.json>")


def parse_args(argv=None):
    from configs.run_config import run_number as default_run
    p = argparse.ArgumentParser(
        description="Build the interactive DRS channel-overlay page from "
                    "histograms produced by an earlier run.")
    p.add_argument("--run", type=int, default=default_run, help="Run number")
    p.add_argument("--channels", default=None, metavar="FILE",
                   help="JSON file mapping labels to channel names. Without "
                        "it every channel found in the file is used, labelled "
                        "by its channel name.")
    p.add_argument("--input", default="drs_profiles.root", metavar="FILE",
                   help="ROOT file under results/root/Run<N>/ to read")
    p.add_argument("--hist-pattern", default="prof_{channel}_blsub_VS_ts_mcp",
                   metavar="PAT",
                   help="Histogram name, with {channel} standing for the "
                        "channel branch name")
    p.add_argument("--output", default="DRS/DRS_Prof_vs_ts_overlay.html",
                   metavar="PATH",
                   help="Output page, relative to results/html/Run<N>/")
    p.add_argument("--xlabel", default="ts")
    p.add_argument("--ylabel", default="Mean DRS Output")
    p.add_argument("--title", default=None)
    return p.parse_args(argv)


def _channel_labels(args, infile):
    """[(channel_name, label), ...] to look for, in a stable order."""
    if args.channels:
        if not os.path.exists(args.channels):
            print(f"ERROR: channel list not found: {args.channels}")
            return None
        with open(args.channels) as f:
            ch_map = json.load(f)
        return [(name, label) for label, name in ch_map.items()]

    # no channel list: take whatever the file holds, in file order
    prefix, suffix = args.hist_pattern.split("{channel}")
    found = []
    for key in infile.GetListOfKeys():
        name = key.GetName()
        if name.startswith(prefix) and name.endswith(suffix):
            channel = name[len(prefix):len(name) - len(suffix)]
            found.append((channel, channel))
    return found


def main(argv=None):
    args = parse_args(argv)
    paths = get_run_paths(args.run)
    root_path = os.path.join(paths["root"], args.input)

    if not os.path.exists(root_path):
        print(f"ERROR: {root_path} not found.")
        print(f"  Run {args.run} has no DRS histograms yet.")
        print(_RUN_FIRST.format(run=args.run))
        return 1

    infile = ROOT.TFile(root_path, "READ")
    if not infile or infile.IsZombie():
        print(f"ERROR: could not open {root_path} (corrupt or truncated).")
        print(_RUN_FIRST.format(run=args.run))
        return 1

    wanted = _channel_labels(args, infile)
    if wanted is None:
        return 1
    if not wanted:
        print(f"ERROR: no channels to plot from {root_path}.")
        print(f"  No histogram matches '{args.hist_pattern}'.")
        print(_RUN_FIRST.format(run=args.run))
        return 1

    entries, missing = [], []
    for channel, label in wanted:
        hist = infile.Get(args.hist_pattern.format(channel=channel))
        if not hist:
            missing.append(label)
            continue
        hist.SetDirectory(0)
        entries.append((hist, label))

    if missing:
        print(f"Warning: {len(missing)} of {len(wanted)} channels have no "
              f"'{args.hist_pattern}' histogram in {args.input}:")
        print("   " + ", ".join(missing[:10])
              + (" ..." if len(missing) > 10 else ""))

    if not entries:
        print(f"ERROR: none of the {len(wanted)} requested channels are in "
              f"{root_path}.")
        print("  The file exists but holds no matching histograms -- it was "
              "probably written for a different channel selection.")
        print(_RUN_FIRST.format(run=args.run))
        return 1

    title = args.title or f"DRS channel overlay — Run {args.run}"
    out = build_channel_overlay(
        entries,
        output_html=os.path.join(paths["html"], args.output),
        xlabel=args.xlabel, ylabel=args.ylabel, title=title,
        intro_text="Tick the channels to superimpose, or click the legend. "
                   f"Curves are {args.hist_pattern.format(channel='<channel>')} "
                   f"from {args.input}.",
        filename=f"DRS_overlay_Run{args.run}")
    infile.Close()
    if not out:
        return 1
    print(f"\n{len(entries)} channels plotted -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
