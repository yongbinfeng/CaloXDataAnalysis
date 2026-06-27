"""Overlay per-channel DRS waveform profiles (ts mcp) between two runs, with ratio.

For every DRS channel, reads the ref+MCP-corrected mean-waveform profile
  prof_{channel}_blsub_VS_ts_mcp
from each run's drs_profiles.root (produced by make_drs_dqm_hists.py /
check_drs_mcp.py) and overlays the runs on one canvas with a ratio panel
(run[i] / run[0]) below each plot — i.e. the "ts mcp" trace from the
DRS_Prof_vs_ts page, compared across runs.

Default: run 1841 vs 1842.

Usage
-----
  python scripts/compare/compare_drs_prof_ts_mcp.py
  python scripts/compare/compare_drs_prof_ts_mcp.py --runs 1841 1842
  python scripts/compare/compare_drs_prof_ts_mcp.py --runs 1841 1842 --labels noFlip Flip
"""

import argparse
import os
import sys
import ROOT

from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.utils import get_channel_var
from channels.channel_map import build_drs_boards
from configs.plot_config import get_drs_prof_plot_ranges
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle

setup_root(n_threads=1, batch_mode=True, load_functions=False)

_COLORS = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kOrange + 1]


def _load_prof(run_number, ch_blsub):
    """Load prof_{ch_blsub}_VS_ts_mcp as a TH1D of the per-bin means (detached).

    The stored object is a TProfile; we ProjectionX() to a plain TH1D so the
    ratio panel divides the means correctly (TProfile::Divide does not).
    """
    path = f"results/root/Run{run_number}/drs_profiles.root"
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        print(f"Warning: could not open {path}")
        return None
    prof = f.Get(f"prof_{ch_blsub}_VS_ts_mcp")
    h = None
    if prof:
        h = prof.ProjectionX(f"{ch_blsub}_Run{run_number}_tsmcp")
        h.SetDirectory(0)
    f.Close()
    return h


def run_comparison(runs, labels, use_jsroot=False):
    out_tag = "_vs_".join(f"Run{r}" for r in runs)
    os.makedirs("results/plots/Compare", exist_ok=True)
    os.makedirs("results/html/Compare", exist_ok=True)

    pm = PlotManager(
        f"results/root/Run{runs[0]}",
        "results/plots/Compare",
        "results/html/Compare",
        run_number=runs[0],
        use_jsroot=use_jsroot,
    )
    pm.set_output_dir(out_tag + "_Prof_ts_mcp")

    colors = _COLORS[:len(runs)]
    style = PlotStyle(dology=False, drawoptions="HIST", mycolors=list(colors),
                      addOverflow=False, addUnderflow=False, legendoptions="L")
    yrlabel = f"/ {labels[0]}"

    drsboards = build_drs_boards(runs[0])

    # Pass 1: load + project every channel (no drawing, so ProjectionX is fast).
    entries = []
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            ch_blsub = chan.get_channel_name(blsub=True)
            present = [(h, l) for h, l in
                       ((_load_prof(r, ch_blsub), l) for r, l in zip(runs, labels))
                       if h]
            if len(present) < 2:
                continue
            entries.append((ch_blsub, get_channel_var(chan), chan,
                            [h for h, _ in present], [l for _, l in present]))

    # Pass 2: draw overlays + ratio panels.
    for ch_blsub, var, chan, hs, ls in entries:
        ymin, ymax = get_drs_prof_plot_ranges(
            subtractMedian=True, is_amplified=chan.is_amplified,
            is6mm=chan.is6mm, is_reference=False, is_cer=chan.isCer,
            run_number=runs[0])
        pm.plot_1d(
            hs, f"DRS_Prof_ts_mcp_{ch_blsub}_{var}",
            "ts (MCP-corrected)", (0, 1024),
            "Mean DRS Output", (ymin, ymax),
            legends=ls, style=style,
            extra_text=var, includeRunNumber=False,
            showratio=True, ratiobase=0,
            yrlabel=yrlabel, yrmin=0.0, yrmax=2.0)

    print(f"Plotted {len(entries)} channels")
    output_html = pm.generate_html(f"{out_tag}_Prof_ts_mcp.html", plots_per_row=6)
    print(f"Output: {output_html}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_labels(runs, user_labels):
    if user_labels:
        if len(user_labels) != len(runs):
            print("Error: --labels count must match --runs count")
            sys.exit(1)
        return [l.replace("_", " ") for l in user_labels]
    return [f"Run{r}" for r in runs]


def _parse_args():
    p = argparse.ArgumentParser(
        description="Overlay per-channel DRS ts-mcp profiles between runs (+ratio)")
    p.add_argument("--runs", nargs="+", type=int, default=[1841, 1842],
                   help="Run numbers to compare (default: 1841 1842)")
    p.add_argument("--labels", nargs="+", type=str,
                   help="Legend labels per run (use _ for spaces)")
    p.add_argument("--no-jsroot", dest="jsroot", action="store_false",
                   help="Save PNG images instead of interactive JSROOT plots")
    p.set_defaults(jsroot=True)
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    if len(args.runs) < 2:
        print("Error: need at least two runs to compare")
        sys.exit(1)
    labels = _build_labels(args.runs, args.labels)
    run_comparison(args.runs, labels, use_jsroot=args.jsroot)
