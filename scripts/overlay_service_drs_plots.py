"""
Overlay service DRS pulse plots across multiple runs.

Reads pre-computed histograms saved by check_service_drs.py
(drs_services_pid.root) and overlays the 1D pulse summary plots.

Usage:
  python scripts/overlay_service_drs_plots.py --runs 1465 1466 1467
  python scripts/overlay_service_drs_plots.py --runs 1465 1466 --labels "20 GeV" "40 GeV"
"""

import argparse
import sys
import ROOT
from utils.data_loader import getRunInfo
from utils.plot_helper import get_run_paths
from utils.root_setup import setup_root
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from configs.plot_config import getServiceDRSProcessedInfoRanges
from configs.selection_config import get_service_drs_cut
from plotting.calox_plot_helper import create_pave_text

setup_root(batch_mode=True)

ALL_DETS = ["HoleVeto", "PSD", "TTUMuonVeto", "Cer474", "Cer519", "Cer537", "KT1", "KT2", "T3", "T4"]

# One color per run; same color used for all plots belonging to the same run
RUN_COLORS = [
    ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2,
    ROOT.kMagenta+1, ROOT.kOrange+1, ROOT.kCyan+2,
    ROOT.kViolet+1, ROOT.kSpring+5,
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Overlay service DRS pulse summary plots across multiple runs"
    )
    parser.add_argument("--runs", nargs="+", type=int, required=True,
                        help="Run numbers to overlay")
    parser.add_argument("--dets", nargs="+", type=str, default=ALL_DETS,
                        help=f"Detectors to include (default: all)")
    parser.add_argument("--labels", nargs="+", type=str,
                        help="Labels for each run (default: Run{N} {type} {energy}GeV)")
    parser.add_argument("--output-run", type=int,
                        help="Run to use for output paths (default: first run)")
    return parser.parse_args()


def build_labels(runs, user_labels):
    if user_labels:
        if len(user_labels) != len(runs):
            print("Error: --labels count must match --runs count")
            sys.exit(1)
        return [l.replace("_", " ") for l in user_labels]
    labels = []
    for r in runs:
        try:
            btype, benergy = getRunInfo(r)
            labels.append(f"Run{r} {btype} {benergy}GeV".replace("_", " "))
        except Exception:
            labels.append(f"Run{r}")
    return labels


def load_hists(run_number, dets):
    """Load 1D pulse histograms from a run's drs_services_pid.root."""
    rootdir = get_run_paths(run_number)["root"]
    infile_name = f"{rootdir}/drs_services_pid.root"
    infile = ROOT.TFile(infile_name, "READ")
    if infile.IsZombie():
        print(f"Error: Could not open {infile_name}")
        return None

    hists = {}
    for det in dets:
        hists[det] = {}
        for metric in ["peak_position", "peak_value", "sum"]:
            h = infile.Get(f"{det}_{metric}")
            if h:
                h.SetDirectory(0)
                hists[det][metric] = h
            else:
                hists[det][metric] = None

    infile.Close()
    return hists


def make_style(colors, dology=False):
    n = len(colors)
    return PlotStyle(
        dology=dology,
        drawoptions=["HIST"] * n,
        mycolors=colors,
        addOverflow=True,
        addUnderflow=True,
        donormalize=True,
        legendPos=[0.20, 0.65, 0.55, 0.90],
    )


def main():
    args = parse_args()
    runs = args.runs
    dets = args.dets
    labels = build_labels(runs, args.labels)
    output_run = args.output_run or runs[0]

    # Load histograms for each run
    loaded = []
    for run, label in zip(runs, labels):
        hists = load_hists(run, dets)
        if hists is None:
            print(f"Skipping Run{run}: file not found")
        else:
            loaded.append((label, hists))

    if not loaded:
        print("No valid runs to plot. Run check_service_drs.py first.")
        sys.exit(1)

    # Assign one color per run; reused across all per-detector plots for that run
    run_colors = {lbl: RUN_COLORS[i % len(RUN_COLORS)] for i, (lbl, _) in enumerate(loaded)}

    paths = get_run_paths(output_run)
    pm = PlotManager(paths["root"], paths["plots"], paths["html"], output_run)
    pm.set_output_dir("drs_service_overlay")

    def gather(det, metric):
        """Return (labels, hists, colors) for entries that have this histogram."""
        entries = [(lbl, h[det][metric]) for lbl, h in loaded if h[det].get(metric)]
        if not entries:
            return [], [], []
        lbls, hs = zip(*entries)
        colors = [run_colors[lbl] for lbl in lbls]
        return list(lbls), list(hs), colors

    for det in dets:
        # Peak position
        lbls, hs, colors = gather(det, "peak_position")
        if hs:
            pm.plot_1d(
                hs, f"{det}_peak_position_overlay",
                f"{det} Peak Position [TS]", (0, 1024),
                ylabel="Frac. of Events", yrange=(0, 10),
                legends=lbls, style=make_style(colors),
                legendoptions=["L"] * len(colors),
                includeRunNumber=False,
            )

        # Peak value
        lbls, hs, colors = gather(det, "peak_value")
        if hs:
            try:
                xmin, xmax = getServiceDRSProcessedInfoRanges(det, "peak_value")
            except Exception:
                xmin, xmax = hs[0].GetXaxis().GetXmin(), hs[0].GetXaxis().GetXmax()
            pm.plot_1d(
                hs, f"{det}_peak_value_overlay",
                f"{det} Peak Value [ADC]", (xmin, xmax),
                ylabel="Frac. of Events", yrange=(0, 10),
                legends=lbls, style=make_style(colors),
                legendoptions=["L"] * len(colors),
                includeRunNumber=False,
            )

        # Sum — log y, normalized fractions are O(1e-3), so ymin = 1e-4
        lbls, hs, colors = gather(det, "sum")
        if hs:
            try:
                xmin, xmax = getServiceDRSProcessedInfoRanges(det, "sum")
            except Exception:
                xmin, xmax = hs[0].GetXaxis().GetXmin(), hs[0].GetXaxis().GetXmax()

            _, _, value_cut, _ = get_service_drs_cut(det)
            cut_line = ROOT.TLine(value_cut, 1e-4, value_cut, 10)
            cut_line.SetLineColor(ROOT.kRed)
            cut_line.SetLineWidth(2)
            cut_line.SetLineStyle(ROOT.kDashed)

            pm.plot_1d(
                hs, f"{det}_sum_overlay",
                f"{det} Sum [ADC]", (xmin, xmax),
                ylabel="Frac. of Events", yrange=(1e-4, 10),
                legends=lbls, style=make_style(colors, dology=True),
                legendoptions=["L"] * len(colors),
                includeRunNumber=False,
                extraToDraw=cut_line,
            )

        pm.add_newline()

    output_html = pm.generate_html(
        "ServiceDRS/overlay.html",
        plots_per_row=3,
        title="Service DRS Overlay",
        intro_text=(
            "Per-detector overlay of peak position, peak value, and sum across runs. "
            "Same color per run across all plots; different color per run."
        ),
    )
    print(f"Output HTML: {output_html}")


if __name__ == "__main__":
    main()
