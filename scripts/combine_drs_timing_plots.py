"""
DRS Timing Plots Script.

Analyzes DRS peak timing relative to MCP for timing calibration.
"""

import json
import ROOT
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle, STYLE_2D_COLZ
from plotting.calox_plot_helper import create_pave_text
from plotting.my_function import LHistos2Hist
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import number_to_string
from variables.drs import calibrateDRSPeakTS

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()
            .apply_hole_veto(flag_only=True)
            )

fersboards = analysis.fersboards
DRSBoards = analysis.drsboards

benergy = analysis.beam_energy
run_number = analysis.run_number
paths = analysis.paths
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]

file_drschannels_bad = "data/drs/badchannels.json"
with open(file_drschannels_bad, "r") as f:
    drschannels_bad = json.load(f)

rdf = analysis.get_particle_analysis("pion")

# Configuration constants
TSmin = -90
TSmax = -10
TSCermin = -70
TSCermax = -50
varsuffix = "DiffRelPeakTS_US"
doDetailedPeakTSPlots = False
doDetailedPeakTS2DPlots = False

# Common plot styles
STYLE_CER_SCI = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 4],
    addOverflow=False,
    addUnderflow=False
)

STYLE_CER_QUARTZ_PLASTIC_SCI = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 6, 4],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90]
)

STYLE_CER_QUARTZ_PLASTIC_SCI_NORM = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 6, 4],
    addOverflow=False,
    addUnderflow=False,
    donormalize=True,
    legendPos=[0.25, 0.75, 0.40, 0.90]
)

STYLE_CER_QUARTZ_PLASTIC = PlotStyle(
    dology=False,
    drawoptions=["hist,C", "hist,C", "hist,C"],
    mycolors=[2, 6],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.30, 0.80, 0.40, 0.90]
)

STYLE_2D_LOG = PlotStyle(
    dology=False,
    dologz=True,
    drawoptions="COLZ",
    addOverflow=False,
    addUnderflow=False,
    zmin=1,
    zmax=1e2
)


def GetMean(hist, useMode=True):
    """Get mean or mode of histogram."""
    if useMode:
        bin_max = hist.GetMaximumBin()
        mean = hist.GetBinCenter(bin_max)
    else:
        mean = hist.GetMean()
    return mean


def makeDRSPeakTSPlots(pm: PlotManager, runs):
    """Generate DRS peak TS plots using PlotManager."""
    pm.reset_plots().set_output_dir("DRSPeakTS_rel_us")

    hists_Cer = []
    hists_Cer_Quartz = []
    hists_Cer_Plastic = []
    hists_Sci = []

    for run in runs:
        infile_name = f"results/root/Run{run}/drspeakts_rel_us_combined.root"
        infile = ROOT.TFile(infile_name, "READ")

        hist_Cer = infile.Get("hist_DRSPeakTS_Cer_Combined")
        hist_Sci = infile.Get("hist_DRSPeakTS_Sci_Combined")
        hist_Cer_Quartz = infile.Get("hist_DRSPeakTS_Cer_Quartz_Combined")
        hist_Cer_Plastic = infile.Get("hist_DRSPeakTS_Cer_Plastic_Combined")

        hists_Cer.append(hist_Cer)
        hists_Sci.append(hist_Sci)
        hists_Cer_Quartz.append(hist_Cer_Quartz)
        hists_Cer_Plastic.append(hist_Cer_Plastic)

    pm.plot_1d(
        hists_Cer_Quartz,
        "DRS_PeakTS_Cer_Quartz",
        "Quartz Peak TS",
        (-75, -55),
        ylabel="Counts",
        yrange=(1, None),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI,
        prepend=False
    )

    pm.plot_1d(
        hists_Cer_Plastic,
        "DRS_PeakTS_Cer_Plastic",
        "Plastic Peak TS",
        (-75, -55),
        ylabel="Counts",
        yrange=(1, None),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI,
        prepend=False
    )

    pm.plot_1d(
        hists_Sci,
        "DRS_PeakTS_Sci",
        "Sci Peak TS",
        (-65, -30),
        ylabel="Counts",
        yrange=(1, None),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI,
        prepend=False
    )

    pm.add_newline()

    output_html = pm.generate_html(
        "DRS/DRSPeakTS_relative_to_MCP.html",
        plots_per_row=4,
        title="DRS Peak Position",
        intro_text="Event-by-event DRS peak (channel-wise shower max) positions with respect to MCP for different towers.\n For now only the central DRS boards are included."
    )

    return output_html


def makeDRSvsTSProfPlots(pm: PlotManager, runs):
    """Generate DRS peak TS plots using PlotManager."""
    pm.reset_plots().set_output_dir("DRS_VS_TS_rel_us")

    hists_Cer = []
    hists_Cer_Quartz = []
    hists_Cer_Plastic = []
    hists_Sci = []

    for run in runs:
        infile_name = f"results/root/Run{run}/drs_vs_ts_calibrated_combined.root"
        infile = ROOT.TFile(infile_name, "READ")

        hist_Cer = infile.Get("prof_DRS_vs_TS_Cer_Combined")
        hist_Sci = infile.Get("prof_DRS_vs_TS_Sci_Combined")
        hist_Cer_Quartz = infile.Get("prof_DRS_vs_TS_Cer_Quartz_Combined")
        hist_Cer_Plastic = infile.Get("prof_DRS_vs_TS_Cer_Plastic_Combined")

        hists_Cer.append(hist_Cer)
        hists_Sci.append(hist_Sci)
        hists_Cer_Quartz.append(hist_Cer_Quartz)
        hists_Cer_Plastic.append(hist_Cer_Plastic)

    pm.plot_1d(
        hists_Cer_Quartz,
        "DRS_VS_TS_Cer_Quartz",
        "Quartz TS",
        (-75, -55),
        ylabel="frac of energy",
        yrange=(0, 0.1),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI_NORM,
        prepend=False,
        includeRunNumber=False
    )

    pm.plot_1d(
        hists_Cer_Plastic,
        "DRS_VS_TS_Cer_Plastic",
        "Plastic TS",
        (-75, -55),
        ylabel="frac of energy",
        yrange=(0, 0.1),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI_NORM,
        prepend=False,
        includeRunNumber=False
    )

    pm.plot_1d(
        hists_Sci,
        "DRS_VS_TS_Sci",
        "Sci TS",
        (-65, -30),
        ylabel="frac of energy",
        yrange=(0, 0.05),
        legends=runs,
        style=STYLE_CER_QUARTZ_PLASTIC_SCI_NORM,
        prepend=False,
        includeRunNumber=False
    )

    pm.add_newline()

    output_html = pm.generate_html(
        "DRS/DRS_VS_TS_relative_to_MCP.html",
        plots_per_row=4,
        title="Profiled DRS ADC vs TS",
        intro_text="Average DRS ADC vs calibrated time slice (TS) for Cer and Sci channels. For each event, the time slice is measured relative to the MCP signal."
    )

    return output_html


def main():
    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        makeDRSPeakTSPlots(pm, [1501, 1507, 1511])
        makeDRSvsTSProfPlots(pm, [1501, 1507, 1511])

        PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
