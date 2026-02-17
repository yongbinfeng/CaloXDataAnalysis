"""
DRS Timing Plots Script.

Analyzes DRS peak timing relative to MCP for timing calibration.
"""

import json
from pdb import pm
import ROOT
from collections import OrderedDict
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

STYLE_CER_QUARTZ_PLASTIC = PlotStyle(
    dology=False,
    drawoptions=["hist,C", "hist,C"],
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


def checkDRSPeakTS(rdf, suffix=""):
    """Book DRS peak timing histograms."""
    h1s_DRSPeakTS = {"Cer": [], "Sci": []}
    h2s_DRSPeakTS_Cer_VS_Sci = []

    channelnames_quartz = []
    channelnames_plastic = []
    channelnames_sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.get_channel_name(blsub=False)
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                h1_DRSPeakTS = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}_{suffix}",
                    f"DRS Peak TS for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
                    TSmax - TSmin, TSmin, TSmax),
                    f"{channelName}_{varsuffix}"
                )
                h1s_DRSPeakTS[var].append(h1_DRSPeakTS)

                if var == "Cer":
                    if chan_DRS.isQuartz:
                        channelnames_quartz.append(
                            f"{channelName}_{varsuffix}")
                    else:
                        channelnames_plastic.append(
                            f"{channelName}_{varsuffix}")
                else:
                    channelnames_sci.append(f"{channelName}_{varsuffix}")

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough good channels found for Board{board_no}, Tower({sTowerX}, {sTowerY})")
                continue

            h2_DRSPeak_Cer_VS_Sci = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}_{suffix}",
                f"DRS Peak TS - CER VS SCI for Board{board_no}, Tower({sTowerX}, {sTowerY}), {suffix};SCI Peak TS;CER Peak TS",
                TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
                f'{channelNames["Sci"]}_{varsuffix}',
                f'{channelNames["Cer"]}_{varsuffix}',
            )
            h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_VS_Sci)

    # Define average TS variables
    rdf = rdf.Define("Cer_Quartz_AvgPeakTS",
                     f"({'+'.join(channelnames_quartz)})/{len(channelnames_quartz)}")
    rdf = rdf.Define("Cer_Plastic_AvgPeakTS",
                     f"({'+'.join(channelnames_plastic)})/{len(channelnames_plastic)}")
    rdf = rdf.Define("Sci_AvgPeakTS",
                     f"({'+'.join(channelnames_sci)})/{len(channelnames_sci)}")

    # Book average histograms
    h1_Cer_Quartz_AvgPeakTS = rdf.Histo1D((
        f"hist_DRSPeakTS_Cer_Quartz_Avg_{suffix}",
        f"DRS Peak TS Avg for Cer Quartz, {suffix};Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Cer_Quartz_AvgPeakTS"
    )
    h1s_DRSPeakTS["Cer"].append(h1_Cer_Quartz_AvgPeakTS)

    h1_Cer_Plastic_AvgPeakTS = rdf.Histo1D((
        f"hist_DRSPeakTS_Cer_Plastic_Avg_{suffix}",
        f"DRS Peak TS Avg for Cer Plastic, {suffix};Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Cer_Plastic_AvgPeakTS"
    )
    h1s_DRSPeakTS["Cer"].append(h1_Cer_Plastic_AvgPeakTS)

    h1_Sci_AvgPeakTS = rdf.Histo1D((
        f"hist_DRSPeakTS_Sci_Avg_{suffix}",
        f"DRS Peak TS Avg for Sci, {suffix};Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Sci_AvgPeakTS"
    )
    h1s_DRSPeakTS["Sci"].append(h1_Sci_AvgPeakTS)

    # Quartz vs Plastic average
    h2_DRSPeak_Cer_Quartz_VS_Plastic = rdf.Histo2D((
        f"hist_DRSPeakTS_Cer_Quartz_VS_Cer_Plastic_Avg_{suffix}",
        f"DRS Peak TS - CER Quartz VS CER Plastic Avg, {suffix};Cer Plastic Peak TS;Cer Quartz Peak TS",
        TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
        "Cer_Plastic_AvgPeakTS",
        "Cer_Quartz_AvgPeakTS",
    )
    h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_Quartz_VS_Plastic)

    return h1s_DRSPeakTS["Cer"], h1s_DRSPeakTS["Sci"], h2s_DRSPeakTS_Cer_VS_Sci


def checkDRSvsCalibrationTS(rdf, suffix=""):
    """Book shower shape vs calibrated time histograms."""
    hprofs_DRS_VS_TS = []
    hprofs_DRS_VS_Z = []
    hists2d_DRS_VS_TS = []
    rdfs_filtered = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.get_channel_name(blsub=False)
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                # rdf = rdf.Define(f"{channelName}_hasSignal",
                #                 f"{channelName}_Sum > 1000.0")

                # rdf_filtered = rdf.Filter(f"{channelName}_hasSignal")
                rdf_filtered = rdf

                hprof_DRS_VS_TS = rdf_filtered.Profile1D((
                    f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}", "", 400, -300, 100),
                    f"{channelName}_AlignedTS", f"{channelName}_blsub")

                hist2d_DRS_VS_TS = rdf_filtered.Histo2D((
                    f"hist2d_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}", "", 400, -300, 100, 100, -200, 1000),
                    f"{channelName}_AlignedTS", f"{channelName}_blsub")
                # hist_DRS_TS = rdf_filtered.Histo1D((
                #    f"hist_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}", "", 400, -300, 100),
                #    f"{channelName}_AlignedTS", f"{channelName}_blsub")

                if var != "Sci":
                    hprof_DRS_VS_Z = rdf_filtered.Profile1D((
                        f"prof_DRS_vs_Z_{var}_{sTowerX}_{sTowerY}_{suffix}", "", 200, -100, 600),
                        f"{channelName}_MeasuredZ", f"{channelName}_blsub")
                    hprofs_DRS_VS_Z.append(hprof_DRS_VS_Z)

                rdfs_filtered.append(rdf_filtered)
                hprofs_DRS_VS_TS.append(hprof_DRS_VS_TS)
                hists2d_DRS_VS_TS.append(hist2d_DRS_VS_TS)
                # hists2d_DRS_VS_TS.append(hist_DRS_TS)

    return hprofs_DRS_VS_TS, hists2d_DRS_VS_TS, hprofs_DRS_VS_Z, rdfs_filtered


def makeDRSPeakTSPlots(suffix=""):
    """Generate DRS peak TS plots using PlotManager."""
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRSPeakTS_rel_us_{suffix}")
    infile_name = f"{rootdir}/drspeakts_rel_us_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")

    hists_Cer = []
    hists_Cer_Quartz = []
    hists_Cer_Plastic = []
    hists_Sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        if board_no > 3:
            continue

        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hists = {}
            channelNos = {}
            isQuartz = False

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist_name = f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}_{suffix}"
                hist = infile.Get(hist_name)

                if var == "Cer" and chan.isQuartz:
                    isQuartz = True

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    hists[var] = None
                    channelNos[var] = "N/A"
                else:
                    hists[var] = hist
                    channelNos[var] = f"({chan.board_no},{chan.group_no},{chan.channel_no})"

            output_name = f"hist_DRSPeakTS_{sTowerX}_{sTowerY}_{suffix}"

            if not hists["Cer"] and not hists["Sci"]:
                print(
                    f"Warning: No histograms found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y})")
                continue

            # Build plot components
            colors = []
            hists_to_draw = []
            labels = []

            pave = create_pave_text(0.20, 0.75, 0.70, 0.90)
            pave.AddText(f"Tower: ({i_tower_x}, {i_tower_y})")

            if hists["Cer"]:
                colors.append(2 if isQuartz else 6)
                hists_to_draw.append(hists["Cer"])
                labels.append("Quartz Cer" if isQuartz else "Plastic Cer")
                pave.AddText(
                    f"Cer: {channelNos['Cer']}, {GetMean(hists['Cer']):.2f} TS")
            if hists["Sci"]:
                colors.append(4)
                hists_to_draw.append(hists["Sci"])
                labels.append("Sci")
                pave.AddText(
                    f"Sci: {channelNos['Sci']}, {GetMean(hists['Sci']):.2f} TS")

            if hists["Cer"] and hists["Sci"]:
                if isQuartz:
                    hists_Cer_Quartz.append(hists["Cer"])
                else:
                    hists_Cer_Plastic.append(hists["Cer"])
                hists_Cer.append(hists["Cer"])
                hists_Sci.append(hists["Sci"])

            style = PlotStyle(
                dology=False,
                drawoptions="HIST",
                mycolors=colors,
                addOverflow=False,
                addUnderflow=False,
                extra_text="Quartz + Sci" if isQuartz else "Plastic + Sci"
            )

            if doDetailedPeakTSPlots:
                pm.plot_1d(
                    hists_to_draw,
                    output_name,
                    "Peak TS",
                    (TSmin, TSmax),
                    ylabel="Counts",
                    yrange=(1, None),
                    legends=labels,
                    style=style,
                    extraToDraw=pave
                )

    # Summary plots
    hist_Cer_Combined = LHistos2Hist(
        hists_Cer, f"hist_DRSPeakTS_Cer_Combined_{suffix}")
    hist_Sci_Combined = LHistos2Hist(
        hists_Sci, f"hist_DRSPeakTS_Sci_Combined_{suffix}")
    hist_Cer_Quartz_Combined = LHistos2Hist(
        hists_Cer_Quartz, f"hist_DRSPeakTS_Cer_Quartz_Combined_{suffix}")
    hist_Cer_Plastic_Combined = LHistos2Hist(
        hists_Cer_Plastic, f"hist_DRSPeakTS_Cer_Plastic_Combined_{suffix}")

    # Combined plot with stats
    pave_combined = create_pave_text(0.50, 0.75, 0.90, 0.90)
    pave_combined.AddText(
        f"Quartz: {GetMean(hist_Cer_Quartz_Combined):.2f} TS")
    pave_combined.AddText(
        f"Plastic: {GetMean(hist_Cer_Plastic_Combined):.2f} TS")
    pave_combined.AddText(f"Sci: {GetMean(hist_Sci_Combined):.2f} TS")

    pm.plot_1d(
        [hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined, hist_Sci_Combined],
        f"DRS_PeakTS_Combined_{suffix}",
        "Peak TS",
        (TSmin, TSmax),
        ylabel="Counts",
        yrange=(1, None),
        legends=["Cer Quartz", "Cer Plastic", "Sci"],
        style=STYLE_CER_QUARTZ_PLASTIC_SCI,
        extraToDraw=pave_combined,
        prepend=True
    )

    pm.plot_1d(
        [hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined],
        f"DRS_PeakTS_Cer_Combined_{suffix}",
        "Peak TS",
        (TSCermin, TSCermax),
        ylabel="Fraction",
        yrange=(0, 0.01),
        legends=["Cer Quartz", "Cer Plastic"],
        style=PlotStyle(
            dology=False,
            drawoptions="HIST",
            mycolors=[2, 6],
            addOverflow=False,
            addUnderflow=False,
            legendPos=[0.30, 0.80, 0.40, 0.90],
            donormalize=True
        ),
        prepend=True
    )

    pm.add_newline()

    output_html = pm.generate_html(
        f"DRS/{suffix}/DRSPeakTS_relative_to_MCP.html",
        plots_per_row=4,
        title=f"DRS Peak Position for {suffix}",
        intro_text="Event-by-event DRS peak (channel-wise shower max) positions with respect to MCP for different towers.\n For now only the central DRS boards are included."
    )

    # Save combined histograms
    outfile_name = f"{rootdir}/drspeakts_rel_us_combined_{suffix}.root"
    outfile = ROOT.TFile(outfile_name, "RECREATE")
    for h in [hist_Cer_Combined, hist_Sci_Combined, hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined]:
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()

    return output_html


def makeDRSPeakTSCerVSSciPlots(suffix=""):
    """Generate DRS peak TS Cer vs Sci correlation plots using PlotManager."""
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRSPeakTSCerVSSci_rel_us_{suffix}")
    infile_name = f"{rootdir}/drspeakts_rel_us_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")

    hists = []
    hists_Quartz = []
    hists_Plastic = []

    # Diagonal reference line
    diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
    diagonal_line.SetLineStyle(2)
    diagonal_line.SetLineWidth(1)
    diagonal_line.SetLineColor(ROOT.kRed)

    for _, DRSBoard in DRSBoards.items():
        if DRSBoard.board_no > 3:
            continue

        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            chan_cer = DRSBoard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            isQuartz = chan_cer and chan_cer.isQuartz

            hist_name = f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}_{suffix}"
            hist = infile.Get(hist_name)
            output_name = f"DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}_{suffix}"
            output_name += "_Quartz" if isQuartz else "_Plastic"

            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue

            hists.append(hist)
            if isQuartz:
                hists_Quartz.append(hist)
            else:
                hists_Plastic.append(hist)

            ytitle = f"Cer ({'Quartz' if isQuartz else 'Plastic'}) Peak TS"

            if doDetailedPeakTS2DPlots:
                pm.plot_2d(
                    hist,
                    output_name,
                    "Sci Peak TS",
                    (TSmin, TSmax),
                    ytitle,
                    (TSmin, TSmax),
                    style=STYLE_2D_LOG,
                    extraToDraw=diagonal_line
                )

    # Summary plots
    hcombined = LHistos2Hist(
        hists, f"hist_DRSPeakTS_Cer_VS_Sci_Combined_{suffix}")
    pm.plot_2d(
        hcombined,
        f"DRS_PeakTS_Cer_VS_Sci_Combined_{suffix}",
        "Sci Peak TS",
        (TSmin, TSmax),
        "Cer Peak TS",
        (TSmin, TSmax),
        style=PlotStyle(dology=False, dologz=True, drawoptions="COLZ",
                        addOverflow=False, zmin=1, zmax=None),
        extraToDraw=diagonal_line,
        prepend=True
    )

    hcombined_Quartz = LHistos2Hist(
        hists_Quartz, f"hist_DRSPeakTS_Cer_VS_Sci_Quartz_Combined_{suffix}")
    pm.plot_2d(
        hcombined_Quartz,
        f"DRS_PeakTS_Cer_Quartz_VS_Sci_Combined_{suffix}",
        "Sci Peak TS",
        (TSmin, TSmax),
        "Cer (Quartz) Peak TS",
        (TSmin, TSmax),
        style=PlotStyle(dology=False, dologz=True, drawoptions="COLZ",
                        addOverflow=False, zmin=1, zmax=None),
        extraToDraw=diagonal_line,
        prepend=True
    )

    hcombined_Plastic = LHistos2Hist(
        hists_Plastic, f"hist_DRSPeakTS_Cer_VS_Sci_Plastic_Combined_{suffix}")
    pm.plot_2d(
        hcombined_Plastic,
        f"DRS_PeakTS_Cer_Plastic_VS_Sci_Combined_{suffix}",
        "Sci Peak TS",
        (TSmin, TSmax),
        "Cer (Plastic) Peak TS",
        (TSmin, TSmax),
        style=PlotStyle(dology=False, dologz=True, drawoptions="COLZ",
                        addOverflow=False, zmin=1, zmax=None),
        extraToDraw=diagonal_line,
        prepend=True
    )

    # Quartz vs Plastic average
    hcombined_Quartz_vs_Plastic = infile.Get(
        f"hist_DRSPeakTS_Cer_Quartz_VS_Cer_Plastic_Avg_{suffix}")
    if hcombined_Quartz_vs_Plastic:
        pm.plot_2d(
            hcombined_Quartz_vs_Plastic,
            f"DRS_PeakTS_Cer_Quartz_VS_Cer_Plastic_Avg_{suffix}",
            "Cer Plastic Peak TS",
            (TSmin, TSmax),
            "Cer Quartz Peak TS",
            (TSmin, TSmax),
            style=PlotStyle(dology=False, dologz=True, drawoptions="COLZ",
                            addOverflow=False, zmin=1, zmax=None),
            extraToDraw=diagonal_line,
            prepend=True
        )

    print(f"quartz: {len(hists_Quartz)}, plastic: {len(hists_Plastic)}")

    return pm.generate_html(f"DRS/{suffix}/DRSPeakTS_Cer_VS_Sci_relative_to_MCP_{suffix}.html", plots_per_row=4)


def makeDRSvsTSProfPlots(suffix: str = ""):
    """Generate profiled DRS ADC vs calibrated TS plots using PlotManager."""
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRS_vs_ts_calibrated_{suffix}")
    infile_name = f"{rootdir}/drs_vs_ts_calibrated_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofs_Cer_Quartz = []
    hprofs_Cer_Plastic = []
    hprofs_Cer = []
    hprofs_Sci = []
    xtitle = "TS"
    ytitle = "Average DRS ADC"

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        pm.add_newline()
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hprofs = {}
            isQuartz = False
            channelNos = {}

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hprof_name = f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}"
                hprof = infile.Get(hprof_name)

                if var == "Cer" and chan.isQuartz:
                    isQuartz = True

                if not hprof:
                    print(
                        f"Warning: Histogram {hprof_name} not found in {infile_name}")
                    hprofs[var] = None
                else:
                    hprofs[var] = hprof.ProjectionX()
                    channelNos[var] = f"({chan.board_no},{chan.group_no},{chan.channel_no})"

            if not hprofs["Cer"] and not hprofs["Sci"]:
                print(
                    f"Warning: No histograms found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y})")
                continue

            # Build plot components
            pave = create_pave_text(0.20, 0.75, 0.70, 0.90)
            pave.AddText(f"Tower: ({i_tower_x}, {i_tower_y})")

            colors = []
            hists_to_draw = []
            labels = []

            if hprofs["Cer"]:
                colors.append(2 if isQuartz else 6)
                hists_to_draw.append(hprofs["Cer"])
                labels.append("Quartz Cer" if isQuartz else "Plastic Cer")
                pave.AddText(f"Cer: {channelNos['Cer']}")
            if hprofs["Sci"]:
                colors.append(4)
                hists_to_draw.append(hprofs["Sci"])
                labels.append("Sci")
                pave.AddText(f"Sci: {channelNos['Sci']}")

            if hprofs["Cer"] and hprofs["Sci"]:
                if isQuartz:
                    hprofs_Cer_Quartz.append(hprofs["Cer"])
                else:
                    hprofs_Cer_Plastic.append(hprofs["Cer"])
                hprofs_Cer.append(hprofs["Cer"])
                hprofs_Sci.append(hprofs["Sci"])

            output_name = f"prof_DRS_vs_TS_{sTowerX}_{sTowerY}_{suffix}"

            pm.plot_1d(
                hists_to_draw,
                output_name,
                xtitle,
                (-200, 100),
                ylabel=ytitle,
                yrange=(-10, None),
                legends=labels,
                style=PlotStyle(
                    dology=False,
                    drawoptions="HIST",
                    mycolors=colors,
                    addOverflow=False,
                    addUnderflow=False
                ),
                extraToDraw=pave
            )

    # Summary plots
    hprof_Cer_Combined = LHistos2Hist(
        hprofs_Cer, f"prof_DRS_vs_TS_Cer_Combined_{suffix}")
    hprof_Sci_Combined = LHistos2Hist(
        hprofs_Sci, f"prof_DRS_vs_TS_Sci_Combined_{suffix}")
    hprof_Cer_Quartz_Combined = LHistos2Hist(
        hprofs_Cer_Quartz, f"prof_DRS_vs_TS_Cer_Quartz_Combined_{suffix}")
    hprof_Cer_Plastic_Combined = LHistos2Hist(
        hprofs_Cer_Plastic, f"prof_DRS_vs_TS_Cer_Plastic_Combined_{suffix}")

    hprof_Cer_Quartz_Combined.GetXaxis().SetRangeUser(-80, -55)

    pm.plot_1d(
        [hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined, hprof_Sci_Combined],
        f"DRS_vs_TS_Cer_Sci_Combined_{suffix}",
        xtitle,
        (-75, 0),
        ylabel=ytitle,
        yrange=(0, None),
        legends=["Cer Quartz", "Cer Plastic", "Sci"],
        style=STYLE_CER_QUARTZ_PLASTIC_SCI,
        prepend=True
    )

    pm.plot_1d(
        [hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined],
        f"DRS_vs_TS_Cer_Combined_{suffix}",
        xtitle,
        (-80, -50),
        ylabel=ytitle,
        yrange=(0, None),
        legends=["Cer Quartz", "Cer Plastic"],
        style=STYLE_CER_QUARTZ_PLASTIC,
        legendoptions=["L", "L"],
        prepend=True
    )

    pm.add_newline()

    output_html = pm.generate_html(
        f"DRS/{suffix}/DRS_VS_TS_relative_to_MCP.html",
        plots_per_row=4,
        title=f"Profiled DRS ADC vs TS for {suffix}",
        intro_text="Average DRS ADC vs calibrated time slice (TS) for Cer and Sci channels. For each event, the time slice is measured relative to the MCP signal."
    )

    # Save combined histograms
    outfile_name = f"{rootdir}/drs_vs_ts_calibrated_combined_{suffix}.root"
    outfile = ROOT.TFile(outfile_name, "RECREATE")
    for h in [hprof_Cer_Combined, hprof_Sci_Combined, hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined]:
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()

    return output_html


def makeDRSvsZProfPlots(suffix: str = ""):
    """Generate profiled DRS ADC vs measured Z plots using PlotManager."""
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRS_vs_z_calibrated_{suffix}")
    infile_name = f"{rootdir}/drs_vs_z_calibrated_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofs_Cer_Quartz = []
    hprofs_Cer_Plastic = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        if board_no > 3:
            continue

        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            isQuartz = False
            for var in ["Cer"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hprof_name = f"prof_DRS_vs_Z_{var}_{sTowerX}_{sTowerY}_{suffix}"
                hprof = infile.Get(hprof_name)

                if var == "Cer" and chan.isQuartz:
                    isQuartz = True

                if not hprof:
                    print(
                        f"Warning: Histogram {hprof_name} not found in {infile_name}")
                    continue

                hprof_proj = hprof.ProjectionX()

                if isQuartz:
                    hprofs_Cer_Quartz.append(hprof_proj)
                else:
                    hprofs_Cer_Plastic.append(hprof_proj)

                output_name = f"prof_DRS_vs_Z_{sTowerX}_{sTowerY}_{suffix}"
                colors = [2, 4] if isQuartz else [6, 4]

                pm.plot_1d(
                    [hprof_proj],
                    output_name,
                    "Measured Z [cm]",
                    (-100, 300),
                    ylabel="DRS ADC",
                    yrange=(-10, hprof_proj.GetMaximum() * 1.5),
                    legends=["Cer"],
                    style=PlotStyle(
                        dology=False,
                        drawoptions="HIST",
                        mycolors=colors,
                        addOverflow=False,
                        addUnderflow=False
                    )
                )

    # Summary plots
    hprof_Cer_Quartz_Combined = LHistos2Hist(
        hprofs_Cer_Quartz, f"prof_DRS_vs_Z_Cer_Quartz_Combined_{suffix}")
    hprof_Cer_Plastic_Combined = LHistos2Hist(
        hprofs_Cer_Plastic, f"prof_DRS_vs_Z_Cer_Plastic_Combined_{suffix}")

    pm.plot_1d(
        [hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined],
        f"DRS_vs_Z_Cer_Combined_{suffix}",
        "Measured Z",
        (-10, 200),
        ylabel="DRS ADC",
        yrange=(0, None),
        legends=["Cer Quartz", "Cer Plastic"],
        style=PlotStyle(
            dology=False,
            drawoptions=["C", "C"],
            mycolors=[2, 6],
            addOverflow=False,
            addUnderflow=False,
            legendPos=[0.55, 0.80, 0.90, 0.90]
        ),
        legendoptions=["P", "P"],
        prepend=True
    )
    pm.add_newline()

    return pm.generate_html(f"DRS/{suffix}/DRS_VS_Z_relative_to_MCP_{suffix}.html", plots_per_row=4)


def makeDRSvsTS2DPlots(suffix: str = ""):
    """Generate 2D DRS ADC vs calibrated TS plots using PlotManager."""
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRS_vs_ts_calibrated_{suffix}")
    infile_name = f"{rootdir}/drs_vs_ts_calibrated_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        if board_no > 3:
            continue

        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                hist2d_name = f"hist2d_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}"
                hist2d = infile.Get(hist2d_name)
                output_name = f"hist2d_DRS_vs_TS_{sTowerX}_{sTowerY}_{var}_{suffix}"

                if not hist2d:
                    print(
                        f"Warning: Histogram {hist2d_name} not found in {infile_name}")
                    pm.add_plot(output_name)  # Placeholder
                    continue

                pm.plot_2d(
                    hist2d,
                    output_name,
                    "Calibrated TS",
                    (-200, 100),
                    "DRS ADC",
                    (-200, 1000),
                    style=PlotStyle(
                        dology=False,
                        dologz=True,
                        drawoptions="COLZ",
                        addOverflow=False,
                        addUnderflow=False,
                        zmin=1,
                        zmax=1e4
                    )
                )

    return pm.generate_html(f"DRS/{suffix}/DRS_VS_TS_relative_to_MCP_2D_{suffix}.html", plots_per_row=4)


# ============================================================================
# Main Execution
# ============================================================================
rdf = analysis.get_rdf()

rdfs = OrderedDict()
rdfs["inclusive"] = rdf
if analysis.beam_type in ["e+", "positron", "e-", "electron", "positrons", "electrons"]:
    rdfs["electron"] = analysis.get_particle_analysis("electron")
    rdfs["pion"] = analysis.get_particle_analysis("pion")
else:
    rdfs["proton"] = analysis.get_particle_analysis("proton")
    rdfs["pion"] = analysis.get_particle_analysis("pion")


def main():
    makeHists = True
    makePlots = True

    booked_results = {}

    if makeHists:
        for cat, rdf in rdfs.items():
            rdf = calibrateDRSPeakTS(rdf, run_number, DRSBoards,
                                     TSminDRS=0, TSmaxDRS=1000, threshold=9.0)

            rdf_prefilterMCP1 = rdf
            map_mcp_channels = get_mcp_channels(run_number)

            condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
            condition += f" && {map_mcp_channels['US'][0]}_PeakTS > 500 && {map_mcp_channels['US'][0]}_PeakTS < 600"
            condition += f" && {map_mcp_channels['US'][0]}_Peak < -300.0"
            rdf_prefilterMCP2 = rdf_prefilterMCP1.Filter(condition,
                                                         "Pre-filter on MCP US channel 0 Peak TS")

            rdf_prefilterMCP2 = rdf_prefilterMCP1.Define(
                "MCP0_DeltaRelPeakTS", f"{map_mcp_channels['DS'][0]}_RelPeakTS - {map_mcp_channels['US'][0]}_RelPeakTS")

            rdf = rdf_prefilterMCP2

            booked_results[cat] = {}
            booked_results[cat]['DRS_PeakTS'], booked_results[cat]['DRS_PeakTS_Sci'], booked_results[cat]['DRS_PeakTS_Cer_VS_Sci'] = checkDRSPeakTS(
                rdf, suffix=cat)
            booked_results[cat]['DRS_VS_TS'], booked_results[cat]['DRS_VS_TS_2D'], booked_results[cat]['DRS_VS_Z'], booked_results[cat]['RDFs_Filtered'] = checkDRSvsCalibrationTS(
                rdf, suffix=cat)

        # collect proxies
        all_proxies = []
        for cat, results in booked_results.items():
            all_proxies.extend(results['DRS_PeakTS'])
            all_proxies.extend(results['DRS_PeakTS_Sci'])
            all_proxies.extend(results['DRS_PeakTS_Cer_VS_Sci'])
            all_proxies.extend(results['DRS_VS_TS'])
            all_proxies.extend(results['DRS_VS_TS_2D'])
            all_proxies.extend(results['DRS_VS_Z'])

        # Trigger the event loop to fill all histograms
        ROOT.RDF.RunGraphs(all_proxies)

        # Save DRS Peak TS histograms
        for cat, hists in booked_results.items():
            outfile_DRSPeakTS = ROOT.TFile(
                f"{rootdir}/drspeakts_rel_us_{cat}.root", "RECREATE")
            for hist in hists['DRS_PeakTS']:
                hist.SetDirectory(outfile_DRSPeakTS)
                hist.Write()
            for hist in hists['DRS_PeakTS_Sci']:
                hist.SetDirectory(outfile_DRSPeakTS)
                hist.Write()
            for hist in hists['DRS_PeakTS_Cer_VS_Sci']:
                hist.SetDirectory(outfile_DRSPeakTS)
                hist.Write()
            outfile_DRSPeakTS.Close()

            # Save DRS vs TS histograms
            outfile_DRS_VS_TS = ROOT.TFile(
                f"{rootdir}/drs_vs_ts_calibrated_{cat}.root", "RECREATE")
            for hprof in hists['DRS_VS_TS']:
                hprof.SetDirectory(outfile_DRS_VS_TS)
                hprof.Write()
            for hist2d in hists['DRS_VS_TS_2D']:
                hist2d.SetDirectory(outfile_DRS_VS_TS)
                hist2d.Write()
            for hprof in hists['DRS_VS_Z']:
                hprof.SetDirectory(outfile_DRS_VS_TS)
                hprof.Write()
            outfile_DRS_VS_TS.Close()

    if makePlots:
        output_htmls = {}
        for cat in rdfs.keys():
            output_htmls[f"DRSPeakTS_{cat}"] = makeDRSPeakTSPlots(suffix=cat)
            output_htmls[f"DRSPeakTSCerVSSci_{cat}"] = makeDRSPeakTSCerVSSciPlots(
                suffix=cat)
            output_htmls[f"DRSvsTSProf_{cat}"] = makeDRSvsTSProfPlots(
                suffix=cat)
            # output_htmls[f"DRSvsTS2D_{cat}"] = makeDRSvsTS2DPlots(suffix=cat)
            # output_htmls[f"DRSvsZProf_{cat}"] = makeDRSvsZProfPlots(suffix=cat)

        for key, html in output_htmls.items():
            print(f"{key}: {html}")


if __name__ == "__main__":
    main()
