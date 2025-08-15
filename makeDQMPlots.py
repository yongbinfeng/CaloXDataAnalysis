import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels
from utils.utils import number2string, round_up_to_1eN
from utils.html_generator import generate_html
from utils.visualization import visualizeFERSBoards
from utils.validateMap import DrawFERSBoards, DrawDRSBoards
from utils.colors import colors
from configs.plotranges import getDRSPlotRanges
from runconfig import runNumber

print("Start running script")
ROOT.gROOT.SetBatch(True)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)
time_reference_channels = buildTimeReferenceChannels(run=runNumber)
hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
hodo_pos_channels = buildHodoPosChannels(run=runNumber)


rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"


def makeConditionsPlots():
    plots = []
    outdir_plots = f"{plotdir}/Conditions_vs_Event"
    infile_name = f"{rootdir}/conditions_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofiles_SipmHV = []
    hprofiles_SipmI = []
    hprofiles_TempDET = []
    hprofiles_TempFPGA = []

    legends = []

    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_SipmHV_name = f"hist_FERS_Board{boardNo}_SipmHV_vs_Event"
        hist_SipmI_name = f"hist_FERS_Board{boardNo}_SipmI_vs_Event"
        hist_TempDET_name = f"hist_FERS_Board{boardNo}_TempDET_vs_Event"
        hist_TempFPGA_name = f"hist_FERS_Board{boardNo}_TempFPGA_vs_Event"

        hist_SipmHV = infile.Get(hist_SipmHV_name)
        hist_SipmI = infile.Get(hist_SipmI_name)
        hist_TempDET = infile.Get(hist_TempDET_name)
        hist_TempFPGA = infile.Get(hist_TempFPGA_name)

        if not (hist_SipmHV and hist_SipmI and hist_TempDET and hist_TempFPGA):
            print(
                f"Warning: Histograms {hist_SipmHV_name}, {hist_SipmI_name}, {hist_TempDET_name}, or {hist_TempFPGA_name} not found in {infile_name}")
            continue

        hprofile_SipmHV = hist_SipmHV.ProfileX(
            f"hprof_FERS_Board{boardNo}_SipmHV_vs_Event")
        hprofiles_SipmHV.append(hprofile_SipmHV)

        hprofile_SipmI = hist_SipmI.ProfileX(
            f"hprof_FERS_Board{boardNo}_SipmI_vs_Event")
        hprofiles_SipmI.append(hprofile_SipmI)

        hprofile_TempDET = hist_TempDET.ProfileX(
            f"hprof_FERS_Board{boardNo}_TempDET_vs_Event")
        hprofiles_TempDET.append(hprofile_TempDET)

        hprofile_TempFPGA = hist_TempFPGA.ProfileX(
            f"hprof_FERS_Board{boardNo}_TempFPGA_vs_Event")
        hprofiles_TempFPGA.append(hprofile_TempFPGA)

        legends.append(str(boardNo))

        nEvents = hist_SipmHV.GetXaxis().GetXmax()
        zmax = nEvents * 10

        extraToDraw = ROOT.TPaveText(0.20, 0.80, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Board: {FERSBoard.boardNo}")

        # output_name = f"Conditions_Board{boardNo}_SipmHV_vs_Event"
        # DrawHistos([hist_SipmHV], f"", 0, nEvents, "Event", 26, 29, "Voltage (V)",
        #           output_name,
        #           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
        #           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        # plots.append(output_name + ".png")

        # output_name = f"Conditions_Board{boardNo}_SipmI_vs_Event"
        # DrawHistos([hist_SipmI], f"", 0, nEvents, "Event", 0.02, 0.2, "Current (mA)",
        #           output_name,
        #           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
        #           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        # plots.append(output_name + ".png")

        # output_name = f"Conditions_Board{boardNo}_TempDET_vs_Event"
        # DrawHistos([hist_TempDET], f"", 0, nEvents, "Event", 10, 30, "Temperature (C)",
        #           output_name,
        #           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
        #           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        # plots.append(output_name + ".png")

        # output_name = f"Conditions_Board{boardNo}_TempFPGA_vs_Event"
        # DrawHistos([hist_TempFPGA], f"", 0, nEvents, "Event", 30, 50, "Temperature (C)",
        #           output_name,
        #           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
        #           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        # plots.append(output_name + ".png")

    # Draw the profiles
    legendPos = [0.3, 0.7, 0.9, 0.9]
    output_name = "Conditions_SipmHV_vs_Event"
    DrawHistos(hprofiles_SipmHV, legends, 0, nEvents, "Event", 26, 29, "Voltage (V)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(0, output_name + ".png")

    output_name = "Conditions_SipmI_vs_Event"
    DrawHistos(hprofiles_SipmI, legends, 0, nEvents, "Event", 0.02, 0.2, "Current (mA)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(1, output_name + ".png")

    output_name = "Conditions_TempDET_vs_Event"
    DrawHistos(hprofiles_TempDET, legends, 0, nEvents, "Event", 15, 30, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(2, output_name + ".png")

    output_name = "Conditions_TempFPGA_vs_Event"
    DrawHistos(hprofiles_TempFPGA, legends, 0, nEvents, "Event", 30, 50, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(3, output_name + ".png")

    output_html = f"{htmldir}/Conditions_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERS1DPlots():
    plots = []

    infile_name = f"{rootdir}fers_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_1D"
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hist_C_name = f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}"
            hist_S_name = f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}"
            hist_C = infile.Get(hist_C_name)
            hist_S = infile.Get(hist_S_name)
            if not hist_C or not hist_S:
                print(
                    f"Warning: Histograms {hist_C_name} or {hist_S_name} not found in {infile_name}")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"Board: {FERSBoard.boardNo}")
            extraToDraw.AddText(f"Tower X: {iTowerX}")
            extraToDraw.AddText(f"Tower Y: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: {FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True).channelNo}")
            extraToDraw.AddText(
                f"Sci Channel: {FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=False).channelNo}")

            output_name = f"Energy_Board{boardNo}_iTowerX{sTowerX}_iTowerY{sTowerY}"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)

            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


def makeFERSStatsPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_Stats"
    # load the json file
    import json
    infile_name = f"{rootdir}/fers_stats.json"
    with open(infile_name, "r") as f:
        stats = json.load(f)

    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100
    valuemaps_HG_mean = {}
    valuemaps_HG_max = {}
    valuemaps_HG_satfreq = {}
    valuemaps_LG_mean = {}
    valuemaps_LG_max = {}
    valuemaps_LG_satfreq = {}

    for channelName, (vmean, vmax, vsatfreq) in stats.items():
        if "energyHG" in channelName:
            valuemaps_HG_mean[channelName] = vmean
            valuemaps_HG_max[channelName] = vmax
            valuemaps_HG_satfreq[channelName] = vsatfreq
        elif "energyLG" in channelName:
            valuemaps_LG_mean[channelName] = vmean
            valuemaps_LG_max[channelName] = vmax
            valuemaps_LG_satfreq[channelName] = vsatfreq

    [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
        FERSBoards, valuemaps_HG_mean, suffix=f"Run{runNumber}_HG_mean", useHG=True)
    [h2_Cer_HG_max, h2_Cer_3mm_HG_max], [h2_Sci_HG_max, h2_Sci_3mm_HG_max] = visualizeFERSBoards(
        FERSBoards, valuemaps_HG_max, suffix=f"Run{runNumber}_HG_max", useHG=True)
    [h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], [h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq] = visualizeFERSBoards(
        FERSBoards, valuemaps_HG_satfreq, suffix=f"Run{runNumber}_HG_satfreq", useHG=True)
    [h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], [h2_Sci_LG_mean, h2_Sci_3mm_LG_mean] = visualizeFERSBoards(
        FERSBoards, valuemaps_LG_mean, suffix=f"Run{runNumber}_LG_mean", useHG=False)
    [h2_Cer_LG_max, h2_Cer_3mm_LG_max], [h2_Sci_LG_max, h2_Sci_3mm_LG_max] = visualizeFERSBoards(
        FERSBoards, valuemaps_LG_max, suffix=f"Run{runNumber}_LG_max", useHG=False)
    [h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], [h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq] = visualizeFERSBoards(
        FERSBoards, valuemaps_LG_satfreq, suffix=f"Run{runNumber}_LG_satfreq", useHG=False)

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_mean"
    DrawHistos([h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_mean, h2_Sci_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_max"
    DrawHistos([h2_Cer_HG_max, h2_Cer_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_max, h2_Sci_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_satfreq"
    DrawHistos([h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_mean"
    DrawHistos([h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_mean, h2_Sci_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_max"
    DrawHistos([h2_Cer_LG_max, h2_Cer_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_max, h2_Sci_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_satfreq"
    DrawHistos([h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_html = f"{htmldir}/FERS_Stats/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# 2D FERS histograms, hg vs lg
def makeFERS2DPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_2D"
    infile_name = f"{rootdir}/fers_all_channels_2D.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_vs_lg"
                hist = infile.Get(hist_name)

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {FERSBoard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_vs_lg"
                DrawHistos([hist], f"", 0, 9000, "HG", 0, 1500, "LG",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True, extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERS_2D/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# FERS output vs event
def trackFERSPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_vs_Event"
    infile_name = f"{rootdir}/fers_all_channels_2D_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {FERSBoard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                nEvents = hist.GetXaxis().GetXmax()

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", 0, nEvents, "Event", 1, 1e5, f"{var} Energy HG",
                           output_name,
                           dology=True, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERS_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# 1D histograms for DRS variables
def makeDRS1DPlots():
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            chan_Cer = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True)
            chan_Sci = DRSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)

            hist_C_name = f"hist_DRS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}"
            hist_S_name = f"hist_DRS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}"
            hist_C = infile.Get(hist_C_name)
            hist_S = infile.Get(hist_S_name)
            if not hist_C or not hist_S:
                print(
                    f"Warning: Histograms {hist_C_name} or {hist_S_name} not found in {infile_name}")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"Board: {DRSBoard.boardNo}")
            extraToDraw.AddText(f"iTowerX: {iTowerX}")
            extraToDraw.AddText(f"iTowerY: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: ({chan_Cer.groupNo}, {chan_Cer.channelNo})")
            extraToDraw.AddText(
                f"Sci Channel: ({chan_Sci.groupNo}, {chan_Sci.channelNo})")

            output_name = f"DRS_Variable_Board{boardNo}_iTowerX{sTowerX}_iTowerY{sTowerY}"
            outdir_plots = f"{plotdir}/DRS_1D"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 1400, 2500, "DRS Output", 1, 1e12, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       legendPos=(0.60, 0.78, 0.90, 0.68),
                       outdir=outdir_plots)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/DRS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


# DRS vs TS
def makeDRS2DPlots(doSubtractMedian=False, doRTS=0):
    suffix = ""
    if doSubtractMedian:
        suffix = "_subtractMedian"
    varTS = "TS"
    if doRTS == 1:
        varTS = "RTSpos"
    elif doRTS == 2:
        varTS = "RTSneg"

    plots = []
    outdir_plots = f"{plotdir}/DRS_vs_{varTS}"
    infile_name = f"{rootdir}/drs_vs_{varTS}.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_Board{boardNo}_{var}_vs_{varTS}_{sTowerX}_{sTowerY}{suffix}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                output_name = f"DRS_{var}_vs_{varTS}_{sTowerX}_{sTowerY}{suffix}"
                plots.append(output_name + ".png")

                ymin_tmp, ymax_tmp = getDRSPlotRanges(
                    subtractMedian=doSubtractMedian, isAmplified=chan.isAmplified)

                # value_mean = hist.GetMean(2)

                extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"B: {DRSBoard.boardNo}, G: {chan.groupNo}, C: {chan.channelNo}")
                extraToDraw.AddText(f"iTowerX: {iTowerX}")
                extraToDraw.AddText(f"iTowerY: {iTowerY}")

                DrawHistos([hist], "", 0, 1024, "Time Slice", ymin_tmp, ymax_tmp, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
    output_html = f"{htmldir}/DRS_vs_{varTS}{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS"
    infile_name = f"{rootdir}/drs_peak_ts.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_Cer = []
    hists_Sci = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hists = {}
            channelNos = {}
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_PeakTS_Board{boardNo}_peakTS_{sTowerX}_{sTowerY}_{var}"
                hist = infile.Get(hist_name)
                output_name = hist_name[5:-4]

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    hists[var] = None
                    channelNos[var] = -1
                else:
                    hists[var] = hist
                    channelNos[var] = chan.channelNo

            if not hists["Cer"] or not hists["Sci"]:
                print(
                    f"Warning: Histograms for Cer or Sci not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY})")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"B: {DRSBoard.boardNo}, G: {chan.groupNo}")
            extraToDraw.AddText(f"Tower X: {iTowerX}")
            extraToDraw.AddText(f"Tower Y: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: {channelNos['Cer']}")
            extraToDraw.AddText(
                f"Sci Channel: {channelNos['Sci']}")

            hists_Cer.append(hists["Cer"])
            hists_Sci.append(hists["Sci"])

            DrawHistos([hists["Cer"], hists["Sci"]], ["Cer", "Sci"], 0, 1000, "Peak TS", 1, None, "Counts",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    # summary plots
    hist_Cer_Sum = ROOT.TH1F("hist_DRS_PeakTS_Cer_Sum",
                             "DRS Peak TS Cer Sum", 1000, 0, 1000)
    hist_Sci_Sum = ROOT.TH1F("hist_DRS_PeakTS_Sci_Sum",
                             "DRS Peak TS Sci Sum", 1000, 0, 1000)
    for hist in hists_Cer:
        if hist:
            hist_Cer_Sum.Add(hist)
    for hist in hists_Sci:
        if hist:
            hist_Sci_Sum.Add(hist)
    DrawHistos([hist_Cer_Sum, hist_Sci_Sum], ["Cer", "Sci"], 0, 1000, "Peak TS", 1, None, "Counts",
               "DRS_PeakTS_Sum",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(0, "DRS_PeakTS_Sum.png")

    output_html = f"{htmldir}/DRSPeakTS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakTS2DPlots():
    plots = []
    hists = []
    outdir_plots = f"{plotdir}/DRSPeakTS2D"
    infile_name = f"{rootdir}/drs_peak_ts.root"
    infile = ROOT.TFile(infile_name, "READ")

    # Create a dashed diagonal line from (0,0) to (1000,1000)
    diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
    diagonal_line.SetLineStyle(2)  # 2 = dashed
    diagonal_line.SetLineWidth(1)
    diagonal_line.SetLineColor(ROOT.kRed)
    extraToDraw = diagonal_line

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hist_name = f"hist_DRSPeak_Cer_vs_Sci_Board{boardNo}_{sTowerX}_{sTowerY}"
            hist = infile.Get(hist_name)
            output_name = hist_name[5:]

            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue

            hists.append(hist)

            DrawHistos([hist], "", 0, 1000, "Cer Peak TS", 0, 1000, f"Sci Peak TS",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)

            plots.append(output_name + ".png")

    # summary plots
    h2 = ROOT.TH2F("hist_DRSPeak_Cer_vs_Sci_Sum",
                   "DRS Peak TS Cer vs Sci Sum", 1000, 0, 1000, 1000, 0, 1000)
    for hist in hists:
        if hist:
            h2.Add(hist)

    output_name = "DRS_PeakTS_Cer_vs_Sci_Sum"
    DrawHistos([h2], "", 0, 1000, "Cer Peak TS", 0, 1000, f"Sci Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/DRSPeakTS2D/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# DRS mean vs event
def trackDRSPlots():
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_2D_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/DRS_vs_Event"
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(f"Board: {DRSBoard.boardNo}")
                extraToDraw.AddText(f"iTowerX: {iTowerX}")
                extraToDraw.AddText(f"iTowerY: {iTowerY}")
                extraToDraw.AddText(f"{var} Group: {chan.groupNo}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                nEvents = hist.GetXaxis().GetXmax()

                output_name = f"DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", 0, nEvents, "Event", 1400, 2300, f"{var} Mean",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/DRS_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# time reference
def compareTimeReferencePlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -2500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/time_reference_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/TimeReference"
    for chan_name in time_reference_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"TimeReference_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/TimeReference{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# trigger
def compareHodoTriggerPlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -1500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/hodo_trigger_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/HodoTrigger"
    for chan_name in hodo_trigger_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"HodoTrigger_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/HodoTrigger{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# hodo position
def compareHodoPosPlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -1500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/hodo_pos_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/HodoPos"
    for board, channels in hodo_pos_channels.items():
        for chan_name in channels:
            hist_name = f"hist_{chan_name}{suffix}"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue
            extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"{chan_name}")
            output_name = f"HodoPos_{chan_name}{suffix}"
            DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                       extraToDraw=extraToDraw,
                       outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/HodoPos{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def checkFERSvsDRSSum():
    """
    Check if the sum of FERS and DRS energies are consistent.
    """
    plots = []
    outdir_plots = f"{plotdir}/FERS_vs_DRS_Sum"
    infile_name = f"{rootdir}/fers_vs_drs.root"
    infile = ROOT.TFile(infile_name, "READ")

    xymax = {
        "Cer": (1000, 4000),
        "Sci": (9000, 9000)
    }
    xymax_LG = {
        "Cer": (1000, 1000),
        "Sci": (9000, 9000)
    }

    # for _, DRSBoard in DRSBoards.items():
    #    boardNo = DRSBoard.boardNo
    #    for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
    #        sTowerX = number2string(iTowerX)
    #        sTowerY = number2string(iTowerY)

    #        for var in ["Cer", "Sci"]:
    #            for gain in ["FERS", "FERSLG"]:
    #                histname = f"hist_{gain}_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
    #                output_name = histname.replace("hist_", "")
    #                plots.append(output_name + ".png")

    #                hist = infile.Get(histname)
    #                if not hist:
    #                    print(
    #                        f"Warning: Histogram {histname} not found in {infile_name}")
    #                    continue

    #                zmax = hist.Integral(0, 10000, 0, 10000)
    #                zmax = round_up_to_1eN(zmax)

    #                tmp = xymax[var] if gain == "FERS" else xymax_LG[var]

    #                DrawHistos([hist], "", 0, tmp[0], "DRS Energy", 0, tmp[1], gain,
    #                           output_name,
    #                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
    #                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")

    # summary plots
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for var in ["Cer", "Sci"]:
            for gain in ["FERS", "FERSLG"]:
                histname = f"hist_{gain}_VS_DRS_Board{boardNo}_{var}_sum"
                output_name = histname.replace("hist_", "")
                # append to plots no matter what
                plots.append(output_name + ".png")

                hist = infile.Get(histname)
                if not hist:
                    print(
                        f"Warning: Histogram {histname} not found in {infile_name}")
                    continue

                zmax = hist.Integral(0, 10000, 0, 10000)
                zmax = round_up_to_1eN(zmax)

                output_name = histname.replace("hist_", "")

                tmp = xymax[var] if gain == "FERS" else xymax_LG[var]
                DrawHistos([hist], "", 0, tmp[0], "DRS Energy", 0, tmp[1], gain,
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")

    output_html = f"{htmldir}/FERS_vs_DRS_Sum/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def checkDRSPeakvsFERS():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeak_vs_FERS"
    infile_name = f"{rootdir}/drs_peak_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if not chan:
                    print(
                        f"Warning: Channel not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY}), Var {var}")
                    continue
                _, ymax = getDRSPlotRanges(
                    subtractMedian=True, isAmplified=chan.isAmplified)
                histname = f"hist_DRSPeak_VS_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                output_name = histname.replace("hist_", "")
                plots.append(output_name + ".png")

                hist = infile.Get(histname)
                if not hist:
                    print(
                        f"Warning: Histogram {histname} not found in {infile_name}")
                    continue

                DrawHistos([hist], "", 0, 9000, "FERS ADC", 0, ymax, "DRS Peak",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")

    output_html = f"{htmldir}/DRSPeak_vs_FERS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


if __name__ == "__main__":
    output_htmls = {}

    output_htmls["conditions plots"] = makeConditionsPlots()

    # validate DRS and FERS boards
    output_htmls["fers mapping"] = DrawFERSBoards(run=runNumber)
    output_htmls["drs mapping"] = DrawDRSBoards(run=runNumber)

    output_htmls["fers 1D"] = makeFERS1DPlots()
    output_htmls["fers stats"] = makeFERSStatsPlots()
    # makeDRS2DPlots()
    output_htmls["drs 2D"] = makeDRS2DPlots(doSubtractMedian=True)
    output_htmls["drs peak ts"] = makeDRSPeakTSPlots()
    output_htmls["drs peak ts 2D"] = makeDRSPeakTS2DPlots()

    # output_htmls["time reference"] = compareTimeReferencePlots(True)
    # output_htmls["hodo trigger"] = compareHodoTriggerPlots(True)
    # output_htmls["hodo pos"] = compareHodoPosPlots(True)

    output_htmls["fers vs drs sum"] = checkFERSvsDRSSum()

    output_htmls["drs peak vs fers"] = checkDRSPeakvsFERS()

    print("\n\n\n")
    print("*" * 30)
    for key, value in output_htmls.items():
        print(f"✅ {key} plots can be viewed at: {value}")

    print("All plots generated successfully.")
