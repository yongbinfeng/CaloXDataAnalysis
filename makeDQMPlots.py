import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels
from utils.utils import number2string
from utils.html_generator import generate_html
from utils.visualization import visualizeFERSBoards
from utils.validateMap import DrawFERSBoards, DrawDRSBoards
from utils.colors import colors
from runconfig import runNumber

print("Start running script")
ROOT.gROOT.SetBatch(True)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)
time_reference_channels = buildTimeReferenceChannels(run=runNumber)
hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
hodo_pos_channels = buildHodoPosChannels(run=runNumber)


rootdir = f"root/Run{runNumber}/"
outdir = f"plots/Run{runNumber}/"


def makeConditionsPlots():
    plots = []
    outdir_plots = outdir + "/Conditions_vs_Event"
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

        output_name = f"Conditions_Board{boardNo}_SipmHV_vs_Event"
        DrawHistos([hist_SipmHV], f"", 0, nEvents, "Event", 26, 29, "Voltage (V)",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

        output_name = f"Conditions_Board{boardNo}_SipmI_vs_Event"
        DrawHistos([hist_SipmI], f"", 0, nEvents, "Event", 0.02, 0.2, "Current (mA)",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

        output_name = f"Conditions_Board{boardNo}_TempDET_vs_Event"
        DrawHistos([hist_TempDET], f"", 0, nEvents, "Event", 10, 30, "Temperature (C)",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

        output_name = f"Conditions_Board{boardNo}_TempFPGA_vs_Event"
        DrawHistos([hist_TempFPGA], f"", 0, nEvents, "Event", 30, 50, "Temperature (C)",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True, extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

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

    output_html = f"html/Run{runNumber}/Conditions_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSEnergySumPlots():
    plots = []
    infile_name = f"{rootdir}/fers_energy_sum.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/FERS_EnergySum"
    hists_CerEnergyHG = []
    hists_SciEnergyHG = []
    hists_CerEnergyLG = []
    hists_SciEnergyLG = []
    legends = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_CerEnergyHG_name = f"hist_FERS_Board{boardNo}_CerEnergyHG"
        hist_SciEnergyHG_name = f"hist_FERS_Board{boardNo}_SciEnergyHG"
        hist_CerEnergyLG_name = f"hist_FERS_Board{boardNo}_CerEnergyLG"
        hist_SciEnergyLG_name = f"hist_FERS_Board{boardNo}_SciEnergyLG"
        hist_CerEnergyHG = infile.Get(hist_CerEnergyHG_name)
        hist_SciEnergyHG = infile.Get(hist_SciEnergyHG_name)
        hist_CerEnergyLG = infile.Get(hist_CerEnergyLG_name)
        hist_SciEnergyLG = infile.Get(hist_SciEnergyLG_name)

        if not (hist_CerEnergyHG and hist_SciEnergyHG and hist_CerEnergyLG and hist_SciEnergyLG):
            print(
                f"Warning: Histograms {hist_CerEnergyHG_name}, {hist_SciEnergyHG_name}, {hist_CerEnergyLG_name}, or {hist_SciEnergyLG_name} not found in {infile_name}")
            continue

        legends.append(str(boardNo))
        hists_CerEnergyHG.append(hist_CerEnergyHG)
        hists_SciEnergyHG.append(hist_SciEnergyHG)
        hists_CerEnergyLG.append(hist_CerEnergyLG)
        hists_SciEnergyLG.append(hist_SciEnergyLG)

    DrawHistos(hists_CerEnergyHG, legends, 0, 20000, "Cer Energy HG", 0, 1000, "Events",
               "FERS_CerEnergyHG",
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.6, 0.90, 0.9])
    plots.append("FERS_CerEnergyHG.png")
    DrawHistos(hists_SciEnergyHG, legends, 0, 20000, "Sci Energy HG", 0, 1000, "Events",
               "FERS_SciEnergyHG",
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.6, 0.90, 0.9])
    plots.append("FERS_SciEnergyHG.png")
    DrawHistos(hists_CerEnergyLG, legends, 0, 20000, "Cer Energy LG", 0, 1000, "Events",
               "FERS_CerEnergyLG",
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.6, 0.90, 0.9])
    plots.append("FERS_CerEnergyLG.png")
    DrawHistos(hists_SciEnergyLG, legends, 0, 20000, "Sci Energy LG", 0, 1000, "Events",
               "FERS_SciEnergyLG",
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.6, 0.90, 0.9])
    plots.append("FERS_SciEnergyLG.png")

    # total energy sum plots
    hist_CerEnergyHG = infile.Get("hist_FERS_CerEnergyHG")
    hist_SciEnergyHG = infile.Get("hist_FERS_SciEnergyHG")
    hist_CerEnergyLG = infile.Get("hist_FERS_CerEnergyLG")
    hist_SciEnergyLG = infile.Get("hist_FERS_SciEnergyLG")
    DrawHistos([hist_CerEnergyHG], "", 0, 2e5, "Cer Energy HG", 0, 1000, "Events",
               "FERS_Total_CerEnergyHG",
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(0, "FERS_Total_CerEnergyHG.png")
    DrawHistos([hist_SciEnergyHG], "", 0, 2e5, "Sci Energy HG", 0, 1000, "Events",
               "FERS_Total_SciEnergyHG",
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(1, "FERS_Total_SciEnergyHG.png")
    DrawHistos([hist_CerEnergyLG], "", 0, 2e5, "Cer Energy LG", 0, 1000, "Events",
               "FERS_Total_CerEnergyLG",
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(2, "FERS_Total_CerEnergyLG.png")
    DrawHistos([hist_SciEnergyLG], "", 0, 2e5, "Sci Energy LG", 0, 1000, "Events",
               "FERS_Total_SciEnergyLG",
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(3, "FERS_Total_SciEnergyLG.png")

    output_html = f"html/Run{runNumber}/FERS_EnergySum/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERS1DPlots():
    plots = []

    infile_name = f"{rootdir}fers_all_channels_1D.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/FERS_1D"
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

    output_html = f"html/Run{runNumber}/FERS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


def makeFERSStatsPlots():
    plots = []
    outdir_plots = outdir + "/FERS_Stats"
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
    valuemaps_LG_mean = {}
    valuemaps_LG_max = {}

    for channelName, (vmean, vmax) in stats.items():
        if "energyHG" in channelName:
            valuemaps_HG_mean[channelName] = vmean
            valuemaps_HG_max[channelName] = vmax
        elif "energyLG" in channelName:
            valuemaps_LG_mean[channelName] = vmean
            valuemaps_LG_max[channelName] = vmax

    [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
        FERSBoards, valuemaps_HG_mean, suffix=f"Run{runNumber}_HG_mean", useHG=True)
    [h2_Cer_HG_max, h2_Cer_3mm_HG_max], [h2_Sci_HG_max, h2_Sci_3mm_HG_max] = visualizeFERSBoards(
        FERSBoards, valuemaps_HG_max, suffix=f"Run{runNumber}_HG_max", useHG=True)
    [h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], [h2_Sci_LG_mean, h2_Sci_3mm_LG_mean] = visualizeFERSBoards(
        FERSBoards, valuemaps_LG_mean, suffix=f"Run{runNumber}_LG_mean", useHG=False)
    [h2_Cer_LG_max, h2_Cer_3mm_LG_max], [h2_Sci_LG_max, h2_Sci_3mm_LG_max] = visualizeFERSBoards(
        FERSBoards, valuemaps_LG_max, suffix=f"Run{runNumber}_LG_max", useHG=False)

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_mean"
    DrawHistos([h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1600)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_mean, h2_Sci_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1600)
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

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_mean"
    DrawHistos([h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1600)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_mean, h2_Sci_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1600)
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

    output_html = f"html/Run{runNumber}/FERS_Stats/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# 2D FERS histograms, hg vs lg
def makeFERS2DPlots():
    plots = []
    outdir_plots = outdir + "/FERS_2D"
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
    output_html = f"html/Run{runNumber}/FERS_2D/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# FERS output vs event
def trackFERSPlots():
    plots = []
    outdir_plots = outdir + "/FERS_vs_Event"
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
    output_html = f"html/Run{runNumber}/FERS_vs_Event/index.html"
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
            outdir_plots = outdir + "/DRS_1D"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 1400, 2500, "DRS Output", 1, 1e12, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       legendPos=(0.60, 0.78, 0.90, 0.68),
                       outdir=outdir_plots)
            plots.append(output_name + ".png")

    output_html = f"html/Run{runNumber}/DRS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


# DRS vs TS
def makeDRS2DPlots(doSubtractMedian=False):
    suffix = ""
    ymin = -50
    ymax = 50
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -20
        ymax = 40
    plots = []
    outdir_plots = outdir + "/DRS_vs_TS"
    infile_name = f"{rootdir}/drs_all_channels_2D.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}{suffix}"
                hist = infile.Get(hist_name)
                output_name = f"DRS_{var}_vs_TS_{sTowerX}_{sTowerY}{suffix}"
                plots.append(output_name + ".png")

                ymax_tmp = ymax
                ymin_tmp = ymin
                if boardNo == 1 and chan.groupNo == 0:
                    if chan.channelNo == 0 or chan.channelNo == 1:
                        ymax_tmp = 300
                        ymin_tmp = -100

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                value_mean = hist.GetMean(2)

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

                DrawHistos([hist], "", 0, 1024, "Time Slice", value_mean + ymin_tmp, value_mean + ymax_tmp, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
    output_html = f"html/Run{runNumber}/DRS_vs_TS{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


# DRS mean vs event
def trackDRSPlots():
    plots = []
    infile_name = f"{rootdir}/drs_all_channels_2D_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/DRS_vs_Event"
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
    output_html = f"html/Run{runNumber}/DRS_vs_Event/index.html"
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
    outdir_plots = outdir + "/TimeReference"
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

    output_html = f"html/Run{runNumber}/TimeReference{suffix}/index.html"
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
    outdir_plots = outdir + "/HodoTrigger"
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

    output_html = f"html/Run{runNumber}/HodoTrigger{suffix}/index.html"
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
    outdir_plots = outdir + "/HodoPos"
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

    output_html = f"html/Run{runNumber}/HodoPos{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


if __name__ == "__main__":
    output_htmls = {}

    output_htmls["conditions plots"] = makeConditionsPlots()

    # validate DRS and FERS boards
    output_htmls["fers mapping"] = DrawFERSBoards(run=runNumber)
    output_htmls["drs mapping"] = DrawDRSBoards(run=runNumber)

    output_htmls["fers 1D"] = makeFERS1DPlots()
    output_htmls["fers energy sum"] = makeFERSEnergySumPlots()
    output_htmls["fers stats"] = makeFERSStatsPlots()
    # makeDRS2DPlots()
    output_htmls["drs 2D"] = makeDRS2DPlots(doSubtractMedian=True)

    output_htmls["time reference"] = compareTimeReferencePlots(True)
    output_htmls["hodo trigger"] = compareHodoTriggerPlots(True)
    output_htmls["hodo pos"] = compareHodoPosPlots(True)

    print("\n\n\n")
    print("*" * 30)
    for key, value in output_htmls.items():
        print(f"✅ {key} plots can be viewed at: {value}")

    print("All plots generated successfully.")
