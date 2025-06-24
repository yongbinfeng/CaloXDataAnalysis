import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTriggerChannels
from utils.utils import number2string
from utils.html_generator import generate_html
from utils.validateMap import DrawFERSBoards, DrawDRSBoards
from runNumber import runNumber

print("Start running script")
ROOT.gROOT.SetBatch(True)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)
trigger_channels = buildTriggerChannels(run=runNumber)


rootdir = f"root/Run{runNumber}/"
outdir = f"plots/Run{runNumber}/"


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
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)

            plots.append(output_name + ".png")

    generate_html(plots, outdir_plots,
                  output_html=f"html/Run{runNumber}/FERS_1D/viewer.html")


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
                           outdir=outdir_plots)
                plots.append(output_name + ".png")
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=f"html/Run{runNumber}/FERS_2D/viewer.html")


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
                           outdir=outdir_plots)
                plots.append(output_name + ".png")
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=f"html/Run{runNumber}/FERS_vs_Event/viewer.html")


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
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, extraToDraw=extraToDraw,
                       legendPos=(0.60, 0.78, 0.90, 0.68),
                       outdir=outdir_plots)
            plots.append(output_name + ".png")

    generate_html(plots, outdir_plots,
                  output_html=f"html/Run{runNumber}/DRS_1D/viewer.html")


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
                output_name = f"DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}{suffix}"
                plots.append(output_name + ".png")

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

                DrawHistos([hist], "", 0, 1024, "Time Slice", value_mean + ymin, value_mean + ymax, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber)
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=f"html/Run{runNumber}/DRS_vs_TS{suffix}/viewer.html")


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
                           outdir=outdir_plots)
                plots.append(output_name + ".png")
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=f"html/Run{runNumber}/DRS_vs_Event/viewer.html")


# trigger
def makeTriggerPlots():
    plots = []
    infile_name = f"{rootdir}/trigger_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = outdir + "/Trigger"
    for chan_name in trigger_channels:
        hist_name = f"hist_{chan_name}"
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
        output_name = f"Trigger_{chan_name}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", 1500, 2200, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots)
        plots.append(output_name + ".png")

        hist_subtractMedian = infile.Get(f"{hist_name}_subtractMedian")
        if hist_subtractMedian:
            extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"{chan_name} (subtract median)")
            output_name = f"Trigger_{chan_name}_subtractMedian"
            DrawHistos([hist_subtractMedian], "", 0, 1024, "Time Slice", -700, 300, "Counts",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                       extraToDraw=extraToDraw,
                       outdir=outdir_plots)
            plots.append(output_name + ".png")

    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=f"html/Run{runNumber}/Trigger/viewer.html")


if __name__ == "__main__":

    # validate DRS and FERS boards
    DrawFERSBoards(run=runNumber)
    DrawDRSBoards(run=runNumber)

    makeFERS1DPlots()
    makeDRS2DPlots()
    makeDRS2DPlots(doSubtractMedian=True)

    print("All plots generated successfully.")
