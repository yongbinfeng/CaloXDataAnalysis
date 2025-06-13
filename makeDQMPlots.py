import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map_new import buildDRSBoards, buildFERSBoards
from utils.utils import number2string
from utils.html_generator import generate_html

print("Start running script")
ROOT.gROOT.SetBatch(True)

runNumber = 583

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)


# 1D FERS histograms
rootdir = f"root/Run{runNumber}/"
infile_name = f"{rootdir}fers_all_channels_1D.root"
infile = ROOT.TFile(infile_name, "READ")

outdir = f"plots/Run{runNumber}/"

plots = []
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
        outdir_plots = outdir + "/FERS_1D"
        DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                   output_name,
                   dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, extraToDraw=extraToDraw,
                   outdir=outdir_plots)

        plots.append(output_name + ".png")

generate_html(plots, outdir_plots, output_html=f"html/FERS_1D/viewer.html")

# 1D histograms for DRS variables
plots = []
infile_name = f"{rootdir}/drs_all_channels_1D.root"
infile = ROOT.TFile(infile_name, "READ")
for _, DRSBoard in DRSBoards.items():
    boardNo = DRSBoard.boardNo
    for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
        sTowerX = number2string(iTowerX)
        sTowerY = number2string(iTowerY)

        chan_Cer = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True)
        chan_Sci = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=False)

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
            f"Cer Channel: {chan_Cer.channelNo}")
        extraToDraw.AddText(
            f"Sci Channel: {chan_Sci.channelNo}")

        output_name = f"DRS_Variable_Board{boardNo}_iTowerX{sTowerX}_iTowerY{sTowerY}"
        outdir_plots = outdir + "/DRS_1D"
        DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 1400, 2500, "DRS Output", 1, 1e12, "Counts",
                   output_name,
                   dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, extraToDraw=extraToDraw,
                   legendPos=(0.60, 0.78, 0.90, 0.68),
                   outdir=outdir_plots)
        plots.append(output_name + ".png")

generate_html(plots, outdir_plots, output_html=f"html/DRS_1D/viewer.html")
