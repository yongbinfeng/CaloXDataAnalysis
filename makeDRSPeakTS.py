import os
import sys
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import buildFERSBoards, buildDRSBoards, getMCPChannels
from utils.dataloader import loadRDF, getRunInfo
from variables.drs import preProcessDRSBoards, calibrateDRSPeakTS
from utils.visualization import visualizeFERSBoards
from utils.html_generator import generate_html
from configs.plotranges import getRangesForFERSEnergySums, getBoardEnergyFitParameters
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection, applyCC1Selection
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos, LHistos2Hist
from utils.timing import auto_timer
from utils.utils import number2string

ROOT.TH1.AddDirectory(False)  # prevents auto-registration in gDirectory


auto_timer("Total Execution Time")

args = get_args()
runNumber = args.run
firstEvent = args.first_event
lastEvent = args.last_event
btype, benergy = getRunInfo(runNumber)

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)
# rdf, rdf_prefilter = vetoMuonCounter(rdf, TSmin=400, TSmax=700, cut=-30)
# rdf, rdf_filterveto = applyUpstreamVeto(rdf, runNumber, applyCut=False)
rdf = applyUpstreamVeto(rdf, runNumber, applyCut=False)

DRSBoards = buildDRSBoards(run=runNumber)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

TSmin = -90
TSmax = -10

varsuffix = "DiffRelPeakTS_US"

if not os.path.exists(rootdir):
    os.makedirs(rootdir)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
if not os.path.exists(htmldir):
    os.makedirs(htmldir)

rdfs = OrderedDict()
rdf = applyPSDSelection(rdf, runNumber, applyCut=False)
rdf = applyCC1Selection(rdf, runNumber, applyCut=False)

rdf_prefilter = rdf
rdf = rdf_prefilter.Filter(
    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 1")

# rdfs["passPSDEle_failCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 0")
# rdfs["failPSDEle_passCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 1")
# rdfs["failPSDEle_failCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 0")


rdf = calibrateDRSPeakTS(rdf, runNumber, DRSBoards,
                         TSminDRS=0, TSmaxDRS=1000, threshold=9.0)

rdf_prefilter2 = rdf
map_mcp_channels = getMCPChannels(runNumber)

condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
condition += f" && {map_mcp_channels['US'][0]}_PeakTS > 500 && {map_mcp_channels['US'][0]}_PeakTS < 600"
rdf = rdf_prefilter2.Filter(condition,
                            "Pre-filter on MCP US channel 0 Peak TS")
# rdf = rdf.Filter(f"{map_mcp_channels['US'][0]}_RelPeakTS > 240 && {map_mcp_channels['US'][0]}_RelPeakTS < 280",
#                 "Pre-filter on MCP US channel 0 Peak TS")


def checkDRSPeakTS():
    h1s_DRSPeakTS = {}
    h1s_DRSPeakTS["Cer"] = []
    h1s_DRSPeakTS["Sci"] = []
    h2s_DRSPeakTS_Cer_VS_Sci = []

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        if boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.GetChannelName(blsub=False)
                channelNames[var] = channelName

                h1_DRSPeakTS = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak TS for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
                    TSmax - TSmin, TSmin, TSmax),
                    f"{channelName}_{varsuffix}"
                )
                h1s_DRSPeakTS[var].append(h1_DRSPeakTS)

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough channels found for Board{boardNo}, Tower({sTowerX}, {sTowerY})")
                continue

            h2_DRSPeak_Cer_VS_Sci = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}",
                f"DRS Peak TS - CER VS SCI for Board{boardNo}, Tower({sTowerX}, {sTowerY});SCI Peak TS;CER Peak TS",
                TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
                f'{channelNames["Sci"]}_{varsuffix}',
                f'{channelNames["Cer"]}_{varsuffix}',
            )
            h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_VS_Sci)
    return h1s_DRSPeakTS["Cer"], h1s_DRSPeakTS["Sci"], h2s_DRSPeakTS_Cer_VS_Sci


hists1d_DRSPeakTS_Cer, hists1d_DRSPeakTS_Sci, hists2d_DRSPeakTS_Cer_VS_Sci = checkDRSPeakTS()
outfile_DRSPeakTS = ROOT.TFile(f"{rootdir}/drspeakts_rel_us.root", "RECREATE")
for hist in hists1d_DRSPeakTS_Cer:
    hist.SetDirectory(outfile_DRSPeakTS)
    hist.Write()
for hist in hists1d_DRSPeakTS_Sci:
    hist.SetDirectory(outfile_DRSPeakTS)
    hist.Write()
for hist in hists2d_DRSPeakTS_Cer_VS_Sci:
    hist.SetDirectory(outfile_DRSPeakTS)
    hist.Write()
outfile_DRSPeakTS.Close()

# make plots


def makeDRSPeakTSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS_rel_us"
    infile_name = f"{rootdir}/drspeakts_rel_us.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_Cer = []
    hists_Cer_Quartz = []
    hists_Cer_Plastic = []
    hists_Sci = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        if boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hists = {}
            channelNos = {}
            isQuartz = False
            colors = [6, 4]
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if var == "Cer" and chan.isQuartz:
                    colors = [2, 4]
                    isQuartz = True

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    hists[var] = None
                    channelNos[var] = "N/A"
                else:
                    hists[var] = hist
                    channelNos[var] = "(" + str(chan.boardNo) + "," + \
                        str(chan.groupNo) + "," + str(chan.channelNo) + ")"

                output_name = f"hist_DRSPeakTS_{sTowerX}_{sTowerY}"

            if not hists["Cer"] or not hists["Sci"]:
                print(
                    f"Warning: Histograms for Cer or Sci not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY})")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.70, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"Tower: ({iTowerX}, {iTowerY})")
            extraToDraw.AddText(
                f"Cer: {channelNos['Cer']}, {hists['Cer'].GetBinCenter(hists['Cer'].GetMaximumBin())} TS")
            extraToDraw.AddText(
                f"Sci: {channelNos['Sci']}, {hists['Sci'].GetBinCenter(hists['Sci'].GetMaximumBin())} TS")
            if isQuartz:
                extraToDraw.AddText("Cer: Quartz")
            else:
                extraToDraw.AddText("Cer: Plastic")

            if isQuartz:
                hists_Cer_Quartz.append(hists["Cer"])
            else:
                hists_Cer_Plastic.append(hists["Cer"])
            hists_Cer.append(hists["Cer"])
            hists_Sci.append(hists["Sci"])

            DrawHistos([hists["Cer"], hists["Sci"]], ["Cer", "Sci"], TSmin, TSmax, "Peak TS", 1, hists["Cer"].GetMaximum() * 1.5, "Counts",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=colors, addOverflow=False, addUnderflow=False, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    # summary plots
    hist_Cer_Combined = LHistos2Hist(hists_Cer, "hist_DRSPeakTS_Cer_Combined")
    hist_Sci_Combined = LHistos2Hist(hists_Sci, "hist_DRSPeakTS_Sci_Combined")
    hist_Cer_Quartz_Combined = LHistos2Hist(
        hists_Cer_Quartz, "hist_DRSPeakTS_Cer_Quartz_Combined")
    hist_Cer_Plastic_Combined = LHistos2Hist(
        hists_Cer_Plastic, "hist_DRSPeakTS_Cer_Plastic_Combined")
    output_name = "DRS_PeakTS_Combined"
    extraToDraw = ROOT.TPaveText(0.50, 0.65, 0.90, 0.90, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(
        f"Quartz: {hist_Cer_Quartz_Combined.GetBinCenter(hist_Cer_Quartz_Combined.GetMaximumBin())} TS")
    extraToDraw.AddText(
        f"Plastic: {hist_Cer_Plastic_Combined.GetBinCenter(hist_Cer_Plastic_Combined.GetMaximumBin())} TS")
    extraToDraw.AddText(
        f"Sci: {hist_Sci_Combined.GetBinCenter(hist_Sci_Combined.GetMaximumBin())} TS")
    DrawHistos([hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined, hist_Sci_Combined], ["Cer Quartz", "Cer Plastic", "Sci"], TSmin, TSmax, "Peak TS", 1, None, "Counts",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 6, 4], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.25, 0.75, 0.40, 0.90], extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    plots.insert(1, "NEWLINE")

    output_html = f"{htmldir}/DRS/DRSPeakTS_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSCerVSSciPlots():
    plots = []
    hists = []
    hists_Quartz = []
    hists_Plastic = []
    outdir_plots = f"{plotdir}/DRSPeakTSCerVSSci_rel_us"
    infile_name = f"{rootdir}/drspeakts_rel_us.root"
    infile = ROOT.TFile(infile_name, "READ")

    # Create a dashed diagonal line from (0,0) to (1000,1000)
    diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
    diagonal_line.SetLineStyle(2)  # 2 = dashed
    diagonal_line.SetLineWidth(1)
    diagonal_line.SetLineColor(ROOT.kRed)
    extraToDraw = diagonal_line

    for _, DRSBoard in DRSBoards.items():
        if DRSBoard.boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            isQuartz = False

            chan_cer = DRSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True)
            if chan_cer and chan_cer.isQuartz:
                isQuartz = True

            hist_name = f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}"
            hist = infile.Get(hist_name)
            output_name = hist_name[5:]
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

            ytitle = "Cer (Quartz)" if isQuartz else "Cer (Plastic)"
            ytitle += " Peak TS"

            DrawHistos([hist], "", TSmin, TSmax, "Sci Peak TS", TSmin, TSmax, ytitle,
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)

            plots.append(output_name + ".png")

    # summary plots
    hcombined = LHistos2Hist(hists, "hist_DRSPeakTS_Cer_VS_Sci_Combined")
    output_name = "DRS_PeakTS_Cer_VS_Sci_Combined"
    DrawHistos([hcombined], "", TSmin, TSmax, "Sci Peak TS", TSmin, TSmax, f"Cer Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    hcombined_Quartz = LHistos2Hist(
        hists_Quartz, "hist_DRSPeakTS_Cer_VS_Sci_Quartz_Combined")
    output_name = "DRS_PeakTS_Cer_Quartz_VS_Sci_Combined"
    DrawHistos([hcombined_Quartz], "", TSmin, TSmax, "Sci Peak TS", TSmin, TSmax, f"Cer (Quartz) Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(1, output_name + ".png")

    hcombined_Plastic = LHistos2Hist(
        hists_Plastic, "hist_DRSPeakTS_Cer_VS_Sci_Plastic_Combined")
    output_name = "DRS_PeakTS_Cer_Plastic_VS_Sci_Combined"
    DrawHistos([hcombined_Plastic], "", TSmin, TSmax, "Sci Peak TS", TSmin, TSmax, f"Cer (Plastic) Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(2, output_name + ".png")
    plots.insert(3, "NEWLINE")

    print("quartz:", len(hists_Quartz), "plastic:", len(hists_Plastic))

    output_html = f"{htmldir}/DRS/DRSPeakTS_Cer_VS_Sci_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


output_html_DRSPeakTS = makeDRSPeakTSPlots()
output_html_DRSPeakTSCerVSSci = makeDRSPeakTSCerVSSciPlots()
print(f"DRS Peak TS plots saved to {output_html_DRSPeakTS}")
print(f"DRS Peak TS Cer VS Sci plots saved to {output_html_DRSPeakTSCerVSSci}")
