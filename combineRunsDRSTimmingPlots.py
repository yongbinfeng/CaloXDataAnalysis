import os
import sys
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import buildFERSBoards, buildDRSBoards, getMCPChannels
from utils.dataloader import loadRDF, getRunInfo
from variables.drs import preProcessDRSBoards, calibrateDRSPeakTS
from utils.html_generator import generate_html
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
jsonFile = args.json_file
btype, benergy = getRunInfo(runNumber)

file_drschannels_bad = "data/drs/badchannels.json"
with open(file_drschannels_bad, "r") as f:
    drschannels_bad = json.load(f)

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions


DRSBoards = buildDRSBoards(run=runNumber)


def collectDRSvsTSProfPlots(runNumber):
    rootdir = f"results/root/Run{runNumber}/"
    infile_name = f"{rootdir}/drs_vs_ts_calibrated.root"
    infile = ROOT.TFile(infile_name, "READ")
    hprofs_Cer_Quartz = []
    hprofs_Cer_Plastic = []
    hprofs_Cer = []
    hprofs_Sci = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        if boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hprofs = {}
            isQuartz = False
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hprof_name = f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}"
                hprof = infile.Get(hprof_name)

                if var == "Cer" and chan.isQuartz:
                    isQuartz = True

                if not hprof:
                    print(
                        f"Warning: Histogram {hprof_name} not found in {infile_name}")
                    hprofs[var] = None
                else:
                    hprofs[var] = hprof.ProjectionX()

            if not hprofs["Cer"] or not hprofs["Sci"]:
                print(
                    f"Warning: Histograms for Cer or Sci not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY})")
                continue

            if isQuartz:
                hprofs_Cer_Quartz.append(hprofs["Cer"])
            else:
                hprofs_Cer_Plastic.append(hprofs["Cer"])
            hprofs_Cer.append(hprofs["Cer"])
            hprofs_Sci.append(hprofs["Sci"])

    hprof_Cer_Quartz_Combined = LHistos2Hist(
        hprofs_Cer_Quartz, "prof_DRS_vs_TS_Cer_Quartz_Combined")
    hprof_Cer_Plastic_Combined = LHistos2Hist(
        hprofs_Cer_Plastic, "prof_DRS_vs_TS_Cer_Plastic_Combined")
    return hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined


def collectPeakTS(hists, colors):
    extraToDraws = []

    for idx, (hist, col) in enumerate(zip(hists, colors)):
        peakBin = hist.GetMaximumBin()
        peakTS = hist.GetXaxis().GetBinCenter(peakBin)

        extraToDraw = ROOT.TPaveText(
            0.60, 0.85 - idx * 0.05, 0.90, 0.90 - idx * 0.05, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.SetTextColor(col)
        extraToDraw.AddText(f"Peak TS: {peakTS:.1f}")
        extraToDraws.append(extraToDraw)
    return extraToDraws


def compareRuns(runLists, xmin_quartz=-80, xmax_quartz=-56, xmin_plastic=-80, xmax_plastic=-48, xmin_sci=-80, xmax_sci=30, ymax_quartz=None, ymax_plastic=None, ymax_sci=None, donormalize=True):
    suffix = "Runs_" + "_".join([str(r) for r in runLists.keys()])
    outdir_plots = f"results/plots/CombineRuns/DRS_vs_TS/"
    htmldir = "results/html/CombineRuns/"
    os.makedirs(outdir_plots, exist_ok=True)
    plots = []

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta,
              ROOT.kCyan+2, ROOT.kOrange+7, ROOT.kGray+2, ROOT.kPink+9]
    hists_quartz = OrderedDict()
    hists_plastic = OrderedDict()
    hists_sci = OrderedDict()
    legends = []
    for runNumber in runLists.keys():
        # hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined = collectDRSvsTSProfPlots(
        #    runNumber)
        indir_root = f"results/root/Run{runNumber}/"
        infile_name = f"{indir_root}/drs_vs_ts_calibrated_combined.root"
        infile = ROOT.TFile(infile_name, "READ")
        hprof_Cer_Quartz_Combined = infile.Get(
            "prof_DRS_vs_TS_Cer_Quartz_Combined")
        hprof_Cer_Plastic_Combined = infile.Get(
            "prof_DRS_vs_TS_Cer_Plastic_Combined")
        hprof_Sci_Combined = infile.Get(
            "prof_DRS_vs_TS_Sci_Combined")
        int_quartz = hprof_Cer_Quartz_Combined.Integral(
            hprof_Cer_Quartz_Combined.FindBin(xmin_quartz), hprof_Cer_Quartz_Combined.FindBin(xmax_quartz))
        int_plastic = hprof_Cer_Plastic_Combined.Integral(
            hprof_Cer_Plastic_Combined.FindBin(xmin_plastic), hprof_Cer_Plastic_Combined.FindBin(xmax_plastic))
        int_sci = hprof_Sci_Combined.Integral(
            hprof_Sci_Combined.FindBin(xmin_sci), hprof_Sci_Combined.FindBin(xmax_sci))
        if donormalize:
            if int_quartz > 0:
                hprof_Cer_Quartz_Combined.Scale(1./int_quartz)
            if int_plastic > 0:
                hprof_Cer_Plastic_Combined.Scale(1./int_plastic)
            if int_sci > 0:
                hprof_Sci_Combined.Scale(1./int_sci)

        hists_quartz[f"{runNumber}"] = hprof_Cer_Quartz_Combined
        hists_plastic[f"{runNumber}"] = hprof_Cer_Plastic_Combined
        hists_sci[f"{runNumber}"] = hprof_Sci_Combined
        legends.append(runLists[runNumber])

    legends_ymin = 0.90 - 0.05*len(legends)

    DrawHistos(hists_quartz.values(), legends, xmin_quartz, xmax_quartz, "Time slice (Quartz)", 0, ymax_quartz, "A.U.",
               f"DRS_vs_TS_Cer_Quartz_{suffix}",
               dology=False, drawoptions=["hist,C"]*len(hists_quartz), mycolors=colors[:len(hists_quartz)], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=None, legendPos=[0.30, legends_ymin, 0.40, 0.90], legendoptions=["L"]*len(hists_quartz), extraToDraw=collectPeakTS(hists_quartz.values(), colors))
    plots.append(f"DRS_vs_TS_Cer_Quartz_{suffix}.png")
    DrawHistos(hists_plastic.values(), legends, xmin_plastic, xmax_plastic, "Time slice (Plastic)", 0, ymax_plastic, "A.U.",
               f"DRS_vs_TS_Cer_Plastic_{suffix}",
               dology=False, drawoptions=["hist,C"]*len(hists_plastic), mycolors=colors[:len(hists_plastic)], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=None, legendPos=[0.30, legends_ymin, 0.40, 0.90], legendoptions=["L"]*len(hists_plastic), extraToDraw=collectPeakTS(hists_plastic.values(), colors))
    plots.append(f"DRS_vs_TS_Cer_Plastic_{suffix}.png")
    DrawHistos(hists_sci.values(), legends, xmin_sci, xmax_sci, "Time slice (Sci)", 0, ymax_sci, "A.U.",
               f"DRS_vs_TS_Sci_{suffix}",
               dology=False, drawoptions=["hist,C"]*len(hists_sci), mycolors=colors[:len(hists_sci)], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=None, legendPos=[0.30, legends_ymin, 0.40, 0.90], legendoptions=["L"]*len(hists_sci), extraToDraw=collectPeakTS(hists_sci.values(), colors))
    plots.append(f"DRS_vs_TS_Sci_{suffix}.png")

    output_html = f"{htmldir}/DRS/DRSPeakTS_relative_to_MCP_{suffix}.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)


if __name__ == "__main__":
    runLists = OrderedDict()
    runLists["1409"] = "e^{+}, 80GeV"
    runLists["1465"] = "#pi^{+}, 80GeV"
    compareRuns(runLists)

    # runLists = OrderedDict()
    # runLists["1416"] = "e^{+}, 30GeV"
    # runLists["1412"] = "e^{+}, 60GeV"
    # runLists["1410"] = "e^{+}, 100GeV"
    # runLists["1411"] = "e^{+}, 120GeV"
    # compareRuns(runLists, xmax_quartz=-55)

    # runLists = OrderedDict()
    # runLists["1442"] = "#pi^{+}, 10GeV"
    # runLists["1441"] = "#pi^{+}, 20GeV"
    # runLists["1439"] = "#pi^{+}, 30GeV"
    # runLists["1437"] = "#pi^{+}, 60GeV"
    # compareRuns(runLists, xmax_quartz=-55)

    # runLists = OrderedDict()
    # runLists["1500"] = "e^{+}, 40GeV, -0.0cm"
    # runLists["1507"] = "e^{+}, 40GeV, -5.0cm"
    # runLists["1511"] = "e^{+}, 40GeV, -10.0cm"
    # compareRuns(runLists, xmax_quartz=-58, ymax_quartz=0.23,
    #            ymax_plastic=0.19, xmax_plastic=-55)
