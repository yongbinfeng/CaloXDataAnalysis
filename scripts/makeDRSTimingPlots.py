import os
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import buildFERSBoards, buildDRSBoards, getMCPChannels
from utils.dataloader import loadRDF, getRunInfo
from variables.drs import preProcessDRSBoards, calibrateDRSPeakTS
from utils.html_generator import generate_html
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection, applyCC1Selection
from utils.parser import get_args
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

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

TSmin = -90
TSmax = -10
TSCermin = -70
TSCermax = -50

varsuffix = "DiffRelPeakTS_US"

if not os.path.exists(rootdir):
    os.makedirs(rootdir)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
if not os.path.exists(htmldir):
    os.makedirs(htmldir)


def GetMean(hist, useMode=True):
    if useMode:
        # peak position of peakTS
        bin_max = hist.GetMaximumBin()
        mean = hist.GetBinCenter(bin_max)
    else:
        # average position of peakTS
        mean = hist.GetMean()
    return mean


def checkDRSPeakTS(rdf):
    h1s_DRSPeakTS = {}
    h1s_DRSPeakTS["Cer"] = []
    h1s_DRSPeakTS["Sci"] = []
    h2s_DRSPeakTS_Cer_VS_Sci = []

    channelnames_quartz = []
    channelnames_plastic = []
    channelnames_sci = []

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
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                h1_DRSPeakTS = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak TS for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
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
                    f"Warning: Not enough good channels found for Board{boardNo}, Tower({sTowerX}, {sTowerY})")
                continue

            h2_DRSPeak_Cer_VS_Sci = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}",
                f"DRS Peak TS - CER VS SCI for Board{boardNo}, Tower({sTowerX}, {sTowerY});SCI Peak TS;CER Peak TS",
                TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
                f'{channelNames["Sci"]}_{varsuffix}',
                f'{channelNames["Cer"]}_{varsuffix}',
            )
            h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_VS_Sci)

    # average of quartz, plastic and sci channels
    rdf = rdf.Define("Cer_Quartz_AvgPeakTS",
                     f"({'+'.join(channelnames_quartz)})/{len(channelnames_quartz)}")
    rdf = rdf.Define("Cer_Plastic_AvgPeakTS",
                     f"({'+'.join(channelnames_plastic)})/{len(channelnames_plastic)}")
    rdf = rdf.Define("Sci_AvgPeakTS",
                     f"({'+'.join(channelnames_sci)})/{len(channelnames_sci)}")
    h1_Cer_Quartz_AvgPeakTS = rdf.Histo1D((
        "hist_DRSPeakTS_Cer_Quartz_Avg",
        f"DRS Peak TS Avg for Cer Quartz;Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Cer_Quartz_AvgPeakTS"
    )
    h1s_DRSPeakTS["Cer"].append(h1_Cer_Quartz_AvgPeakTS)
    h1_Cer_Plastic_AvgPeakTS = rdf.Histo1D((
        "hist_DRSPeakTS_Cer_Plastic_Avg",
        f"DRS Peak TS Avg for Cer Plastic;Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Cer_Plastic_AvgPeakTS"
    )
    h1s_DRSPeakTS["Cer"].append(h1_Cer_Plastic_AvgPeakTS)
    h1_Sci_AvgPeakTS = rdf.Histo1D((
        "hist_DRSPeakTS_Sci_Avg",
        f"DRS Peak TS Avg for Sci;Peak TS;Counts",
        TSmax - TSmin, TSmin, TSmax),
        "Sci_AvgPeakTS"
    )
    h1s_DRSPeakTS["Sci"].append(h1_Sci_AvgPeakTS)
    # avg of quartz vs avg of plastic
    h2_DRSPeak_Cer_Quartz_VS_Plastic = rdf.Histo2D((
        "hist_DRSPeakTS_Cer_Quartz_VS_Cer_Plastic_Avg",
        f"DRS Peak TS - CER Quartz VS CER Plastic Avg;Cer Plastic Peak TS;Cer Quartz Peak TS",
        TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
        "Cer_Plastic_AvgPeakTS",
        "Cer_Quartz_AvgPeakTS",
    )
    h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_Quartz_VS_Plastic)
    return h1s_DRSPeakTS["Cer"], h1s_DRSPeakTS["Sci"], h2s_DRSPeakTS_Cer_VS_Sci


def checkDRSvsCalibrationTS(rdf):
    """
    shower shapes vs z (calibrated time)
    """
    hprofs_DRS_VS_TS = []
    hprofs_DRS_VS_Z = []
    hists2d_DRS_VS_TS = []
    rdfs_filtered = []

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
                if channelName in drschannels_bad:
                    print(
                        f"Warning: DRS Channel {channelName} is marked as bad channel, skipping...")
                    continue
                channelNames[var] = channelName

                rdf = rdf.Define(f"{channelName}_hasSignal",
                                 f"{channelName}_Sum > 100.0")

                rdf_filtered = rdf.Filter(f"{channelName}_hasSignal")

                hprof_DRS_VS_TS = rdf_filtered.Profile1D((
                    f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}", "", 400, -300, 100), f"{channelName}_AlignedTS", f"{channelName}_blsub")

                hist2d_DRS_VS_TS = rdf_filtered.Histo2D((
                    f"hist2d_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}", "", 400, -300, 100, 100, -200, 1000), f"{channelName}_AlignedTS", f"{channelName}_blsub")

                if var != "Sci":
                    # map aligned TS to real distance in z (cm)
                    # only for Cer (Quartz and Plastic) channels for now
                    hprof_DRS_VS_Z = rdf_filtered.Profile1D((
                        f"prof_DRS_vs_Z_{var}_{sTowerX}_{sTowerY}", "", 200, -100, 600), f"{channelName}_MeasuredZ", f"{channelName}_blsub")
                    hprofs_DRS_VS_Z.append(hprof_DRS_VS_Z)

                rdfs_filtered.append(rdf_filtered)

                hprofs_DRS_VS_TS.append(hprof_DRS_VS_TS)
                hists2d_DRS_VS_TS.append(hist2d_DRS_VS_TS)
    return hprofs_DRS_VS_TS, hists2d_DRS_VS_TS, hprofs_DRS_VS_Z, rdfs_filtered


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

            extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.70, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"Tower: ({iTowerX}, {iTowerY})")
            extraToDraw.AddText(
                f"Cer: {channelNos['Cer']}, {GetMean(hists['Cer']):.2f} TS")
            extraToDraw.AddText(
                f"Sci: {channelNos['Sci']}, {GetMean(hists['Sci']):.2f} TS")
            # if isQuartz:
            #    extraToDraw.AddText("Cer: Quartz")
            # else:
            #    extraToDraw.AddText("Cer: Plastic")

            if isQuartz:
                hists_Cer_Quartz.append(hists["Cer"])
            else:
                hists_Cer_Plastic.append(hists["Cer"])
            hists_Cer.append(hists["Cer"])
            hists_Sci.append(hists["Sci"])

            DrawHistos([hists["Cer"], hists["Sci"]], ["Cer", "Sci"], TSmin, TSmax, "Peak TS", 1, hists["Cer"].GetMaximum() * 1.5, "Counts",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=colors, addOverflow=False, addUnderflow=False, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber, extraText="Quartz + Sci" if isQuartz else "Plastic + Sci")
            plots.append(output_name + ".png")

    # summary plots
    hist_Cer_Combined = LHistos2Hist(hists_Cer, "hist_DRSPeakTS_Cer_Combined")
    hist_Sci_Combined = LHistos2Hist(hists_Sci, "hist_DRSPeakTS_Sci_Combined")
    hist_Cer_Quartz_Combined = LHistos2Hist(
        hists_Cer_Quartz, "hist_DRSPeakTS_Cer_Quartz_Combined")
    hist_Cer_Plastic_Combined = LHistos2Hist(
        hists_Cer_Plastic, "hist_DRSPeakTS_Cer_Plastic_Combined")
    output_name = "DRS_PeakTS_Combined"
    extraToDraw = ROOT.TPaveText(0.50, 0.75, 0.90, 0.90, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(
        f"Quartz: {GetMean(hist_Cer_Quartz_Combined):.2f} TS")
    extraToDraw.AddText(
        f"Plastic: {GetMean(hist_Cer_Plastic_Combined):.2f} TS")
    extraToDraw.AddText(
        f"Sci: {GetMean(hist_Sci_Combined):.2f} TS")
    DrawHistos([hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined, hist_Sci_Combined], ["Cer Quartz", "Cer Plastic", "Sci"], TSmin, TSmax, "Peak TS", 1, None, "Counts",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 6, 4], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.25, 0.75, 0.40, 0.90], extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    output_name = "DRS_PeakTS_Cer_Combined"
    DrawHistos([hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined], ["Cer Quartz", "Cer Plastic"], TSCermin, TSCermax, "Peak TS", 0, 0.15, "Fraction",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 6], addOverflow=False, addUnderflow=False, donormalize=True,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.30, 0.80, 0.40, 0.90])
    plots.insert(1, output_name + ".png")

    plots.insert(2, "NEWLINE")

    output_html = f"{htmldir}/DRS/DRSPeakTS_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)

    outfile_name = f"{rootdir}/drspeakts_rel_us_combined.root"
    outfile = ROOT.TFile(outfile_name, "RECREATE")
    for h in [hist_Cer_Combined, hist_Sci_Combined, hist_Cer_Quartz_Combined, hist_Cer_Plastic_Combined]:
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()
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

    # average of quartz vs average of plastic
    hcombined_Quartz_vs_Plastic = infile.Get(
        "hist_DRSPeakTS_Cer_Quartz_VS_Cer_Plastic_Avg")
    output_name = "DRS_PeakTS_Cer_Quartz_VS_Cer_Plastic_Avg"
    DrawHistos([hcombined_Quartz_vs_Plastic], "", TSmin, TSmax, "Cer Plastic Peak TS", TSmin, TSmax, f"Cer Quartz Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(3, output_name + ".png")

    print("quartz:", len(hists_Quartz), "plastic:", len(hists_Plastic))

    output_html = f"{htmldir}/DRS/DRSPeakTS_Cer_VS_Sci_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)

    return output_html


def makeDRSvsTSProfPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS_vs_ts_calibrated"
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
            colors = [6, 4]
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hprof_name = f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}"
                hprof = infile.Get(hprof_name)

                if var == "Cer" and chan.isQuartz:
                    colors = [2, 4]
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

            output_name = f"prof_DRS_vs_TS_{sTowerX}_{sTowerY}"

            DrawHistos([hprofs["Cer"], hprofs["Sci"]], ["Cer", "Sci"], -200, 100, "Time slice", -10, hprofs["Cer"].GetMaximum() * 1.5, "DRS ADC",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=colors, addOverflow=False, addUnderflow=False,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    # summary plots
    hprof_Cer_Combined = LHistos2Hist(
        hprofs_Cer, "prof_DRS_vs_TS_Cer_Combined")
    hprof_Sci_Combined = LHistos2Hist(
        hprofs_Sci, "prof_DRS_vs_TS_Sci_Combined")
    hprof_Cer_Quartz_Combined = LHistos2Hist(
        hprofs_Cer_Quartz, "prof_DRS_vs_TS_Cer_Quartz_Combined")
    hprof_Cer_Plastic_Combined = LHistos2Hist(
        hprofs_Cer_Plastic, "prof_DRS_vs_TS_Cer_Plastic_Combined")
    output_name = "DRS_vs_TS_Cer_Sci_Combined"
    hprof_Cer_Quartz_Combined.GetXaxis().SetRangeUser(-80, -55)
    DrawHistos([hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined, hprof_Sci_Combined], ["Cer Quartz", "Cer Plastic", "Sci"], -75, 0, "Time slice", 0, None, "DRS ADC",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 6, 4], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.25, 0.75, 0.40, 0.90])
    plots.insert(0, output_name + ".png")
    output_name = "DRS_vs_TS_Cer_Combined"
    DrawHistos([hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined], ["Cer Quartz", "Cer Plastic"], -80, -50, "Time slice", 0, None, "ADC",
               output_name,
               dology=False, drawoptions=["hist,C", "hist,C"], mycolors=[2, 6], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.30, 0.80, 0.40, 0.90], legendoptions=["L", "L"])
    plots.insert(1, output_name + ".png")

    plots.insert(2, "NEWLINE")

    output_html = f"{htmldir}/DRS/DRS_VS_TS_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)

    # save combined histograms to root file
    outfile_name = f"{rootdir}/drs_vs_ts_calibrated_combined.root"
    outfile = ROOT.TFile(outfile_name, "RECREATE")
    for h in [hprof_Cer_Combined, hprof_Sci_Combined, hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined]:
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()

    return output_html


def makeDRSvsZProfPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS_vs_ts_calibrated"
    infile_name = f"{rootdir}/drs_vs_ts_calibrated.root"
    infile = ROOT.TFile(infile_name, "READ")
    hprofs_Cer_Quartz = []
    hprofs_Cer_Plastic = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        if boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hprofs = {}
            isQuartz = False
            colors = [6, 4]
            for var in ["Cer"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hprof_name = f"prof_DRS_vs_Z_{var}_{sTowerX}_{sTowerY}"
                hprof = infile.Get(hprof_name)

                if var == "Cer" and chan.isQuartz:
                    colors = [2, 4]
                    isQuartz = True

                if not hprof:
                    print(
                        f"Warning: Histogram {hprof_name} not found in {infile_name}")
                    hprofs[var] = None
                else:
                    hprofs[var] = hprof.ProjectionX()

            if not hprofs["Cer"]:
                print(
                    f"Warning: Histograms for Cer not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY})")
                continue

            if isQuartz:
                hprofs_Cer_Quartz.append(hprofs["Cer"])
            else:
                hprofs_Cer_Plastic.append(hprofs["Cer"])

            output_name = f"prof_DRS_vs_Z_{sTowerX}_{sTowerY}"

            DrawHistos([hprofs["Cer"]], ["Cer"], -100, 300, "Measured Z [cm]", -10, hprofs["Cer"].GetMaximum() * 1.5, "DRS ADC",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=colors, addOverflow=False, addUnderflow=False,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    # summary plots
    hprof_Cer_Quartz_Combined = LHistos2Hist(
        hprofs_Cer_Quartz, "prof_DRS_vs_Z_Cer_Quartz_Combined")
    hprof_Cer_Plastic_Combined = LHistos2Hist(
        hprofs_Cer_Plastic, "prof_DRS_vs_Z_Cer_Plastic_Combined")
    output_name = "DRS_vs_Z_Cer_Combined"
    DrawHistos([hprof_Cer_Quartz_Combined, hprof_Cer_Plastic_Combined], ["Cer Quartz", "Cer Plastic"], -10, 200, "Measured Z", 0, None, "DRS ADC",
               output_name,
               dology=False, drawoptions=["C", "C"], mycolors=[2, 6], addOverflow=False, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.55, 0.80, 0.90, 0.90], legendoptions=["P", "P"])
    plots.insert(0, output_name + ".png")
    plots.insert(1, "NEWLINE")

    output_html = f"{htmldir}/DRS/DRS_VS_Z_relative_to_MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSvsTS2DPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS_vs_ts_calibrated"
    infile_name = f"{rootdir}/drs_vs_ts_calibrated.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        if boardNo > 3:
            continue
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist2d_name = f"hist2d_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}"
                hist2d = infile.Get(hist2d_name)

                output_name = f"hist2d_DRS_vs_TS_{sTowerX}_{sTowerY}_{var}"
                plots.append(output_name + ".png")

                if not hist2d:
                    print(
                        f"Warning: Histogram {hist2d_name} not found in {infile_name}")
                    continue

                DrawHistos([hist2d], "", -200, 100, "Calibrated TS", -200, 1000, "DRS ADC",
                           output_name,
                           dology=False, addOverflow=False, addUnderflow=False, doth2=True, drawoptions="COLZ", zmin=1, zmax=1e4, dologz=True,
                           outdir=outdir_plots, runNumber=runNumber)

    output_html = f"{htmldir}/DRS/DRS_VS_TS_relative_to_MCP_2D.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def main():
    makeHists = True
    makePlots = True

    if makeHists:
        rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent, jsonFile)
        rdf = preProcessDRSBoards(rdf, runNumber=runNumber)
        rdf, rdf_prefilter = vetoMuonCounter(
            rdf, runNumber, TSmin=200, TSmax=700, cut=-100)
        rdf = applyUpstreamVeto(rdf, runNumber, applyCut=False)

        # rdfs = OrderedDict()
        # rdf, _ = applyPSDSelection(rdf, runNumber, applyCut=True)
        # rdf, _ = applyCC1Selection(rdf, runNumber, applyCut=True)

        rdf = calibrateDRSPeakTS(rdf, runNumber, DRSBoards,
                                 TSminDRS=0, TSmaxDRS=1000, threshold=9.0)

        rdf_prefilterMCP1 = rdf
        map_mcp_channels = getMCPChannels(runNumber)

        condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
        condition += f" && {map_mcp_channels['US'][0]}_PeakTS > 500 && {map_mcp_channels['US'][0]}_PeakTS < 600"
        condition += f" && {map_mcp_channels['US'][0]}_Peak < -300.0"
        rdf_prefilterMCP2 = rdf_prefilterMCP1.Filter(condition,
                                                     "Pre-filter on MCP US channel 0 Peak TS")

        rdf_prefilterMCP2 = rdf_prefilterMCP1.Define(
            "MCP0_DeltaRelPeakTS", f"{map_mcp_channels['DS'][0]}_RelPeakTS - {map_mcp_channels['US'][0]}_RelPeakTS")

        # requirement on Downstream MCP and also the delta between US and DS
        # so that MCP timing is reliable (this is probably too tight; can be relaxed later)
        # condition = f"{map_mcp_channels['DS'][0]}_RelPeakTS > -350 && {map_mcp_channels['DS'][0]}_RelPeakTS < -100"
        # condition += f" && {map_mcp_channels['DS'][0]}_PeakTS > 500 && {map_mcp_channels['DS'][0]}_PeakTS < 600"
        # condition += f" && {map_mcp_channels['DS'][0]}_Peak < -300.0"
        # condition += f" && MCP0_DeltaRelPeakTS > 0 && MCP0_DeltaRelPeakTS < 6"
        # rdf = rdf_prefilterMCP2.Filter(condition,
        #                               "Filter on MCP DS channel 0 Peak TS and DeltaRelPeakTS")
        rdf = rdf_prefilterMCP2
        # rdf_old = rdf
        # event_ns = list(rdf_old.Take["unsigned int"]("event_n").GetValue())
        # print("event numbers in this run:", event_ns)
        # rdf = rdf_old.Filter(
        #    f"event_n == {event_ns[0]}", "Select first event only for DRS timing checks")

        hists1d_DRSPeakTS_Cer, hists1d_DRSPeakTS_Sci, hists2d_DRSPeakTS_Cer_VS_Sci = checkDRSPeakTS(
            rdf)
        hprofs_DRS_VS_TS, hists2d_DRS_VS_TS, hprofs_DRS_VS_Z, rdfs_filtered = checkDRSvsCalibrationTS(
            rdf)

        outfile_DRSPeakTS = ROOT.TFile(
            f"{rootdir}/drspeakts_rel_us.root", "RECREATE")
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

        outfile_DRS_VS_TS = ROOT.TFile(
            f"{rootdir}/drs_vs_ts_calibrated.root", "RECREATE")
        for hprof in hprofs_DRS_VS_TS:
            hprof.SetDirectory(outfile_DRS_VS_TS)
            hprof.Write()
        for hist2d in hists2d_DRS_VS_TS:
            hist2d.SetDirectory(outfile_DRS_VS_TS)
            hist2d.Write()
        for hprof in hprofs_DRS_VS_Z:
            hprof.SetDirectory(outfile_DRS_VS_TS)
            hprof.Write()
        outfile_DRS_VS_TS.Close()

    if makePlots:
        output_html_DRSPeakTS = makeDRSPeakTSPlots()
        output_html_DRSPeakTSCerVSSci = makeDRSPeakTSCerVSSciPlots()
        output_html_DRS_VS_TS_Prof = makeDRSvsTSProfPlots()
        # output_html_DRS_VS_TS_2D = makeDRSvsTS2DPlots()
        output_html_DRS_VS_Z_Prof = makeDRSvsZProfPlots()
        print(f"DRS Peak TS plots saved to {output_html_DRSPeakTS}")
        print(
            f"DRS Peak TS Cer VS Sci plots saved to {output_html_DRSPeakTSCerVSSci}")
        print(
            f"DRS VS TS profiled plots saved to {output_html_DRS_VS_TS_Prof}")
        # print(f"DRS VS TS 2D plots saved to {output_html_DRS_VS_TS_2D}")
        print(f"DRS VS Z profiled plots saved to {output_html_DRS_VS_Z_Prof}")


if __name__ == "__main__":
    main()
