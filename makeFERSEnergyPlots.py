import os
import sys
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import buildFERSBoards
from utils.dataloader import loadRDF, getRunInfo
from variables.fers import getFERSEnergySum, vectorizeFERS, calibrateFERSChannels, subtractFERSPedestal, getFERSEnergyWeightedCenter, addFERSPosXY, mixFERSHGLG
from variables.drs import preProcessDRSBoards
from utils.visualization import visualizeFERSBoards
from utils.html_generator import generate_html
from utils.fitter import eventFit
from utils.colors import colors
from configs.plotranges import getRangesForFERSEnergySums, getBoardEnergyFitParameters
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection, applyCC1Selection
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos, LHistos2Hist
from utils.timing import auto_timer

ROOT.TH1.AddDirectory(False)  # prevents auto-registration in gDirectory


auto_timer("Total Execution Time")

args = get_args()
runNumber = args.run
firstEvent = args.first_event
lastEvent = args.last_event
btype, benergy = getRunInfo(runNumber)

# HE = (runNumber >= 1200)
HE = (benergy >= 50)  # GeV
# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

# load the gains and pedestals from SiPM fits
# file_gains = f"results/root/Run{runNumber}/valuemaps_gain.json"
# file_pedestals = f"data/fers/FERS_pedestals_run1259.json"
file_pedestals_HG = f"results/root/Run{runNumber}/fers_pedestals_hg.json"
file_pedestals_LG = f"results/root/Run{runNumber}/fers_pedestals_lg.json"

file_calibrations = "data/fers/FERS_response.json"
file_HG2LG = "data/fers/FERS_HG2LG.json"
file_deadchannels = f"data/fers/deadchannels.json"

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)
# rdf, rdf_prefilter = vetoMuonCounter(rdf, TSmin=400, TSmax=700, cut=-30)
# rdf, rdf_filterveto = applyUpstreamVeto(rdf, runNumber, applyCut=False)
rdf = applyUpstreamVeto(rdf, runNumber, applyCut=False)

fersboards = buildFERSBoards(run=runNumber)

rdf = vectorizeFERS(rdf, fersboards)
# subtract pedestals
rdf = subtractFERSPedestal(
    rdf, fersboards, gain="HG", file_pedestals=file_pedestals_HG)
rdf = subtractFERSPedestal(
    rdf, fersboards, gain="LG", file_pedestals=file_pedestals_LG)
# mix HG and LG
# rdf = mixFERSHGLG(
#    rdf, fersboards, file_HG2LG=file_HG2LG)
# calibrate Mix gain
# rdf = calibrateFERSChannels(
#    rdf, fersboards, file_calibrations=file_calibrations, gain="Mix", file_deadchannels=file_deadchannels)

# GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]
GainCalibs = [("HG", False), ("LG", False)]

# calculate energy sums
for gain, calib in GainCalibs:
    rdf = getFERSEnergySum(rdf, fersboards, pdsub=True,
                           calib=calib, gain=gain)
    rdf = getFERSEnergyWeightedCenter(
        rdf, fersboards, pdsub=True, calib=calib, gain=gain)

rdf = addFERSPosXY(rdf, fersboards)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

if not os.path.exists(rootdir):
    os.makedirs(rootdir)
if not os.path.exists(plotdir):
    os.makedirs(plotdir)
if not os.path.exists(htmldir):
    os.makedirs(htmldir)

# study PSD and CC1 selections
rdfs = OrderedDict()
rdf = applyPSDSelection(rdf, runNumber, applyCut=False)
rdf = applyCC1Selection(rdf, runNumber, applyCut=False)

rdfs["inc"] = rdf
rdfs["passHaloVeto"] = rdf.Filter("pass_upstream_veto == 1")
rdfs["passHaloVeto_passPSDEle_passCC1Ele"] = rdf.Filter(
    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 1 && pass_upstream_veto == 1")
# rdfs["passPSDEle_failCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 0")
# rdfs["failPSDEle_passCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 1")
# rdfs["failPSDEle_failCC1Ele"] = rdf.Filter(
#    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 0")


def makeFERSEnergySumHists(rdf=rdf, suffix=""):
    hists_FERS_EnergySum = []
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        for cat in ["cer", "sci"]:
            # per-board sum
            for fersboard in fersboards.values():
                varname = fersboard.GetEnergySumName(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                hist = rdf.Histo1D((
                    f"hist_{varname}_{suffix}",
                    f"hist_{varname}_{suffix}",
                    500, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"]),
                    varname
                )
                hists_FERS_EnergySum.append(hist)
            # per-event sum
            varname = fersboards.GetEnergySumName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
            hist = rdf.Histo1D((
                f"hist_{varname}_{suffix}",
                f"hist_{varname}_{suffix}",
                500, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                varname
            )
            hists_FERS_EnergySum.append(hist)

    return hists_FERS_EnergySum


def makeFERSCervsSciHists(rdf=rdf, suffix=""):
    hists_FERS_Cer_vs_Sci = []
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        # per-board Cer vs Sci
        for fersboard in fersboards.values():
            var_cer = fersboard.GetEnergySumName(
                gain=gain, isCer=True, pdsub=True, calib=calib)
            var_sci = fersboard.GetEnergySumName(
                gain=gain, isCer=False, pdsub=True, calib=calib)
            hist_Cer_vs_Sci = rdf.Histo2D((
                f"hist_{var_cer}_VS_{var_sci}_{suffix}",
                f"hist_{var_cer}_VS_{var_sci}_{suffix}",
                500, config["xmin_board"][f"{gain}_sci"], config["xmax_board"][f"{gain}_sci"],
                500, config["xmin_board"][f"{gain}_cer"], config["xmax_board"][f"{gain}_cer"]),
                var_sci,
                var_cer
            )
            hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci)

        # per-event Cer vs Sci
        var_cer = fersboards.GetEnergySumName(
            gain=gain, isCer=True, pdsub=True, calib=calib)
        var_sci = fersboards.GetEnergySumName(
            gain=gain, isCer=False, pdsub=True, calib=calib)
        hist_Cer_vs_Sci_Event = rdf.Histo2D((
            f"hist_{var_cer}_VS_{var_sci}_{suffix}",
            f"hist_{var_cer}_VS_{var_sci}_{suffix}",
            500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
            500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
            var_sci,
            var_cer
        )
        hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci_Event)

    return hists_FERS_Cer_vs_Sci


def makeFERSEnergyWeightedCenterHists(rdf=rdf, suffix=""):
    hists_FERS_EnergyWeightedCenter = []
    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            # per-event
            varname_X = fersboards.GetEnergyWeightedCenterName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.GetEnergyWeightedCenterName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)

            histX = rdf.Histo1D((
                f"hist_{varname_X}_{suffix}",
                f"hist_{varname_X}_{suffix}",
                300, -15, 15),
                varname_X
            )
            hists_FERS_EnergyWeightedCenter.append(histX)
            histY = rdf.Histo1D((
                f"hist_{varname_Y}_{suffix}",
                f"hist_{varname_Y}_{suffix}",
                300, -15, 15),
                varname_Y
            )
            hists_FERS_EnergyWeightedCenter.append(histY)
            hist2D = rdf.Histo2D((
                f"hist_{varname_X}_VS_{varname_Y}_{suffix}",
                f"hist_{varname_X}_VS_{varname_Y}_{suffix}",
                300, -15, 15,
                300, -15, 15),
                varname_X,
                varname_Y
            )
            hists_FERS_EnergyWeightedCenter.append(hist2D)

    return hists_FERS_EnergyWeightedCenter


def makeFERSShowerShapeHists(rdf=rdf, suffix=""):
    hists_X = []
    hists_Y = []
    hists_R = []
    hists_Y_VS_X = []

    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            hists_tmp_X = []
            hists_tmp_Y = []
            hists_tmp_R = []
            hists_tmp_Y_VS_X = []
            for fersboard in fersboards.values():
                for channel in fersboard.GetListOfChannels(isCer=(cat == "cer")):
                    channelName = channel.GetChannelName(
                        gain=gain, pdsub=True, calib=calib)
                    hist = rdf.Histo1D((
                        f"hist_RealX_{channelName}_{suffix}",
                        f"hist_RealX_{channelName}_{suffix}",
                        100, -20, 20),
                        channel.GetRealPosName(isX=True), channelName
                    )
                    hists_tmp_X.append(hist)
                    hist = rdf.Histo1D((
                        f"hist_RealY_{channelName}_{suffix}",
                        f"hist_RealY_{channelName}_{suffix}",
                        100, -20, 20),
                        channel.GetRealPosName(isX=False), channelName
                    )
                    hists_tmp_Y.append(hist)
                    # calculate radius
                    rdf = rdf.Define(
                        "RealR_" + channelName, f"std::sqrt(std::pow({channel.GetRealPosName(isX=True)}, 2) + std::pow({channel.GetRealPosName(isX=False)}, 2))")
                    hist = rdf.Histo1D((
                        f"hist_RealR_{channelName}_{suffix}",
                        f"hist_RealR_{channelName}_{suffix}",
                        25, 0, 25),
                        "RealR_" + channelName, channelName
                    )
                    hists_tmp_R.append(hist)
                    hist = rdf.Histo2D((
                        f"hist_RealY_VS_RealX_{channelName}_{suffix}",
                        f"hist_RealY_VS_RealX_{channelName}_{suffix}",
                        100, -20, 20,
                        100, -20, 20),
                        channel.GetRealPosName(isX=True), channel.GetRealPosName(
                            isX=False), channelName
                    )
                    hists_tmp_Y_VS_X.append(hist)

            hists_X.append(
                (hists_tmp_X, f"hist_RealX_{gain}_{cat}_{suffix}"))
            hists_Y.append(
                (hists_tmp_Y, f"hist_RealY_{gain}_{cat}_{suffix}"))
            hists_R.append(
                (hists_tmp_R, f"hist_RealR_{gain}_{cat}_{suffix}"))
            hists_Y_VS_X.append(
                (hists_tmp_Y_VS_X, f"hist_RealY_VS_RealX_{gain}_{cat}_{suffix}"))

    nEvents = rdf.Count()

    return hists_X, hists_Y,  hists_R, hists_Y_VS_X, nEvents


def collectFERSStats(rdf):
    stats = {}
    # mean
    # and how frequent the saturation value is reached
    for fersboard in fersboards.values():
        for iTowerX, iTowerY in fersboard.GetListOfTowers():
            channel_cer = fersboard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            channel_sci = fersboard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)

            for gain, calib in GainCalibs:
                varname_cer = channel_cer.GetChannelName(
                    gain=gain, pdsub=True, calib=calib)
                varname_sci = channel_sci.GetChannelName(
                    gain=gain, pdsub=True, calib=calib)
                # rdf = rdf.Define(
                #    f"{varname_cer}_over_{varname_sci}", f"{varname_cer} / ({varname_sci} + 1e-6)")

                stats[channel_cer.GetChannelName(gain=gain)] = rdf.Mean(
                    f"{varname_cer}")
                stats[channel_sci.GetChannelName(gain=gain)] = rdf.Mean(
                    f"{varname_sci}")

    return stats


def makeFERSEnergySumPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergySum_{suffix}"

    # per-board energy sum plots
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        for cat in ["cer", "sci"]:
            hists = []
            legends = []
            for fersboard in fersboards.values():
                boardNo = fersboard.boardNo
                varname = fersboard.GetEnergySumName(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                hist_name = f"hist_{varname}_{suffix}"
                hist = infile.Get(hist_name)
                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} for board {boardNo} not found in {infile_name}")
                    continue
                legends.append(str(boardNo))
                hists.append(hist)

            output_name = f"FERS_Boards_{gain}_{cat}{suffix}"
            DrawHistos(hists, legends, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"], f"{cat.capitalize()} {gain} {config[f'title_{gain}']}", 1, None, "Events",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
            plots.append(output_name + ".png")

    # per-event energy sum plot ranges
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        for cat in ["cer", "sci"]:
            varname = fersboards.GetEnergySumName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
            hist_name = f"hist_{varname}_{suffix}"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue
            output_name = f"FERS_Total_{gain}_{cat}{suffix}"
            DrawHistos([hist], "", config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"], f"{cat.capitalize()} {gain} {config[f'title_{gain}']}", 1, None, "Events",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/FERS_EnergySum{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=6,
                  output_html=output_html)
    return output_html


def makeFERSCerVsSciPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci_{suffix}"

    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        for fersboard in fersboards.values():
            boardNo = fersboard.boardNo
            var_cer = fersboard.GetEnergySumName(
                gain=gain, isCer=True, pdsub=True, calib=calib)
            var_sci = fersboard.GetEnergySumName(
                gain=gain, isCer=False, pdsub=True, calib=calib)
            hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
            hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
            if not hist_Cer_vs_Sci:
                print(
                    f"Warning: Histogram {hist_Cer_vs_Sci_name} for board {boardNo} not found in {infile_name}")
                continue

            output_name = f"FERS_Board{boardNo}_Cer_VS_Sci_{gain}{suffix}"
            DrawHistos([hist_Cer_vs_Sci], "", config["xmin_board"][f"{gain}_sci"], config["xmax_board"][f"{gain}_sci"], f"Sci {gain} {config[f'title_{gain}']}", config["xmin_board"][f"{gain}_cer"], config["xmax_board"][f"{gain}_cer"], f"Cer {gain} {config[f'title_{gain}']}",
                       output_name,
                       dology=False, drawoptions=["colz"],
                       outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
            plots.append(output_name + ".png")

    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, runNumber=runNumber)
        var_cer = fersboards.GetEnergySumName(
            gain=gain, isCer=True, pdsub=True, calib=calib)
        var_sci = fersboards.GetEnergySumName(
            gain=gain, isCer=False, pdsub=True, calib=calib)
        hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue
        output_name = f"FERS_Total_Cer_VS_Sci_{gain}{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], f"Sci {gain} {config[f'title_{gain}']}", config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], f"Cer {gain} {config[f'title_{gain}']}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/FERS_Cer_VS_Sci{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeFERSEnergyWeightedCenterPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_weighted_center_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergyWeightedCenter_{suffix}"

    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            varname_X = fersboards.GetEnergyWeightedCenterName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.GetEnergyWeightedCenterName(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)

            extraToDrawBase = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
            extraToDrawBase.SetTextAlign(11)
            extraToDrawBase.SetFillColorAlpha(0, 0)
            extraToDrawBase.SetBorderSize(0)
            extraToDrawBase.SetTextFont(42)
            extraToDrawBase.SetTextSize(0.04)

            histX_name = f"hist_{varname_X}_{suffix}"
            histX = infile.Get(histX_name)
            if histX:
                valcenter = histX.GetMean()
                rms = histX.GetRMS()
                extraToDraw = extraToDrawBase.Clone()
                extraToDraw.AddText(f"Center = {valcenter:.2f} +/- {rms:.2f}")
                output_name = f"FERS_Total_{gain}_{cat}_EWC_X{suffix}"
                DrawHistos([histX], "", -15, 15, f"{cat.capitalize()} {gain} EWC X [cm]", 1, None, "Events",
                           output_name,
                           dology=False, drawoptions="HIST", mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, runNumber=runNumber, extraToDraw=extraToDraw)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {histX_name} not found in {infile_name}")

            histY_name = f"hist_{varname_Y}_{suffix}"
            histY = infile.Get(histY_name)
            if histY:
                valcenter = histY.GetMean()
                rms = histY.GetRMS()
                extraToDraw = extraToDrawBase.Clone()
                extraToDraw.AddText(f"Center = {valcenter:.2f} +/- {rms:.2f}")
                output_name = f"FERS_Total_{gain}_{cat}_EWC_Y{suffix}"
                DrawHistos([histY], "", -15, 15, f"{cat.capitalize()} {gain} EWC Y [cm]", 1, None, "Events",
                           output_name,
                           dology=False, drawoptions="HIST", mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, runNumber=runNumber, extraToDraw=extraToDraw)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {histY_name} not found in {infile_name}")

            hist2D_name = f"hist_{varname_X}_VS_{varname_Y}_{suffix}"
            hist2D = infile.Get(hist2D_name)
            if hist2D:
                output_name = f"FERS_Total_{gain}_{cat}_EWC_X_vs_Y{suffix}"
                DrawHistos([hist2D], "", -15, 15, f"{cat.capitalize()} {gain} EWC X [cm]", -15, 15, f"{cat.capitalize()} {gain} EWC Y [cm]",
                           output_name,
                           dology=False, drawoptions=["colz"], addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {hist2D_name} not found in {infile_name}")

    output_html = f"{htmldir}/FERS_EnergyWeightedCenter{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeFERSShowerShapePlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_shower_shape_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_ShowerShape_{suffix}"

    for gain, calib in GainCalibs:
        hists_R = []
        for cat in ["cer", "sci"]:
            hist_X_name = f"hist_RealX_{gain}_{cat}_{suffix}"
            hist_X = infile.Get(hist_X_name)
            if hist_X:
                output_name = f"FERS_ShowerShape_RealX_{gain}_{cat}_{suffix}"
                DrawHistos([hist_X], "", -20, 20, f"X {cat.capitalize()} {gain} [cm]", 1e-4, 1, "Frac. of Energy",
                           output_name, drawoptions="HIST",
                           dology=True, mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, runNumber=runNumber, donormalize=True)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {hist_X_name} not found in {infile_name}")

            hist_Y_name = f"hist_RealY_{gain}_{cat}_{suffix}"
            hist_Y = infile.Get(hist_Y_name)
            if hist_Y:
                output_name = f"FERS_ShowerShape_RealY_{gain}_{cat}_{suffix}"
                DrawHistos([hist_Y], "", -20, 20, f"Y {cat.capitalize()} {gain} [cm]", 1e-4, 1, "Frac. of Energy",
                           output_name, drawoptions="HIST",
                           dology=True, mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, runNumber=runNumber, donormalize=True)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {hist_Y_name} not found in {infile_name}")

            hist_R_name = f"hist_RealR_{gain}_{cat}_{suffix}"
            hist_R = infile.Get(hist_R_name)
            if hist_R:
                hists_R.append(hist_R)
            else:
                print(
                    f"Warning: Histogram {hist_R_name} not found in {infile_name}")

            # hist_Y_VS_X_name = f"hist_RealY_VS_RealX_{gain}_{cat}_{suffix}"
            # hist_Y_VS_X = infile.Get(hist_Y_VS_X_name)
            # if hist_Y_VS_X:
            #    output_name = f"FERS_ShowerShape_RealY_VS_RealX_{gain}_{cat}_{suffix}"
            #    DrawHistos([hist_Y_VS_X], "", -20, 20, f"X {cat.capitalize()} {gain} [cm]", -20, 20, f"Y {cat.capitalize()} {gain} [cm]",
            #               output_name,
            #               dology=False, drawoptions=["colz"], addOverflow=True, addUnderflow=True, dologz=True,
            #               outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1e-4, zmax=1, donormalize=True, zlabel="Frac. of Energy")
            #    plots.append(output_name + ".png")
            # else:
            #    print(
            #        f"Warning: Histogram {hist_Y_VS_X_name} not found in {infile_name}")

        output_name = f"FERS_ShowerShape_RealR_{gain}_{suffix}"
        ymax = max([h.GetMaximum() for h in hists_R]) * 1.2 if hists_R else 1
        DrawHistos(hists_R, ["Cer", "Sci"], 0, 25, f"R [cm]", 1e-4, 1, "Frac. of Energy",
                   output_name,
                   dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
                   outdir=outdir_plots, runNumber=runNumber, donormalize=True)
        plots.append(output_name + ".png")

        if len(hists_R) == 2:
            # plot ratio of Cer/Sci
            hists_R[0].Scale(1.0 / hists_R[0].Integral(0,
                             hists_R[0].GetNbinsX() + 1))
            hists_R[1].Scale(1.0 / hists_R[1].Integral(0,
                             hists_R[1].GetNbinsX() + 1))
            hist_ratio = hists_R[1].Clone()
            hist_ratio.Divide(hists_R[0])
            output_name = f"FERS_ShowerShape_RealR_Cer_over_Sci_{gain}_{suffix}"
            DrawHistos([hist_ratio], "", 0, 25, f"R [cm]", 0, 1.5, "Cer/Sci",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=[1], addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_ShowerShape{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=6,
                  output_html=output_html)
    return output_html


def makeFERSStatsPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_Stats_Cer_ovs_Sci"
    # load the json file
    import json
    infile_name = f"{rootdir}/fers_stats_cer_ovs_sci.json"
    with open(infile_name, "r") as f:
        stats = json.load(f)

    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100

    for gain, calib in GainCalibs:
        [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm] = visualizeFERSBoards(
            fersboards, stats, suffix=f"Run{runNumber}_{gain}_mean", gain=gain)

        h2_Cer_Over_Sci = h2_Cer.Clone(f"h2_Cer_over_Sci_{gain}")
        h2_Cer_Over_Sci.Divide(h2_Sci)
        h2_Cer_3mm_Over_Sci = h2_Cer_3mm.Clone(f"h2_Cer_3mm_over_Sci_{gain}")
        h2_Cer_3mm_Over_Sci.Divide(h2_Sci_3mm)

        output_name = f"FERS_Boards_Run{runNumber}_Stats_Cer_{gain}"
        DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=None)
        plots.append(output_name + ".png")
        output_name = f"FERS_Boards_Run{runNumber}_Stats_Sci_{gain}"
        DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=None)
        plots.append(output_name + ".png")
        output_name = f"FERS_Boards_Run{runNumber}_Stats_Cer_over_Sci_{gain}"
        DrawHistos([h2_Cer_Over_Sci, h2_Cer_3mm_Over_Sci], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer / Sci", runNumber=runNumber, zmin=0, zmax=1.5)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_Stats_Cer_over_Sci/index.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeBoardFits():
    suffix = "subtracted_calibd"
    filename = f"{rootdir}/fers_energy_sum_{suffix}.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run makeFERSEnergyPlots.py with makeHists=True first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    plots = []
    outdir = f"{plotdir}/boardfits"
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        hCer = ifile.Get(f"hist_FERS_Board{boardNo}_CerEnergyHG_{suffix}")
        hSci = ifile.Get(f"hist_FERS_Board{boardNo}_SciEnergyHG_{suffix}")

        args_cer = getBoardEnergyFitParameters(
            runNumber, is3mm=fersboard.Is3mm(), isCer=True)
        args_sci = getBoardEnergyFitParameters(
            runNumber, is3mm=fersboard.Is3mm(), isCer=False)

        output_name = eventFit(hCer, f"Run{runNumber}_Board{boardNo}_CerHG",
                               outdir=outdir, xlabel="Cer # p.e.",
                               **args_cer)
        plots.append(output_name)
        output_name = eventFit(hSci, f"Run{runNumber}_Board{boardNo}_SciHG",
                               outdir=outdir, xlabel="Sci # p.e.",
                               **args_sci)
        plots.append(output_name)

    output_html = f"{htmldir}/boardfits/index.html"
    generate_html(plots, outdir, plots_per_row=2,
                  output_html=output_html)
    print(f"Generated HTML file: {output_html}")
    return output_html


if __name__ == "__main__":
    makeHists = True
    makePlots = True
    outputs_html = {}

    for cat, rdf in rdfs.items():
        if makeHists:
            hists_raw = makeFERSEnergySumHists(rdf=rdf, suffix=cat)

            hists_cer_vs_sci_raw = makeFERSCervsSciHists(rdf=rdf, suffix=cat)

            hists_energy_weighted_center = makeFERSEnergyWeightedCenterHists(
                rdf=rdf, suffix=cat)

            hists_shower_shapes_X, hists_shower_shapes_Y, hists_shower_shapes_R, hists_shower_shapes_Y_VS_X, nEvents = makeFERSShowerShapeHists(
                rdf=rdf, suffix=cat)

            # stats = collectFERSStats(rdf=rdf)

            # import json
            # stats_results = {}
            # for channelName, mean in stats.items():
            #     stats_results[channelName] = mean.GetValue()
            # with open(f"{rootdir}/fers_stats_cer_ovs_sci.json", "w") as f:
            #     json.dump(stats_results, f, indent=4)

            # save histograms to ROOT files
            outfile_name = f"{rootdir}/fers_energy_sum_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_raw:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved raw histograms to {outfile_name}")

            # save cer vs sci histograms

            outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_cer_vs_sci_raw:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved CER vs SCI histograms to {outfile_name}")

            # save energy weighted center histograms
            outfile_name = f"{rootdir}/fers_energy_weighted_center_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_energy_weighted_center:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved energy weighted center histograms to {outfile_name}")

            nEvts = nEvents.GetValue()
            # save shower shape histograms
            # process the shower shape histograms to merge per-channel histograms
            hists_shower_shape = []
            for hists_X, hname_X in hists_shower_shapes_X:
                hists_X_new = [h.GetValue() for h in hists_X]
                hist_combined = LHistos2Hist(hists_X_new, hname_X)
                # normalize to number of events
                hist_combined.Scale(1.0 / nEvts)
                hists_shower_shape.append(hist_combined)
            for hists_Y, hname_Y in hists_shower_shapes_Y:
                hists_Y_new = [h.GetValue() for h in hists_Y]
                hist_combined = LHistos2Hist(hists_Y_new, hname_Y)
                hist_combined.Scale(1.0 / nEvts)
                hists_shower_shape.append(hist_combined)
            for hists_R, hname_R in hists_shower_shapes_R:
                hists_R_new = [h.GetValue() for h in hists_R]
                hist_combined = LHistos2Hist(hists_R_new, hname_R)
                hist_combined.Scale(1.0 / nEvts)
                hists_shower_shape.append(hist_combined)
            for hists_Y_VS_X, hname_Y_VS_X in hists_shower_shapes_Y_VS_X:
                hists_Y_VS_X_new = [h.GetValue() for h in hists_Y_VS_X]
                hist_combined = LHistos2Hist(hists_Y_VS_X_new, hname_Y_VS_X)
                hist_combined.Scale(1.0 / nEvts)
                hists_shower_shape.append(hist_combined)
            outfile_name = f"{rootdir}/fers_shower_shape_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_shower_shape:
                    # Python won't try to delete
                    ROOT.SetOwnership(hist, False)
                    histC = hist.Clone()
                    histC.SetDirectory(outfile)
                    histC.Write()
            print(f"Saved shower shape histograms to {outfile_name}")

        # make plots
        if makePlots:
            outputs_html[f"raw_{cat}"] = makeFERSEnergySumPlots(suffix=cat)

            outputs_html[f"cer_vs_sci_raw_{cat}"] = makeFERSCerVsSciPlots(
                suffix=cat)

            outputs_html[f"energy_weighted_center_{cat}"] = makeFERSEnergyWeightedCenterPlots(
                suffix=cat)

            outputs_html[f"shower_shape_{cat}"] = makeFERSShowerShapePlots(
                suffix=cat)

            # outputs_html[f"stats_cer_ovs_sci"] = makeFERSStatsPlots()

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")
