import os
import sys
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import buildFERSBoards
from utils.dataloader import loadRDF, getRunInfo
from variables.fers import getFERSEnergySum, vectorizeFERS, calibrateFERSChannels, subtractFERSPedestal
from variables.drs import preProcessDRSBoards
from utils.html_generator import generate_html
from utils.fitter import eventFit
from utils.colors import colors
from configs.plotranges import getRangesForFERSEnergySums, getBoardEnergyFitParameters
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection, applyCC1Selection
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.timing import auto_timer

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
# file_pedestals = f"results/root/Run{runNumber}/valuemaps_pedestal.json"
file_pedestals_HG = f"results/root/Run{runNumber}/fers_pedestals_hg.json"
file_pedestals_LG = f"results/root/Run{runNumber}/fers_pedestals_lg.json"

pedestals_HG = json.load(open(file_pedestals_HG))
pedestals_LG = json.load(open(file_pedestals_LG))

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)
rdf, rdf_prefilter = vetoMuonCounter(rdf, TSmin=400, TSmax=700, cut=-30)

rdf, rdf_filterveto = applyUpstreamVeto(rdf, runNumber)
# rdf, rdf_psd = PSDSelection(rdf, runNumber, isHadron=False)

fersboards = buildFERSBoards(run=runNumber)

rdf = vectorizeFERS(rdf, fersboards)
# define energy sums with different configurations
# rdf = calibrateFERSChannels(
#    rdf, fersboards, file_gains=file_gains, file_pedestals=file_pedestals)
rdf = subtractFERSPedestal(
    rdf, fersboards, pedestals_HG, pedestalsLG=pedestals_LG)
rdf = getFERSEnergySum(
    rdf, fersboards, pdsub=True, calib=False)
# rdf = getFERSEnergySum(
#    rdf, fersboards, pdsub=True, calib=False, clip=False)
# rdf = getFERSEnergySum(
#    rdf, fersboards, pdsub=True, calib=True, clip=False)
# rdf = getFERSEnergySum(
#    rdf, fersboards, pdsub=True, calib=True, clip=True)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

# study PSD and CC1 selections
rdfs = OrderedDict()
rdf = applyPSDSelection(rdf, runNumber, applyCut=False)
rdf = applyCC1Selection(rdf, runNumber, applyCut=False)

rdfs["inc"] = rdf
rdfs["passPSDEle_passCC1Ele"] = rdf.Filter(
    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 1")
rdfs["passPSDEle_failCC1Ele"] = rdf.Filter(
    "pass_PSDEle_selection == 1 && pass_CC1Ele_selection == 0")
rdfs["failPSDEle_passCC1Ele"] = rdf.Filter(
    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 1")
rdfs["failPSDEle_failCC1Ele"] = rdf.Filter(
    "pass_PSDEle_selection == 0 && pass_CC1Ele_selection == 0")


def makeFERSEnergySumHists(rdf=rdf, pdsub=False, calib=False, clip=False, suffix=""):
    suffix_type, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _, xmin_LG_board, xmax_LG_board, xmax_LG_board_cer, xmin_LG_total, xmax_LG_total, xmax_LG_total_cer, _ = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)
    suffix += suffix_type

    hists_FERS_EnergySum = []
    for fersboard in fersboards.values():
        varname = fersboard.GetEnergySumName(
            useHG=True, isCer=True, pdsub=pdsub, calib=calib)
        hist_CerEnergyHG_Board = rdf.Histo1D((
            f"hist_{varname}",
            f"hist_{varname}",
            500, xmin_board, xmax_board_cer),
            varname
        )
        varname = fersboard.GetEnergySumName(
            useHG=False, isCer=True, pdsub=pdsub, calib=calib)
        hist_CerEnergyLG_Board = rdf.Histo1D((
            f"hist_{varname}",
            f"hist_{varname}",
            500, xmin_LG_board, xmax_LG_board_cer),
            varname
        )
        varname = fersboard.GetEnergySumName(
            useHG=True, isCer=False, pdsub=pdsub, calib=calib)
        hist_SciEnergyHG_Board = rdf.Histo1D((
            f"hist_{varname}",
            f"hist_{varname}",
            500, xmin_board, xmax_board),
            varname
        )
        varname = fersboard.GetEnergySumName(
            useHG=False, isCer=False, pdsub=pdsub, calib=calib)
        hist_SciEnergyLG_Board = rdf.Histo1D((
            f"hist_{varname}",
            f"hist_{varname}",
            500, xmin_LG_board, xmax_LG_board),
            varname
        )
        hists_FERS_EnergySum.append(hist_CerEnergyHG_Board)
        hists_FERS_EnergySum.append(hist_CerEnergyLG_Board)
        hists_FERS_EnergySum.append(hist_SciEnergyHG_Board)
        hists_FERS_EnergySum.append(hist_SciEnergyLG_Board)

    # per-event energy sum
    varname = fersboards.GetEnergySumName(
        useHG=True, isCer=True, pdsub=pdsub, calib=calib)
    hist_CerEnergyHG = rdf.Histo1D((
        f"hist_{varname}",
        f"hist_{varname}",
        500, xmin_total, xmax_total_cer),
        varname
    )
    varname = fersboards.GetEnergySumName(
        useHG=False, isCer=True, pdsub=pdsub, calib=calib)
    hist_CerEnergyLG = rdf.Histo1D((
        f"hist_{varname}",
        f"hist_{varname}",
        500, xmin_LG_total, xmax_LG_total_cer),
        varname
    )
    varname = fersboards.GetEnergySumName(
        useHG=True, isCer=False, pdsub=pdsub, calib=calib)
    hist_SciEnergyHG = rdf.Histo1D((
        f"hist_{varname}",
        f"hist_{varname}",
        500, xmin_total, xmax_total),
        varname
    )
    varname = fersboards.GetEnergySumName(
        useHG=False, isCer=False, pdsub=pdsub, calib=calib)
    hist_SciEnergyLG = rdf.Histo1D((
        f"hist_{varname}",
        f"hist_{varname}",
        500, xmin_LG_total, xmax_LG_total),
        varname
    )
    hists_FERS_EnergySum.append(hist_CerEnergyHG)
    hists_FERS_EnergySum.append(hist_CerEnergyLG)
    hists_FERS_EnergySum.append(hist_SciEnergyHG)
    hists_FERS_EnergySum.append(hist_SciEnergyLG)

    return hists_FERS_EnergySum


def makeFERSCervsSciHists(rdf=rdf, pdsub=False, calib=False, clip=False, suffix=""):
    suffix_type, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _, xmin_LG_board, xmax_LG_board, xmax_LG_board_cer, xmin_LG_total, xmax_LG_total, xmax_LG_total_cer, _ = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)
    suffix += suffix_type

    hists_FERS_Cer_vs_Sci = []
    for fersboard in fersboards.values():
        var_cer = fersboard.GetEnergySumName(
            useHG=True, isCer=True, pdsub=pdsub, calib=calib)
        var_sci = fersboard.GetEnergySumName(
            useHG=True, isCer=False, pdsub=pdsub, calib=calib)
        hist_Cer_vs_Sci = rdf.Histo2D((
            f"hist_{var_cer}_vs_{var_sci}",
            f"hist_{var_cer}_vs_{var_sci}",
            500, xmin_board, xmax_board,
            500, xmin_board, xmax_board_cer),
            var_sci,
            var_cer
        )
        hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci)

        var_cer = fersboard.GetEnergySumName(
            useHG=False, isCer=True, pdsub=pdsub, calib=calib)
        var_sci = fersboard.GetEnergySumName(
            useHG=False, isCer=False, pdsub=pdsub, calib=calib)
        hist_Cer_vs_Sci_LG = rdf.Histo2D((
            f"hist_{var_cer}_vs_{var_sci}",
            f"hist_{var_cer}_vs_{var_sci}",
            500, xmin_LG_board, xmax_LG_board,
            500, xmin_LG_board, xmax_LG_board_cer),
            var_sci,
            var_cer
        )
        hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci_LG)

    # total CER vs SCI energy plot
    var_cer = fersboards.GetEnergySumName(
        useHG=True, isCer=True, pdsub=pdsub, calib=calib)
    var_sci = fersboards.GetEnergySumName(
        useHG=True, isCer=False, pdsub=pdsub, calib=calib)
    hist_Cer_vs_Sci_Total = rdf.Histo2D((
        f"hist_{var_cer}_vs_{var_sci}",
        f"hist_{var_cer}_vs_{var_sci}",
        500, xmin_total, xmax_total,
        500, xmin_total, xmax_total_cer),
        var_sci,
        var_cer
    )
    var_cer = fersboards.GetEnergySumName(
        useHG=False, isCer=True, pdsub=pdsub, calib=calib)
    var_sci = fersboards.GetEnergySumName(
        useHG=False, isCer=False, pdsub=pdsub, calib=calib)
    hist_Cer_vs_Sci_LG_Total = rdf.Histo2D((
        f"hist_{var_cer}_vs_{var_sci}",
        f"hist_{var_cer}_vs_{var_sci}",
        500, xmin_LG_total, xmax_LG_total,
        500, xmin_LG_total, xmax_LG_total_cer),
        var_sci,
        var_cer
    )
    hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci_Total)
    hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci_LG_Total)

    return hists_FERS_Cer_vs_Sci


def makeFERSEnergySumPlots(pdsub=False, calib=False, clip=False, suffix=""):
    suffix_type, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, xtitle, xmin_LG_board, xmax_LG_board, xmax_LG_board_cer, xmin_LG_total, xmax_LG_total, xmax_LG_total_cer, xtitle_LG = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)
    suffix += suffix_type

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergySum_{suffix}"
    hists_CerEnergyHG = []
    hists_SciEnergyHG = []
    hists_CerEnergyLG = []
    hists_SciEnergyLG = []
    legends = []
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        varname = fersboard.GetEnergySumName(
            useHG=True, isCer=True, pdsub=pdsub, calib=calib)
        hist_CerEnergyHG = infile.Get(f"hist_{varname}")
        varname = fersboard.GetEnergySumName(
            useHG=False, isCer=True, pdsub=pdsub, calib=calib)
        hist_CerEnergyLG = infile.Get(f"hist_{varname}")
        varname = fersboard.GetEnergySumName(
            useHG=True, isCer=False, pdsub=pdsub, calib=calib)
        hist_SciEnergyHG = infile.Get(f"hist_{varname}")
        varname = fersboard.GetEnergySumName(
            useHG=False, isCer=False, pdsub=pdsub, calib=calib)
        hist_SciEnergyLG = infile.Get(f"hist_{varname}")

        if not (hist_CerEnergyHG and hist_SciEnergyHG and hist_CerEnergyLG and hist_SciEnergyLG):
            # if not (hist_CerEnergyHG and hist_SciEnergyHG):
            print(
                f"Warning: Histogram(s) for board {boardNo} not found in {infile_name}")
            continue

        legends.append(str(boardNo))
        hists_CerEnergyHG.append(hist_CerEnergyHG)
        hists_SciEnergyHG.append(hist_SciEnergyHG)
        hists_CerEnergyLG.append(hist_CerEnergyLG)
        hists_SciEnergyLG.append(hist_SciEnergyLG)

    output_name = "FERS_CerEnergyHG" + suffix
    DrawHistos(hists_CerEnergyHG, legends, xmin_board, xmax_board_cer, f"Cer {xtitle}", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_SciEnergyHG" + suffix
    DrawHistos(hists_SciEnergyHG, legends, xmin_board, xmax_board, f"Sci {xtitle}", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_CerEnergyLG" + suffix
    DrawHistos(hists_CerEnergyLG, legends, xmin_LG_board, xmax_LG_board_cer, f"Cer {xtitle_LG}", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_SciEnergyLG" + suffix
    DrawHistos(hists_SciEnergyLG, legends, xmin_LG_board, xmax_LG_board, f"Sci {xtitle_LG}", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")

    # total energy sum plots
    varname = fersboards.GetEnergySumName(
        useHG=True, isCer=True, pdsub=pdsub, calib=calib)
    hist_CerEnergyHG = infile.Get(f"hist_{varname}")
    varname = fersboards.GetEnergySumName(
        useHG=False, isCer=True, pdsub=pdsub, calib=calib)
    hist_CerEnergyLG = infile.Get(f"hist_{varname}")
    varname = fersboards.GetEnergySumName(
        useHG=True, isCer=False, pdsub=pdsub, calib=calib)
    hist_SciEnergyHG = infile.Get(f"hist_{varname}")
    varname = fersboards.GetEnergySumName(
        useHG=False, isCer=False, pdsub=pdsub, calib=calib)
    hist_SciEnergyLG = infile.Get(f"hist_{varname}")

    ymax = None  # let plotter decide the ymax
    output_name = "FERS_Total_CerEnergyHG" + suffix
    DrawHistos([hist_CerEnergyHG], "", xmin_total, xmax_total_cer, f"Cer {xtitle}", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(0, output_name + ".png")
    output_name = "FERS_Total_SciEnergyHG" + suffix
    DrawHistos([hist_SciEnergyHG], "", xmin_total, xmax_total, f"Sci {xtitle}", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(1, output_name + ".png")
    output_name = "FERS_Total_CerEnergyLG" + suffix
    DrawHistos([hist_CerEnergyLG], "", xmin_LG_total, xmax_LG_total_cer, f"Cer {xtitle_LG}", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(2, output_name + ".png")
    output_name = "FERS_Total_SciEnergyLG" + suffix
    DrawHistos([hist_SciEnergyLG], "", xmin_LG_total, xmax_LG_total, f"Sci {xtitle_LG}", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(3, output_name + ".png")

    output_html = f"{htmldir}/FERS_EnergySum{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSCerVsSciPlots(pdsub=False, calib=False, clip=False, suffix=""):
    suffix_type, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, xtitle, xmin_LG_board, xmax_LG_board, xmax_LG_board_cer, xmin_LG_total, xmax_LG_total, xmax_LG_total_cer, xtitle_LG = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)
    suffix += suffix_type

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci_{suffix}"
    # for fersboard in fersboards.values():
    #    boardNo = fersboard.boardNo
    #    hist_Cer_vs_Sci_name = f"hist_FERS_Board{boardNo}_Cer_vs_Sci_HG{suffix}"
    #    hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
    #    if not hist_Cer_vs_Sci:
    #        print(
    #            f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
    #        continue

    #    output_name = f"FERS_Board{boardNo}_Cer_vs_Sci{suffix}"
    #    DrawHistos([hist_Cer_vs_Sci], "", xmin_board, xmax_board, f"Sci {xtitle}", xmin_board, xmax_board_cer, f"Cer {xtitle}",
    #               output_name,
    #               dology=False, drawoptions=["colz"],
    #               outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
    #    plots.append(output_name + ".png")

    #    # LG
    #    hist_Cer_vs_Sci_LG_name = f"hist_FERS_Board{boardNo}_Cer_vs_Sci_LG{suffix}"
    #    hist_Cer_vs_Sci_LG = infile.Get(hist_Cer_vs_Sci_LG_name)
    #    if not hist_Cer_vs_Sci_LG:
    #        print(
    #            f"Warning: Histogram {hist_Cer_vs_Sci_LG_name} not found in {infile_name}")
    #        continue

    #    output_name = f"FERS_Board{boardNo}_Cer_vs_Sci_LG{suffix}"
    #    DrawHistos([hist_Cer_vs_Sci_LG], "", xmin_LG_board, xmax_LG_board, f"Sci {xtitle_LG}", xmin_LG_board, xmax_LG_board_cer, f"Cer {xtitle_LG}",
    #               output_name,
    #               dology=False, drawoptions=["colz"],
    #               outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
    #    plots.append(output_name + ".png")

    # total CER vs SCI energy plot
    var_cer = fersboards.GetEnergySumName(
        useHG=True, isCer=True, pdsub=pdsub, calib=calib)
    var_sci = fersboards.GetEnergySumName(
        useHG=True, isCer=False, pdsub=pdsub, calib=calib)
    hist_Cer_vs_Sci_Total = infile.Get(f"hist_{var_cer}_vs_{var_sci}")
    if not hist_Cer_vs_Sci_Total:
        print(
            f"Warning: Histogram hist_{var_cer}_vs_{var_sci} not found in {infile_name}")
    else:
        output_name = "FERS_Total_Cer_vs_Sci" + suffix
        DrawHistos([hist_Cer_vs_Sci_Total], "", xmin_total, xmax_total, f"Sci {xtitle}", xmin_total, xmax_total_cer, f"Cer {xtitle}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_Cer_vs_Sci{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)

    # LG
    var_cer = fersboards.GetEnergySumName(
        useHG=False, isCer=True, pdsub=pdsub, calib=calib)
    var_sci = fersboards.GetEnergySumName(
        useHG=False, isCer=False, pdsub=pdsub, calib=calib)
    hist_Cer_vs_Sci_LG_Total = infile.Get(f"hist_{var_cer}_vs_{var_sci}")
    if not hist_Cer_vs_Sci_LG_Total:
        print(
            f"Warning: Histogram hist_{var_cer}_vs_{var_sci}_LG{suffix} not found in {infile_name}")
    else:
        output_name = "FERS_Total_Cer_vs_Sci_LG" + suffix
        DrawHistos([hist_Cer_vs_Sci_LG_Total], "", xmin_LG_total, xmax_LG_total, f"Sci {xtitle_LG}", xmin_LG_total, xmax_LG_total_cer, f"Cer {xtitle_LG}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_Cer_vs_Sci{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
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
            hists_raw = makeFERSEnergySumHists(
                rdf=rdf, pdsub=True, calib=False, clip=False, suffix=cat)
            # hists_subtracted = makeFERSEnergySumHists(
            #    pdsub=True, calib=False, clip=False)
            # hists_subtracted_calibd = makeFERSEnergySumHists(
            #    pdsub=True, calib=True, clip=False)
            # hists_subtracted_calibd_clipped = makeFERSEnergySumHists(
            #    pdsub=True, calib=True, clip=True)

            hists_cer_vs_sci_raw = makeFERSCervsSciHists(rdf=rdf,
                                                         pdsub=True, calib=False, clip=False, suffix=cat)
            # hists_cer_vs_sci_subtracted = makeFERSCervsSciHists(
            #    pdsub=True, calib=False, clip=False)
            # hists_cer_vs_sci_subtracted_calibd = makeFERSCervsSciHists(
            #    pdsub=True, calib=True, clip=False)
            # hists_cer_vs_sci_subtracted_calibd_clipped = makeFERSCervsSciHists(
            #    pdsub=True, calib=True, clip=True)

            # save histograms to ROOT files
            outfile_name = f"{rootdir}/fers_energy_sum_{cat}_pdsub.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_raw:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved raw histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_subtracted.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_subtracted:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(f"Saved subtracted histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_subtracted_calibd.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_subtracted_calibd:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(f"Saved subtracted and calibd histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_subtracted_calibd_clipped.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_subtracted_calibd_clipped:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(
            #    f"Saved subtracted, calibd, and clipped histograms to {outfile_name}")

            # save cer vs sci histograms

            outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{cat}_pdsub.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_cer_vs_sci_raw:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved CER vs SCI histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_subtracted.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_cer_vs_sci_subtracted:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(f"Saved subtracted CER vs SCI histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_subtracted_calibd.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_cer_vs_sci_subtracted_calibd:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(
            #    f"Saved subtracted and calibd CER vs SCI histograms to {outfile_name}")

            # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_subtracted_calibd_clipped.root"
            # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
            #    for hist in hists_cer_vs_sci_subtracted_calibd_clipped:
            #        hist.SetDirectory(outfile)
            #        hist.Write()
            # print(
            #    f"Saved subtracted, calibd, and clipped CER vs SCI histograms to {outfile_name}")

        # make plots
        if makePlots:
            outputs_html[f"raw_{cat}"] = makeFERSEnergySumPlots(
                pdsub=True, calib=False, suffix=cat)
            # outputs_html["subtracted"] = makeFERSEnergySumPlots(
            #    pdsub=True, calib=False)
            # outputs_html["subtracted_calibd"] = makeFERSEnergySumPlots(
            #     pdsub=True, calib=True)
            # outputs_html["subtracted_calibd_clipped"] = makeFERSEnergySumPlots(
            #    pdsub=True, calib=True, clip=True)

            outputs_html[f"cer_vs_sci_raw_{cat}"] = makeFERSCerVsSciPlots(
                pdsub=True, calib=False, suffix=cat)
            # outputs_html["cer_vs_sci_subtracted"] = makeFERSCerVsSciPlots(
            #    pdsub=True, calib=False)
            # outputs_html["cer_vs_sci_subtracted_calibd"] = makeFERSCerVsSciPlots(
            #     pdsub=True, calib=True)
            # outputs_html["cer_vs_sci_subtracted_calibd_clipped"] = makeFERSCerVsSciPlots(
            #    pdsub=True, calib=True, clip=True)

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")
