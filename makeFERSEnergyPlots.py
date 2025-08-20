import os
import sys
import ROOT
from utils.channel_map import buildFERSBoards
from utils.utils import loadRDF, calculateEnergySumFERS, vectorizeFERS, calibrateFERSChannels, preProcessDRSBoards
from utils.html_generator import generate_html
from utils.fitter import eventFit
from utils.colors import colors
from configs.plotranges import getRangesForFERSEnergySums, getBoardEnergyFitParameters, getEventEnergyFitParameters
from selections.selections import vetoMuonCounter, applyUpstreamVeto, PSDSelection
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos


print("Start running makeFERSEnergyPlots.py")


args = get_args()
runNumber = args.run
firstEvent = args.first_event
lastEvent = args.last_event

HE = (runNumber >= 1200)
# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

# load the gains and pedestals from SiPM fits
file_gains = f"results/root/Run{runNumber}/valuemaps_gain.json"
file_pedestals = f"results/root/Run{runNumber}/valuemaps_pedestal.json"

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)
rdf, rdf_prefilter = vetoMuonCounter(rdf, TSmin=400, TSmax=700, cut=-30)

rdf, rdf_filterveto = applyUpstreamVeto(rdf, runNumber)
rdf, rdf_psd = PSDSelection(rdf, runNumber, isHadron=True)

FERSBoards = buildFERSBoards(run=runNumber)

rdf = vectorizeFERS(rdf, FERSBoards)
# define energy sums with different configurations
# rdf = calibrateFERSChannels(
#    rdf, FERSBoards, file_gains=file_gains, file_pedestals=file_pedestals)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=False, calibrate=False, clip=False)
# rdf = calculateEnergySumFERS(
#    rdf, FERSBoards, subtractPedestal=True, calibrate=False, clip=False)
# rdf = calculateEnergySumFERS(
#    rdf, FERSBoards, subtractPedestal=True, calibrate=True, clip=False)
# rdf = calculateEnergySumFERS(
#    rdf, FERSBoards, subtractPedestal=True, calibrate=True, clip=True)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

random_board_comparisons = [
    ("Board0_CerEnergyHG", "Board1_SciEnergyHG"),
    ("Board1_CerEnergyHG", "Board3_SciEnergyHG"),
    ("Board4_CerEnergyHG", "Board5_SciEnergyHG"),
    ("Board4_CerEnergyHG", "Board5_CerEnergyHG"),
    ("Board4_SciEnergyHG", "Board8_CerEnergyHG"),
    ("Board10_CerEnergyHG", "Board12_SciEnergyHG"),
    ("Board12_CerEnergyHG", "Board13_SciEnergyHG"),
]

random_channel_comparisons = []
for chan in range(64):
    random_channel_comparisons.append(
        ("Board4_energyHG_0", f"Board4_energyHG_{chan}"))
for chan in range(32):
    random_channel_comparisons.append(
        ("Board4_energyHG_32", f"Board5_energyHG_{chan}"))
for chan in range(32):
    random_channel_comparisons.append(
        ("Board3_energyHG_1", f"Board3_energyHG_{chan}"))


def makeFERSEnergySumHists(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _ = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    hists_FERS_EnergySum = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_CerEnergyHG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_CerEnergyHG{suffix}",
            f"FERS Board {boardNo} - CER Energy HG;CER Energy HG;Counts",
            500, xmin_board, xmax_board_cer),
            f"FERS_Board{boardNo}_CerEnergyHG{suffix}"
        )
        # hist_CerEnergyLG_Board = rdf.Histo1D((
        #    f"hist_FERS_Board{boardNo}_CerEnergyLG{suffix}",
        #    f"FERS Board {boardNo} - CER Energy LG;CER Energy LG;Counts",
        #    500, 0, 40000),
        #    f"FERS_Board{boardNo}_CerEnergyLG{suffix}"
        # )
        hist_SciEnergyHG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_SciEnergyHG{suffix}",
            f"FERS Board {boardNo} - SCI Energy HG;SCI Energy HG;Counts",
            500, xmin_board, xmax_board),
            f"FERS_Board{boardNo}_SciEnergyHG{suffix}"
        )
        # hist_SciEnergyLG_Board = rdf.Histo1D((
        #    f"hist_FERS_Board{boardNo}_SciEnergyLG{suffix}",
        #    f"FERS Board {boardNo} - SCI Energy LG;SCI Energy LG;Counts",
        #    500, 0, 40000),
        #    f"FERS_Board{boardNo}_SciEnergyLG{suffix}"
        # )
        hists_FERS_EnergySum.append(hist_CerEnergyHG_Board)
        # hists_FERS_EnergySum.append(hist_CerEnergyLG_Board)
        hists_FERS_EnergySum.append(hist_SciEnergyHG_Board)
        # hists_FERS_EnergySum.append(hist_SciEnergyLG_Board)

    # per-event energy sum
    hist_CerEnergyHG = rdf.Histo1D((
        f"hist_FERS_CerEnergyHG{suffix}",
        "FERS - CER Energy HG;CER Energy HG;Counts",
        500, xmin_total, xmax_total_cer),
        f"FERS_CerEnergyHG{suffix}"
    )
    # hist_CerEnergyLG = rdf.Histo1D((
    #     f"hist_FERS_CerEnergyLG{suffix}",
    #     "FERS - CER Energy LG;CER Energy LG;Counts",
    #     500, 0, 2e5),
    #     f"FERS_CerEnergyLG{suffix}"
    # )
    hist_SciEnergyHG = rdf.Histo1D((
        f"hist_FERS_SciEnergyHG{suffix}",
        "FERS - SCI Energy HG;SCI Energy HG;Counts",
        500, xmin_total, xmax_total),
        f"FERS_SciEnergyHG{suffix}"
    )
    # hist_SciEnergyLG = rdf.Histo1D((
    #     f"hist_FERS_SciEnergyLG{suffix}",
    #     "FERS - SCI Energy LG;SCI Energy LG;Counts",
    #     500, 0, 2e5),
    #     f"FERS_SciEnergyLG{suffix}"
    # )
    hists_FERS_EnergySum.append(hist_CerEnergyHG)
    # hists_FERS_EnergySum.append(hist_CerEnergyLG)
    hists_FERS_EnergySum.append(hist_SciEnergyHG)
    # hists_FERS_EnergySum.append(hist_SciEnergyLG)

    return hists_FERS_EnergySum


def makeFERSCervsSciHists(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _ = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    hists_FERS_Cer_vs_Sci = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_Cer_vs_Sci = rdf.Histo2D((
            f"hist_FERS_Board{boardNo}_Cer_vs_Sci{suffix}",
            f"FERS Board {boardNo} - CER vs SCI Energy",
            500, xmin_board, xmax_board,
            500, xmin_board, xmax_board_cer),
            f"FERS_Board{boardNo}_SciEnergyHG{suffix}",
            f"FERS_Board{boardNo}_CerEnergyHG{suffix}"
        )
        hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci)

    # total CER vs SCI energy plot
    hist_Cer_vs_Sci_Total = rdf.Histo2D((
        f"hist_FERS_Cer_vs_Sci{suffix}",
        "FERS - CER vs SCI Energy;CER Energy HG;SCI Energy HG;Counts",
        500, xmin_total, xmax_total,
        500, xmin_total, xmax_total_cer),
        f"FERS_SciEnergyHG{suffix}",
        f"FERS_CerEnergyHG{suffix}"
    )
    hists_FERS_Cer_vs_Sci.append(hist_Cer_vs_Sci_Total)

    return hists_FERS_Cer_vs_Sci


def makeFERSCervsSciRandomHists(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _ = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    hists_FERS_Cer_vs_Sci_random = []
    for key, value in random_board_comparisons:
        hist_Cer_vs_Sci = rdf.Histo2D((
            f"hist_FERS_{key}_vs_{value}{suffix}",
            f"FERS {key} vs {value}",
            500, xmin_board, xmax_board,
            500, xmin_board, xmax_board),
            f"FERS_{key}{suffix}",
            f"FERS_{value}{suffix}"
        )
        hists_FERS_Cer_vs_Sci_random.append(hist_Cer_vs_Sci)

    return hists_FERS_Cer_vs_Sci_random


def makeFERSChannelComparisonHists(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, _ = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    xmin_channel = xmin_board / 32
    xmax_channel = xmax_board / 32

    hists_FERS_ChannelComparison = []
    for key, value in random_channel_comparisons:
        hist_Cer_vs_Sci = rdf.Histo2D((
            f"hist_FERS_{key}_vs_{value}{suffix}",
            f"FERS {key} vs {value}",
            500, xmin_channel, xmax_channel,
            500, xmin_channel, xmax_channel),
            f"FERS_{key}{suffix}",
            f"FERS_{value}{suffix}"
        )
        hists_FERS_ChannelComparison.append(hist_Cer_vs_Sci)

    return hists_FERS_ChannelComparison


def makeFERSEnergySumPlots(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, xtitle = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergySum{suffix}"
    hists_CerEnergyHG = []
    hists_SciEnergyHG = []
    # hists_CerEnergyLG = []
    # hists_SciEnergyLG = []
    legends = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_CerEnergyHG_name = f"hist_FERS_Board{boardNo}_CerEnergyHG{suffix}"
        hist_SciEnergyHG_name = f"hist_FERS_Board{boardNo}_SciEnergyHG{suffix}"
        # hist_CerEnergyLG_name = f"hist_FERS_Board{boardNo}_CerEnergyLG{suffix}"
        # hist_SciEnergyLG_name = f"hist_FERS_Board{boardNo}_SciEnergyLG{suffix}"
        hist_CerEnergyHG = infile.Get(hist_CerEnergyHG_name)
        hist_SciEnergyHG = infile.Get(hist_SciEnergyHG_name)
        # hist_CerEnergyLG = infile.Get(hist_CerEnergyLG_name)
        # hist_SciEnergyLG = infile.Get(hist_SciEnergyLG_name)

        # if not (hist_CerEnergyHG and hist_SciEnergyHG and hist_CerEnergyLG and hist_SciEnergyLG):
        if not (hist_CerEnergyHG and hist_SciEnergyHG):
            print(
                f"Warning: Histograms {hist_CerEnergyHG_name}, {hist_SciEnergyHG_name} not found in {infile_name}")
            continue

        legends.append(str(boardNo))
        hists_CerEnergyHG.append(hist_CerEnergyHG)
        hists_SciEnergyHG.append(hist_SciEnergyHG)
        # hists_CerEnergyLG.append(hist_CerEnergyLG)
        # hists_SciEnergyLG.append(hist_SciEnergyLG)

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
    # output_name = "FERS_CerEnergyLG" + suffix
    # DrawHistos(hists_CerEnergyLG, legends, 0, xmax_total_cer, "Cer Energy LG", 0, None, "Events",
    #           output_name,
    #           dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
    #           outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.7, 0.90, 0.9])
    # plots.append(output_name + ".png")
    # output_name = "FERS_SciEnergyLG" + suffix
    # DrawHistos(hists_SciEnergyLG, legends, 0, xmax_total, "Sci Energy LG", 0, None, "Events",
    #           output_name,
    #           dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
    #           outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.7, 0.90, 0.9])
    # plots.append(output_name + ".png")

    # total energy sum plots
    hist_CerEnergyHG = infile.Get(f"hist_FERS_CerEnergyHG{suffix}")
    hist_SciEnergyHG = infile.Get(f"hist_FERS_SciEnergyHG{suffix}")
    # hist_CerEnergyLG = infile.Get(f"hist_FERS_CerEnergyLG{suffix}")
    # hist_SciEnergyLG = infile.Get(f"hist_FERS_SciEnergyLG{suffix}")
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
    # output_name = "FERS_Total_CerEnergyLG" + suffix
    # DrawHistos([hist_CerEnergyLG], "", 0, xmax_total_cer, "Cer Energy LG", 0, ymax, "Events",
    #           output_name,
    #           dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
    #           outdir=outdir_plots, runNumber=runNumber)
    # plots.insert(2, output_name + ".png")
    # output_name = "FERS_Total_SciEnergyLG" + suffix
    # DrawHistos([hist_SciEnergyLG], "", 0, xmax_total, "Sci Energy LG", 0, ymax, "Events",
    #           output_name,
    #           dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
    #           outdir=outdir_plots, runNumber=runNumber)
    # plots.insert(3, output_name + ".png")

    output_html = f"{htmldir}/FERS_EnergySum{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSCerVsSciPlots(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, xmin_total, xmax_total, xmax_total_cer, xtitle = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci{suffix}"
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_Cer_vs_Sci_name = f"hist_FERS_Board{boardNo}_Cer_vs_Sci{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue

        output_name = f"FERS_Board{boardNo}_Cer_vs_Sci{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", xmin_board, xmax_board, f"Sci {xtitle}", xmin_board, xmax_board_cer, f"Cer {xtitle}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.append(output_name + ".png")

    # total CER vs SCI energy plot
    hist_Cer_vs_Sci_Total = infile.Get(f"hist_FERS_Cer_vs_Sci{suffix}")
    if not hist_Cer_vs_Sci_Total:
        print(
            f"Warning: Histogram hist_FERS_Cer_vs_Sci{suffix} not found in {infile_name}")
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
    return output_html


def makeFERSCerVsSciRandomPlots(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, _, _, _, xtitle = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_random{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci_Random{suffix}"
    for key, value in random_board_comparisons:
        hist_Cer_vs_Sci_name = f"hist_FERS_{key}_vs_{value}{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue

        xtitle_temp = key.split(
            "_")[0] + (" Cer " if "Cer" in key else " Sci ") + xtitle
        ytitle_temp = value.split(
            "_")[0] + (" Cer " if "Cer" in value else " Sci ") + xtitle

        output_name = f"FERS_{key}_vs_{value}{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", xmin_board, xmax_board, xtitle_temp, xmin_board, xmax_board, ytitle_temp,
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_Cer_vs_Sci_Random{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSChannelComparisonPlots(subtractPedestal=False, calibrate=False, clip=False):
    suffix, xmin_board, xmax_board, xmax_board_cer, _, _, _, xtitle = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip, HE=HE)
    xmin_channel = xmin_board / 32
    xmax_channel = xmax_board / 32

    plots = []
    infile_name = f"{rootdir}/fers_energy_channel_comparison{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_ChannelComparison{suffix}"
    for key, value in random_channel_comparisons:
        hist_Cer_vs_Sci_name = f"hist_FERS_{key}_vs_{value}{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue

        xtitle_temp = key.split("_")[0] + " " + \
            key.split("_")[2] + " " + xtitle
        ytitle_temp = value.split(
            "_")[0] + " " + value.split("_")[2] + " " + xtitle

        output_name = f"FERS_{key}_vs_{value}{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", xmin_channel, xmax_channel, xtitle_temp, xmin_channel, xmax_channel, ytitle_temp,
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_ChannelComparison{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeBoardFits():
    suffix = "subtracted_calibrated"
    filename = f"{rootdir}/fers_energy_sum_{suffix}.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run makeFERSEnergyPlots.py with makeHists=True first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    plots = []
    outdir = f"{plotdir}/boardfits"
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hCer = ifile.Get(f"hist_FERS_Board{boardNo}_CerEnergyHG_{suffix}")
        hSci = ifile.Get(f"hist_FERS_Board{boardNo}_SciEnergyHG_{suffix}")

        args_cer = getBoardEnergyFitParameters(
            runNumber, is3mm=FERSBoard.Is3mm(), isCer=True)
        args_sci = getBoardEnergyFitParameters(
            runNumber, is3mm=FERSBoard.Is3mm(), isCer=False)

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

    if makeHists:
        hists_raw = makeFERSEnergySumHists(
            subtractPedestal=False, calibrate=False, clip=False)
        # hists_subtracted = makeFERSEnergySumHists(
        #    subtractPedestal=True, calibrate=False, clip=False)
        # hists_subtracted_calibrated = makeFERSEnergySumHists(
        #    subtractPedestal=True, calibrate=True, clip=False)
        # hists_subtracted_calibrated_clipped = makeFERSEnergySumHists(
        #    subtractPedestal=True, calibrate=True, clip=True)

        hists_cer_vs_sci_raw = makeFERSCervsSciHists(
            subtractPedestal=False, calibrate=False, clip=False)
        # hists_cer_vs_sci_subtracted = makeFERSCervsSciHists(
        #    subtractPedestal=True, calibrate=False, clip=False)
        # hists_cer_vs_sci_subtracted_calibrated = makeFERSCervsSciHists(
        #    subtractPedestal=True, calibrate=True, clip=False)
        # hists_cer_vs_sci_subtracted_calibrated_clipped = makeFERSCervsSciHists(
        #    subtractPedestal=True, calibrate=True, clip=True)

        # hists_cer_vs_sci_random_raw = makeFERSCervsSciRandomHists(
        #    subtractPedestal=False, calibrate=False, clip=False)
        # hists_cer_vs_sci_random_subtracted = makeFERSCervsSciRandomHists(
        #    subtractPedestal=True, calibrate=False, clip=False)
        # hists_cer_vs_sci_random_subtracted_calibrated = makeFERSCervsSciRandomHists(
        #    subtractPedestal=True, calibrate=True, clip=False)
        # hists_cer_vs_sci_random_subtracted_calibrated_clipped = makeFERSCervsSciRandomHists(
        #    subtractPedestal=True, calibrate=True, clip=True)

        # hists_channel_comparison = makeFERSChannelComparisonHists(
        #    subtractPedestal=False, calibrate=False)
        # hists_channel_comparison_subtracted = makeFERSChannelComparisonHists(
        #    subtractPedestal=True, calibrate=False)
        # hists_channel_comparison_subtracted_calibrated = makeFERSChannelComparisonHists(
        #    subtractPedestal=True, calibrate=True)
        # hists_channel_comparison_subtracted_calibrated_clipped = makeFERSChannelComparisonHists(
        #    subtractPedestal=True, calibrate=True, clip=True)

        # save histograms to ROOT files
        outfile_name = f"{rootdir}/fers_energy_sum.root"
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

        # outfile_name = f"{rootdir}/fers_energy_sum_subtracted_calibrated.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_subtracted_calibrated:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(f"Saved subtracted and calibrated histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_sum_subtracted_calibrated_clipped.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_subtracted_calibrated_clipped:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted, calibrated, and clipped histograms to {outfile_name}")

        # save cer vs sci histograms

        outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci.root"
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

        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_subtracted_calibrated.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_subtracted_calibrated:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted and calibrated CER vs SCI histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_subtracted_calibrated_clipped.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_subtracted_calibrated_clipped:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted, calibrated, and clipped CER vs SCI histograms to {outfile_name}")

        # save random CER vs SCI histograms
        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_random.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_random_raw:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(f"Saved random CER vs SCI histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_random_subtracted.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_random_subtracted:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(f"Saved random CER vs SCI histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_random_subtracted_calibrated.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_random_subtracted_calibrated:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved random CER vs SCI subtracted and calibrated histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_random_subtracted_calibrated_clipped.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_cer_vs_sci_random_subtracted_calibrated_clipped:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved random CER vs SCI subtracted, calibrated, and clipped histograms to {outfile_name}")

        # save channel comparison histograms
        # outfile_name = f"{rootdir}/fers_energy_channel_comparison.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_channel_comparison:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(f"Saved channel comparison histograms to {outfile_name}")
        # outfile_name = f"{rootdir}/fers_energy_channel_comparison_subtracted.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_channel_comparison_subtracted:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted channel comparison histograms to {outfile_name}")
        # outfile_name = f"{rootdir}/fers_energy_channel_comparison_subtracted_calibrated.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_channel_comparison_subtracted_calibrated:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted and calibrated channel comparison histograms to {outfile_name}")

        # outfile_name = f"{rootdir}/fers_energy_channel_comparison_subtracted_calibrated_clipped.root"
        # with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        #    for hist in hists_channel_comparison_subtracted_calibrated_clipped:
        #        hist.SetDirectory(outfile)
        #        hist.Write()
        # print(
        #    f"Saved subtracted, calibrated, and clipped channel comparison histograms to {outfile_name}")

    # make plots
    if makePlots:
        outputs_html["raw"] = makeFERSEnergySumPlots(
            subtractPedestal=False, calibrate=False)
        # outputs_html["subtracted"] = makeFERSEnergySumPlots(
        #    subtractPedestal=True, calibrate=False)
        # outputs_html["subtracted_calibrated"] = makeFERSEnergySumPlots(
        #     subtractPedestal=True, calibrate=True)
        # outputs_html["subtracted_calibrated_clipped"] = makeFERSEnergySumPlots(
        #    subtractPedestal=True, calibrate=True, clip=True)

        outputs_html["cer_vs_sci_raw"] = makeFERSCerVsSciPlots(
            subtractPedestal=False, calibrate=False)
        # outputs_html["cer_vs_sci_subtracted"] = makeFERSCerVsSciPlots(
        #    subtractPedestal=True, calibrate=False)
        # outputs_html["cer_vs_sci_subtracted_calibrated"] = makeFERSCerVsSciPlots(
        #     subtractPedestal=True, calibrate=True)
        # outputs_html["cer_vs_sci_subtracted_calibrated_clipped"] = makeFERSCerVsSciPlots(
        #    subtractPedestal=True, calibrate=True, clip=True)

        # outputs_html["cer_vs_sci_random_random"] = makeFERSCerVsSciRandomPlots(
        #    subtractPedestal=False, calibrate=False)
        # outputs_html["cer_vs_sci_random_subtracted"] = makeFERSCerVsSciRandomPlots(
        # subtractPedestal=True, calibrate=False)
        # outputs_html["cer_vs_sci_random_subtracted_calibrated"] = makeFERSCerVsSciRandomPlots(
        #    subtractPedestal=True, calibrate=True)
        # outputs_html["cer_vs_sci_random_subtracted_calibrated_clipped"] = makeFERSCerVsSciRandomPlots(
        #    subtractPedestal=True, calibrate=True, clip=True)

        # outputs_html["channel_comparison"] = makeFERSChannelComparisonPlots(
        #    subtractPedestal=False, calibrate=False)
        # outputs_html["channel_comparison_subtracted"] = makeFERSChannelComparisonPlots(
        # subtractPedestal=True, calibrate=False)
        # outputs_html["channel_comparison_subtracted_calibrated"] = makeFERSChannelComparisonPlots(
        #    subtractPedestal=True, calibrate=True)
        # outputs_html["channel_comparison_subtracted_calibrated_clipped"] = makeFERSChannelComparisonPlots(
        #    subtractPedestal=True, calibrate=True, clip=True)

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")
