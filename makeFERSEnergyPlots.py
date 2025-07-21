import os
import sys
import ROOT
from utils.channel_map import buildFERSBoards
from utils.utils import number2string, filterPrefireEvents, loadRDF, calculateEnergySumFERS, vectorizeFERS, calibrateFERSChannels
from utils.html_generator import generate_html
from utils.fitter import eventFit
from utils.colors import colors
from runconfig import runNumber, firstEvent, lastEvent
import time
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos


print("Start running prepareDQMPlots.py")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

debugDRS = False


file_gains = f"results/root/Run{runNumber}/valuemaps_gain.json"
file_pedestals = f"results/root/Run{runNumber}/valuemaps_pedestal.json"

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf, rdf_prefilter = filterPrefireEvents(rdf, runNumber)

FERSBoards = buildFERSBoards(run=runNumber)

# Get total number of entries
n_entries = rdf.Count().GetValue()
nEvents = int(n_entries)
nbins_Event = min(max(int(nEvents / 100), 1), 500)
print(f"Total number of events to process: {nEvents} in run {runNumber}")

rdf = vectorizeFERS(rdf, FERSBoards)
rdf = calibrateFERSChannels(
    rdf, FERSBoards, file_gains=file_gains, file_pedestals=file_pedestals)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=False, calibrate=False)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=True, calibrate=False)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=True, calibrate=True)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"


def makeFERSEnergySumPlots(subtractPedestal=False, calibrate=False):
    suffix = ""
    if subtractPedestal:
        suffix = "_subtracted"
    if calibrate:
        suffix += "_calibrated"
    xmax_board = 15000
    xmax_total = 200000
    if subtractPedestal:
        suffix = "_subtracted"
        xmax_board = 10000
        xmax_total = 5e4
    if calibrate:
        suffix += "_calibrated"
        xmax_board = 200
        xmax_total = 1e3
    hists_FERS_EnergySum = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_CerEnergyHG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_CerEnergyHG{suffix}",
            f"FERS Board {boardNo} - CER Energy HG;CER Energy HG;Counts",
            500, 0, xmax_board),
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
            500, 0, xmax_board),
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
        500, 0, xmax_total),
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
        500, 0, xmax_total),
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


def plotFERSEnergySumPlots(subtractPedestal=False, calibrate=False):
    suffix = ""
    xmax_board = 15000
    xmax_total = 200000
    if subtractPedestal:
        suffix = "_subtracted"
        xmax_board = 10000
        xmax_total = 5e4
    if calibrate:
        suffix += "_calibrated"
        xmax_board = 200
        xmax_total = 1e3
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
    DrawHistos(hists_CerEnergyHG, legends, 0, xmax_board, "Cer Energy HG", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.7, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_SciEnergyHG" + suffix
    DrawHistos(hists_SciEnergyHG, legends, 0, xmax_board, "Sci Energy HG", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=[0.30, 0.7, 0.90, 0.9])
    plots.append(output_name + ".png")
    # output_name = "FERS_CerEnergyLG" + suffix
    # DrawHistos(hists_CerEnergyLG, legends, 0, xmax_total, "Cer Energy LG", 0, None, "Events",
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
    DrawHistos([hist_CerEnergyHG], "", 0, xmax_total, "Cer Energy HG", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(0, output_name + ".png")
    output_name = "FERS_Total_SciEnergyHG" + suffix
    DrawHistos([hist_SciEnergyHG], "", 0, xmax_total, "Sci Energy HG", 1, ymax, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(1, output_name + ".png")
    # output_name = "FERS_Total_CerEnergyLG" + suffix
    # DrawHistos([hist_CerEnergyLG], "", 0,  xmax_total, "Cer Energy LG", 0, ymax, "Events",
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


def makeEventFits():
    suffix = "subtracted_calibrated"
    filename = f"{rootdir}/fers_energy_sum_{suffix}.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run prepareDQMPlots.py first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    hCer = ifile.Get(f"hist_FERS_CerEnergyHG_{suffix}")
    hSci = ifile.Get(f"hist_FERS_SciEnergyHG_{suffix}")

    plots = []
    outdir = f"{plotdir}/energyfits"
    output_name = eventFit(hCer, f"Run{runNumber}_CerHG",
                           outdir=outdir, addMIP=True, addHE=False, xlabel="Cer # p.e.",
                           xmin=0, xmax=1000,
                           xfitmin=0, xfitmax=250,
                           xgausmean=105, xgausmin=50, xgausmax=120,
                           wgausmean=40, wgausmin=20, wgausmax=50,
                           xmipmean=170, xmipmin=120, xmipmax=200,
                           wmipmean=60, wmipmin=30, wmipmax=80,
                           runNumber=runNumber)
    plots.append(output_name)
    output_name = eventFit(hSci, f"Run{runNumber}_SciHG",
                           outdir=outdir, addMIP=True, addHE=False, xlabel="Sci # p.e.",
                           xmin=0, xmax=1000,
                           xfitmin=0, xfitmax=450,
                           xgausmean=110, xgausmin=50, xgausmax=150,
                           wgausmean=40, wgausmin=20, wgausmax=60,
                           xmipmean=300, xmipmin=200, xmipmax=350,
                           wmipmean=60, wmipmin=30, wmipmax=80,
                           runNumber=runNumber)
    plots.append(output_name)
    output_html = f"{htmldir}/energyfits/index.html"
    generate_html(plots, outdir, plots_per_row=2,
                  output_html=output_html)
    print(f"Generated HTML file: {output_html}")
    return output_html


if __name__ == "__main__":
    hists_raw = makeFERSEnergySumPlots(subtractPedestal=False, calibrate=False)
    hists_subtracted = makeFERSEnergySumPlots(
        subtractPedestal=True, calibrate=False)
    hists_subtractged_calibrated = makeFERSEnergySumPlots(
        subtractPedestal=True, calibrate=True)

    # save histograms to ROOT files
    outfile_name = f"{rootdir}/fers_energy_sum.root"
    with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        for hist in hists_raw:
            hist.SetDirectory(outfile)
            hist.Write()
    print(f"Saved raw histograms to {outfile_name}")

    outfile_name = f"{rootdir}/fers_energy_sum_subtracted.root"
    with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        for hist in hists_subtracted:
            hist.SetDirectory(outfile)
            hist.Write()
    print(f"Saved subtracted histograms to {outfile_name}")

    outfile_name = f"{rootdir}/fers_energy_sum_subtracted_calibrated.root"
    with ROOT.TFile(outfile_name, "RECREATE") as outfile:
        for hist in hists_subtractged_calibrated:
            hist.SetDirectory(outfile)
            hist.Write()
    print(f"Saved subtracted and calibrated histograms to {outfile_name}")

    # make plots
    outputs_html = {}
    outputs_html["raw"] = plotFERSEnergySumPlots(
        subtractPedestal=False, calibrate=False)
    outputs_html["subtracted"] = plotFERSEnergySumPlots(
        subtractPedestal=True, calibrate=False)
    outputs_html["subtracted_calibrated"] = plotFERSEnergySumPlots(
        subtractPedestal=True, calibrate=True)

    # run event fits
    outputs_html["event_fits"] = makeEventFits()

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")
