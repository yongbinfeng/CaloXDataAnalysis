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
#    rdf, fersboards, pdsub=True, calib=True, clip=False)

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
    config = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)

    hists_FERS_EnergySum = []
    for gain in ["HG", "LG"]:
        for cat in ["cer", "sci"]:
            # per-board sum
            for fersboard in fersboards.values():
                varname = fersboard.GetEnergySumName(
                    useHG=(gain == "HG"), isCer=(cat == "cer"), pdsub=pdsub, calib=calib)
                hist = rdf.Histo1D((
                    f"hist_{varname}_{suffix}",
                    f"hist_{varname}_{suffix}",
                    500, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"]),
                    varname
                )
                hists_FERS_EnergySum.append(hist)
            # per-event sum
            varname = fersboards.GetEnergySumName(
                useHG=(gain == "HG"), isCer=(cat == "cer"), pdsub=pdsub, calib=calib)
            hist = rdf.Histo1D((
                f"hist_{varname}_{suffix}",
                f"hist_{varname}_{suffix}",
                500, config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"]),
                varname
            )
            hists_FERS_EnergySum.append(hist)

    return hists_FERS_EnergySum


def makeFERSCervsSciHists(rdf=rdf, pdsub=False, calib=False, clip=False, suffix=""):
    config = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)

    hists_FERS_Cer_vs_Sci = []
    for gain in ["HG", "LG"]:
        # per-board Cer vs Sci
        for fersboard in fersboards.values():
            var_cer = fersboard.GetEnergySumName(
                useHG=(gain == "HG"), isCer=True, pdsub=pdsub, calib=calib)
            var_sci = fersboard.GetEnergySumName(
                useHG=(gain == "HG"), isCer=False, pdsub=pdsub, calib=calib)
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
            useHG=(gain == "HG"), isCer=True, pdsub=pdsub, calib=calib)
        var_sci = fersboards.GetEnergySumName(
            useHG=(gain == "HG"), isCer=False, pdsub=pdsub, calib=calib)
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


def makeFERSEnergySumPlots(pdsub=False, calib=False, clip=False, suffix=""):
    config = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergySum_{suffix}"

    # per-board energy sum plots
    for gain in ["HG", "LG"]:
        for cat in ["cer", "sci"]:
            hists = []
            legends = []
            for fersboard in fersboards.values():
                boardNo = fersboard.boardNo
                varname = fersboard.GetEnergySumName(
                    useHG=(gain == "HG"), isCer=(cat == "cer"), pdsub=pdsub, calib=calib)
                hist_name = f"hist_{varname}_{suffix}"
                hist = infile.Get(hist_name)
                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} for board {boardNo} not found in {infile_name}")
                    continue
                legends.append(str(boardNo))
                hists.append(hist)

            output_name = f"FERS_Boards_{gain}_{cat}{suffix}"
            DrawHistos(hists, legends, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"], f"{cat.capitalize()} {gain} {config['title']}", 1, None, "Events",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
            plots.append(output_name + ".png")

    # per-event energy sum plot ranges
    for gain in ["HG", "LG"]:
        for cat in ["cer", "sci"]:
            varname = fersboards.GetEnergySumName(
                useHG=(gain == "HG"), isCer=(cat == "cer"), pdsub=pdsub, calib=calib)
            hist_name = f"hist_{varname}_{suffix}"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue
            output_name = f"FERS_Total_{gain}_{cat}{suffix}"
            DrawHistos([hist], "", config["xmin_total"][f"{gain}_{cat}"], config["xmax_total"][f"{gain}_{cat}"], f"{cat.capitalize()} {gain} {config['title']}", 1, None, "Events",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/FERS_EnergySum{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSCerVsSciPlots(pdsub=False, calib=False, clip=False, suffix=""):
    config = getRangesForFERSEnergySums(
        pdsub=pdsub, calib=calib, clip=clip, HE=HE)

    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci_{suffix}"

    for gain in ["HG", "LG"]:
        for fersboard in fersboards.values():
            boardNo = fersboard.boardNo
            var_cer = fersboard.GetEnergySumName(
                useHG=(gain == "HG"), isCer=True, pdsub=pdsub, calib=calib)
            var_sci = fersboard.GetEnergySumName(
                useHG=(gain == "HG"), isCer=False, pdsub=pdsub, calib=calib)
            hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
            hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
            if not hist_Cer_vs_Sci:
                print(
                    f"Warning: Histogram {hist_Cer_vs_Sci_name} for board {boardNo} not found in {infile_name}")
                continue

            output_name = f"FERS_Board{boardNo}_Cer_VS_Sci_{gain}{suffix}"
            DrawHistos([hist_Cer_vs_Sci], "", config["xmin_board"][f"{gain}_sci"], config["xmax_board"][f"{gain}_sci"], f"Sci {gain} {config['title']}", config["xmin_board"][f"{gain}_cer"], config["xmax_board"][f"{gain}_cer"], f"Cer {gain} {config['title']}",
                       output_name,
                       dology=False, drawoptions=["colz"],
                       outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
            plots.append(output_name + ".png")

    for gain in ["HG", "LG"]:
        var_cer = fersboards.GetEnergySumName(
            useHG=(gain == "HG"), isCer=True, pdsub=pdsub, calib=calib)
        var_sci = fersboards.GetEnergySumName(
            useHG=(gain == "HG"), isCer=False, pdsub=pdsub, calib=calib)
        hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue
        output_name = f"FERS_Total_Cer_VS_Sci_{gain}{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], f"Sci {gain} {config['title']}", config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], f"Cer {gain} {config['title']}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, runNumber=runNumber, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_Cer_VS_Sci{suffix}/index.html"
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

            hists_cer_vs_sci_raw = makeFERSCervsSciHists(rdf=rdf,
                                                         pdsub=True, calib=False, clip=False, suffix=cat)
            # hists_cer_vs_sci_subtracted = makeFERSCervsSciHists(
            #    pdsub=True, calib=False, clip=False)
            # hists_cer_vs_sci_subtracted_calibd = makeFERSCervsSciHists(
            #    pdsub=True, calib=True, clip=False)

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

        # make plots
        if makePlots:
            outputs_html[f"raw_{cat}"] = makeFERSEnergySumPlots(
                pdsub=True, calib=False, suffix=cat)
            # outputs_html["subtracted"] = makeFERSEnergySumPlots(
            #    pdsub=True, calib=False)
            # outputs_html["subtracted_calibd"] = makeFERSEnergySumPlots(
            #     pdsub=True, calib=True)

            outputs_html[f"cer_vs_sci_raw_{cat}"] = makeFERSCerVsSciPlots(
                pdsub=True, calib=False, suffix=cat)
            # outputs_html["cer_vs_sci_subtracted"] = makeFERSCerVsSciPlots(
            #    pdsub=True, calib=False)
            # outputs_html["cer_vs_sci_subtracted_calibd"] = makeFERSCerVsSciPlots(
            #     pdsub=True, calib=True)

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")
