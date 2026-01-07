import os
from collections import OrderedDict
import ROOT
from plotting.my_function import DrawHistos, LHistos2Hist
from configs.plot_ranges import getRangesForFERSEnergySums, getBoardEnergyFitParameters
from core.analysis_manager import CaloXAnalysisManager
from utils.colors import colors
from utils.fitter import eventFit
from utils.html_generator import generate_html
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.visualization import visualizeFERSBoards

auto_timer("Total Execution Time")

args = get_args()
setup_root(n_threads=10, batch_mode=True, load_functions=True)

analysis = (CaloXAnalysisManager(args)
            .prepare()                   # Baseline and vectorization
            .calibrate_fers()            # Pedestals, mixing, and response
            .apply_selections())         # Muon and hole vetoes

GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]

# calculate energy sums
for gain, calib in GainCalibs:
    analysis = analysis.define_physics_variables(
        gain=gain, calib=calib, pdsub=True)

rdf = analysis.get_rdf()  # Get the final RDF after all transformations
fersboards = analysis.fersboards

benergy = analysis.beam_energy
run_number = analysis.run_number
paths = analysis.paths
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]

doPerBoardPlots = False
HE = (benergy >= 50)  # GeV


def makeFERSEnergySumHists(rdf=rdf, suffix=""):
    hists_FERS_EnergySum = []
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
        for cat in ["cer", "sci"]:
            # per-board sum
            for fersboard in fersboards.values():
                varname = fersboard.get_energy_sum_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                hist = rdf.Histo1D((
                    f"hist_{varname}_{suffix}",
                    f"hist_{varname}_{suffix}",
                    500, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"]),
                    varname
                )
                hists_FERS_EnergySum.append(hist)
            # per-event sum
            varname = fersboards.get_energy_sum_name(
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
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
        # per-board Cer vs Sci
        for fersboard in fersboards.values():
            var_cer = fersboard.get_energy_sum_name(
                gain=gain, isCer=True, pdsub=True, calib=calib)
            var_sci = fersboard.get_energy_sum_name(
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
        var_cer = fersboards.get_energy_sum_name(
            gain=gain, isCer=True, pdsub=True, calib=calib)
        var_sci = fersboards.get_energy_sum_name(
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


def makeFERSDRHists(rdf=rdf, suffix=""):
    calib = True
    gain = "Mix"
    config = getRangesForFERSEnergySums(
        pdsub=True, calib=True, clip=False, HE=HE, run_number=run_number)
    varname_Cer = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    varname_Sci = fersboards.get_energy_sum_name(
        gain=gain, isCer=False, pdsub=True, calib=calib)
    hists_DR = []
    # leakage
    varname = varname_Cer + "_OuterRing"
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 500, 0, 20.0), varname)
    hists_DR.append(hist)
    varname = varname_Sci + "_OuterRing"
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 500, 0, 20.0), varname)
    hists_DR.append(hist)
    # leakage + central
    varname = varname_Cer + "_LeakCorr"
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]), varname)
    hists_DR.append(hist)
    varname = varname_Sci + "_LeakCorr"
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]), varname)
    hists_DR.append(hist)
    # per-event C/S ratio
    varname = "COverS"
    hist = rdf.Histo1D(
        (f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 100, 0, 2.0), varname)
    hists_DR.append(hist)
    # per-event fEM
    varname = "fEM"
    hist = rdf.Histo1D(
        (f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 100, 0, 2.0), varname)
    hists_DR.append(hist)
    # sum Cer + Sci
    varname = varname_Cer.replace("Cer", "CerSci")
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}",
        f"hist_{varname}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]*2.5),
        varname
    )
    hists_DR.append(hist)
    # dual readout
    for method in ["", "_method2", "_method3"]:
        varname = varname_Cer.replace("Cer", "DR") + method
        hist = rdf.Histo1D((
            f"hist_{varname}_{suffix}",
            f"hist_{varname}_{suffix}",
            500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]),
            varname
        )
        hists_DR.append(hist)

    # S vs Sum Cer + Sci
    var_sum = varname_Cer.replace("Cer", "CerSci")
    hist_Sci_vs_CerSci = rdf.Histo2D((
        f"hist_{varname_Sci}_VS_{var_sum}_{suffix}",
        f"hist_{varname_Sci}_VS_{var_sum}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"] * 2.5,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_sum,
        varname_Sci
    )
    hists_DR.append(hist_Sci_vs_CerSci)
    hist_Cer_vs_CerSci = rdf.Histo2D((
        f"hist_{varname_Cer}_VS_{var_sum}_{suffix}",
        f"hist_{varname_Cer}_VS_{var_sum}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"] * 2.5,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_sum,
        varname_Cer
    )
    hists_DR.append(hist_Cer_vs_CerSci)
    # S vs fEM
    var_fem = "fEM"
    hist_Sci_vs_fEM = rdf.Histo2D((
        f"hist_{varname_Sci}_VS_{var_fem}_{suffix}",
        f"hist_{varname_Sci}_VS_{var_fem}_{suffix}",
        500, 0., 2.0,
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]),
        var_fem,
        varname_Sci
    )
    hists_DR.append(hist_Sci_vs_fEM)
    # C vs fEM
    hist_Cer_vs_fEM = rdf.Histo2D((
        f"hist_{varname_Cer}_VS_{var_fem}_{suffix}",
        f"hist_{varname_Cer}_VS_{var_fem}_{suffix}",
        500, 0., 2.0,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_fem,
        varname_Cer
    )
    hists_DR.append(hist_Cer_vs_fEM)
    # C/S vs DR/DR_method2
    var_dr_base = varname_Cer.replace("Cer", "DR")
    for dr_method in ["", "_method2", "_method3"]:
        var_dr = var_dr_base + dr_method
        for cat, var_y, y_range_key in [("Cer", varname_Cer, "cer"), ("Sci", varname_Sci, "sci")]:
            hist = rdf.Histo2D((
                f"hist_{var_y}_VS_{var_dr}_{suffix}",
                f"hist_{var_y}_VS_{var_dr}_{suffix}",
                500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                500, config["xmin_total"][f"{gain}_{y_range_key}"], config["xmax_total"][f"{gain}_{y_range_key}"]),
                var_dr,
                var_y
            )
            hists_DR.append(hist)

    return hists_DR


def makeFERSEnergyWeightedCenterHists(rdf=rdf, suffix=""):
    hists_FERS_EnergyWeightedCenter = []
    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            # per-event
            varname_X = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.get_energy_weighted_center_name(
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
                for channel in fersboard.get_list_of_channels(isCer=(cat == "cer")):
                    channelName = channel.get_channel_name(
                        gain=gain, pdsub=True, calib=calib)
                    hist = rdf.Histo1D((
                        f"hist_RealX_{channelName}_{suffix}",
                        f"hist_RealX_{channelName}_{suffix}",
                        100, -20, 20),
                        channel.get_real_pos_name(isX=True), channelName
                    )
                    hists_tmp_X.append(hist)
                    hist = rdf.Histo1D((
                        f"hist_RealY_{channelName}_{suffix}",
                        f"hist_RealY_{channelName}_{suffix}",
                        100, -20, 20),
                        channel.get_real_pos_name(isX=False), channelName
                    )
                    hists_tmp_Y.append(hist)
                    # calculate radius
                    rdf = rdf.Define(
                        "RealR_" + channelName, f"std::sqrt(std::pow({channel.get_real_pos_name(isX=True)}, 2) + std::pow({channel.get_real_pos_name(isX=False)}, 2))")
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
                        channel.get_real_pos_name(isX=True), channel.get_real_pos_name(
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

    n_events = rdf.Count()

    return hists_X, hists_Y,  hists_R, hists_Y_VS_X, n_events


def collectFERSStats(rdf):
    stats = {}
    # mean
    # and how frequent the saturation value is reached
    for fersboard in fersboards.values():
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            channel_cer = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            channel_sci = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=False)

            for gain, calib in GainCalibs:
                varname_cer = channel_cer.get_channel_name(
                    gain=gain, pdsub=True, calib=calib)
                varname_sci = channel_sci.get_channel_name(
                    gain=gain, pdsub=True, calib=calib)
                # rdf = rdf.Define(
                #    f"{varname_cer}_over_{varname_sci}", f"{varname_cer} / ({varname_sci} + 1e-6)")

                stats[channel_cer.get_channel_name(gain=gain)] = rdf.Mean(
                    f"{varname_cer}")
                stats[channel_sci.get_channel_name(gain=gain)] = rdf.Mean(
                    f"{varname_sci}")

    return stats


def makeFERSEnergySumPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergySum_{suffix}"

    # per-board energy sum plots
    if doPerBoardPlots:
        for gain, calib in GainCalibs:
            config = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
            for cat in ["cer", "sci"]:
                hists = []
                legends = []
                for fersboard in fersboards.values():
                    board_no = fersboard.board_no
                    varname = fersboard.get_energy_sum_name(
                        gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                    hist_name = f"hist_{varname}_{suffix}"
                    hist = infile.Get(hist_name)
                    if not hist:
                        print(
                            f"Warning: Histogram {hist_name} for board {board_no} not found in {infile_name}")
                        continue
                    legends.append(str(board_no))
                    hists.append(hist)

                output_name = f"FERS_Boards_{gain}_{cat}{suffix}"
                DrawHistos(hists, legends, config["xmin_board"][f"{gain}_{cat}"], config["xmax_board"][f"{gain}_{cat}"], f"{cat.capitalize()} {gain} {config[f'title_{gain}']}", 1, None, "Events",
                           output_name,
                           dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
                           outdir=outdir_plots, run_number=run_number, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
                plots.append(output_name + ".png")

    # per-event energy sum plot ranges
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
        for cat in ["cer", "sci"]:
            varname = fersboards.get_energy_sum_name(
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
                       dology=False, drawoptions="HIST", mycolors=[2] if cat == "cer" else [4], addOverflow=True, addUnderflow=True,
                       outdir=outdir_plots, run_number=run_number)
            plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/FERS/EnergySum_{suffix}.html"
    generate_html(plots, outdir_plots, plots_per_row=6,
                  output_html=output_html)
    return output_html


def makeFERSCerVsSciPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_sum_cer_vs_sci_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_Cer_vs_Sci_{suffix}"

    if doPerBoardPlots:
        for gain, calib in GainCalibs:
            config = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
            for fersboard in fersboards.values():
                board_no = fersboard.board_no
                var_cer = fersboard.get_energy_sum_name(
                    gain=gain, isCer=True, pdsub=True, calib=calib)
                var_sci = fersboard.get_energy_sum_name(
                    gain=gain, isCer=False, pdsub=True, calib=calib)
                hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
                hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
                if not hist_Cer_vs_Sci:
                    print(
                        f"Warning: Histogram {hist_Cer_vs_Sci_name} for board {board_no} not found in {infile_name}")
                    continue

                output_name = f"FERS_Board{board_no}_Cer_VS_Sci_{gain}{suffix}"
                DrawHistos([hist_Cer_vs_Sci], "", config["xmin_board"][f"{gain}_sci"], config["xmax_board"][f"{gain}_sci"], f"Sci {gain} {config[f'title_{gain}']}", config["xmin_board"][f"{gain}_cer"], config["xmax_board"][f"{gain}_cer"], f"Cer {gain} {config[f'title_{gain}']}",
                           output_name,
                           dology=False, drawoptions=["colz"],
                           outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
                plots.append(output_name + ".png")

    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
        var_cer = fersboards.get_energy_sum_name(
            gain=gain, isCer=True, pdsub=True, calib=calib)
        var_sci = fersboards.get_energy_sum_name(
            gain=gain, isCer=False, pdsub=True, calib=calib)
        hist_Cer_vs_Sci_name = f"hist_{var_cer}_VS_{var_sci}_{suffix}"
        hist_Cer_vs_Sci = infile.Get(hist_Cer_vs_Sci_name)
        if not hist_Cer_vs_Sci:
            print(
                f"Warning: Histogram {hist_Cer_vs_Sci_name} not found in {infile_name}")
            continue
        extraObjs = []
        if gain == "Mix":
            # make a TF1 for Cer = Sci
            f11 = ROOT.TF1("f11", "x", 0, 120)
            f12 = ROOT.TF1("f12", "0.5 * x", 0, 120)
            line_x = ROOT.TLine(80, 0, 80, 120)
            line_y = ROOT.TLine(0, 80, 120, 80)

            extraObjs = [f11, f12, line_x, line_y]
            for line in extraObjs:
                line.SetLineColor(ROOT.kRed)
                line.SetLineStyle(2)
                line.SetLineWidth(2)
        output_name = f"FERS_Total_Cer_VS_Sci_{gain}{suffix}"
        DrawHistos([hist_Cer_vs_Sci], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], f"Sci {gain} {config[f'title_{gain}']}", config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], f"Cer {gain} {config[f'title_{gain}']}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   extraToDraw=extraObjs,
                   outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
        plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/FERS/EnergySum_Cer_VS_Sci_{suffix}.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeFERSDRPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_DR_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_DR_{suffix}"

    gain = "Mix"
    calib = True
    config = getRangesForFERSEnergySums(
        pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number)
    varname_Cer = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    varname_Sci = fersboards.get_energy_sum_name(
        gain=gain, isCer=False, pdsub=True, calib=calib)

    # leakage
    hists = []
    for cat, var in [("Cer", varname_Cer), ("Sci", varname_Sci)]:
        varname_leak = var + "_OuterRing"
        hist_name = f"hist_{varname_leak}_{suffix}"
        hist = infile.Get(hist_name)
        hists.append(hist)
    output_name = f"FERS_Total_{gain}_OuterRing_{suffix}"
    DrawHistos(hists, ["Cer", "Sci"], 0, 20.0, f"Leakage Energy", 0, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number)
    plots.append(output_name + ".png")

    # leakage + central
    hists = []
    for cat, var in [("Cer", varname_Cer), ("Sci", varname_Sci)]:
        varname_leak = var + "_LeakCorr"
        hist_name = f"hist_{varname_leak}_{suffix}"
        hist = infile.Get(hist_name)
        hists.append(hist)
    output_name = f"FERS_Total_{gain}_LeakCorr_{suffix}"
    DrawHistos(hists, ["Cer", "Sci"], config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], f"Leakage Corrected Energy {gain} {config[f'title_{gain}']}", 1, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number)
    plots.append(output_name + ".png")
    plots.append("NEWLINE")

    varname = "COverS"
    hist_name = f"hist_{varname}_{suffix}"
    hist = infile.Get(hist_name)
    output_name = f"FERS_Total_{gain}_CerOverSci_{suffix}"
    DrawHistos([hist], "", 0, 2, "C/S", 0, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number)
    plots.append(output_name + ".png")

    varname = "fEM"
    hist_name = f"hist_{varname}_{suffix}"
    hist = infile.Get(hist_name)
    output_name = f"FERS_Total_{gain}_fEM_{suffix}"
    DrawHistos([hist], "", 0, 1.5, "f_{EM}", 0, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number)
    plots.append(output_name + ".png")

    # sum Cer + Sci
    varname = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    varname = varname.replace("Cer", "CerSci")
    hist_name = f"hist_{varname}_{suffix}"
    hist = infile.Get(hist_name)
    output_name = f"FERS_Total_{gain}_CerSci_{suffix}"
    DrawHistos([hist], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"] * 2.5, f"Cer+Sci {gain} {config[f'title_{gain}']}", 1, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[6], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number)
    plots.append(output_name + ".png")

    # dual readout
    varname = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    for method in ["", "_method2", "_method3"]:
        varname_dr = varname.replace("Cer", "DR") + method
        hist_name = f"hist_{varname_dr}_{suffix}"
        hist = infile.Get(hist_name)
        output_name = f"FERS_Total_{gain}_DR{method}_{suffix}"
        DrawHistos([hist], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], f"DR{method} {gain} {config[f'title_{gain}']}", 1, None, "Events",
                   output_name,
                   dology=False, drawoptions="HIST", mycolors=[1 if method == "" else (7 if method == "_method2" else 8)], addOverflow=True, addUnderflow=True,
                   outdir=outdir_plots, run_number=run_number)
        plots.append(output_name + ".png")

    # plot sci, cer, and dr energy in the same plot, and fit the energy to get resolution
    varname_cer = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    varname_sci = fersboards.get_energy_sum_name(
        gain=gain, isCer=False, pdsub=True, calib=calib)
    varname_dr = varname_cer.replace("Cer", "DR")

    infile_ferssum = ROOT.TFile(
        f"{rootdir}/fers_energy_sum_{suffix}.root", "READ")
    hist_cer = infile_ferssum.Get(f"hist_{varname_cer}_{suffix}")
    hist_sci = infile_ferssum.Get(f"hist_{varname_sci}_{suffix}")
    hist_dr = infile.Get(f"hist_{varname_dr}_{suffix}")
    hist_dr_method2 = infile.Get(f"hist_{varname_dr}_method2_{suffix}")
    hist_dr_method3 = infile.Get(f"hist_{varname_dr}_method3_{suffix}")
    extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    # fit with gaus, plot the parameters
    tf1s = []
    for cat, hist in [("Cer", hist_cer), ("Sci", hist_sci), ("DR", hist_dr), ("DR method2", hist_dr_method2), ("DR method3", hist_dr_method3)]:
        if hist:
            fit_result = hist.Fit("gaus", "S")
            if fit_result.IsValid():
                fit_func = hist.GetFunction("gaus")
                mean = fit_func.GetParameter(1)
                sigma = fit_func.GetParameter(2)
                chi2 = fit_func.GetChisquare()
                ndf = fit_func.GetNDF()
                extraToDraw.AddText(
                    f"{cat}: #mu#pm#sigma = {mean:.1f}#pm{sigma:.1f}, #chi^{{2}}/ndf = {chi2:.1f}/{ndf}")

                fit_func.SetLineColor(
                    2 if cat == "Cer" else (4 if cat == "Sci" else 1))
                tf1s.append(fit_func)

    output_name = f"FERS_Total_{gain}_Energy_{suffix}"
    DrawHistos([hist_cer, hist_sci, hist_dr, hist_dr_method2, hist_dr_method3], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"], "Energy [GeV]", 0, None, "Events",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[2, 4, 1, 7, 8], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, extraToDraw=[extraToDraw] + tf1s)
    plots.append(output_name + ".png")

    plots.append("NEWLINE")

    # Draw Sci vs CerSci
    var_sum = varname_Cer.replace("Cer", "CerSci")
    hist_Sci_vs_CerSci_name = f"hist_{varname_Sci}_VS_{var_sum}_{suffix}"
    hist_Sci_vs_CerSci = infile.Get(hist_Sci_vs_CerSci_name)
    output_name = f"FERS_Total_Sci_VS_CerSci_{suffix}"
    DrawHistos([hist_Sci_vs_CerSci], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]*2.5, f"Sci+Cer {gain} {config[f'title_{gain}']}", config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], f"Sci {gain} {config[f'title_{gain}']}",
               output_name,
               dology=False, drawoptions=["colz"],
               outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
    plots.append(output_name + ".png")

    # Draw Cer vs CerSci
    hist_Cer_vs_CerSci_name = f"hist_{varname_Cer}_VS_{var_sum}_{suffix}"
    hist_Cer_vs_CerSci = infile.Get(hist_Cer_vs_CerSci_name)
    output_name = f"FERS_Total_Cer_VS_CerSci_{suffix}"
    DrawHistos([hist_Cer_vs_CerSci], "", config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]*2.5, f"Sci+Cer {gain} {config[f'title_{gain}']}", config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"], f"Cer {gain} {config[f'title_{gain}']}",
               output_name,
               dology=False, drawoptions=["colz"],
               outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
    plots.append(output_name + ".png")

    # Draw Sci vs fEM
    var_fem = "fEM"
    for cat in ["Sci", "Cer"]:
        varname = varname_Sci if cat == "Sci" else varname_Cer
        hist_varname_vs_fEM_name = f"hist_{varname}_VS_{var_fem}_{suffix}"
        hist_varname_vs_fEM = infile.Get(hist_varname_vs_fEM_name)
        hist_prof = hist_varname_vs_fEM.ProfileX()
        hist_prof.Fit("pol1")
        hist_prof.SetLineColor(ROOT.kRed)
        hist_prof.SetMarkerColor(ROOT.kRed)
        fit_result = hist_prof.Fit("pol1", "S", "", 0.4, 0.8)
        extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        if fit_result.IsValid():
            fit_func = hist_prof.GetFunction("pol1")
            p0 = fit_func.GetParameter(0)
            p1 = fit_func.GetParameter(1)
            fit_func.SetLineColor(ROOT.kRed)
            fit_func.SetLineWidth(2)
            fit_func.SetLineStyle(2)
            extraToDraw.AddText(f"{cat} = {p0:.1f} + {p1:.1f} * f_{{EM}}")

        output_name = f"FERS_Total_{cat}_VS_fEM_{suffix}"
        DrawHistos([hist_varname_vs_fEM], "", 0, 1.5, "f_{EM}", config["xmin_total"][f"{gain}_{cat.lower()}"], config["xmax_total"][f"{gain}_{cat.lower()}"], f"{cat} {gain} {config[f'title_{gain}']}",
                   output_name,
                   dology=False, drawoptions=["colz"],
                   outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True, extraToDraw=[extraToDraw, fit_func])
        plots.append(output_name + ".png")

    # Draw Cer/Sci vs DR/DR_method2
    var_dr_base = varname_Cer.replace("Cer", "DR")
    for dr_method, dr_label in [("", "DR"), ("_method2", "DR method2"), ("_method3", "DR method3")]:
        var_dr = var_dr_base + dr_method
        for cat, var_y, y_range_key in [("Cer", varname_Cer, "cer"), ("Sci", varname_Sci, "sci")]:
            hist_name = f"hist_{var_y}_VS_{var_dr}_{suffix}"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue

            output_name = f"FERS_Total_{cat}_VS_{var_dr.split('_')[-1]}_{suffix}"
            DrawHistos([hist], "",
                       config["xmin_total"][f"{gain}_sci"], config["xmax_total"][
                           f"{gain}_sci"], f"{dr_label} {gain} {config[f'title_{gain}']}",
                       config["xmin_total"][f"{gain}_{y_range_key}"], config["xmax_total"][
                           f"{gain}_{y_range_key}"], f"{cat} {gain} {config[f'title_{gain}']}",
                       output_name,
                       dology=False, drawoptions=["colz"],
                       outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None, addOverflow=True, addUnderflow=True)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS/EnergySum_DR_{suffix}.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSEnergyWeightedCenterPlots(suffix=""):
    plots = []
    infile_name = f"{rootdir}/fers_energy_weighted_center_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_EnergyWeightedCenter_{suffix}"

    for gain, calib in GainCalibs:
        for cat in ["cer", "sci"]:
            varname_X = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.get_energy_weighted_center_name(
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
                           outdir=outdir_plots, run_number=run_number, extraToDraw=extraToDraw)
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
                           outdir=outdir_plots, run_number=run_number, extraToDraw=extraToDraw)
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
                           outdir=outdir_plots, run_number=run_number, doth2=True, zmin=1, zmax=None)
                plots.append(output_name + ".png")
            else:
                print(
                    f"Warning: Histogram {hist2D_name} not found in {infile_name}")

    output_html = f"{htmldir}/FERS/EnergyWeightedCenter_{suffix}.html"
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
                           outdir=outdir_plots, run_number=run_number, donormalize=True)
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
                           outdir=outdir_plots, run_number=run_number, donormalize=True)
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
                   outdir=outdir_plots, run_number=run_number, donormalize=True)
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
                       outdir=outdir_plots, run_number=run_number)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS/ShowerShape_{suffix}.html"
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
            fersboards, stats, suffix=f"Run{run_number}_{gain}_mean", gain=gain)

        h2_Cer_Over_Sci = h2_Cer.Clone(f"h2_Cer_over_Sci_{gain}")
        h2_Cer_Over_Sci.Divide(h2_Sci)
        h2_Cer_3mm_Over_Sci = h2_Cer_3mm.Clone(f"h2_Cer_3mm_over_Sci_{gain}")
        h2_Cer_3mm_Over_Sci.Divide(h2_Sci_3mm)

        output_name = f"FERS_Boards_Run{run_number}_Stats_Cer_{gain}"
        DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=None)
        plots.append(output_name + ".png")
        output_name = f"FERS_Boards_Run{run_number}_Stats_Sci_{gain}"
        DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=None)
        plots.append(output_name + ".png")
        output_name = f"FERS_Boards_Run{run_number}_Stats_Cer_over_Sci_{gain}"
        DrawHistos([h2_Cer_Over_Sci, h2_Cer_3mm_Over_Sci], "", xmin, xmax, "iX", ymin,
                   ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
                   outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer / Sci", run_number=run_number, zmin=0, zmax=1.5)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS/Stats_Cer_over_Sci.html"
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
        board_no = fersboard.board_no
        hCer = ifile.Get(f"hist_FERS_Board{board_no}_CerEnergyHG_{suffix}")
        hSci = ifile.Get(f"hist_FERS_Board{board_no}_SciEnergyHG_{suffix}")

        args_cer = getBoardEnergyFitParameters(
            run_number, is3mm=fersboard.is_3mm_size(), isCer=True)
        args_sci = getBoardEnergyFitParameters(
            run_number, is3mm=fersboard.is_3mm_size(), isCer=False)

        output_name = eventFit(hCer, f"Run{run_number}_Board{board_no}_CerHG",
                               outdir=outdir, xlabel="Cer # p.e.",
                               **args_cer)
        plots.append(output_name)
        output_name = eventFit(hSci, f"Run{run_number}_Board{board_no}_SciHG",
                               outdir=outdir, xlabel="Sci # p.e.",
                               **args_sci)
        plots.append(output_name)

    output_html = f"{htmldir}/boardfits/index.html"
    generate_html(plots, outdir, plots_per_row=2,
                  output_html=output_html)
    print(f"Generated HTML file: {output_html}")
    return output_html


rdfs = OrderedDict()
rdfs["passHaloVeto"] = rdf


def main():
    makeHists = True
    makePlots = True
    outputs_html = {}

    for cat, rdf in rdfs.items():
        if makeHists:
            hists_raw = makeFERSEnergySumHists(rdf=rdf, suffix=cat)

            hists_cer_vs_sci_raw = makeFERSCervsSciHists(rdf=rdf, suffix=cat)

            hists_dr = makeFERSDRHists(rdf=rdf, suffix=cat)

            hists_energy_weighted_center = makeFERSEnergyWeightedCenterHists(
                rdf=rdf, suffix=cat)

            hists_shower_shapes_X, hists_shower_shapes_Y, hists_shower_shapes_R, hists_shower_shapes_Y_VS_X, n_events = makeFERSShowerShapeHists(
                rdf=rdf, suffix=cat)

            # stats = collectFERSStats(rdf=rdf)

            # import json
            # stats_results = {}
            # for channelName, mean in stats.items():
            #     stats_results[channelName] = mean.GetValue()
            # with open(f"{rootdir}/fers_stats_cer_ovs_sci.json", "w") as f:
            #     json.dump(stats_results, f, indent=4)

            analysis.sel_mgr.print_cutflow()

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

            # save dr histograms
            outfile_name = f"{rootdir}/fers_DR_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_dr:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved DR histograms to {outfile_name}")

            # save energy weighted center histograms
            outfile_name = f"{rootdir}/fers_energy_weighted_center_{cat}.root"
            with ROOT.TFile(outfile_name, "RECREATE") as outfile:
                for hist in hists_energy_weighted_center:
                    hist.SetDirectory(outfile)
                    hist.Write()
            print(f"Saved energy weighted center histograms to {outfile_name}")

            nEvts = n_events.GetValue()
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

            outputs_html[f"dr_{cat}"] = makeFERSDRPlots(suffix=cat)

            outputs_html[f"energy_weighted_center_{cat}"] = makeFERSEnergyWeightedCenterPlots(
                suffix=cat)

            outputs_html[f"shower_shape_{cat}"] = makeFERSShowerShapePlots(
                suffix=cat)

            # outputs_html[f"stats_cer_ovs_sci"] = makeFERSStatsPlots()

    print("Generated HTML files:")
    for key, html in outputs_html.items():
        print(f"{key}: {html}")


if __name__ == "__main__":
    main()
