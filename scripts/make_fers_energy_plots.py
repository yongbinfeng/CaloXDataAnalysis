"""
FERS Energy Plot Generation Script.

This script generates FERS energy sum plots, Cer vs Sci correlations,
dual readout plots, energy weighted center plots, and shower shape plots.
"""

import os
from collections import OrderedDict
import ROOT
from plotting.my_function import DrawHistos, LHistos2Hist
from configs.plot_config import getRangesForFERSEnergySums, getBoardEnergyFitParameters, get_ttu_hodo_ranges
from core.analysis_manager import CaloXAnalysisManager
from utils.colors import colors
from utils.fitter import eventFit
from utils.html_generator import generate_html
from utils.parser import get_args
from configs.plot_style import PlotStyle, STYLE_CER, STYLE_SCI, STYLE_CER_SCI, STYLE_2D_COLZ
from plotting.calox_plot_helper import create_pave_text
from core.plot_manager import PlotManager
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.visualization import visualizeFERSBoards

auto_timer("Total Execution Time")

args = get_args()
setup_root(n_threads=10, batch_mode=True, load_functions=True)

analysis = (CaloXAnalysisManager(args)
            .prepare()
            .calibrate_fers()
            .apply_hole_veto(flag_only=True)
            )

GainCalibs = [("HG", False), ("LG", False), ("Mix", True)]

# Calculate energy sums
for gain, calib in GainCalibs:
    analysis = analysis.define_physics_variables(
        gain=gain, calib=calib, pdsub=True)

fersboards = analysis.fersboards

benergy = analysis.beam_energy
run_number = analysis.run_number
paths = analysis.paths
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]

doPerBoardPlots = False
HE = (benergy >= 50)  # GeV

# Common styles
STYLE_BOARD_MULTI = PlotStyle(
    dology=True,
    drawoptions="HIST",
    mycolors=colors,
    legendNCols=5,
    legendPos=[0.20, 0.75, 0.90, 0.9]
)

STYLE_2D_LOG = PlotStyle(
    dology=False,
    dologz=False,
    drawoptions=["colz"],
    addOverflow=True,
    addUnderflow=True,
    zmin=1,
    zmax=None
)


def makeFERSEnergySumHists(rdf, suffix=""):
    """Book FERS energy sum histograms."""
    hists_FERS_EnergySum = []
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        for cat in ["cer", "sci"]:
            # Per-board sum
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
            # Per-event sum
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


def makeFERSCervsSciHists(rdf, suffix=""):
    """Book FERS Cer vs Sci 2D histograms."""
    hists_FERS_Cer_vs_Sci = []
    for gain, calib in GainCalibs:
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
        # Per-board Cer vs Sci
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

        # Per-event Cer vs Sci
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


def makeFERSDRHists(rdf, suffix=""):
    """Book dual readout histograms."""
    calib = True
    gain = "Mix"
    config = getRangesForFERSEnergySums(
        pdsub=True, calib=True, clip=False, HE=HE, run_number=run_number, beam_energy=benergy)
    varname_Cer = fersboards.get_energy_sum_name(
        gain=gain, isCer=True, pdsub=True, calib=calib)
    varname_Sci = fersboards.get_energy_sum_name(
        gain=gain, isCer=False, pdsub=True, calib=calib)
    hists_DR = []

    # Leakage
    for var in [varname_Cer, varname_Sci]:
        varname = var + "_OuterRing"
        hist = rdf.Histo1D((
            f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}", 500, 0, 20.0), varname)
        hists_DR.append(hist)

    # Leakage corrected
    for var in [varname_Cer, varname_Sci]:
        varname = var + "_LeakCorr"
        hist = rdf.Histo1D((
            f"hist_{varname}_{suffix}", f"hist_{varname}_{suffix}",
            500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]), varname)
        hists_DR.append(hist)

    # C/S ratio
    hist = rdf.Histo1D(
        (f"hist_COverS_{suffix}", f"hist_COverS_{suffix}", 100, 0, 2.0), "COverS")
    hists_DR.append(hist)

    # fEM
    hist = rdf.Histo1D(
        (f"hist_fEM_{suffix}", f"hist_fEM_{suffix}", 100, 0, 2.0), "fEM")
    hists_DR.append(hist)

    # Sum Cer + Sci
    varname = varname_Cer.replace("Cer", "CerSci")
    hist = rdf.Histo1D((
        f"hist_{varname}_{suffix}",
        f"hist_{varname}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]*2.5),
        varname
    )
    hists_DR.append(hist)

    # Dual readout methods
    for method in ["", "_method2", "_method3"]:
        varname = varname_Cer.replace("Cer", "DR") + method
        hist = rdf.Histo1D((
            f"hist_{varname}_{suffix}",
            f"hist_{varname}_{suffix}",
            500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]),
            varname
        )
        hists_DR.append(hist)

    # 2D correlations
    var_sum = varname_Cer.replace("Cer", "CerSci")

    # S vs CerSci
    hist_Sci_vs_CerSci = rdf.Histo2D((
        f"hist_{varname_Sci}_VS_{var_sum}_{suffix}",
        f"hist_{varname_Sci}_VS_{var_sum}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"] * 2.5,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_sum, varname_Sci
    )
    hists_DR.append(hist_Sci_vs_CerSci)

    # C vs CerSci
    hist_Cer_vs_CerSci = rdf.Histo2D((
        f"hist_{varname_Cer}_VS_{var_sum}_{suffix}",
        f"hist_{varname_Cer}_VS_{var_sum}_{suffix}",
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"] * 2.5,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_sum, varname_Cer
    )
    hists_DR.append(hist_Cer_vs_CerSci)

    # S vs fEM
    var_fem = "fEM"
    hist_Sci_vs_fEM = rdf.Histo2D((
        f"hist_{varname_Sci}_VS_{var_fem}_{suffix}",
        f"hist_{varname_Sci}_VS_{var_fem}_{suffix}",
        500, 0., 2.0,
        500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"]),
        var_fem, varname_Sci
    )
    hists_DR.append(hist_Sci_vs_fEM)

    # C vs fEM
    hist_Cer_vs_fEM = rdf.Histo2D((
        f"hist_{varname_Cer}_VS_{var_fem}_{suffix}",
        f"hist_{varname_Cer}_VS_{var_fem}_{suffix}",
        500, 0., 2.0,
        500, config["xmin_total"][f"{gain}_cer"], config["xmax_total"][f"{gain}_cer"]),
        var_fem, varname_Cer
    )
    hists_DR.append(hist_Cer_vs_fEM)

    # C/S vs DR methods
    var_dr_base = varname_Cer.replace("Cer", "DR")
    for dr_method in ["", "_method2", "_method3"]:
        var_dr = var_dr_base + dr_method
        for cat, var_y, y_range_key in [("Cer", varname_Cer, "cer"), ("Sci", varname_Sci, "sci")]:
            hist = rdf.Histo2D((
                f"hist_{var_y}_VS_{var_dr}_{suffix}",
                f"hist_{var_y}_VS_{var_dr}_{suffix}",
                500, config["xmin_total"][f"{gain}_sci"], config["xmax_total"][f"{gain}_sci"],
                500, config["xmin_total"][f"{gain}_{y_range_key}"], config["xmax_total"][f"{gain}_{y_range_key}"]),
                var_dr, var_y
            )
            hists_DR.append(hist)

    return hists_DR


def makeFERSEnergyWeightedCenterHists(rdf, suffix=""):
    """Book energy weighted center histograms."""
    hists_FERS_EnergyWeightedCenter = []
    for gain, calib in GainCalibs:
        if gain != "Mix":
            continue
        for cat in ["cer", "sci"]:
            varname_X = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)

            histX = rdf.Histo1D((
                f"hist_{varname_X}_{suffix}", f"hist_{varname_X}_{suffix}",
                300, -15, 15), varname_X)
            hists_FERS_EnergyWeightedCenter.append(histX)

            histY = rdf.Histo1D((
                f"hist_{varname_Y}_{suffix}", f"hist_{varname_Y}_{suffix}",
                300, -15, 15), varname_Y)
            hists_FERS_EnergyWeightedCenter.append(histY)

            hist2D = rdf.Histo2D((
                f"hist_{varname_Y}_VS_{varname_X}_{suffix}",
                f"hist_{varname_Y}_VS_{varname_X}_{suffix}",
                300, -15, 15, 300, -15, 15),
                varname_X, varname_Y)
            hists_FERS_EnergyWeightedCenter.append(hist2D)

            # average energy with respect to EWC
            hprof2D_energy = rdf.Profile2D((
                f"hprof_{varname_Y}_VS_{varname_X}_WithEnergy_{suffix}",
                f"hprof_{varname_Y}_VS_{varname_X}_WithEnergy_{suffix}",
                300, -15, 15, 300, -15, 15),
                varname_X, varname_Y,
                fersboards.get_energy_sum_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
            )
            hists_FERS_EnergyWeightedCenter.append(hprof2D_energy)

    return hists_FERS_EnergyWeightedCenter


def makeFERSEWCvsHodoHists(rdf, suffix=""):
    """Book EWC vs Hodoscope histograms."""
    hists_EWC_vs_Hodo = []
    hodo_min, hodo_max, hodo_nbins = get_ttu_hodo_ranges()

    for gain, calib in GainCalibs:
        if gain != "Mix":
            continue
        for cat in ["cer", "sci"]:
            varname_X = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
            varname_Y = fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)
            energy_var = fersboards.get_energy_sum_name(
                gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)

            # EWC X vs Hodo X
            hists_EWC_vs_Hodo.append(rdf.Histo2D((
                f"hist_{varname_X}_VS_HodoX_{suffix}",
                f"hist_{varname_X}_VS_HodoX_{suffix}",
                hodo_nbins, hodo_min, hodo_max, 300, -15, 15),
                "TTU_Hodo_X", varname_X))

            # EWC Y vs Hodo Y
            hists_EWC_vs_Hodo.append(rdf.Histo2D((
                f"hist_{varname_Y}_VS_HodoY_{suffix}",
                f"hist_{varname_Y}_VS_HodoY_{suffix}",
                hodo_nbins, hodo_min, hodo_max, 300, -15, 15),
                "TTU_Hodo_Y", varname_Y))

            # average energy vs EWC X and Hodo X
            hists_EWC_vs_Hodo.append(rdf.Profile2D((
                f"hprof_{varname_X}_VS_HodoX_WithEnergy_{suffix}",
                f"hprof_{varname_X}_VS_HodoX_WithEnergy_{suffix}",
                hodo_nbins, hodo_min, hodo_max, 300, -15, 15),
                "TTU_Hodo_X", varname_X, energy_var))

            # average energy vs EWC Y and Hodo Y
            hists_EWC_vs_Hodo.append(rdf.Profile2D((
                f"hprof_{varname_Y}_VS_HodoY_WithEnergy_{suffix}",
                f"hprof_{varname_Y}_VS_HodoY_WithEnergy_{suffix}",
                hodo_nbins, hodo_min, hodo_max, 300, -15, 15),
                "TTU_Hodo_Y", varname_Y, energy_var))

            # average energy vs Hodo Y and Hodo X
            hists_EWC_vs_Hodo.append(rdf.Profile2D((
                f"hprof_HodoY_VS_HodoX_WithEnergy_{suffix}",
                f"hprof_HodoY_VS_HodoX_WithEnergy_{suffix}",
                hodo_nbins, hodo_min, hodo_max, hodo_nbins, hodo_min, hodo_max),
                "TTU_Hodo_X", "TTU_Hodo_Y", energy_var))

    return hists_EWC_vs_Hodo


def makeFERSShowerShapeHists(rdf, suffix=""):
    """Book shower shape histograms."""
    hists_X = []
    hists_Y = []
    hists_R = []
    hists_Y_VS_X = []

    for gain, calib in GainCalibs:
        if gain != "Mix":
            continue
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
                        channel.get_real_pos_name(isX=True), channelName)
                    hists_tmp_X.append(hist)

                    hist = rdf.Histo1D((
                        f"hist_RealY_{channelName}_{suffix}",
                        f"hist_RealY_{channelName}_{suffix}",
                        100, -20, 20),
                        channel.get_real_pos_name(isX=False), channelName)
                    hists_tmp_Y.append(hist)

                    # Calculate radius
                    rdf = rdf.Define(
                        "RealR_" + channelName,
                        f"std::sqrt(std::pow({channel.get_real_pos_name(isX=True)}, 2) + "
                        f"std::pow({channel.get_real_pos_name(isX=False)}, 2))")
                    hist = rdf.Histo1D((
                        f"hist_RealR_{channelName}_{suffix}",
                        f"hist_RealR_{channelName}_{suffix}",
                        25, 0, 25),
                        "RealR_" + channelName, channelName)
                    hists_tmp_R.append(hist)

                    hist = rdf.Histo2D((
                        f"hist_RealY_VS_RealX_{channelName}_{suffix}",
                        f"hist_RealY_VS_RealX_{channelName}_{suffix}",
                        100, -20, 20, 100, -20, 20),
                        channel.get_real_pos_name(isX=True),
                        channel.get_real_pos_name(isX=False), channelName)
                    hists_tmp_Y_VS_X.append(hist)

            hists_X.append((hists_tmp_X, f"hist_RealX_{gain}_{cat}_{suffix}"))
            hists_Y.append((hists_tmp_Y, f"hist_RealY_{gain}_{cat}_{suffix}"))
            hists_R.append((hists_tmp_R, f"hist_RealR_{gain}_{cat}_{suffix}"))
            hists_Y_VS_X.append(
                (hists_tmp_Y_VS_X, f"hist_RealY_VS_RealX_{gain}_{cat}_{suffix}"))

    n_events = rdf.Count()
    return hists_X, hists_Y, hists_R, hists_Y_VS_X, n_events


# ============================================================================
# Plotting Functions
# ============================================================================

def makeFERSEnergySumPlots(suffix=""):
    """Plot FERS energy sum distributions."""
    filename = f"fers_energy_sum_{suffix}.root"

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_EnergySum_{suffix}")

        # Per-board plots
        if doPerBoardPlots:
            for gain, calib in GainCalibs:
                config = getRangesForFERSEnergySums(
                    pdsub=True, calib=calib, clip=False, HE=HE,
                    run_number=run_number, beam_energy=benergy)

                for cat in ["cer", "sci"]:
                    hists, legends = [], []
                    for fersboard in fersboards.values():
                        varname = fersboard.get_energy_sum_name(
                            gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                        hist = pm.get_histogram_by_pattern(
                            filename, varname, suffix, required=False)
                        if hist:
                            hists.append(hist)
                            legends.append(str(fersboard.board_no))

                    if hists:
                        pm.plot_1d(
                            hists,
                            f"FERS_Boards_{gain}_{cat}{suffix}",
                            f"{cat.capitalize()} {gain} {config[f'title_{gain}']}",
                            (config["xmin_board"][f"{gain}_{cat}"],
                             config["xmax_board"][f"{gain}_{cat}"]),
                            yrange=(1, None),
                            legends=legends,
                            style=STYLE_BOARD_MULTI
                        )

        # Per-event plots
        for gain, calib in GainCalibs:
            config = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE,
                run_number=run_number, beam_energy=benergy)

            for cat in ["cer", "sci"]:
                varname = fersboards.get_energy_sum_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib)
                hist = pm.get_histogram_by_pattern(
                    filename, varname, suffix, required=False)

                if hist:
                    style = STYLE_CER if cat == "cer" else STYLE_SCI
                    pm.plot_1d(
                        hist,
                        f"FERS_Total_{gain}_{cat}{suffix}",
                        f"{cat.capitalize()} {gain} {config[f'title_{gain}']}",
                        (config["xmin_total"][f"{gain}_{cat}"],
                         config["xmax_total"][f"{gain}_{cat}"]),
                        style=style,
                        prepend=True
                    )

        return pm.generate_html(f"FERS/{suffix}/EnergySum.html", plots_per_row=6)


def makeFERSCerVsSciPlots(suffix=""):
    """Plot FERS Cer vs Sci 2D correlations."""
    filename = f"fers_energy_sum_cer_vs_sci_{suffix}.root"

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_Cer_vs_Sci_{suffix}")

        # Per-board plots
        if doPerBoardPlots:
            for gain, calib in GainCalibs:
                config = getRangesForFERSEnergySums(
                    pdsub=True, calib=calib, clip=False, HE=HE,
                    run_number=run_number, beam_energy=benergy)

                for fersboard in fersboards.values():
                    board_no = fersboard.board_no
                    var_cer = fersboard.get_energy_sum_name(
                        gain=gain, isCer=True, pdsub=True, calib=calib)
                    var_sci = fersboard.get_energy_sum_name(
                        gain=gain, isCer=False, pdsub=True, calib=calib)

                    hist = pm.get_histogram(
                        filename, f"hist_{var_cer}_VS_{var_sci}_{suffix}", required=False)
                    if hist:
                        pm.plot_2d(
                            hist,
                            f"FERS_Board{board_no}_Cer_VS_Sci_{gain}{suffix}",
                            f"Sci {gain} {config[f'title_{gain}']}",
                            (config["xmin_board"][f"{gain}_sci"],
                             config["xmax_board"][f"{gain}_sci"]),
                            f"Cer {gain} {config[f'title_{gain}']}",
                            (config["xmin_board"][f"{gain}_cer"],
                             config["xmax_board"][f"{gain}_cer"]),
                            style=STYLE_2D_LOG
                        )

        # Per-event plots
        for gain, calib in GainCalibs:
            config = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE,
                run_number=run_number, beam_energy=benergy)

            var_cer = fersboards.get_energy_sum_name(
                gain=gain, isCer=True, pdsub=True, calib=calib)
            var_sci = fersboards.get_energy_sum_name(
                gain=gain, isCer=False, pdsub=True, calib=calib)

            hist = pm.get_histogram(
                filename, f"hist_{var_cer}_VS_{var_sci}_{suffix}", required=False)
            if not hist:
                continue

            extraObjs = []
            if gain == "Mix":
                f11 = ROOT.TF1("f11", "x", 0, 120)
                f12 = ROOT.TF1("f12", "0.5 * x", 0, 120)
                line_x = ROOT.TLine(benergy, 0, benergy, 120)
                line_y = ROOT.TLine(0, benergy, 120, benergy)
                extraObjs = [f11, f12, line_x, line_y]
                for line in extraObjs:
                    line.SetLineColor(ROOT.kRed)
                    line.SetLineStyle(2)
                    line.SetLineWidth(2)

            pm.plot_2d(
                hist,
                f"FERS_Total_Cer_VS_Sci_{gain}{suffix}",
                f"Sci {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"]),
                f"Cer {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_cer"],
                 config["xmax_total"][f"{gain}_cer"]),
                style=STYLE_2D_LOG,
                extraToDraw=extraObjs if extraObjs else None,
                prepend=True
            )

        return pm.generate_html(f"FERS/{suffix}/EnergySum_Cer_VS_Sci.html", plots_per_row=3)


def makeFERSDRPlots(suffix=""):
    """Plot dual readout plots."""
    filename = f"fers_DR_{suffix}.root"
    filename_sum = f"fers_energy_sum_{suffix}.root"

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_DR_{suffix}")

        gain = "Mix"
        calib = True
        config = getRangesForFERSEnergySums(
            pdsub=True, calib=calib, clip=False, HE=HE,
            run_number=run_number, beam_energy=benergy)
        varname_Cer = fersboards.get_energy_sum_name(
            gain=gain, isCer=True, pdsub=True, calib=calib)
        varname_Sci = fersboards.get_energy_sum_name(
            gain=gain, isCer=False, pdsub=True, calib=calib)

        # Leakage plots
        hists_leak = []
        for var in [varname_Cer, varname_Sci]:
            hist = pm.get_histogram(
                filename, f"hist_{var}_OuterRing_{suffix}", required=False)
            if hist:
                hists_leak.append(hist)

        if hists_leak:
            pm.plot_1d(
                hists_leak,
                f"FERS_Total_{gain}_OuterRing_{suffix}",
                "Leakage Energy",
                (0, 20.0),
                legends=["Cer", "Sci"],
                style=STYLE_CER_SCI
            )

        # Leakage corrected
        hists_leak_corr = []
        for var in [varname_Cer, varname_Sci]:
            hist = pm.get_histogram(
                filename, f"hist_{var}_LeakCorr_{suffix}", required=False)
            if hist:
                hists_leak_corr.append(hist)

        if hists_leak_corr:
            pm.plot_1d(
                hists_leak_corr,
                f"FERS_Total_{gain}_LeakCorr_{suffix}",
                f"Leakage Corrected Energy {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"]),
                yrange=(1, None),
                legends=["Cer", "Sci"],
                style=STYLE_CER_SCI
            )
        pm.add_newline()

        # C/S ratio
        hist = pm.get_histogram(
            filename, f"hist_COverS_{suffix}", required=False)
        if hist:
            pm.plot_1d(hist, f"FERS_Total_{gain}_CerOverSci_{suffix}", "C/S", (0, 2),
                       style=STYLE_CER)

        # fEM
        hist = pm.get_histogram(filename, f"hist_fEM_{suffix}", required=False)
        if hist:
            pm.plot_1d(hist, f"FERS_Total_{gain}_fEM_{suffix}", "f_{EM}", (0, 1.5),
                       style=STYLE_SCI)

        # Sum Cer + Sci
        varname_cersc = varname_Cer.replace("Cer", "CerSci")
        hist = pm.get_histogram(
            filename, f"hist_{varname_cersc}_{suffix}", required=False)
        if hist:
            pm.plot_1d(
                hist,
                f"FERS_Total_{gain}_CerSci_{suffix}",
                f"Cer+Sci {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"] * 2.5),
                yrange=(1, None),
                style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[6])
            )

        # Dual readout methods
        dr_colors = {"": 1, "_method2": 7, "_method3": 8}
        for method, color in dr_colors.items():
            varname_dr = varname_Cer.replace("Cer", "DR") + method
            hist = pm.get_histogram(
                filename, f"hist_{varname_dr}_{suffix}", required=False)
            if hist:
                pm.plot_1d(
                    hist,
                    f"FERS_Total_{gain}_DR{method}_{suffix}",
                    f"DR{method} {gain} {config[f'title_{gain}']}",
                    (config["xmin_total"][f"{gain}_sci"],
                     config["xmax_total"][f"{gain}_sci"]),
                    yrange=(1, None),
                    style=PlotStyle(
                        dology=False, drawoptions="HIST", mycolors=[color])
                )

        # Combined energy plot with fits
        infile_ferssum = pm._get_file(filename_sum)
        hist_cer = infile_ferssum.Get(f"hist_{varname_Cer}_{suffix}")
        hist_sci = infile_ferssum.Get(f"hist_{varname_Sci}_{suffix}")

        infile_dr = pm._get_file(filename)
        varname_dr = varname_Cer.replace("Cer", "DR")
        hist_dr = infile_dr.Get(f"hist_{varname_dr}_{suffix}")
        hist_dr_method2 = infile_dr.Get(f"hist_{varname_dr}_method2_{suffix}")
        hist_dr_method3 = infile_dr.Get(f"hist_{varname_dr}_method3_{suffix}")

        hists_energy = [hist_cer, hist_sci, hist_dr,
                        hist_dr_method2, hist_dr_method3]
        hists_energy = [h for h in hists_energy if h]

        if hists_energy:
            pave = create_pave_text(0.20, 0.65, 0.60, 0.90)
            tf1s = []

            for cat, hist in [("Cer", hist_cer), ("Sci", hist_sci), ("DR", hist_dr),
                              ("DR method2", hist_dr_method2), ("DR method3", hist_dr_method3)]:
                if hist and hist.Integral() > 0:
                    fit_result = hist.Fit("gaus", "S")
                    if fit_result.IsValid():
                        fit_func = hist.GetFunction("gaus")
                        mean = fit_func.GetParameter(1)
                        sigma = fit_func.GetParameter(2)
                        chi2 = fit_func.GetChisquare()
                        ndf = fit_func.GetNDF()
                        pave.AddText(
                            f"{cat}: #mu#pm#sigma = {mean:.1f}#pm{sigma:.1f}, #chi^{{2}}/ndf = {chi2:.1f}/{ndf}")
                        fit_func.SetLineColor(
                            2 if cat == "Cer" else (4 if cat == "Sci" else 1))
                        tf1s.append(fit_func)

            pm.plot_1d(
                [hist_cer, hist_sci, hist_dr, hist_dr_method2, hist_dr_method3],
                f"FERS_Total_{gain}_Energy_{suffix}",
                "Energy [GeV]",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"]),
                style=PlotStyle(dology=False, drawoptions="HIST",
                                mycolors=[2, 4, 1, 7, 8]),
                extraToDraw=[pave] + tf1s,
                prepend=True
            )

        pm.add_newline()

        # 2D correlation plots
        var_sum = varname_Cer.replace("Cer", "CerSci")

        # Sci vs CerSci
        hist = pm.get_histogram(
            filename, f"hist_{varname_Sci}_VS_{var_sum}_{suffix}", required=False)
        if hist:
            pm.plot_2d(
                hist,
                f"FERS_Total_Sci_VS_CerSci_{suffix}",
                f"Sci+Cer {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"] * 2.5),
                f"Sci {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_cer"],
                 config["xmax_total"][f"{gain}_cer"]),
                style=STYLE_2D_LOG
            )

        # Cer vs CerSci
        hist = pm.get_histogram(
            filename, f"hist_{varname_Cer}_VS_{var_sum}_{suffix}", required=False)
        if hist:
            pm.plot_2d(
                hist,
                f"FERS_Total_Cer_VS_CerSci_{suffix}",
                f"Sci+Cer {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_sci"],
                 config["xmax_total"][f"{gain}_sci"] * 2.5),
                f"Cer {gain} {config[f'title_{gain}']}",
                (config["xmin_total"][f"{gain}_cer"],
                 config["xmax_total"][f"{gain}_cer"]),
                style=STYLE_2D_LOG
            )

        # Sci/Cer vs fEM
        var_fem = "fEM"
        for cat, varname, y_key in [("Sci", varname_Sci, "sci"), ("Cer", varname_Cer, "cer")]:
            hist = pm.get_histogram(
                filename, f"hist_{varname}_VS_{var_fem}_{suffix}", required=False)
            if hist and hist.GetEntries() > 0:
                hist_prof = hist.ProfileX()
                hist_prof.SetLineColor(ROOT.kRed)
                hist_prof.SetMarkerColor(ROOT.kRed)
                fit_result = hist_prof.Fit("pol1", "S", "", 0.4, 0.8)

                pave = create_pave_text(0.20, 0.85, 0.60, 0.90)
                fit_func = None
                if fit_result and fit_result.Get() and fit_result.IsValid():
                    fit_func = hist_prof.GetFunction("pol1")
                    p0 = fit_func.GetParameter(0)
                    p1 = fit_func.GetParameter(1)
                    fit_func.SetLineColor(ROOT.kRed)
                    fit_func.SetLineWidth(2)
                    fit_func.SetLineStyle(2)
                    pave.AddText(f"{cat} = {p0:.1f} + {p1:.1f} * f_{{EM}}")

                extra = [pave]
                if fit_func:
                    extra.append(fit_func)

                pm.plot_2d(
                    hist,
                    f"FERS_Total_{cat}_VS_fEM_{suffix}",
                    "f_{EM}",
                    (0, 1.5),
                    f"{cat} {gain} {config[f'title_{gain}']}",
                    (config["xmin_total"][f"{gain}_{y_key}"],
                     config["xmax_total"][f"{gain}_{y_key}"]),
                    style=STYLE_2D_LOG,
                    extraToDraw=extra
                )

        # Cer/Sci vs DR methods
        var_dr_base = varname_Cer.replace("Cer", "DR")
        for dr_method, dr_label in [("", "DR"), ("_method2", "DR method2"), ("_method3", "DR method3")]:
            var_dr = var_dr_base + dr_method
            for cat, var_y, y_key in [("Cer", varname_Cer, "cer"), ("Sci", varname_Sci, "sci")]:
                hist = pm.get_histogram(
                    filename, f"hist_{var_y}_VS_{var_dr}_{suffix}", required=False)
                if hist:
                    pm.plot_2d(
                        hist,
                        f"FERS_Total_{cat}_VS_{var_dr.split('_')[-1]}_{suffix}",
                        f"{dr_label} {gain} {config[f'title_{gain}']}",
                        (config["xmin_total"][f"{gain}_sci"],
                         config["xmax_total"][f"{gain}_sci"]),
                        f"{cat} {gain} {config[f'title_{gain}']}",
                        (config["xmin_total"][f"{gain}_{y_key}"],
                         config["xmax_total"][f"{gain}_{y_key}"]),
                        style=STYLE_2D_LOG
                    )

        return pm.generate_html(f"FERS/{suffix}/EnergySum_DR.html", plots_per_row=4)


def makeFERSEnergyWeightedCenterPlots(suffix=""):
    """Plot energy weighted center distributions."""
    filename = f"fers_energy_weighted_center_{suffix}.root"

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_EnergyWeightedCenter_{suffix}")

        for gain, calib in GainCalibs:
            if gain != "Mix":
                continue

            for cat in ["cer", "sci"]:
                varname_X = fersboards.get_energy_weighted_center_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
                varname_Y = fersboards.get_energy_weighted_center_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)
                style = STYLE_CER if cat == "cer" else STYLE_SCI

                # X and Y 1D plots with stats
                for axis, varname in [("X", varname_X), ("Y", varname_Y)]:
                    hist = pm.get_histogram_by_pattern(
                        filename, varname, suffix, required=False)
                    if hist:
                        pave = create_pave_text(0.20, 0.85, 0.60, 0.90)
                        pave.AddText(
                            f"Center = {hist.GetMean():.2f} +/- {hist.GetRMS():.2f}")

                        pm.plot_1d(
                            hist,
                            f"FERS_Total_{gain}_{cat}_EWC_{axis}{suffix}",
                            f"{cat.capitalize()} {gain} EWC {axis} [cm]",
                            (-15, 15),
                            style=style,
                            extraToDraw=pave
                        )

                # 2D Y vs X
                hist2D = pm.get_histogram(
                    filename, f"hist_{varname_Y}_VS_{varname_X}_{suffix}", required=False)
                if hist2D:
                    pm.plot_2d(
                        hist2D,
                        f"FERS_Total_{gain}_{cat}_EWC_Y_vs_X{suffix}",
                        f"{cat.capitalize()} {gain} EWC X [cm]",
                        (-15, 15),
                        f"{cat.capitalize()} {gain} EWC Y [cm]",
                        (-15, 15),
                        style=PlotStyle(
                            drawoptions=["colz"], addOverflow=False, addUnderflow=False, zmin=1)
                    )

                # 2D with energy
                hprof2D_energy = pm.get_histogram(
                    filename, f"hprof_{varname_Y}_VS_{varname_X}_WithEnergy_{suffix}", required=False)
                if hprof2D_energy and hist2D:
                    zmin = 0.7 * benergy if cat == "sci" else 0.5 * benergy
                    zmax = 1.2 * benergy if cat == "sci" else 1.1 * benergy

                    pm.plot_2d(
                        hprof2D_energy,
                        f"FERS_Total_{gain}_{cat}_EWC_Y_vs_X_WithEnergy{suffix}",
                        f"{cat.capitalize()} {gain} EWC X [cm]",
                        (-15, 15),
                        f"{cat.capitalize()} {gain} EWC Y [cm]",
                        (-15, 15),
                        style=PlotStyle(drawoptions=["colz"], addOverflow=False, addUnderflow=False,
                                        zmin=zmin, zmax=zmax, zlabel=f"Avg Energy {cat} {gain}")
                    )

        return pm.generate_html(
            f"FERS/{suffix}/EnergyWeightedCenter.html",
            plots_per_row=4,
            title="FERS Energy Weighted Center"
        )


def makeFERSEWCvsHodoPlots(suffix=""):
    """Plot EWC vs Hodoscope correlations."""
    filename = f"fers_ewc_vs_hodo_{suffix}.root"
    hodo_min, hodo_max, hodo_nbins = get_ttu_hodo_ranges()

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_EWC_vs_Hodo_{suffix}")

        for gain, calib in GainCalibs:
            if gain != "Mix":
                continue

            for cat in ["cer", "sci"]:
                varname_X = fersboards.get_energy_weighted_center_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=True)
                varname_Y = fersboards.get_energy_weighted_center_name(
                    gain=gain, isCer=(cat == "cer"), pdsub=True, calib=calib, isX=False)

                zmin = 0.7 * benergy if cat == "sci" else 0.5 * benergy
                zmax = 1.2 * benergy if cat == "sci" else 1.1 * benergy

                # EWC vs Hodo
                for axis, varname, hodo_var in [("X", varname_X, "HodoX"), ("Y", varname_Y, "HodoY")]:
                    hist = pm.get_histogram(
                        filename, f"hist_{varname}_VS_{hodo_var}_{suffix}", required=False)
                    if hist:
                        pm.plot_2d(
                            hist,
                            f"FERS_Total_{gain}_{cat}_EWC_{axis}_vs_{hodo_var}{suffix}",
                            f"Hodo {axis} [cm]",
                            (hodo_min, hodo_max),
                            f"{cat.capitalize()} {gain} EWC {axis} [cm]",
                            (-15, 15),
                            style=PlotStyle(
                                drawoptions=["colz"], addOverflow=False, addUnderflow=False, zmin=1)
                        )

                    # profiled energy
                    hprof_energy = pm.get_histogram(
                        filename, f"hprof_{varname}_VS_{hodo_var}_WithEnergy_{suffix}", required=False)
                    if hprof_energy:
                        pm.plot_2d(
                            hprof_energy,
                            f"FERS_Total_{gain}_{cat}_EWC_{axis}_vs_{hodo_var}_WithEnergy{suffix}",
                            f"Hodo {axis} [cm]",
                            (hodo_min, hodo_max),
                            f"{cat.capitalize()} {gain} EWC {axis} [cm]",
                            (-15, 15),
                            style=PlotStyle(drawoptions=["colz"], addOverflow=False, addUnderflow=False,
                                            zmin=zmin, zmax=zmax, zlabel=f"Avg Energy {cat} {gain}")
                        )

                # Hodo Y vs X
                # hist_hodo = pm.get_histogram(
                #    filename, f"hist_HodoY_VS_HodoX_{suffix}", required=False)
                # if hist_hodo:
                #    pm.plot_2d(
                #        hist_hodo,
                #        f"FERS_Total_{gain}_{cat}_HodoY_vs_HodoX{suffix}",
                #        "Hodo X [cm]",
                #        (hodo_min, hodo_max),
                #        "Hodo Y [cm]",
                #        (hodo_min, hodo_max),
                #        style=PlotStyle(
                #            drawoptions=["colz"], addOverflow=False, addUnderflow=False, zmin=1)
                #    )

                hprof_hodo_energy = pm.get_histogram(
                    filename, f"hprof_HodoY_VS_HodoX_WithEnergy_{suffix}", required=False)
                if hprof_hodo_energy:
                    pm.plot_2d(
                        hprof_hodo_energy,
                        f"FERS_Total_{gain}_{cat}_HodoY_vs_HodoX_WithEnergy{suffix}",
                        "Hodo X [cm]",
                        (hodo_min, hodo_max),
                        "Hodo Y [cm]",
                        (hodo_min, hodo_max),
                        style=PlotStyle(drawoptions=["colz"], addOverflow=False, addUnderflow=False,
                                        zmin=zmin, zmax=zmax, zlabel=f"Avg Energy {cat} {gain}")
                    )

        intro_text = """Hodo positions are from TTU hodoscope reconstructed positions, vs FERS energy weighted center positions.
* The right plot is the beam position using hodoscope, and the average energy as a function of position."""

        return pm.generate_html(
            f"FERS/{suffix}/EWC_vs_Hodo.html",
            plots_per_row=5,
            title="TTU Hodoscope and vs FERS EWC",
            intro_text=intro_text
        )


def makeFERSShowerShapePlots(suffix=""):
    """Plot shower shape distributions."""
    filename = f"fers_shower_shape_{suffix}.root"

    with PlotManager(rootdir, plotdir, htmldir, run_number) as pm:
        pm.set_output_dir(f"FERS_ShowerShape_{suffix}")

        for gain, calib in GainCalibs:
            if gain != "Mix":
                continue

            hists_R = []
            for cat in ["cer", "sci"]:
                style = PlotStyle(dology=True, drawoptions="HIST",
                                  mycolors=[2] if cat == "cer" else [4], donormalize=True)

                # X distribution
                hist_X = pm.get_histogram(
                    filename, f"hist_RealX_{gain}_{cat}_{suffix}", required=False)
                if hist_X:
                    pm.plot_1d(
                        hist_X,
                        f"FERS_ShowerShape_RealX_{gain}_{cat}_{suffix}",
                        f"X {cat.capitalize()} {gain} [cm]",
                        (-20, 20),
                        ylabel="Frac. of Energy",
                        yrange=(1e-4, 1),
                        style=style
                    )

                # Y distribution
                hist_Y = pm.get_histogram(
                    filename, f"hist_RealY_{gain}_{cat}_{suffix}", required=False)
                if hist_Y:
                    pm.plot_1d(
                        hist_Y,
                        f"FERS_ShowerShape_RealY_{gain}_{cat}_{suffix}",
                        f"Y {cat.capitalize()} {gain} [cm]",
                        (-20, 20),
                        ylabel="Frac. of Energy",
                        yrange=(1e-4, 1),
                        style=style
                    )

                # R distribution
                hist_R = pm.get_histogram(
                    filename, f"hist_RealR_{gain}_{cat}_{suffix}", required=False)
                if hist_R:
                    hists_R.append(hist_R)

            # Combined R plot
            if hists_R:
                pm.plot_1d(
                    hists_R,
                    f"FERS_ShowerShape_RealR_{gain}_{suffix}",
                    "R [cm]",
                    (0, 25),
                    ylabel="Frac. of Energy",
                    yrange=(1e-4, 1),
                    legends=["Cer", "Sci"],
                    style=PlotStyle(dology=True, drawoptions="HIST", mycolors=[
                                    2, 4], donormalize=True)
                )

                # Ratio plot
                if len(hists_R) == 2:
                    hists_R[0].Scale(
                        1.0 / (hists_R[0].Integral(0, hists_R[0].GetNbinsX() + 1) + 1e-6))
                    hists_R[1].Scale(
                        1.0 / (hists_R[1].Integral(0, hists_R[1].GetNbinsX() + 1) + 1e-6))
                    hist_ratio = hists_R[1].Clone()
                    hist_ratio.Divide(hists_R[0])

                    pm.plot_1d(
                        hist_ratio,
                        f"FERS_ShowerShape_RealR_Cer_over_Sci_{gain}_{suffix}",
                        "R [cm]",
                        (0, 25),
                        ylabel="Cer/Sci",
                        yrange=(0, 1.5),
                        style=PlotStyle(
                            dology=False, drawoptions="HIST", mycolors=[1])
                    )

        return pm.generate_html(f"FERS/{suffix}/ShowerShape.html", plots_per_row=6)


def makeFERSStatsPlots():
    """Plot FERS Cer/Sci ratio statistics."""
    import json

    outdir_plots = f"{plotdir}/FERS_Stats_Cer_ovs_Sci"
    infile_name = f"{rootdir}/fers_stats_cer_ovs_sci.json"

    if not os.path.exists(infile_name):
        print(f"Warning: {infile_name} not found, skipping makeFERSStatsPlots")
        return None

    with open(infile_name, "r") as f:
        stats = json.load(f)

    plots = []
    xmax, xmin = 14, -14
    ymax, ymin = 10, -10
    W_ref, H_ref = 1000, 1100

    for gain, calib in GainCalibs:
        [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm] = visualizeFERSBoards(
            fersboards, stats, suffix=f"Run{run_number}_{gain}_mean", gain=gain)

        h2_Cer_Over_Sci = h2_Cer.Clone(f"h2_Cer_over_Sci_{gain}")
        h2_Cer_Over_Sci.Divide(h2_Sci)
        h2_Cer_3mm_Over_Sci = h2_Cer_3mm.Clone(f"h2_Cer_3mm_over_Sci_{gain}")
        h2_Cer_3mm_Over_Sci.Divide(h2_Sci_3mm)

        for name, hists, extra, zmax in [
            (f"Stats_Cer_{gain}", [h2_Cer, h2_Cer_3mm], "Cer", None),
            (f"Stats_Sci_{gain}", [h2_Sci, h2_Sci_3mm], "Sci", None),
            (f"Stats_Cer_over_Sci_{gain}", [
             h2_Cer_Over_Sci, h2_Cer_3mm_Over_Sci], "Cer / Sci", 1.5),
        ]:
            output_name = f"FERS_Boards_Run{run_number}_{name}"
            DrawHistos(hists, "", xmin, xmax, "iX", ymin, ymax, "iY",
                       output_name, dology=False, drawoptions=["col,text", "col,text"],
                       outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref,
                       extra_text=extra, run_number=run_number, zmin=0, zmax=zmax)
            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS/Stats_Cer_over_Sci.html"
    generate_html(plots, outdir_plots, plots_per_row=3,
                  output_html=output_html)
    return output_html


def makeBoardFits():
    """Fit board energy distributions."""
    suffix = "subtracted_calibd"
    filename = f"{rootdir}/fers_energy_sum_{suffix}.root"

    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run with makeHists=True first.")
        return None

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
                               outdir=outdir, xlabel="Cer # p.e.", **args_cer)
        plots.append(output_name)

        output_name = eventFit(hSci, f"Run{run_number}_Board{board_no}_SciHG",
                               outdir=outdir, xlabel="Sci # p.e.", **args_sci)
        plots.append(output_name)

    output_html = f"{htmldir}/boardfits/index.html"
    generate_html(plots, outdir, plots_per_row=2, output_html=output_html)
    print(f"Generated HTML file: {output_html}")
    return output_html


# ============================================================================
# Main Execution
# ============================================================================

rdf = analysis.get_rdf()

rdfs = OrderedDict()
rdfs["inclusive"] = rdf
if analysis.beam_type in ["e+", "positron", "e-", "electron", "positrons", "electrons"]:
    rdfs["electron"] = analysis.get_particle_analysis("electron")
    rdfs["pion"] = analysis.get_particle_analysis("pion")
else:
    rdfs["proton"] = analysis.get_particle_analysis("proton")
    rdfs["pion"] = analysis.get_particle_analysis("pion")
    rdfs["muon"] = analysis.get_particle_analysis("muon")


def main():
    makeHists = True
    makePlots = True
    outputs_html = {}

    # 1. Booking Phase
    booked_results = {}

    if makeHists:
        for cat, rdf in rdfs.items():
            booked_results[cat] = {
                "energy_sum": makeFERSEnergySumHists(rdf, suffix=cat),
                "cer_vs_sci": makeFERSCervsSciHists(rdf, suffix=cat),
                "dr": makeFERSDRHists(rdf, suffix=cat),
                "weighted_center": makeFERSEnergyWeightedCenterHists(rdf, suffix=cat),
                "ewc_vs_hodo": makeFERSEWCvsHodoHists(rdf, suffix=cat),
                "shower_shape_data": makeFERSShowerShapeHists(rdf, suffix=cat)
            }

        # Collect all proxies
        all_proxies = []
        for cat in booked_results:
            res = booked_results[cat]
            all_proxies.extend(res["energy_sum"])
            all_proxies.extend(res["cer_vs_sci"])
            all_proxies.extend(res["dr"])
            all_proxies.extend(res["weighted_center"])
            all_proxies.extend(res["ewc_vs_hodo"])

            ss_X, ss_Y, ss_R, ss_YX, n_evt_proxy = res["shower_shape_data"]
            for hists_showers in [ss_X, ss_Y, ss_R, ss_YX]:
                for hists_proxies, h_combined_name in hists_showers:
                    all_proxies.extend(hists_proxies)
            all_proxies.append(n_evt_proxy)

        # Execute all at once
        ROOT.RDF.RunGraphs(all_proxies)

        # Save histograms
        for cat in rdfs.keys():
            res = booked_results[cat]

            # Energy sum
            with ROOT.TFile(f"{rootdir}/fers_energy_sum_{cat}.root", "RECREATE") as f:
                for h in res["energy_sum"]:
                    h.Write()

            # Cer vs Sci
            with ROOT.TFile(f"{rootdir}/fers_energy_sum_cer_vs_sci_{cat}.root", "RECREATE") as f:
                for hist in res["cer_vs_sci"]:
                    hist.SetDirectory(f)
                    hist.Write()
            print(f"Saved CER vs SCI histograms for {cat}")

            # DR
            with ROOT.TFile(f"{rootdir}/fers_DR_{cat}.root", "RECREATE") as f:
                for hist in res["dr"]:
                    hist.SetDirectory(f)
                    hist.Write()
            print(f"Saved DR histograms for {cat}")

            # Energy weighted center
            with ROOT.TFile(f"{rootdir}/fers_energy_weighted_center_{cat}.root", "RECREATE") as f:
                for hist in res["weighted_center"]:
                    hist.SetDirectory(f)
                    hist.Write()
            print(f"Saved Energy Weighted Center histograms for {cat}")

            # EWC vs Hodo
            with ROOT.TFile(f"{rootdir}/fers_ewc_vs_hodo_{cat}.root", "RECREATE") as f:
                for hist in res["ewc_vs_hodo"]:
                    hist.SetDirectory(f)
                    hist.Write()
            print(f"Saved EWC vs Hodo histograms for {cat}")

            # Shower shapes
            ss_X, ss_Y, ss_R, ss_YX, n_events_proxy = res["shower_shape_data"]
            nEvts = n_events_proxy.GetValue()

            hists_shower_shape = []
            for hists_showers in [ss_X, ss_Y, ss_R, ss_YX]:
                for hists_proxies, h_combined_name in hists_showers:
                    hists_list = [h.GetValue() for h in hists_proxies]
                    hist_combined = LHistos2Hist(hists_list, h_combined_name)
                    hist_combined.Scale(1.0 / (nEvts + 1e-6))
                    hists_shower_shape.append(hist_combined)

            with ROOT.TFile(f"{rootdir}/fers_shower_shape_{cat}.root", "RECREATE") as f:
                for hist in hists_shower_shape:
                    ROOT.SetOwnership(hist, False)
                    hist.Write()

        # Make plots
        if makePlots:
            for cat in rdfs.keys():
                outputs_html[f"raw_{cat}"] = makeFERSEnergySumPlots(suffix=cat)
                outputs_html[f"cer_vs_sci_raw_{cat}"] = makeFERSCerVsSciPlots(
                    suffix=cat)
                outputs_html[f"dr_{cat}"] = makeFERSDRPlots(suffix=cat)
                outputs_html[f"energy_weighted_center_{cat}"] = makeFERSEnergyWeightedCenterPlots(
                    suffix=cat)
                outputs_html[f"shower_shape_{cat}"] = makeFERSShowerShapePlots(
                    suffix=cat)
                outputs_html[f"ewc_vs_hodo_{cat}"] = makeFERSEWCvsHodoPlots(
                    suffix=cat)

    print("\nGenerated HTML files:")
    for key, html in outputs_html.items():
        print(f"  {key}: {html}")


if __name__ == "__main__":
    main()
