"""
DQM (Data Quality Monitoring) Plot Generation Script.
This script generates various DQM plots for FERS and DRS readout channels.
"""

import ROOT
from channels.channel_map import (build_drs_boards, build_fers_boards,
                                  build_time_reference_channels,
                                  get_mcp_channels, get_service_drs_channels)
from channels.validate_map import DrawDRSBoards, DrawFERSBoards
from plotting.my_function import LHistos2Hist
from configs.plot_config import get_drs_plot_ranges, get_drs_prof_plot_ranges, get_service_drs_plot_ranges
from utils.colors import colors
from utils.parser import get_args
from utils.plot_helper import get_run_paths
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle, STYLE_CER_SCI
from plotting.calox_plot_helper import BoardPlotHelper, create_pave_text, create_board_info_pave
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import number_to_string, round_up_to_1eN, get_channel_var
from utils.visualization import visualizeFERSBoards

auto_timer("Total Execution Time")

setup_root(n_threads=1, batch_mode=True, load_functions=False)

do_detailed_plots = True

args = get_args()
run_number = args.run

DRSBoards = build_drs_boards(run=run_number)
fersboards = build_fers_boards(run=run_number)
time_reference_channels = build_time_reference_channels(run=run_number)
service_drs_channels = get_service_drs_channels(run=run_number)
mcp_channels = get_mcp_channels(run=run_number)

paths = get_run_paths(run_number)

# Common styles
STYLE_BOARD_MULTI = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=colors,
    legendNCols=5,
    legendPos=[0.2, 0.15, 0.9, 0.35]
)

STYLE_2D_LOG = PlotStyle(
    dology=False,
    dologz=True,
    drawoptions="COLZ",
    zmin=1,
    zmax=1e4
)

STYLE_1D = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 6, 4],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90]
)

STYLE_1D_LOG = PlotStyle(
    dology=True,
    drawoptions="HIST",
    mycolors=[2, 6, 4],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90]
)


def make_conditions_plots():
    """Plot FERS board conditions (voltage, current, temperature) vs event."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("Conditions_VS_Event")

        infile = pm._get_file("conditions_vs_event.root")

        # Collect profiles from all boards
        profiles = {
            "SipmHV": [], "SipmI": [], "TempDET": [], "TempFPGA": []
        }
        legends = []

        for fersboard in fersboards.values():
            board_no = fersboard.board_no

            hprof_sipm_hv = infile.Get(
                f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event")
            hprof_sipm_i = infile.Get(
                f"hprof_{fersboard.get_sipm_i_name()}_VS_Event")
            hprof_temp_det = infile.Get(
                f"hprof_{fersboard.get_temp_det_name()}_VS_Event")
            hprof_temp_fpga = infile.Get(
                f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event")

            if not all([hprof_sipm_hv, hprof_sipm_i, hprof_temp_det, hprof_temp_fpga]):
                print(f"Warning: Some profiles not found for board {board_no}")
                continue

            profiles["SipmHV"].append(hprof_sipm_hv)
            profiles["SipmI"].append(hprof_sipm_i)
            profiles["TempDET"].append(hprof_temp_det)
            profiles["TempFPGA"].append(hprof_temp_fpga)
            legends.append(str(board_no))

        if not profiles["SipmHV"]:
            print("Warning: No condition profiles found")
            return None

        n_events = profiles["SipmHV"][0].GetXaxis().GetXmax()

        # Plot all conditions
        plot_configs = [
            ("SipmHV", "Voltage (V)", (26, 30)),
            ("SipmI", "Current (mA)", (0.0, 0.25)),
            ("TempDET", "Temperature (C)", (14, 40)),
            ("TempFPGA", "Temperature (C)", (32, 50)),
        ]

        for var, ylabel, yrange in plot_configs:
            pm.plot_1d(
                profiles[var],
                f"Conditions_{var}_VS_Event",
                "Event",
                (0, n_events),
                ylabel=ylabel,
                yrange=yrange,
                legends=legends,
                style=STYLE_BOARD_MULTI,
                prepend=True
            )

        return pm.generate_html("Conditions/conditions_vs_event.html", plots_per_row=4)


def make_fers_sum_plots():
    """Plot FERS energy sums vs event."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_EnergySum_VS_Event")

        infile = pm._get_file("fers_energysum_vs_event.root")

        # Collect per-board profiles
        profiles = {
            "Cer_HG": [], "Sci_HG": [], "Cer_LG": [], "Sci_LG": []
        }
        legends = []

        for fersboard in fersboards.values():
            board_no = fersboard.board_no

            for cat in ["Cer", "Sci"]:
                for gain in ["HG", "LG"]:
                    key = f"{cat}_{gain}"
                    hprof_name = f"hprof_{fersboard.get_energy_sum_name(gain=gain, isCer=(cat == 'Cer'))}_VS_Event"
                    hprof = infile.Get(hprof_name)
                    if hprof:
                        profiles[key].append(hprof)

            legends.append(str(board_no))

        if not profiles["Cer_HG"]:
            print("Warning: No FERS sum profiles found")
            return None

        n_events = profiles["Cer_HG"][0].GetXaxis().GetXmax()

        # Per-board plots
        board_plots = [
            ("Cer_HG", "Cer FERS Sum HG (ADC)", 5e4),
            ("Sci_HG", "Sci FERS Sum HG (ADC)", 1.6e5),
            ("Cer_LG", "Cer FERS Sum LG (ADC)", 1.4e4),
            ("Sci_LG", "Sci FERS Sum LG (ADC)", 5e4),
        ]

        for key, ylabel, _ in board_plots:
            cat, gain = key.split("_")
            pm.plot_1d(
                profiles[key],
                f"FERS_{cat}_{gain}_EnergySum_VS_Event",
                "Event",
                (0, n_events),
                ylabel=ylabel,
                yrange=(0, None),
                legends=legends,
                style=STYLE_BOARD_MULTI,
                prepend=True
            )

        # Total sum plots
        for gain in ["HG", "LG"]:
            hprof_cer = infile.Get(
                f"hprof_{fersboards.get_energy_sum_name(gain=gain, isCer=True)}_VS_Event")
            hprof_sci = infile.Get(
                f"hprof_{fersboards.get_energy_sum_name(gain=gain, isCer=False)}_VS_Event")

            if hprof_cer and hprof_sci:
                pm.plot_1d(
                    [hprof_cer, hprof_sci],
                    f"FERS_Total_{gain}_EnergySum_VS_Event",
                    "Event",
                    (0, n_events),
                    ylabel=f"FERS Total Sum {gain} (ADC)",
                    yrange=(0, None),
                    legends=["Cer", "Sci"],
                    style=STYLE_CER_SCI,
                    legendPos=[0.6, 0.85, 0.9, 0.9],
                    legendNCols=2,
                    prepend=True
                )

        return pm.generate_html("Conditions/FERS_energysum_vs_event.html", plots_per_row=6)


def make_fers_1d_plots():
    """Plot 1D FERS channel energy distributions."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_1D")

        infile = pm._get_file("fers_all_channels_1d.root")

        for fersboard in fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                hist_c = infile.Get(
                    f"hist_FERS_Board{board_no}_Cer_{s_tower_x}_{s_tower_y}")
                hist_s = infile.Get(
                    f"hist_FERS_Board{board_no}_Sci_{s_tower_x}_{s_tower_y}")

                if not hist_c or not hist_s:
                    print(
                        f"Warning: Histograms not found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y})")
                    continue

                pave = create_board_info_pave(
                    board_no, i_tower_x, i_tower_y,
                    channel_info={
                        "Cer": fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=True).channel_no,
                        "Sci": fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=False).channel_no
                    }
                )

                pm.plot_1d(
                    [hist_c, hist_s],
                    f"Energy_Board{board_no}_iTowerX{s_tower_x}_iTowerY{s_tower_y}",
                    "Energy HG",
                    (0, 1000),
                    ylabel="Counts",
                    yrange=(1, 1e5),
                    legends=["Cer", "Sci"],
                    style=PlotStyle(
                        dology=True, drawoptions="HIST", mycolors=[2, 4]),
                    extraToDraw=pave
                )

        return pm.generate_html("FERS/ChannelADC.html")


def make_fers_stats_plots(include_pedestals=False):
    """Plot FERS channel statistics as 2D board maps."""
    import json

    output_htmls = []

    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_Stats")

        # Load stats
        with open(f"{paths['root']}/fers_stats.json", "r") as f:
            stats = json.load(f)

        pedestals_hg = {}
        pedestals_lg = {}
        if include_pedestals:
            with open(f"{paths['root']}/fers_pedestals_hg.json", "r") as f:
                pedestals_hg = json.load(f)
            with open(f"{paths['root']}/fers_pedestals_lg.json", "r") as f:
                pedestals_lg = json.load(f)

        # Organize stats by type
        valuemaps = {}
        for gain in ["HG", "LG"]:
            for stat in ["mean", "max", "satfreq", "pedestal"]:
                valuemaps[f"{gain}_{stat}"] = {}

        for channel_name, (vmean, vmax, v_satfreq) in stats.items():
            gain = "HG" if "energyHG" in channel_name else "LG"
            valuemaps[f"{gain}_mean"][channel_name] = vmean
            valuemaps[f"{gain}_max"][channel_name] = vmax
            valuemaps[f"{gain}_satfreq"][channel_name] = v_satfreq
            if gain == "HG":
                valuemaps[f"{gain}_pedestal"][channel_name] = pedestals_hg.get(
                    channel_name, 0.)
            else:
                valuemaps[f"{gain}_pedestal"][channel_name] = pedestals_lg.get(
                    channel_name, 0.)

        # Create visualizations
        board_hists = {}
        for gain in ["HG", "LG"]:
            for stat in ["mean", "max", "satfreq", "pedestal"]:
                key = f"{gain}_{stat}"
                cer_hists, sci_hists = visualizeFERSBoards(
                    fersboards, valuemaps[key], suffix=f"Run{run_number}_{key}", gain=gain
                )
                board_hists[f"{key}_Cer"] = cer_hists
                board_hists[f"{key}_Sci"] = sci_hists

        helper = BoardPlotHelper(pm)

        # Plot mean values
        pm.reset_plots()
        for gain in ["HG", "LG"]:
            helper.plot_cer_sci_pair(
                board_hists[f"{gain}_mean_Cer"],
                board_hists[f"{gain}_mean_Sci"],
                f"FERS_Boards_Run{run_number}_Stats_{gain}_mean",
                zmin=0, zmax=8000
            )
        output_htmls.append(pm.generate_html(
            "FERS/Stat/Channel_Mean.html", plots_per_row=2))

        # Plot max values
        pm.reset_plots()
        for gain in ["HG", "LG"]:
            helper.plot_cer_sci_pair(
                board_hists[f"{gain}_max_Cer"],
                board_hists[f"{gain}_max_Sci"],
                f"FERS_Boards_Run{run_number}_Stats_{gain}_max",
                zmin=0, zmax=8000
            )
        output_htmls.append(pm.generate_html(
            "FERS/Stat/Channel_Max.html", plots_per_row=2))

        # Plot saturation frequency
        pm.reset_plots()
        for gain in ["HG", "LG"]:
            helper.plot_cer_sci_pair(
                board_hists[f"{gain}_satfreq_Cer"],
                board_hists[f"{gain}_satfreq_Sci"],
                f"FERS_Boards_Run{run_number}_Stats_{gain}_satfreq",
                zmin=0, zmax=1, nTextDigits=2
            )
        output_htmls.append(pm.generate_html(
            "FERS/Stat/Channel_SatFreq.html", plots_per_row=2))

        # Plot pedestals
        pm.reset_plots()
        for gain in ["HG", "LG"]:
            helper.plot_cer_sci_pair(
                board_hists[f"{gain}_pedestal_Cer"],
                board_hists[f"{gain}_pedestal_Sci"],
                f"FERS_Boards_Run{run_number}_Stats_{gain}_pedestal",
                zmin=100, zmax=300, nTextDigits=0
            )
        output_htmls.append(pm.generate_html(
            "FERS/Stat/Channel_Pedestal.html", plots_per_row=2))

    return output_htmls


def make_fers_max_value_plots():
    """Plot FERS max value distributions."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_MaxValues")

        infile = pm._get_file("fers_max_values.root")
        xmin, xmax = 7500, 8500

        # Collect board histograms
        board_hists = {"Cer_HG": [], "Sci_HG": [], "Cer_LG": [], "Sci_LG": []}
        legends = []

        for fersboard in fersboards.values():
            board_no = fersboard.board_no
            for cat in ["Cer", "Sci"]:
                for gain in ["HG", "LG"]:
                    hist = infile.Get(
                        f'hist_{fersboard.get_energy_max_name(gain=gain, isCer=(cat == "Cer"))}')
                    if hist:
                        board_hists[f"{cat}_{gain}"].append(hist)
            legends.append(str(board_no))

        # Per-board plots
        board_style = PlotStyle(
            dology=True, drawoptions="HIST", mycolors=colors,
            legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9]
        )

        for cat in ["Cer", "Sci"]:
            for gain in ["HG", "LG"]:
                key = f"{cat}_{gain}"
                pm.plot_1d(
                    board_hists[key],
                    f"FERS_Boards_{cat}Energy{gain}_max",
                    f"{gain} {cat} Max (Board)",
                    (xmin, xmax),
                    yrange=(1, None),
                    legends=legends,
                    style=board_style
                )

        # Total plots with saturation fractions
        v_sat = 8000
        for gain in ["HG", "LG"]:
            hist_cer = infile.Get(
                f'hist_{fersboards.get_energy_max_name(gain=gain, isCer=True)}')
            hist_sci = infile.Get(
                f'hist_{fersboards.get_energy_max_name(gain=gain, isCer=False)}')

            if hist_cer and hist_sci:
                frac_sci = hist_sci.Integral(hist_sci.FindBin(
                    v_sat), 100000) / (hist_sci.Integral(0, 100000) + 1e-6)
                frac_cer = hist_cer.Integral(hist_cer.FindBin(
                    v_sat), 100000) / (hist_cer.Integral(0, 100000) + 1e-6)

                pave = create_pave_text(0.20, 0.63, 0.90, 0.72)
                pave.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
                pave.AddText(f"Sat Frac Cer : {frac_cer:.3f}")

                pm.plot_1d(
                    [hist_cer, hist_sci],
                    f"FERS_Energy{gain}_max",
                    f"{gain} Max (All Boards)",
                    (xmin, xmax),
                    yrange=(1, None),
                    legends=["Cer", "Sci"],
                    style=PlotStyle(
                        dology=True, drawoptions="HIST", mycolors=[2, 4]),
                    legendPos=[0.30, 0.75, 0.50, 0.9],
                    extraToDraw=pave
                )

        return pm.generate_html("FERS/Stat/Channel_Max_1D.html", plots_per_row=2)


def make_fers_2d_plots():
    """Plot 2D FERS HG vs LG correlations."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_2D")

        infile = pm._get_file("fers_all_channels_2d.root")

        for fersboard in fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                for var in ["Cer", "Sci"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist_name = f"hist_FERS_Board{board_no}_{var}_{s_tower_x}_{s_tower_y}_hg_VS_lg"
                    hist = infile.Get(hist_name)

                    if not hist:
                        continue

                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={var: chan.channel_no},
                        position=(0.20, 0.70, 0.60, 0.90)
                    )

                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_tower_x}_{s_tower_y}_hg_VS_lg",
                        "HG", (0, 9000),
                        "LG", (0, 1500),
                        style=STYLE_2D_LOG,
                        extraToDraw=pave
                    )

        return pm.generate_html("FERS/LG_vs_HG.html", plots_per_row=4)


def track_fers_plots():
    """Plot FERS output vs event number."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("FERS_VS_Event")

        infile = pm._get_file("fers_all_channels_2D_VS_event.root")

        for fersboard in fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                for var in ["Cer", "Sci"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist_name = f"hist_FERS_Board{board_no}_{var}_VS_Event_{s_tower_x}_{s_tower_y}"
                    hist = infile.Get(hist_name)

                    if not hist:
                        print(f"Warning: Histogram {hist_name} not found")
                        continue

                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={var: chan.channel_no},
                        position=(0.20, 0.70, 0.60, 0.90)
                    )

                    n_events = hist.GetXaxis().GetXmax()

                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_tower_x}_{s_tower_y}_VS_Event",
                        "Event", (0, n_events),
                        f"{var} Energy HG", (1, 1e5),
                        style=PlotStyle(dology=True, dologz=True,
                                        drawoptions="COLZ", zmin=1, zmax=1e4),
                        extraToDraw=pave
                    )

        return pm.generate_html("FERS_VS_Event/index.html", plots_per_row=4)


def make_drs_vs_ts_plots(use_calibrated_ts=False, use_cfd_aligned_ts=False):
    """Plot DRS output vs time slice."""
    suffix = "Calib" if use_calibrated_ts else "CFD" if use_cfd_aligned_ts else "Raw"
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRS_VS_TS")

        infile = pm._get_file("drs_vs_ts.root")

        for _, DRSBoard in DRSBoards.items():
            for chan in DRSBoard:
                channelName = chan.get_channel_name(blsub=True)
                if use_calibrated_ts:
                    hist_name = f"hist_{channelName}_VS_AlignedTS"
                elif use_cfd_aligned_ts:
                    hist_name = f"hist_{channelName}_VS_CFDAlignedTS"
                else:
                    hist_name = f"hist_{channelName}_VS_TS"
                hist = infile.Get(hist_name)

                if not hist:
                    print(f"Warning: Histogram {hist_name} not found")
                    continue

                ymin_tmp, ymax_tmp = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm, is_reference=chan.is_reference
                )

                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")

                var = get_channel_var(chan)

                xaxis = "Aligned TS" if use_calibrated_ts else "CFD Aligned TS" if use_cfd_aligned_ts else "TS"

                pm.plot_2d(
                    hist,
                    f"DRS_ADC_VS_TS_{channelName}_{suffix}_{var}",
                    xaxis, (0, 1024),
                    "DRS Output", (ymin_tmp, ymax_tmp),
                    style=STYLE_2D_LOG,
                    extraToDraw=pave,
                    extra_text=var
                )

        return pm.generate_html(f"DRS/DRS_vs_TS_{suffix}.html", plots_per_row=9)


def make_drs_prof_vs_ts_plots():
    """Plot DRS output profile vs time slice."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRS_Prof_VS_TS")

        infile = pm._get_file("drs_vs_ts.root")

        for _, DRSBoard in DRSBoards.items():
            for chan in DRSBoard:
                channelName = chan.get_channel_name(blsub=True)
                hist_name = f"prof_{channelName}_VS_AlignedTS"
                hist_alignedTS = infile.Get(hist_name).ProjectionX()
                hist_name = f"prof_{channelName}_VS_TS"
                hist = infile.Get(hist_name).ProjectionX()
                hist_name = f"prof_{channelName}_VS_CFDAlignedTS"
                hist_cfdAlignedTS = infile.Get(hist_name).ProjectionX()

                if not hist:
                    print(f"Warning: Histogram {hist_name} not found")
                    continue

                pave = create_pave_text(0.20, 0.85, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")

                xaxis = "TS"

                ymin_tmp, ymax_tmp = get_drs_prof_plot_ranges(
                    subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm, is_reference=chan.is_reference
                )

                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")

                var = get_channel_var(chan)

                pm.plot_1d(
                    [hist, hist_alignedTS, hist_cfdAlignedTS],
                    f"DRS_ADC_Prof_VS_TS_{channelName}_{var}",
                    xaxis, (0, 1024),
                    "Mean DRS Output", (ymin_tmp, ymax_tmp),
                    legends=["Raw TS", "Aligned TS", "CFD Aligned TS"],
                    legendPos=[0.55, 0.75, 0.90, 0.90],
                    style=STYLE_1D,
                    extraToDraw=pave,
                    extra_text=var
                )

        return pm.generate_html(f"DRS/DRS_Prof_vs_TS.html", plots_per_row=9)


def make_drs_sum_plots():
    # plot DRS sum distributions
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRS_Stats")

        infile = pm._get_file("drs_stats.root")

        for _, DRSBoard in DRSBoards.items():
            for chan in DRSBoard:
                channelSumName = chan.get_channel_sum_name()
                hist_name = f"hist_{channelSumName}"
                hist = infile.Get(hist_name)

                channelName_blsub = chan.get_channel_name(blsub=True)
                hist_name_cfd = f"hist_{channelName_blsub}_cfd"
                hist_cfd = infile.Get(hist_name_cfd)

                if not hist or not hist_cfd:
                    print(
                        f"Warning: Histogram {hist_name} or {hist_name_cfd} not found")
                    continue

                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")

                var = get_channel_var(chan)

                pm.plot_1d(
                    [hist, hist_cfd],
                    f"DRS_Sum_{channelSumName}_{var}",
                    "DRS Sum", (0, 6000),
                    "Counts", (1, None),
                    legends=["Sum", "CFD Integral"],
                    legendPos=[0.55, 0.80, 0.90, 0.90],
                    style=STYLE_1D_LOG,
                    extraToDraw=pave,
                    extra_text=var
                )
                # pm.plot_1d(
                #    hist_cfd,
                #    f"DRS_Sum_CFD_{channelSumName}",
                #    "DRS CFD Integral", (0, 5000),
                #    "Counts", (1, None),
                #    style=STYLE_1D_LOG,
                #    extraToDraw=pave
                # )

        return pm.generate_html("DRS/DRS_Sum.html", plots_per_row=9)


def make_drs_peak_plots():
    # plot DRS peak distributions
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRS_Stats")

        infile = pm._get_file("drs_stats.root")

        for _, DRSBoard in DRSBoards.items():
            for chan in DRSBoard:
                channelPeakName = chan.get_channel_peak_name()
                hist_name = f"hist_{channelPeakName}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(f"Warning: Histogram {hist_name} not found")
                    continue

                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                var = get_channel_var(chan)

                pm.plot_1d(
                    hist,
                    f"DRS_Peak_{channelPeakName}_{var}",
                    "DRS Peak", (0, 800),
                    "Counts", (1, None),
                    style=STYLE_1D_LOG,
                    extraToDraw=pave,
                    extra_text=var
                )

        return pm.generate_html("DRS/DRS_Peak.html", plots_per_row=9)


def make_drs_time_plots():
    # plot DRS peak time distributions
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRS_Stats")

        infile = pm._get_file("drs_stats.root")

        for _, DRSBoard in DRSBoards.items():
            for chan in DRSBoard:
                channelName_blsub = chan.get_channel_name(blsub=True)
                hist_name = f"hist_{channelName_blsub}_RelPeakTS"
                hist = infile.Get(hist_name)

                hist_name_cfd = f"hist_{channelName_blsub}_RelCFDTS"
                hist_cfd = infile.Get(hist_name_cfd)

                if not hist or not hist_cfd:
                    print(
                        f"Warning: Histogram {hist_name} or {hist_name_cfd} not found")
                    continue

                pave = create_pave_text(0.20, 0.85, 0.60, 0.90)
                pave.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")

                pm.plot_1d(
                    [hist, hist_cfd],
                    f"DRS_Time_{channelName_blsub}",
                    "Pulse TS", (0, 1024),
                    "Counts", (1, None),
                    legends=["Peak", "0.2 CFD"],
                    legendPos=[0.55, 0.80, 0.90, 0.90],
                    style=STYLE_1D_LOG,
                    extraToDraw=pave
                )

        return pm.generate_html("DRS/DRS_Time.html", plots_per_row=9)


def make_drs_peak_ts_plots():
    """Plot DRS peak time slice distributions."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRSPeakTS")

        infile = pm._get_file("drspeakts.root")

        hists_cer = []
        hists_sci = []

        for _, DRSBoard in DRSBoards.items():
            board_no = DRSBoard.board_no
            for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                hists = {}
                channel_nos = {}
                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist_name = f"hist_DRSPeakTS_{var}_{s_tower_x}_{s_tower_y}"
                    hist = infile.Get(hist_name)
                    hists[var] = hist
                    channel_nos[var] = chan.channel_no if chan else -1

                if not hists["Cer"] or not hists["Sci"]:
                    print(
                        f"Warning: Histograms not found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y})")
                    continue

                hists_cer.append(hists["Cer"])
                hists_sci.append(hists["Sci"])

                if do_detailed_plots:
                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={
                            "Cer": channel_nos["Cer"], "Sci": channel_nos["Sci"]}
                    )
                    # Add group info
                    pave_lines = list(pave.GetListOfLines())
                    if chan:
                        pave.AddText(f"Group: {chan.group_no}")

                    pm.plot_1d(
                        [hists["Cer"], hists["Sci"]],
                        f"hist_DRSPeakTS_{s_tower_x}_{s_tower_y}",
                        "Peak TS", (400, 600),
                        yrange=(1, None),
                        legends=["Cer", "Sci"],
                        style=STYLE_CER_SCI,
                        extraToDraw=pave
                    )

        # Summary plot
        if hists_cer and hists_sci:
            hist_cer_combined = LHistos2Hist(
                hists_cer, "hist_DRSPeakTS_Cer_Combined")
            hist_sci_combined = LHistos2Hist(
                hists_sci, "hist_DRSPeakTS_Sci_Combined")

            pm.plot_1d(
                [hist_cer_combined, hist_sci_combined],
                "DRS_PeakTS_Combined",
                "Peak TS", (400, 600),
                yrange=(1, None),
                legends=["Cer", "Sci"],
                style=STYLE_CER_SCI,
                prepend=True
            )

        return pm.generate_html("DRS/DRS_PeakTS.html", plots_per_row=4)


def make_drs_peak_ts_cer_vs_sci_plots():
    """Plot DRS peak time slice Cer vs Sci correlations."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRSPeakTSCerVSSci")

        infile = pm._get_file("drspeakts.root")

        # Diagonal line for reference
        diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
        diagonal_line.SetLineStyle(2)
        diagonal_line.SetLineWidth(1)
        diagonal_line.SetLineColor(ROOT.kRed)

        hists = []

        for _, DRSBoard in DRSBoards.items():
            for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                hist_name = f"hist_DRSPeakTS_Cer_VS_Sci_{s_tower_x}_{s_tower_y}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(f"Warning: Histogram {hist_name} not found")
                    continue

                hists.append(hist)

                if do_detailed_plots:
                    pm.plot_2d(
                        hist,
                        f"DRSPeakTS_Cer_VS_Sci_{s_tower_x}_{s_tower_y}",
                        "Sci Peak TS", (400, 600),
                        "Cer Peak TS", (400, 600),
                        style=PlotStyle(
                            dologz=True, drawoptions="COLZ", zmin=1, zmax=1e2, addOverflow=False),
                        extraToDraw=diagonal_line
                    )

        # Summary plot
        if hists:
            hcombined = LHistos2Hist(
                hists, "hist_DRSPeakTS_Cer_VS_Sci_Combined")
            pm.plot_2d(
                hcombined,
                "DRS_PeakTS_Cer_VS_Sci_Combined",
                "Sci Peak TS", (400, 600),
                "Cer Peak TS", (400, 600),
                style=PlotStyle(dologz=True, drawoptions="COLZ",
                                zmin=1, zmax=1e2, addOverflow=False),
                extraToDraw=diagonal_line,
                prepend=True
            )

        return pm.generate_html("DRS/DRS_PeakTS_Cer_VS_Sci.html", plots_per_row=4)


def compare_time_reference_plots(do_subtract_median=False):
    """Plot time reference channel distributions."""
    suffix = "_subtractMedian" if do_subtract_median else ""
    ymin, ymax = (-2500, 500) if do_subtract_median else (500, 2500)

    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("TimeReference")

        infile = pm._get_file("time_reference_channels.root")

        for chan_name in time_reference_channels:
            hist_name = f"hist_{chan_name}{suffix}"
            hist = infile.Get(hist_name)

            if not hist:
                print(f"Warning: Histogram {hist_name} not found")
                continue

            pave = create_pave_text(0.20, 0.70, 0.60, 0.90)
            pave.AddText(f"{chan_name}")

            pm.plot_2d(
                hist,
                f"TimeReference_{chan_name}{suffix}",
                "Time Slice", (0, 1024),
                "Counts", (ymin, ymax),
                style=STYLE_2D_LOG,
                extraToDraw=pave
            )

        return pm.generate_html(f"TimeReference{suffix}/index.html", plots_per_row=2)


def compare_service_drs_plots():
    """Plot service DRS channel distributions."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("ServiceDRS")

        infile = pm._get_file("service_drs_channels.root")

        for det_name, chan_name in service_drs_channels.items():
            hist_name = f"hist_{chan_name}_blsub"
            hist = infile.Get(hist_name)

            if not hist:
                print(f"Warning: Histogram {hist_name} not found")
                continue

            pave = create_pave_text(0.60, 0.20, 0.90, 0.30)
            pave.AddText(f"{det_name}")

            ymin, ymax = get_service_drs_plot_ranges(
                chan_name, subtractMedian=True)

            pm.plot_2d(
                hist,
                f"ServiceDRS_{chan_name}",
                "Time Slice", (0, 1024),
                "ADC", (ymin, ymax),
                style=STYLE_2D_LOG,
                extraToDraw=pave
            )

        return pm.generate_html("ServiceDRS/detectors.html", plots_per_row=2)


def compare_mcp_plots():
    """Plot MCP channel distributions."""
    ymin, ymax = -1500, 500

    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("MCP")

        infile = pm._get_file("mcp_channels.root")

        for det_name, channels in mcp_channels.items():
            for chan_name in channels:
                hist_name = f"hist_{chan_name}_blsub"
                hist = infile.Get(hist_name)

                if not hist:
                    print(f"Warning: Histogram {hist_name} not found")
                    continue

                pave = create_pave_text(0.60, 0.20, 0.90, 0.30)
                pave.AddText(f"{det_name}")

                pm.plot_2d(
                    hist,
                    f"MCP_{chan_name}",
                    "Time Slice", (0, 1024),
                    "Counts", (ymin, ymax),
                    style=STYLE_2D_LOG,
                    extraToDraw=pave
                )

        return pm.generate_html("ServiceDRS/MCP.html", plots_per_row=4)


def make_drs_sum_vs_fers_plots():
    """Plot DRS sum vs FERS correlations."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRSSum_VS_FERS")

        infile = pm._get_file("drssum_vs_fers.root")

        xymax = {"Cer": (20000, 8500), "Sci": (30000, 8500)}
        xymax_LG = {"Cer": (20000, 2000), "Sci": (30000, 4000)}

        hists = []

        for _, DRSBoard in DRSBoards.items():
            for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                for var in ["Cer", "Sci"]:
                    for gain in ["FERS", "FERSLG"]:
                        histname = f"hist_DRSSum_VS_{gain}_{var}_{s_tower_x}_{s_tower_y}"
                        hist = infile.Get(histname)

                        if not hist:
                            print(f"Warning: Histogram {histname} not found")
                            continue

                        hists.append(hist)

                        if do_detailed_plots:
                            zmax = round_up_to_1eN(
                                hist.Integral(0, 10000, 0, 10000))
                            tmp = xymax[var] if gain == "FERS" else xymax_LG[var]

                            pm.plot_2d(
                                hist,
                                histname.replace("hist_", ""),
                                gain, (0, tmp[1]),
                                "DRSSum", (0, tmp[0]),
                                style=PlotStyle(
                                    dologz=True, drawoptions="COLZ", zmin=1, zmax=zmax),
                                extra_text=var
                            )

        # Summary plots
        for var in ["Cer", "Sci"]:
            for gain in ["FERS", "FERSLG"]:
                histname = f"hist_DRSSum_VS_{gain}_{var}_Combined"
                hists_to_combine = [
                    h for h in hists if f"_{gain}_{var}_" in h.GetName()]

                if hists_to_combine:
                    hist_combined = LHistos2Hist(hists_to_combine, histname)
                    print(
                        f"Combining histograms for {histname}, count: {len(hists_to_combine)}")

                    zmax = round_up_to_1eN(
                        hist_combined.Integral(0, 10000, 0, 10000))
                    tmp = xymax[var] if gain == "FERS" else xymax_LG[var]

                    pm.plot_2d(
                        hist_combined,
                        histname.replace("hist_", ""),
                        gain, (0, tmp[1]),
                        "DRSSum", (0, tmp[0]),
                        style=PlotStyle(
                            dologz=True, drawoptions="COLZ", zmin=1, zmax=zmax),
                        extra_text=var,
                        prepend=True
                    )

        return pm.generate_html("DRS_VS_FERS/DRSSum_vs_FERS.html", plots_per_row=4)


def make_drs_peak_vs_fers_plots():
    """Plot DRS peak vs FERS correlations."""
    with PlotManager(paths["root"], paths["plots"], paths["html"], run_number) as pm:
        pm.set_output_dir("DRSPeak_VS_FERS")

        infile = pm._get_file("drspeak_vs_fers.root")

        for _, DRSBoard in DRSBoards.items():
            board_no = DRSBoard.board_no
            for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
                s_tower_x = number_to_string(i_tower_x)
                s_tower_y = number_to_string(i_tower_y)

                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if not chan:
                        print(
                            f"Warning: Channel not found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y}), Var {var}")
                        continue

                    _, ymax = get_drs_plot_ranges(
                        subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm
                    )

                    histname = f"hist_DRSPeak_VS_FERS_{var}_{s_tower_x}_{s_tower_y}"
                    hist = infile.Get(histname)

                    if not hist:
                        print(f"Warning: Histogram {histname} not found")
                        continue

                    pm.plot_2d(
                        hist,
                        histname.replace("hist_", ""),
                        "FERS ADC", (0, 9000),
                        "DRS Peak", (0, ymax),
                        style=PlotStyle(
                            dologz=True, drawoptions="COLZ", zmin=1, zmax=None),
                        extra_text=var
                    )

        return pm.generate_html("DRS_VS_FERS/DRSPeak_vs_FERS.html", plots_per_row=4)


def main():
    print("Starting DQM Plot Generation...\n")

    # Task Registry: (Label, Function)
    plot_tasks = [
        ("Conditions", make_conditions_plots),
        ("FERS Sum", make_fers_sum_plots),
        ("FERS Mapping", lambda: DrawFERSBoards(run=run_number)),
        ("DRS Mapping", lambda: DrawDRSBoards(run=run_number)),
        ("FERS Stats", lambda: make_fers_stats_plots(include_pedestals=True)),
        # ("DRS Peak TS", make_drs_peak_ts_plots),
        # ("DRS Peak TS Cer vs Sci", make_drs_peak_ts_cer_vs_sci_plots),
        # ("Service DRS", compare_service_drs_plots),
        # ("MCP", compare_mcp_plots),
        ("FERS Max Values", make_fers_max_value_plots),
        # ("DRS Sum vs FERS", make_drs_sum_vs_fers_plots),
    ]

    if do_detailed_plots:
        plot_tasks.extend([
            # ("FERS 1D", make_fers_1d_plots),
            ("DRS vs TS", make_drs_vs_ts_plots),
            ("DRS vs TS Calib", lambda: make_drs_vs_ts_plots(use_calibrated_ts=True)),
            ("DRS vs TS CFD", lambda: make_drs_vs_ts_plots(use_cfd_aligned_ts=True)),
            ("DRS Prof vs TS", make_drs_prof_vs_ts_plots),
            ("DRS Sum", make_drs_sum_plots),
            ("DRS Peak", make_drs_peak_plots),
            ("DRS Time", make_drs_time_plots)
            # ("DRS Peak vs FERS", make_drs_peak_vs_fers_plots),
        ])

    output_htmls = {}
    for label, func in plot_tasks:
        print(f"Generating {label} plots...")
        output_htmls[label] = func()

    print("\n" + "*" * 30)
    print("Plot Generation Summary:")
    for label, url in output_htmls.items():
        if isinstance(url, str):
            print(f"✅ {label} plots: {url}")
        elif isinstance(url, list):
            print(f"✅ {label} plots:")
            for u in url:
                print(f"   - {u}")
        else:
            print(f"✅ {label} plots: (output type: {type(url)})")


if __name__ == "__main__":
    main()
