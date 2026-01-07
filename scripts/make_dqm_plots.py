import ROOT
from channels.channel_map import (build_drs_boards, build_fers_boards,
                                  build_hodo_trigger_channels, build_hodo_pos_channels,
                                  build_time_reference_channels, get_downstream_muon_channel,
                                  get_mcp_channels, get_service_drs_channels,
                                  get_upstream_veto_channel)
from channels.validate_map import DrawDRSBoards, DrawFERSBoards
from plotting.my_function import DrawHistos, LHistos2Hist
from configs.plot_ranges import get_drs_plot_ranges, get_service_drs_plot_ranges
from utils.colors import colors
from utils.html_generator import generate_html
from utils.parser import get_args
from utils.plot_helper import get_run_paths
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import number_to_string, round_up_to_1eN
from utils.visualization import visualizeFERSBoards
auto_timer("Total Execution Time")

setup_root(n_threads=1, batch_mode=True, load_functions=False)

args = get_args()
run_number = args.run

DRSBoards = build_drs_boards(run=run_number)
fersboards = build_fers_boards(run=run_number)
time_reference_channels = build_time_reference_channels(run=run_number)
hodo_trigger_channels = build_hodo_trigger_channels(run=run_number)
hodo_pos_channels = build_hodo_pos_channels(run=run_number)
upstream_veto_channel = get_upstream_veto_channel(run=run_number)
downstream_muon_channel = get_downstream_muon_channel(run=run_number)
service_drs_channels = get_service_drs_channels(run=run_number)
mcp_channels = get_mcp_channels(run=run_number)

paths = get_run_paths(run_number)


def makeConditionsPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/Conditions_VS_Event"
    infile_name = f"{paths['root']}/conditions_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofiles_SipmHV = []
    hprofiles_SipmI = []
    hprofiles_TempDET = []
    hprofiles_TempFPGA = []

    legends = []

    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        hprof_SipmHV_name = f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event"
        hprof_SipmI_name = f"hprof_{fersboard.get_sipm_i_name()}_VS_Event"
        hprof_TempDET_name = f"hprof_{fersboard.get_temp_det_name()}_VS_Event"
        hprof_TempFPGA_name = f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event"

        hprof_SipmHV = infile.Get(hprof_SipmHV_name)
        hprof_SipmI = infile.Get(hprof_SipmI_name)
        hprof_TempDET = infile.Get(hprof_TempDET_name)
        hprof_TempFPGA = infile.Get(hprof_TempFPGA_name)

        if not (hprof_SipmHV and hprof_SipmI and hprof_TempDET and hprof_TempFPGA):
            print(
                f"Warning: Profiles {hprof_SipmHV_name}, {hprof_SipmI_name}, {hprof_TempDET_name}, or {hprof_TempFPGA_name} not found in {infile_name}")
            continue

        hprofiles_SipmHV.append(hprof_SipmHV)
        hprofiles_SipmI.append(hprof_SipmI)
        hprofiles_TempDET.append(hprof_TempDET)
        hprofiles_TempFPGA.append(hprof_TempFPGA)

        legends.append(str(board_no))

    # Draw the profiles
    n_events = hprofiles_SipmHV[0].GetXaxis().GetXmax()
    legendPos = [0.3, 0.7, 0.9, 0.9]
    output_name = "Conditions_SipmHV_VS_Event"
    DrawHistos(hprofiles_SipmHV, legends, 0, n_events, "Event", 26, 30, "Voltage (V)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(0, output_name + ".png")

    output_name = "Conditions_SipmI_VS_Event"
    DrawHistos(hprofiles_SipmI, legends, 0, n_events, "Event", 0.0, 0.3, "Current (mA)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(1, output_name + ".png")

    output_name = "Conditions_TempDET_VS_Event"
    DrawHistos(hprofiles_TempDET, legends, 0, n_events, "Event", 14, 35, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(2, output_name + ".png")

    output_name = "Conditions_TempFPGA_VS_Event"
    DrawHistos(hprofiles_TempFPGA, legends, 0, n_events, "Event", 32, 50, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(3, output_name + ".png")

    output_html = f"{paths['html']}/Conditions/conditions_vs_event.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSSumPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/FERS_EnergySum_VS_Event"
    infile_name = f"{paths['root']}/fers_energysum_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofiles_Cer_HG_sum = []
    hprofiles_Sci_HG_sum = []
    hprofiles_Cer_LG_sum = []
    hprofiles_Sci_LG_sum = []

    legends = []

    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        hprof_Cer_HG_sum_name = f"hprof_{fersboard.get_energy_sum_name(gain='HG', isCer=True)}_VS_Event"
        hprof_Sci_HG_sum_name = f"hprof_{fersboard.get_energy_sum_name(gain='HG', isCer=False)}_VS_Event"
        hprof_Cer_LG_sum_name = f"hprof_{fersboard.get_energy_sum_name(gain='LG', isCer=True)}_VS_Event"
        hprof_Sci_LG_sum_name = f"hprof_{fersboard.get_energy_sum_name(gain='LG', isCer=False)}_VS_Event"

        hprof_Cer_HG_sum = infile.Get(hprof_Cer_HG_sum_name)
        hprof_Sci_HG_sum = infile.Get(hprof_Sci_HG_sum_name)
        hprof_Cer_LG_sum = infile.Get(hprof_Cer_LG_sum_name)
        hprof_Sci_LG_sum = infile.Get(hprof_Sci_LG_sum_name)

        if not (hprof_Cer_HG_sum and hprof_Sci_HG_sum and hprof_Cer_LG_sum and hprof_Sci_LG_sum):
            print(
                f"Warning: Profiles {hprof_Cer_HG_sum_name}, {hprof_Sci_HG_sum_name}, {hprof_Cer_LG_sum_name}, or {hprof_Sci_LG_sum_name} not found in {infile_name}")
            continue

        hprofiles_Cer_HG_sum.append(hprof_Cer_HG_sum)
        hprofiles_Sci_HG_sum.append(hprof_Sci_HG_sum)
        hprofiles_Cer_LG_sum.append(hprof_Cer_LG_sum)
        hprofiles_Sci_LG_sum.append(hprof_Sci_LG_sum)

        legends.append(str(board_no))

    n_events = hprofiles_Cer_HG_sum[0].GetXaxis().GetXmax()
    legendPos = [0.3, 0.7, 0.9, 0.9]
    output_name = "FERS_Cer_HG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Cer_HG_sum, legends, 0, n_events, "Event", 0, 5e4, "Cer FERS Sum HG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(0, output_name + ".png")

    output_name = "FERS_Sci_HG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Sci_HG_sum, legends, 0, n_events, "Event", 0, 1.6e5, "Sci FERS Sum HG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(1, output_name + ".png")

    output_name = "FERS_Cer_LG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Cer_LG_sum, legends, 0, n_events, "Event", 0, 1.4e4, "Cer FERS Sum LG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(2, output_name + ".png")

    output_name = "FERS_Sci_LG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Sci_LG_sum, legends, 0, n_events, "Event", 0, 5e4, "Sci FERS Sum LG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=3, legendPos=legendPos)
    plots.insert(3, output_name + ".png")

    # total sum
    hprof_Cer_HG_sum = infile.Get(
        f"hprof_{fersboards.get_energy_sum_name(gain='HG', isCer=True)}_VS_Event")
    hprof_Sci_HG_sum = infile.Get(
        f"hprof_{fersboards.get_energy_sum_name(gain='HG', isCer=False)}_VS_Event")
    hprof_Cer_LG_sum = infile.Get(
        f"hprof_{fersboards.get_energy_sum_name(gain='LG', isCer=True)}_VS_Event")
    hprof_Sci_LG_sum = infile.Get(
        f"hprof_{fersboards.get_energy_sum_name(gain='LG', isCer=False)}_VS_Event")

    DrawHistos([hprof_Cer_HG_sum, hprof_Sci_HG_sum], ["Cer", "Sci"], 0, n_events, "Event", 0, 7e5, "FERS Total Sum HG (ADC)",
               "FERS_Total_HG_EnergySum_VS_Event",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=2, legendPos=[0.6, 0.85, 0.9, 0.9])
    plots.insert(0, "FERS_Total_HG_EnergySum_VS_Event.png")
    DrawHistos([hprof_Cer_LG_sum, hprof_Sci_LG_sum], ["Cer", "Sci"], 0, n_events, "Event", 0, 2e5, "FERS Total Sum LG (ADC)",
               "FERS_Total_LG_EnergySum_VS_Event",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=2, legendPos=[0.6, 0.85, 0.9, 0.9])
    plots.insert(1, "FERS_Total_LG_EnergySum_VS_Event.png")

    output_html = f"{paths['html']}/Conditions/FERS_energysum_vs_event.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html, plots_per_row=6)

    return output_html


def makeFERS1DPlots():
    plots = []

    infile_name = f"{paths['root']}/fers_all_channels_1d.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{paths['plots']}/FERS_1D"
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hist_C_name = f"hist_FERS_Board{board_no}_Cer_{sTowerX}_{sTowerY}"
            hist_S_name = f"hist_FERS_Board{board_no}_Sci_{sTowerX}_{sTowerY}"
            hist_C = infile.Get(hist_C_name)
            hist_S = infile.Get(hist_S_name)
            if not hist_C or not hist_S:
                print(
                    f"Warning: Histograms {hist_C_name} or {hist_S_name} not found in {infile_name}")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"Board: {fersboard.board_no}")
            extraToDraw.AddText(f"Tower X: {i_tower_x}")
            extraToDraw.AddText(f"Tower Y: {i_tower_y}")
            extraToDraw.AddText(
                f"Cer Channel: {fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=True).channel_no}")
            extraToDraw.AddText(
                f"Sci Channel: {fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=False).channel_no}")

            output_name = f"Energy_Board{board_no}_iTowerX{sTowerX}_iTowerY{sTowerY}"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       outdir=outdir_plots, run_number=run_number)

            plots.append(output_name + ".png")

    output_html = f"{paths['html']}/FERS/ChannelADC.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


def makeFERSStatsPlots(includePedestals=False):
    output_htmls = []
    outdir_plots = f"{paths['plots']}/FERS_Stats"
    # load the json file
    import json
    infile_name = f"{paths['root']}/fers_stats.json"
    with open(infile_name, "r") as f:
        stats = json.load(f)

    if includePedestals:
        infile_name_HG = f"{paths['root']}/fers_pedestals_hg.json"
        infile_name_LG = f"{paths['root']}/fers_pedestals_lg.json"
        with open(infile_name_HG, "r") as f:
            pedestals_HG = json.load(f)
        with open(infile_name_LG, "r") as f:
            pedestals_LG = json.load(f)

    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100
    valuemaps_HG_mean = {}
    valuemaps_HG_max = {}
    valuemaps_HG_satfreq = {}
    valuemaps_HG_pedestal = {}
    valuemaps_LG_mean = {}
    valuemaps_LG_max = {}
    valuemaps_LG_satfreq = {}
    valuemaps_LG_pedestal = {}

    for channelName, (vmean, vmax, VSatfreq) in stats.items():
        if "energyHG" in channelName:
            valuemaps_HG_mean[channelName] = vmean
            valuemaps_HG_max[channelName] = vmax
            valuemaps_HG_satfreq[channelName] = VSatfreq
            valuemaps_HG_pedestal[channelName] = pedestals_HG.get(
                channelName, None) if includePedestals else 0.
        elif "energyLG" in channelName:
            valuemaps_LG_mean[channelName] = vmean
            valuemaps_LG_max[channelName] = vmax
            valuemaps_LG_satfreq[channelName] = VSatfreq
            valuemaps_LG_pedestal[channelName] = pedestals_LG.get(
                channelName, None) if includePedestals else 0.

    [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
        fersboards, valuemaps_HG_mean, suffix=f"Run{run_number}_HG_mean", gain="HG")
    [h2_Cer_HG_max, h2_Cer_3mm_HG_max], [h2_Sci_HG_max, h2_Sci_3mm_HG_max] = visualizeFERSBoards(
        fersboards, valuemaps_HG_max, suffix=f"Run{run_number}_HG_max", gain="HG")
    [h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], [h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq] = visualizeFERSBoards(
        fersboards, valuemaps_HG_satfreq, suffix=f"Run{run_number}_HG_satfreq", gain="HG")
    [h2_Cer_HG_pedestal, h2_Cer_3mm_HG_pedestal], [h2_Sci_HG_pedestal, h2_Sci_3mm_HG_pedestal] = visualizeFERSBoards(
        fersboards, valuemaps_HG_pedestal, suffix=f"Run{run_number}_HG_pedestal", gain="HG")
    [h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], [h2_Sci_LG_mean, h2_Sci_3mm_LG_mean] = visualizeFERSBoards(
        fersboards, valuemaps_LG_mean, suffix=f"Run{run_number}_LG_mean", gain="LG")
    [h2_Cer_LG_max, h2_Cer_3mm_LG_max], [h2_Sci_LG_max, h2_Sci_3mm_LG_max] = visualizeFERSBoards(
        fersboards, valuemaps_LG_max, suffix=f"Run{run_number}_LG_max", gain="LG")
    [h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], [h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq] = visualizeFERSBoards(
        fersboards, valuemaps_LG_satfreq, suffix=f"Run{run_number}_LG_satfreq", gain="LG")
    [h2_Cer_LG_pedestal, h2_Cer_3mm_LG_pedestal], [h2_Sci_LG_pedestal, h2_Sci_3mm_LG_pedestal] = visualizeFERSBoards(
        fersboards, valuemaps_LG_pedestal, suffix=f"Run{run_number}_LG_pedestal", gain="LG")

    plots = []
    output_name = f"FERS_Boards_Run{run_number}_Stats_HG_mean"
    DrawHistos([h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_mean, h2_Sci_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{run_number}_Stats_LG_mean"
    DrawHistos([h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_mean, h2_Sci_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_html = f"{paths['html']}/FERS/Channel_Mean.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    output_htmls.append(output_html)

    plots = []
    output_name = f"FERS_Boards_Run{run_number}_Stats_HG_max"
    DrawHistos([h2_Cer_HG_max, h2_Cer_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_max, h2_Sci_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{run_number}_Stats_LG_max"
    DrawHistos([h2_Cer_LG_max, h2_Cer_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_max, h2_Sci_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_html = f"{paths['html']}/FERS/Channel_Max.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    output_htmls.append(output_html)

    plots = []
    output_name = f"FERS_Boards_Run{run_number}_Stats_HG_satfreq"
    DrawHistos([h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{run_number}_Stats_LG_satfreq"
    DrawHistos([h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_html = f"{paths['html']}/FERS/Channel_SatFreq.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    output_htmls.append(output_html)

    plots = []
    output_name = f"FERS_Boards_Run{run_number}_Stats_HG_pedestal"
    DrawHistos([h2_Cer_HG_pedestal, h2_Cer_3mm_HG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_pedestal, h2_Sci_3mm_HG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{run_number}_Stats_LG_pedestal"
    DrawHistos([h2_Cer_LG_pedestal, h2_Cer_3mm_LG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Cer", run_number=run_number, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_pedestal, h2_Sci_3mm_LG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extra_text="Sci", run_number=run_number, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Sci.png")

    output_html = f"{paths['html']}/FERS/Channel_Pedestal.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    output_htmls.append(output_html)
    return output_htmls


def makeFERSMaxValuePlots():
    plots = []
    xmin = 7500
    xmax = 8500
    outdir_plots = f"{paths['plots']}/FERS_MaxValues"
    infile_name = f"{paths['root']}/fers_max_values.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_board_cer_HG_max = []
    hists_board_cer_LG_max = []
    hists_board_sci_HG_max = []
    hists_board_sci_LG_max = []
    legends = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        hist_board_cer_HG_max = infile.Get(
            f'hist_{fersboard.get_energy_max_name(gain="HG", isCer=True)}')
        hist_board_sci_HG_max = infile.Get(
            f'hist_{fersboard.get_energy_max_name(gain="HG", isCer=False)}')
        hist_board_cer_LG_max = infile.Get(
            f'hist_{fersboard.get_energy_max_name(gain="LG", isCer=True)}')
        hist_board_sci_LG_max = infile.Get(
            f'hist_{fersboard.get_energy_max_name(gain="LG", isCer=False)}')
        hists_board_cer_HG_max.append(hist_board_cer_HG_max)
        hists_board_sci_HG_max.append(hist_board_sci_HG_max)
        hists_board_cer_LG_max.append(hist_board_cer_LG_max)
        hists_board_sci_LG_max.append(hist_board_sci_LG_max)
        legends.append(str(board_no))

    hist_cer_HG_max = infile.Get(
        f'hist_{fersboards.get_energy_max_name(gain="HG", isCer=True)}')
    hist_sci_HG_max = infile.Get(
        f'hist_{fersboards.get_energy_max_name(gain="HG", isCer=False)}')
    hist_cer_LG_max = infile.Get(
        f'hist_{fersboards.get_energy_max_name(gain="LG", isCer=True)}')
    hist_sci_LG_max = infile.Get(
        f'hist_{fersboards.get_energy_max_name(gain="LG", isCer=False)}')

    output_name = "FERS_Boards_CerEnergyHG_max"
    DrawHistos(hists_board_cer_HG_max, legends, xmin, xmax, f"HG Cer Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_SciEnergyHG_max"
    DrawHistos(hists_board_sci_HG_max, legends, xmin, xmax, f"HG Sci Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_CerEnergyLG_max"
    DrawHistos(hists_board_cer_LG_max, legends, xmin, xmax, f"LG Cer Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_SciEnergyLG_max"
    DrawHistos(hists_board_sci_LG_max, legends, xmin, xmax, f"LG Sci Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")

    output_name = "FERS_EnergyHG_max"
    VSat = 8000
    frac_sci = hist_sci_HG_max.Integral(hist_sci_HG_max.FindBin(VSat), 100000) / \
        hist_sci_HG_max.Integral(0, 100000)
    frac_cer = hist_cer_HG_max.Integral(hist_cer_HG_max.FindBin(VSat), 100000) / \
        hist_cer_HG_max.Integral(0, 100000)
    extraToDraw = ROOT.TPaveText(0.20, 0.63, 0.90, 0.72, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
    extraToDraw.AddText(f"Sat Frac Cer : {frac_cer:.3f}")
    DrawHistos([hist_cer_HG_max, hist_sci_HG_max],
               ["Cer",
                   "Sci"], xmin, xmax, f"HG Max (All Boards)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendPos=[0.30, 0.75, 0.50, 0.9], extraToDraw=extraToDraw)
    plots.append(output_name + ".png")
    output_name = "FERS_EnergyLG_max"
    frac_sci = hist_sci_LG_max.Integral(hist_sci_LG_max.FindBin(VSat), 100000) / \
        hist_sci_LG_max.Integral(0, 100000)
    frac_cer = hist_cer_LG_max.Integral(hist_cer_LG_max.FindBin(VSat), 100000) / \
        hist_cer_LG_max.Integral(0, 100000)
    extraToDraw = ROOT.TPaveText(0.20, 0.55, 0.90, 0.63, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
    extraToDraw.AddText(f"Sat Frac Cer : {frac_cer:.3f}")
    DrawHistos([hist_cer_LG_max, hist_sci_LG_max],
               ["Cer",
                   "Sci"], xmin, xmax, f"LG Max (All Boards)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, run_number=run_number, legendPos=[0.30, 0.75, 0.50, 0.9], extraToDraw=extraToDraw)
    plots.append(output_name + ".png")

    output_html = f"{paths['html']}/FERS/channelmax.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html

# 2D FERS histograms, hg VS lg


def makeFERS2DPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/FERS_2D"
    infile_name = f"{paths['root']}/fers_all_channels_2d.root"
    infile = ROOT.TFile(infile_name, "READ")
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{board_no}_{var}_{sTowerX}_{sTowerY}_hg_VS_lg"
                hist = infile.Get(hist_name)

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {fersboard.board_no}")
                extraToDraw.AddText(f"Tower X: {i_tower_x}")
                extraToDraw.AddText(f"Tower Y: {i_tower_y}")
                extraToDraw.AddText(f"{var} Channel: {chan.channel_no}")

                output_name = f"FERS_Board{board_no}_{var}_{sTowerX}_{sTowerY}_hg_VS_lg"
                DrawHistos([hist], f"", 0, 9000, "HG", 0, 1500, "LG",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True, extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, run_number=run_number)
                plots.append(output_name + ".png")
    output_html = f"{paths['html']}/FERS/LG_vs_HG.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# FERS output VS event
def trackFERSPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/FERS_VS_Event"
    infile_name = f"{paths['root']}/fers_all_channels_2D_VS_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{board_no}_{var}_VS_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {fersboard.board_no}")
                extraToDraw.AddText(f"Tower X: {i_tower_x}")
                extraToDraw.AddText(f"Tower Y: {i_tower_y}")
                extraToDraw.AddText(f"{var} Channel: {chan.channel_no}")

                n_events = hist.GetXaxis().GetXmax()

                output_name = f"FERS_Board{board_no}_{var}_{sTowerX}_{sTowerY}_VS_Event"
                DrawHistos([hist], "", 0, n_events, "Event", 1, 1e5, f"{var} Energy HG",
                           output_name,
                           dology=True, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, run_number=run_number)
                plots.append(output_name + ".png")
    output_html = f"{paths['html']}/FERS_VS_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# DRS VS TS
def makeDRSVSTSPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/DRS_VS_TS"
    infile_name = f"{paths['root']}/drs_vs_ts.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan is None:
                    print(
                        f"Warning: No channel found for Board {DRSBoard.board_no}, Tower ({i_tower_x}, {i_tower_y}), Var {var}")
                    continue
                channelName = chan.get_channel_name(blsub=True)
                hist_name = f"hist_{channelName}_VS_TS"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                output_name = f"DRS_{var}_VS_TS_{sTowerX}_{sTowerY}"
                plots.append(output_name + ".png")

                ymin_tmp, ymax_tmp = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm)

                extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"B: {DRSBoard.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                extraToDraw.AddText(f"i_tower_x: {i_tower_x}")
                extraToDraw.AddText(f"i_tower_y: {i_tower_y}")

                DrawHistos([hist], "", 0, 1024, "Time Slice", ymin_tmp, ymax_tmp, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extra_text=var, run_number=run_number, addOverflow=True)
    output_html = f"{paths['html']}/DRS/DRS_vs_TS.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/DRSPeakTS"
    infile_name = f"{paths['root']}/drspeakts.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_Cer = []
    hists_Sci = []
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hists = {}
            channelNos = {}
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist_name = f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    hists[var] = None
                    channelNos[var] = -1
                else:
                    hists[var] = hist
                    channelNos[var] = chan.channel_no

                output_name = f"hist_DRSPeakTS_{sTowerX}_{sTowerY}"

            if not hists["Cer"] or not hists["Sci"]:
                print(
                    f"Warning: Histograms for Cer or Sci not found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y})")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"B: {DRSBoard.board_no}, G: {chan.group_no}")
            extraToDraw.AddText(f"Tower X: {i_tower_x}")
            extraToDraw.AddText(f"Tower Y: {i_tower_y}")
            extraToDraw.AddText(
                f"Cer Channel: {channelNos['Cer']}")
            extraToDraw.AddText(
                f"Sci Channel: {channelNos['Sci']}")

            hists_Cer.append(hists["Cer"])
            hists_Sci.append(hists["Sci"])

            DrawHistos([hists["Cer"], hists["Sci"]], ["Cer", "Sci"], 400, 600, "Peak TS", 1, None, "Counts",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False, extraToDraw=extraToDraw,
                       outdir=outdir_plots, run_number=run_number)
            plots.append(output_name + ".png")

    # summary plots
    hist_Cer_Combined = LHistos2Hist(hists_Cer, "hist_DRSPeakTS_Cer_Combined")
    hist_Sci_Combined = LHistos2Hist(hists_Sci, "hist_DRSPeakTS_Sci_Combined")
    DrawHistos([hist_Cer_Combined, hist_Sci_Combined], ["Cer", "Sci"], 400, 600, "Peak TS", 1, None, "Counts",
               "DRS_PeakTS_Combined",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False,
               outdir=outdir_plots, run_number=run_number)
    plots.insert(0, "DRS_PeakTS_Combined.png")

    output_html = f"{paths['html']}/DRS/DRS_PeakTS.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSCerVSSciPlots():
    plots = []
    hists = []
    outdir_plots = f"{paths['plots']}/DRSPeakTSCerVSSci"
    infile_name = f"{paths['root']}/drspeakts.root"
    infile = ROOT.TFile(infile_name, "READ")

    # Create a dashed diagonal line from (0,0) to (1000,1000)
    diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
    diagonal_line.SetLineStyle(2)  # 2 = dashed
    diagonal_line.SetLineWidth(1)
    diagonal_line.SetLineColor(ROOT.kRed)
    extraToDraw = diagonal_line

    for _, DRSBoard in DRSBoards.items():
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hist_name = f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}"
            hist = infile.Get(hist_name)
            output_name = hist_name[5:]

            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue

            hists.append(hist)

            DrawHistos([hist], "", 400, 600, "Sci Peak TS", 400, 600, f"Cer Peak TS",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, run_number=run_number, extraToDraw=extraToDraw)

            plots.append(output_name + ".png")

    # summary plots
    hcombined = LHistos2Hist(hists, "hist_DRSPeakTS_Cer_VS_Sci_Combined")
    output_name = "DRS_PeakTS_Cer_VS_Sci_Combined"
    DrawHistos([hcombined], "", 400, 600, "Sci Peak TS", 400, 600, f"Cer Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
               outdir=outdir_plots, addOverflow=False, run_number=run_number, extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    output_html = f"{paths['html']}/DRS/DRS_PeakTS_Cer_VS_Sci.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# time reference
def compareTimeReferencePlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -2500
        ymax = 500
    plots = []
    infile_name = f"{paths['root']}/time_reference_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{paths['plots']}/TimeReference"
    for chan_name in time_reference_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"TimeReference_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, run_number=run_number)
        plots.append(output_name + ".png")

    output_html = f"{paths['html']}/TimeReference{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def compareServiceDRSPlots():
    ymin = 500
    ymax = 2500
    plots = []
    infile_name = f"{paths['root']}/service_drs_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{paths['plots']}/ServiceDRS"

    for det_name, chan_name in service_drs_channels.items():
        hist_name = f"hist_{chan_name}_blsub"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.60, 0.20, 0.90, 0.30, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{det_name}")
        output_name = f"ServiceDRS_{chan_name}"

        ymin, ymax = get_service_drs_plot_ranges(
            chan_name, subtractMedian=True)
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, run_number=run_number)
        plots.append(output_name + ".png")

    output_html = f"{paths['html']}/ServiceDRS/detectors.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)

    return output_html


def compareMCPPlots():
    ymin = -1500
    ymax = 500
    plots = []
    infile_name = f"{paths['root']}/mcp_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{paths['plots']}/MCP"

    for det_name, channels in mcp_channels.items():
        for chan_name in channels:
            hist_name = f"hist_{chan_name}_blsub"
            hist = infile.Get(hist_name)
            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue
            extraToDraw = ROOT.TPaveText(0.60, 0.20, 0.90, 0.30, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(f"{det_name}")
            output_name = f"MCP_{chan_name}"

            DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                       extraToDraw=extraToDraw,
                       outdir=outdir_plots, addOverflow=True, run_number=run_number)
            plots.append(output_name + ".png")

    output_html = f"{paths['html']}/ServiceDRS/MCP.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSSumVSFERSPlots():
    """
    Check if the sum of FERS and DRS energies are consistent.
    """
    plots = []
    outdir_plots = f"{paths['plots']}/DRSSum_VS_FERS"
    infile_name = f"{paths['root']}/drssum_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")

    xymax = {
        "Cer": (20000, 8500),
        "Sci": (30000, 8500)
    }
    xymax_LG = {
        "Cer": (20000, 2000),
        "Sci": (30000, 4000)
    }

    hists = []
    for _, DRSBoard in DRSBoards.items():
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                for gain in ["FERS", "FERSLG"]:
                    histname = f"hist_DRSSum_VS_{gain}_{var}_{sTowerX}_{sTowerY}"
                    hist = infile.Get(histname)
                    if not hist:
                        print(
                            f"Warning: Histogram {histname} not found in {infile_name}")
                        continue

                    zmax = hist.Integral(0, 10000, 0, 10000)
                    zmax = round_up_to_1eN(zmax)

                    tmp = xymax[var] if gain == "FERS" else xymax_LG[var]

                    hists.append(hist)

                    output_name = histname.replace("hist_", "")
                    DrawHistos([hist], "", 0, tmp[1], gain, 0, tmp[0], "DRSSum",
                               output_name,
                               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
                               outdir=outdir_plots, addOverflow=True, run_number=run_number, extra_text=f"{var}")
                    plots.append(output_name + ".png")

    # summary plots
    for var in ["Cer", "Sci"]:
        for gain in ["FERS", "FERSLG"]:
            histname = f"hist_DRSSum_VS_{gain}_{var}_Combined"
            hists_to_combine = [
                h for h in hists if f"_{gain}_{var}_" in h.GetName()]
            hist_combined = LHistos2Hist(hists_to_combine, histname)
            print("Combining histograms for", histname, " is ",
                  hist_combined, " len ", len(hists_to_combine))
            output_name = histname.replace("hist_", "")

            zmax = hist_combined.Integral(0, 10000, 0, 10000)
            zmax = round_up_to_1eN(zmax)

            output_name = histname.replace("hist_", "")

            tmp = xymax[var] if gain == "FERS" else xymax_LG[var]
            DrawHistos([hist_combined], "", 0, tmp[1], gain, 0, tmp[0], "DRSSum",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
                       outdir=outdir_plots, addOverflow=True, run_number=run_number, extra_text=f"{var}")
            plots.insert(0, output_name + ".png")

    output_html = f"{paths['html']}/DRS_VS_FERS/DRSSum_vs_FERS.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakVSFERSPlots():
    plots = []
    outdir_plots = f"{paths['plots']}/DRSPeak_VS_FERS"
    infile_name = f"{paths['root']}/drspeak_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if not chan:
                    print(
                        f"Warning: Channel not found for Board {board_no}, Tower ({i_tower_x}, {i_tower_y}), Var {var}")
                    continue
                _, ymax = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm)
                histname = f"hist_DRSPeak_VS_FERS_{var}_{sTowerX}_{sTowerY}"
                output_name = histname.replace("hist_", "")
                plots.append(output_name + ".png")

                hist = infile.Get(histname)
                if not hist:
                    print(
                        f"Warning: Histogram {histname} not found in {infile_name}")
                    continue

                DrawHistos([hist], "", 0, 9000, "FERS ADC", 0, ymax, "DRS Peak",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
                           outdir=outdir_plots, addOverflow=True, run_number=run_number, extra_text=f"{var}")

    output_html = f"{paths['html']}/DRS_VS_FERS/DRSPeak_vs_FERS.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def main():
    print("Starting DQM Plot Generation...\n")
    # Task Registry: (Label, Function)
    plot_tasks = [
        ("Conditions", makeConditionsPlots),
        ("FERS Sum", makeFERSSumPlots),
        ("FERS Mapping", lambda: DrawFERSBoards(run=run_number)),
        ("DRS Mapping", lambda: DrawDRSBoards(run=run_number)),
        ("FERS 1D", makeFERS1DPlots),
        ("FERS Stats", lambda: makeFERSStatsPlots(includePedestals=True)),
        ("DRS vs TS", makeDRSVSTSPlots),
        ("DRS Peak TS", makeDRSPeakTSPlots),
        ("DRS Peak TS Cer vs Sci", makeDRSPeakTSCerVSSciPlots),
        ("Service DRS", compareServiceDRSPlots),
        ("MCP", compareMCPPlots),
        # ("Time Reference", lambda: compareTimeReferencePlots(True)),
        ("FERS Max Values", makeFERSMaxValuePlots),
        ("DRS Sum vs FERS", makeDRSSumVSFERSPlots),
        ("DRS Peak vs FERS", makeDRSPeakVSFERSPlots),
    ]

    output_htmls = {}
    for label, func in plot_tasks:
        print(f"Generating {label} plots...")
        output_htmls[label] = func()

    print("\n" + "*"*30)
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
