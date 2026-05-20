"""
Plot-making functions for each CaloX sequence.

Each public function has the signature:
    plot_<name>(ctx: CaloXAnalysisManager) -> html_path | list[html_path] | None
"""

from configs.selection_config import get_service_drs_cut
from configs.plot_config import get_service_drs_processed_info_ranges
import ROOT
from channels.channel_map import (
    get_mcp_channels, get_pid_channels,
)
from channels.validate_map import DrawDRSBoards, DrawFERSBoards
from plotting.my_function import LHistos2Hist
from configs.plot_config import get_drs_plot_ranges, get_drs_cfd_finebins_range, get_drs_prof_plot_ranges
from utils.colors import colors
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle, STYLE_CER_SCI
from plotting.calox_plot_helper import BoardPlotHelper, create_pave_text, create_board_info_pave
from utils.utils import number_to_string, round_up_to_1eN, get_channel_var, get_hist_mpv
from variables.drs import get_ts_arr_name
from utils.visualization import visualizeFERSBoards, visualizeDRSBoards


# ---------------------------------------------------------------------------
# Shared styles (module-level, cheap to define)
# ---------------------------------------------------------------------------

_STYLE_BOARD_MULTI = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=colors,
    legendNCols=5,
    legendPos=[0.2, 0.15, 0.9, 0.35],
)
_STYLE_2D_LOG = PlotStyle(
    dology=False,
    dologz=True,
    drawoptions="COLZ",
    zmin=1,
    zmax=1e4,
)
_STYLE_1D = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 6, 4, 8],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90],
)
_STYLE_1D_LOG = PlotStyle(
    dology=True,
    drawoptions="HIST",
    mycolors=[2, 6, 4, 8],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90],
)


def _pm(ctx):
    """Convenience: create a PlotManager bound to ctx paths."""
    return PlotManager(
        ctx.paths["root"], ctx.paths["plots"], ctx.paths["html"], ctx.run_number,
        selection_text=getattr(ctx, "selection_summary", ""))


# ---------------------------------------------------------------------------
# FERS: conditions
# ---------------------------------------------------------------------------

def plot_monitor_conditions(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("Conditions_VS_Event")
        infile = pm._get_file("conditions_vs_event.root")

        profiles = {"SipmHV": [], "SipmI": [], "TempDET": [], "TempFPGA": []}
        legends = []

        for fersboard in ctx.fersboards.values():
            board_no = fersboard.board_no
            hp_hv = infile.Get(
                f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event")
            hp_i = infile.Get(f"hprof_{fersboard.get_sipm_i_name()}_VS_Event")
            hp_tdet = infile.Get(
                f"hprof_{fersboard.get_temp_det_name()}_VS_Event")
            hp_tfpga = infile.Get(
                f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event")

            if not all([hp_hv, hp_i, hp_tdet, hp_tfpga]):
                print(
                    f"Warning: Some condition profiles not found for board {board_no}")
                continue
            profiles["SipmHV"].append(hp_hv)
            profiles["SipmI"].append(hp_i)
            profiles["TempDET"].append(hp_tdet)
            profiles["TempFPGA"].append(hp_tfpga)
            legends.append(str(board_no))

        if not profiles["SipmHV"]:
            print("Warning: No condition profiles found")
            return None

        n_events = profiles["SipmHV"][0].GetXaxis().GetXmax()
        for var, ylabel, yrange in [
            ("SipmHV",  "Voltage (V)",       (26, 30)),
            ("SipmI",   "Current (mA)",       (0.0, 0.25)),
            ("TempDET", "Temperature (C)",    (14, 40)),
            ("TempFPGA", "Temperature (C)",    (32, 50)),
        ]:
            pm.plot_1d(
                profiles[var],
                f"Conditions_{var}_VS_Event",
                "Event", (0, n_events),
                ylabel=ylabel, yrange=yrange,
                legends=legends, style=_STYLE_BOARD_MULTI, prepend=True)

        return pm.generate_html("Conditions/conditions_vs_event.html", plots_per_row=4)


# ---------------------------------------------------------------------------
# FERS: energy sum
# ---------------------------------------------------------------------------

def plot_fers_energy_sum(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_EnergySum_VS_Event")
        infile = pm._get_file("fers_energysum_vs_event.root")

        profiles = {"Cer_HG": [], "Sci_HG": [], "Cer_LG": [], "Sci_LG": []}
        legends = []

        for fersboard in ctx.fersboards.values():
            for cat in ["Cer", "Sci"]:
                for gain in ["HG", "LG"]:
                    hname = (f"hprof_{fersboard.get_energy_sum_name(gain=gain, isCer=(cat == 'Cer'))}"
                             "_VS_Event")
                    hp = infile.Get(hname)
                    if hp:
                        profiles[f"{cat}_{gain}"].append(hp)
            legends.append(str(fersboard.board_no))

        if not profiles["Cer_HG"]:
            print("Warning: No FERS sum profiles found")
            return None

        n_events = profiles["Cer_HG"][0].GetXaxis().GetXmax()

        for cat, gain, ylabel in [
            ("Cer", "HG", "Cer FERS Sum HG (ADC)"),
            ("Sci", "HG", "Sci FERS Sum HG (ADC)"),
            ("Cer", "LG", "Cer FERS Sum LG (ADC)"),
            ("Sci", "LG", "Sci FERS Sum LG (ADC)"),
        ]:
            pm.plot_1d(
                profiles[f"{cat}_{gain}"],
                f"FERS_{cat}_{gain}_EnergySum_VS_Event",
                "Event", (0, n_events),
                ylabel=ylabel, yrange=(0, None),
                legends=legends, style=_STYLE_BOARD_MULTI, prepend=True)

        for gain in ["HG", "LG"]:
            hp_cer = infile.Get(
                f"hprof_{ctx.fersboards.get_energy_sum_name(gain=gain, isCer=True)}_VS_Event")
            hp_sci = infile.Get(
                f"hprof_{ctx.fersboards.get_energy_sum_name(gain=gain, isCer=False)}_VS_Event")
            if hp_cer and hp_sci:
                pm.plot_1d(
                    [hp_cer, hp_sci],
                    f"FERS_Total_{gain}_EnergySum_VS_Event",
                    "Event", (0, n_events),
                    ylabel=f"FERS Total Sum {gain} (ADC)", yrange=(0, None),
                    legends=["Cer", "Sci"], style=STYLE_CER_SCI,
                    legendPos=[0.6, 0.85, 0.9, 0.9], legendNCols=2, prepend=True)

        return pm.generate_html("Conditions/FERS_energysum_vs_event.html", plots_per_row=6)


# ---------------------------------------------------------------------------
# FERS: channel mapping
# ---------------------------------------------------------------------------

def plot_fers_mapping(ctx):
    DrawFERSBoards(run=ctx.run_number)
    DrawDRSBoards(run=ctx.run_number)


# ---------------------------------------------------------------------------
# FERS: per-channel 1D distributions
# ---------------------------------------------------------------------------

def plot_fers_channels(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_1D")
        infile = pm._get_file("fers_all_channels_1d.root")

        for fersboard in ctx.fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                hist_c = infile.Get(
                    f"hist_FERS_Board{board_no}_Cer_{s_x}_{s_y}")
                hist_s = infile.Get(
                    f"hist_FERS_Board{board_no}_Sci_{s_x}_{s_y}")
                if not hist_c or not hist_s:
                    print(
                        f"Warning: Hists not found Board{board_no} Tower({i_tower_x},{i_tower_y})")
                    continue
                pave = create_board_info_pave(
                    board_no, i_tower_x, i_tower_y,
                    channel_info={
                        "Cer": fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=True).channel_no,
                        "Sci": fersboard.get_channel_by_tower(i_tower_x, i_tower_y, isCer=False).channel_no,
                    })
                pm.plot_1d(
                    [hist_c, hist_s],
                    f"Energy_Board{board_no}_iTowerX{s_x}_iTowerY{s_y}",
                    "Energy HG", (0, 1000),
                    ylabel="Counts", yrange=(1, 1e5),
                    legends=["Cer", "Sci"],
                    style=PlotStyle(
                        dology=True, drawoptions="HIST", mycolors=[2, 4]),
                    extraToDraw=pave)

        return pm.generate_html("FERS/ChannelADC.html")


# ---------------------------------------------------------------------------
# FERS: board-level statistics (mean, max, saturation, pedestal)
# ---------------------------------------------------------------------------

def plot_fers_stats(ctx):
    import json

    output_htmls = []
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_Stats")

        with open(f"{ctx.paths['root']}/fers_stats.json") as f:
            stats = json.load(f)

        pedestals_hg, pedestals_lg = {}, {}
        if ctx.include_pedestals:
            try:
                with open(f"{ctx.paths['root']}/fers_pedestals_hg.json") as f:
                    pedestals_hg = json.load(f)
                with open(f"{ctx.paths['root']}/fers_pedestals_lg.json") as f:
                    pedestals_lg = json.load(f)
            except FileNotFoundError:
                print("Warning: Pedestal JSON files not found; skipping pedestals")

        valuemaps = {f"{gain}_{stat}": {}
                     for gain in ["HG", "LG"]
                     for stat in ["mean", "max", "satfreq", "pedestal"]}

        for channel_name, (vmean, vmax, v_satfreq) in stats.items():
            gain = "HG" if "energyHG" in channel_name else "LG"
            valuemaps[f"{gain}_mean"][channel_name] = vmean
            valuemaps[f"{gain}_max"][channel_name] = vmax
            valuemaps[f"{gain}_satfreq"][channel_name] = v_satfreq
            peds = pedestals_hg if gain == "HG" else pedestals_lg
            valuemaps[f"{gain}_pedestal"][channel_name] = peds.get(
                channel_name, 0.0)

        board_hists = {}
        for gain in ["HG", "LG"]:
            for stat in ["mean", "max", "satfreq", "pedestal"]:
                key = f"{gain}_{stat}"
                cer_hists, sci_hists = visualizeFERSBoards(
                    ctx.fersboards, valuemaps[key],
                    suffix=f"Run{ctx.run_number}_{key}", gain=gain)
                board_hists[f"{key}_Cer"] = cer_hists
                board_hists[f"{key}_Sci"] = sci_hists

        helper = BoardPlotHelper(pm)
        plot_configs = [
            ("mean",    "FERS/Stat/Channel_Mean.html",    0,    8000, 0),
            ("max",     "FERS/Stat/Channel_Max.html",     0,    8000, 0),
            ("satfreq", "FERS/Stat/Channel_SatFreq.html", 0,    1,    2),
            ("pedestal", "FERS/Stat/Channel_Pedestal.html", 100, 300,  0),
        ]
        for stat, html_path, zmin, zmax, digits in plot_configs:
            pm.reset_plots()
            for gain in ["HG", "LG"]:
                helper.plot_cer_sci_pair(
                    board_hists[f"{gain}_{stat}_Cer"],
                    board_hists[f"{gain}_{stat}_Sci"],
                    f"FERS_Boards_Run{ctx.run_number}_Stats_{gain}_{stat}",
                    zmin=zmin, zmax=zmax, nTextDigits=digits)
            output_htmls.append(pm.generate_html(html_path, plots_per_row=2))

    return output_htmls


# ---------------------------------------------------------------------------
# FERS: max-value distributions
# ---------------------------------------------------------------------------

def plot_fers_max(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_MaxValues")
        infile = pm._get_file("fers_max_values.root")
        xmin, xmax = 7500, 8500
        board_hists = {"Cer_HG": [], "Sci_HG": [], "Cer_LG": [], "Sci_LG": []}
        legends = []

        for fersboard in ctx.fersboards.values():
            for cat, isCer in [("Cer", True), ("Sci", False)]:
                for gain in ["HG", "LG"]:
                    hist = infile.Get(
                        f'hist_{fersboard.get_energy_max_name(gain=gain, isCer=isCer)}')
                    if hist:
                        board_hists[f"{cat}_{gain}"].append(hist)
            legends.append(str(fersboard.board_no))

        board_style = PlotStyle(
            dology=True, drawoptions="HIST", mycolors=colors,
            legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])

        for cat, isCer in [("Cer", True), ("Sci", False)]:
            for gain in ["HG", "LG"]:
                pm.plot_1d(
                    board_hists[f"{cat}_{gain}"],
                    f"FERS_Boards_{cat}Energy{gain}_max",
                    f"{gain} {cat} Max (Board)", (xmin, xmax),
                    yrange=(1, None), legends=legends, style=board_style)

        v_sat = 8000
        for gain in ["HG", "LG"]:
            hist_cer = infile.Get(
                f'hist_{ctx.fersboards.get_energy_max_name(gain=gain, isCer=True)}')
            hist_sci = infile.Get(
                f'hist_{ctx.fersboards.get_energy_max_name(gain=gain, isCer=False)}')
            if hist_cer and hist_sci:
                frac_sci = hist_sci.Integral(hist_sci.FindBin(v_sat), 100000) \
                    / (hist_sci.Integral(0, 100000) + 1e-6)
                frac_cer = hist_cer.Integral(hist_cer.FindBin(v_sat), 100000) \
                    / (hist_cer.Integral(0, 100000) + 1e-6)
                pave = create_pave_text(0.20, 0.63, 0.90, 0.72)
                pave.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
                pave.AddText(f"Sat Frac Cer : {frac_cer:.3f}")
                pm.plot_1d(
                    [hist_cer, hist_sci],
                    f"FERS_Energy{gain}_max",
                    f"{gain} Max (All Boards)", (xmin, xmax),
                    yrange=(1, None), legends=["Cer", "Sci"],
                    style=PlotStyle(
                        dology=True, drawoptions="HIST", mycolors=[2, 4]),
                    legendPos=[0.30, 0.75, 0.50, 0.9], extraToDraw=pave)

        return pm.generate_html("FERS/Stat/Channel_Max_1D.html", plots_per_row=2)


# ---------------------------------------------------------------------------
# FERS: 2D correlations
# ---------------------------------------------------------------------------

def plot_fers_2d(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_2D")
        infile = pm._get_file("fers_all_channels_2d.root")
        for fersboard in ctx.fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist = infile.Get(
                        f"hist_FERS_Board{board_no}_{var}_{s_x}_{s_y}_hg_VS_lg")
                    if not hist:
                        continue
                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={var: chan.channel_no},
                        position=(0.20, 0.70, 0.60, 0.90))
                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_x}_{s_y}_hg_VS_lg",
                        "HG", (0, 9000), "LG", (0, 1500),
                        style=_STYLE_2D_LOG, extraToDraw=pave)
        return pm.generate_html("FERS/LG_vs_HG.html", plots_per_row=4)


# ---------------------------------------------------------------------------
# FERS: vs-event 2D tracking
# ---------------------------------------------------------------------------

def plot_fers_track(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_VS_Event")
        infile = pm._get_file("fers_all_channels_2D_VS_event.root")
        for fersboard in ctx.fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hname = f"hist_FERS_Board{board_no}_{var}_VS_Event_{s_x}_{s_y}"
                    hist = infile.Get(hname)
                    if not hist:
                        print(f"Warning: {hname} not found")
                        continue
                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={var: chan.channel_no},
                        position=(0.20, 0.70, 0.60, 0.90))
                    n_events = hist.GetXaxis().GetXmax()
                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_x}_{s_y}_VS_Event",
                        "Event", (0, n_events),
                        f"{var} Energy HG", (1, 1e5),
                        style=PlotStyle(dology=True, dologz=True,
                                        drawoptions="COLZ", zmin=1, zmax=1e4),
                        extraToDraw=pave)
        return pm.generate_html("FERS_VS_Event/index.html", plots_per_row=4)


# ---------------------------------------------------------------------------
# DRS: waveforms vs time-slice
# ---------------------------------------------------------------------------

def _plot_drs_channel_vs_ts(pm, infile, channel_name, mode, ymin, ymax, pave,
                            plot_name, extra_text=None):
    suffix = f"_{mode}" if mode != "raw" else ""
    ch_blsub = f"{channel_name}_blsub"
    hist = infile.Get(f"hist_{ch_blsub}_VS_ts{suffix}")
    if not hist:
        print(f"Warning: hist_{ch_blsub}_VS_ts{suffix} not found")
        return
    ts_col = "ts" if mode == "raw" else get_ts_arr_name(
        channel_name, use_mcp=(mode == "mcp"))
    pm.plot_2d(
        hist, plot_name, ts_col, (0, 1024), "DRS Output", (ymin, ymax),
        style=_STYLE_2D_LOG, extraToDraw=pave, extra_text=extra_text)


def plot_drs_waveforms(ctx):
    """Plot DRS waveform 2D hists for raw / ref / mcp TS correction modes
    plus service-DRS profiles. Returns list of HTML paths."""
    output_htmls = []
    map_mcp_channels = get_mcp_channels(ctx.run_number)
    list_mcp_channels = list(map_mcp_channels.values())

    #for mode in ["raw", "ref", "mcp"]:
    for mode in ["mcp"]:
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_VS_TS")
            infile = pm._get_file("drs_vs_ts.root")
            for _, board in ctx.drsboards.items():
                for chan in board:
                    ch = chan.get_channel_name(blsub=False)
                    ch_blsub = chan.get_channel_name(blsub=True)
                    ymin, ymax = get_drs_plot_ranges(
                        subtractMedian=True,
                        is_amplified=chan.is_amplified,
                        is6mm=chan.is6mm,
                        is_reference=chan.is_reference)
                    pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                    pave.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    if not chan.is_reference:
                        pave.AddText(
                            f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                    var = get_channel_var(chan)
                    _plot_drs_channel_vs_ts(
                        pm, infile, ch, mode, ymin, ymax, pave,
                        f"DRS_ADC_VS_ts_{mode}_{ch_blsub}_{var}",
                        extra_text=var)
            output_htmls.append(
                pm.generate_html(f"DRS/DRS_vs_ts_{mode}.html", plots_per_row=9))

    # Service DRS (MCP channels) raw mode
    with _pm(ctx) as pm:
        pm.set_output_dir("Service_DRS_VS_TS")
        infile = pm._get_file("drs_vs_ts.root")
        for channel_name in list_mcp_channels:
            ymin, ymax = get_drs_plot_ranges(subtractMedian=True, is_mcp=True)
            pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
            pave.AddText(channel_name)
            _plot_drs_channel_vs_ts(
                pm, infile, channel_name, "raw", ymin, ymax, pave,
                f"Service_DRS_ADC_VS_ts_raw_{channel_name}_blsub",
                extra_text="Service")
        output_htmls.append(
            pm.generate_html("DRS/Service_DRS_vs_ts_raw.html", plots_per_row=4))

    # Profile plots
    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_Prof_VS_TS")
        infile = pm._get_file("drs_vs_ts.root")
        for _, board in ctx.drsboards.items():
            for chan in board:
                ch_blsub = chan.get_channel_name(blsub=True)
                hist = infile.Get(f"prof_{ch_blsub}_VS_ts")
                hist_ref = infile.Get(f"prof_{ch_blsub}_VS_ts_ref")
                hist_mcp = infile.Get(f"prof_{ch_blsub}_VS_ts_mcp")
                if not hist:
                    print(f"Warning: prof_{ch_blsub}_VS_ts not found")
                    continue
                for h in [hist, hist_ref, hist_mcp]:
                    if h:
                        h = h.ProjectionX()
                ymin_tmp, ymax_tmp = get_drs_prof_plot_ranges(
                    subtractMedian=True,
                    is_amplified=chan.is_amplified,
                    is6mm=chan.is6mm,
                    is_reference=chan.is_reference,
                    is_cer=chan.isCer)
                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                var = get_channel_var(chan)
                pm.plot_1d(
                    [h for h in [hist, hist_ref, hist_mcp] if h],
                    f"DRS_ADC_Prof_VS_ts_{ch_blsub}_{var}",
                    "TS", (0, 1024),
                    "Mean DRS Output", (ymin_tmp, ymax_tmp),
                    legends=["raw ts", "ts_ref", "ts_mcp"],
                    legendPos=[0.55, 0.70, 0.90, 0.90],
                    style=_STYLE_1D, extraToDraw=pave, extra_text=var)
        output_htmls.append(
            pm.generate_html("DRS/DRS_Prof_vs_TS.html", plots_per_row=9))

    return output_htmls


# ---------------------------------------------------------------------------
# DRS: pulse statistics (peak, energy, timing)
# ---------------------------------------------------------------------------

def plot_drs_stats(ctx):
    """Plot DRS peak, energy, and timing distributions. Returns list of HTML paths."""
    output_htmls = []

    # -- peak values --
    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_Stats")
        infile = pm._get_file("drs_stats.root")
        for _, board in ctx.drsboards.items():
            for chan in board:
                ch = chan.get_channel_name(blsub=False)
                hist = infile.Get(f"hist_{ch}_peak_value")
                if not hist:
                    print(f"Warning: hist_{ch}_peak_value not found")
                    continue
                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                var = get_channel_var(chan)
                pm.plot_1d(
                    hist, f"DRS_PeakValue_{ch}_{var}",
                    "DRS peak value", (0, 800), "Counts", (0.9, None),
                    style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)
        output_htmls.append(pm.generate_html(
            "DRS/DRS_Peak.html", plots_per_row=9))

    # -- energy (CFD integral) --
    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_Stats")
        infile = pm._get_file("drs_stats.root")
        for _, board in ctx.drsboards.items():
            for chan in board:
                ch = chan.get_channel_name(blsub=False)
                hist = infile.Get(f"hist_{ch}_energy")
                if not hist:
                    print(f"Warning: hist_{ch}_energy not found")
                    continue
                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                if not chan.is_reference:
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                var = get_channel_var(chan)
                pm.plot_1d(
                    hist, f"DRS_Energy_{ch}_{var}",
                    "CFD Energy", (0, 6000), "Counts", (0.9, None),
                    style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)
        output_htmls.append(pm.generate_html(
            "DRS/DRS_Sum.html", plots_per_row=9))

    # -- timing --
    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_Stats")
        infile = pm._get_file("drs_stats.root")
        for _, board in ctx.drsboards.items():
            for chan in board:
                if chan.is_reference:
                    continue
                ch = chan.get_channel_name(blsub=False)
                hp_ref = infile.Get(f"hist_{ch}_TS_peak_ref")
                hc_ref = infile.Get(f"hist_{ch}_TS_cfd_ref")
                hp_mcp = infile.Get(f"hist_{ch}_TS_peak_mcp")
                hc_mcp = infile.Get(f"hist_{ch}_TS_cfd_mcp")
                hc_fine = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
                if not hp_ref or not hc_ref:
                    print(f"Warning: TS histograms for {ch} not found")
                    continue
                var = get_channel_var(chan)
                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(
                    f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                pave.AddText(f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                pm.plot_1d(
                    [hp_ref, hc_ref, hp_mcp, hc_mcp],
                    f"DRS_Time_{ch}_{var}",
                    "Pulse TS", (0, 1024), "Counts", (0.9, None),
                    legends=["TS_peak_ref", "TS_cfd_ref",
                             "TS_peak_mcp", "TS_cfd_mcp"],
                    legendPos=[0.55, 0.70, 0.90, 0.90],
                    style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)

        output_htmls.append(pm.generate_html(
            "DRS/DRS_Time.html", plots_per_row=9))

        # fine-binned CFD MCP
        pm.reset_plots()
        for _, board in ctx.drsboards.items():
            for chan in board:
                if chan.is_reference:
                    continue
                ch = chan.get_channel_name(blsub=False)
                hist = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
                if not hist:
                    continue
                var = get_channel_var(chan)
                pave = create_pave_text(0.15, 0.68, 0.70, 0.90)
                pave.AddText(
                    f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                pave.AddText(f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                if hist.GetEntries() >= 5:
                    mpv, mpv_err = get_hist_mpv(hist)
                    pave.AddText(f"MPV: {mpv:.1f} #pm {mpv_err:.1f} TS")
                pm.plot_1d(
                    hist, f"DRS_Time_FineBins_{ch}_{var}",
                    "CFD TS (MCP-corrected)", get_drs_cfd_finebins_range(chan.isCer),
                    "Counts", (0.9, None),
                    style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)
        output_htmls.append(pm.generate_html(
            "DRS/DRS_Time_FineBins.html", plots_per_row=9))

    return output_htmls


# ---------------------------------------------------------------------------
# DRS: CFD MPV board map (one value per channel)
# ---------------------------------------------------------------------------

def plot_drs_cfd_mpv(ctx):
    """Read the MPV of the fine-binned TS_cfd_mcp histogram for every channel
    and display the result as a 2D DRS board map."""
    infile = ctx.hbook.open_file("drs_stats.root")

    mpv_map = {}
    for _, board in ctx.drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            ch = chan.get_channel_name(blsub=False)
            hist = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
            if not hist:
                continue
            if hist.GetEntries() < 5:
                mpv_map[ch] = 0.0
            else:
                mpv, _ = get_hist_mpv(hist)
                mpv_map[ch] = mpv - 400

    cer_hists, sci_hists = visualizeDRSBoards(
        ctx.drsboards, valuemaps=mpv_map,
        suffix=f"MPV_cfd_Run{ctx.run_number}")

    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_CFD_MPV")
        helper = BoardPlotHelper(pm)
        helper.plot_cer_sci_pair(
            cer_hists, sci_hists,
            f"DRS_CFD_MPV_Run{ctx.run_number}",
            zmin=0, zmax=100, nTextDigits=1)
        return pm.generate_html(
            "DRS/DRS_CFD_MPV.html", plots_per_row=2,
            title="DRS CFD MPV Board Map",
            intro_text=(
                "Per-channel most probable value (MPV) of the fine-binned "
                "MCP-corrected CFD time slice, displayed on the DRS board map. "
                "Values are shifted by -400 TS for readability. "
                "Channels with fewer than 5 entries are shown as 0."
            ))



# ---------------------------------------------------------------------------
# DRS: peak time-slice
# ---------------------------------------------------------------------------

def plot_drs_peak_ts(ctx):
    output_htmls = []
    with _pm(ctx) as pm:
        pm.set_output_dir("DRSPeakTS")
        infile = pm._get_file("drspeakts.root")
        hists_cer, hists_sci = [], []

        for _, board in ctx.drsboards.items():
            board_no = board.board_no
            for i_tower_x, i_tower_y in board.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                hists, channel_nos = {}, {}
                for var in ["Cer", "Sci"]:
                    chan = board.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist = infile.Get(f"hist_DRSPeakTS_{var}_{s_x}_{s_y}")
                    hists[var] = hist
                    channel_nos[var] = chan.channel_no if chan else -1
                if not hists.get("Cer") or not hists.get("Sci"):
                    continue
                hists_cer.append(hists["Cer"])
                hists_sci.append(hists["Sci"])
                if ctx.do_detailed_plots:
                    pave = create_board_info_pave(
                        board_no, i_tower_x, i_tower_y,
                        channel_info={"Cer": channel_nos["Cer"],
                                      "Sci": channel_nos["Sci"]})
                    pm.plot_1d(
                        [hists["Cer"], hists["Sci"]],
                        f"hist_DRSPeakTS_{s_x}_{s_y}",
                        "Peak TS", (400, 600),
                        yrange=(1, None), legends=["Cer", "Sci"],
                        style=STYLE_CER_SCI, extraToDraw=pave)

        if hists_cer and hists_sci:
            pm.plot_1d(
                [LHistos2Hist(hists_cer, "hist_DRSPeakTS_Cer_Combined"),
                 LHistos2Hist(hists_sci, "hist_DRSPeakTS_Sci_Combined")],
                "DRS_PeakTS_Combined",
                "Peak TS", (400, 600),
                yrange=(1, None), legends=["Cer", "Sci"],
                style=STYLE_CER_SCI, prepend=True)
        output_htmls.append(pm.generate_html(
            "DRS/DRS_PeakTS.html", plots_per_row=4))

    # Cer vs Sci 2D
    with _pm(ctx) as pm:
        pm.set_output_dir("DRSPeakTSCerVSSci")
        infile = pm._get_file("drspeakts.root")
        diagonal = ROOT.TLine(0, 0, 1000, 1000)
        diagonal.SetLineStyle(2)
        diagonal.SetLineWidth(1)
        diagonal.SetLineColor(ROOT.kRed)
        hists = []
        for _, board in ctx.drsboards.items():
            for i_tower_x, i_tower_y in board.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                hist = infile.Get(f"hist_DRSPeakTS_Cer_VS_Sci_{s_x}_{s_y}")
                if not hist:
                    continue
                hists.append(hist)
                if ctx.do_detailed_plots:
                    pm.plot_2d(
                        hist,
                        f"DRSPeakTS_Cer_VS_Sci_{s_x}_{s_y}",
                        "Sci Peak TS", (400, 600), "Cer Peak TS", (400, 600),
                        style=PlotStyle(dologz=True, drawoptions="COLZ",
                                        zmin=1, zmax=1e2, addOverflow=False),
                        extraToDraw=diagonal)
        if hists:
            pm.plot_2d(
                LHistos2Hist(hists, "hist_DRSPeakTS_Cer_VS_Sci_Combined"),
                "DRS_PeakTS_Cer_VS_Sci_Combined",
                "Sci Peak TS", (400, 600), "Cer Peak TS", (400, 600),
                style=PlotStyle(dologz=True, drawoptions="COLZ",
                                zmin=1, zmax=1e2, addOverflow=False),
                extraToDraw=diagonal, prepend=True)
        output_htmls.append(
            pm.generate_html("DRS/DRS_PeakTS_Cer_VS_Sci.html", plots_per_row=4))

    return output_htmls


# ---------------------------------------------------------------------------
# DRS ↔ FERS correlations
# ---------------------------------------------------------------------------

def plot_drs_sum_vs_fers(ctx):
    xymax = {"Cer": (20000, 8500), "Sci": (30000, 8500)}
    xymax_LG = {"Cer": (20000, 2000), "Sci": (30000, 4000)}

    with _pm(ctx) as pm:
        pm.set_output_dir("DRSSum_VS_FERS")
        infile = pm._get_file("drssum_vs_fers.root")
        hists = []
        for _, board in ctx.drsboards.items():
            for i_tower_x, i_tower_y in board.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    for gain_tag in ["FERS", "FERSLG"]:
                        hname = f"hist_DRSSum_VS_{gain_tag}_{var}_{s_x}_{s_y}"
                        hist = infile.Get(hname)
                        if not hist:
                            print(f"Warning: {hname} not found")
                            continue
                        hists.append(hist)
                        if ctx.do_detailed_plots:
                            zmax = round_up_to_1eN(
                                hist.Integral(0, 10000, 0, 10000))
                            tmp = xymax[var] if gain_tag == "FERS" else xymax_LG[var]
                            pm.plot_2d(
                                hist, hname.replace("hist_", ""),
                                gain_tag, (0, tmp[1]), "DRSSum", (0, tmp[0]),
                                style=PlotStyle(dologz=True, drawoptions="COLZ",
                                                zmin=1, zmax=zmax),
                                extra_text=var)

        for var in ["Cer", "Sci"]:
            for gain_tag in ["FERS", "FERSLG"]:
                subset = [
                    h for h in hists if f"_{gain_tag}_{var}_" in h.GetName()]
                if subset:
                    combined = LHistos2Hist(
                        subset, f"hist_DRSSum_VS_{gain_tag}_{var}_Combined")
                    zmax = round_up_to_1eN(
                        combined.Integral(0, 10000, 0, 10000))
                    tmp = xymax[var] if gain_tag == "FERS" else xymax_LG[var]
                    pm.plot_2d(
                        combined, f"DRSSum_VS_{gain_tag}_{var}_Combined",
                        gain_tag, (0, tmp[1]), "DRSSum", (0, tmp[0]),
                        style=PlotStyle(dologz=True, drawoptions="COLZ",
                                        zmin=1, zmax=zmax),
                        extra_text=var, prepend=True)

        return pm.generate_html("DRS_VS_FERS/DRSSum_vs_FERS.html", plots_per_row=4)


def plot_drs_peak_vs_fers(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("DRSPeak_VS_FERS")
        infile = pm._get_file("drspeak_vs_fers.root")
        for _, board in ctx.drsboards.items():
            board_no = board.board_no
            for i_tower_x, i_tower_y in board.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    chan = board.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if not chan:
                        continue
                    _, ymax = get_drs_plot_ranges(
                        subtractMedian=True,
                        is_amplified=chan.is_amplified,
                        is6mm=chan.is6mm)
                    hname = f"hist_DRSPeak_VS_FERS_{var}_{s_x}_{s_y}"
                    hist = infile.Get(hname)
                    if not hist:
                        print(f"Warning: {hname} not found")
                        continue
                    pm.plot_2d(
                        hist, hname.replace("hist_", ""),
                        "FERS ADC", (0, 9000), "DRS Peak", (0, ymax),
                        style=PlotStyle(dologz=True, drawoptions="COLZ",
                                        zmin=1, zmax=None),
                        extra_text=var)
        return pm.generate_html("DRS_VS_FERS/DRSPeak_vs_FERS.html", plots_per_row=4)


# ===========================================================================
# Service DRS sequences  (ctx: CaloXAnalysisManager)
# ===========================================================================


_STYLE_SVC_1D = PlotStyle(dology=False, drawoptions="HIST", mycolors=[1])
_STYLE_SVC_1D_LOG = PlotStyle(dology=True,  drawoptions="HIST", mycolors=[1])
_STYLE_SVC_2D_LOG = PlotStyle(dology=False, dologz=True, drawoptions="COLz",
                              zmin=1, zmax=None)


def _plot_pulse_distributions(channels, infile, pm, suffix):
    """Per-detector pulse-shape distributions (1D + profiles)."""
    for det in channels.keys():
        wf_ymin, wf_ymax = get_service_drs_processed_info_ranges(
            det, "waveform")
        pv_xmin, pv_xmax = get_service_drs_processed_info_ranges(
            det, "peak_value")
        en_xmin, en_xmax = get_service_drs_processed_info_ranges(det, "sum")
        _, _, _, _, value_cut, cut_method = get_service_drs_cut(det)

        hist = infile.Get(f"{det}_ADC_vs_TS")
        if hist:
            pm.plot_2d(hist, f"{det}_ADC_vs_TS", "Time Slice", (0, 1024),
                       "ADC Counts", (wf_ymin, wf_ymax), style=_STYLE_SVC_2D_LOG)
        prof = infile.Get(f"{det}_ADC_vs_TS_prof")
        if prof:
            prof_px = prof.ProjectionX()
            pm.plot_1d(prof_px, f"{det}_ADC_vs_TS_prof", "Time Slice", (0, 1024),
                       "Mean ADC Counts", (wf_ymin / 10.0, wf_ymax / 5.0),
                       style=_STYLE_SVC_1D)

        hist = infile.Get(f"{det}_cfd_ts")
        if hist:
            pm.plot_1d(hist, f"{det}_cfd_ts", "CFD TS", (0, 1024),
                       style=_STYLE_SVC_1D)

        hist = infile.Get(f"{det}_peak_value_vs_cfd_ts")
        if hist:
            pm.plot_2d(hist, f"{det}_peak_value_vs_cfd_ts",
                       "CFD TS", (0, 1024), "Peak Value", (pv_xmin, pv_xmax),
                       style=_STYLE_SVC_2D_LOG)

        hist_pv = infile.Get(f"{det}_peak_value")
        if hist_pv:
            if cut_method == "PeakValue":
                nhad = hist_pv.Integral(hist_pv.FindBin(value_cut), 10000)
                nele = hist_pv.Integral(0, hist_pv.FindBin(value_cut))
                ntot = nhad + nele + 1e-6
                pave = create_pave_text(0.23, 0.75, 0.55, 0.85)
                pave.SetFillColor(0)
                pave.AddText(
                    f"N (peak > {value_cut:.2g}): {nhad:.0f} ({nhad/ntot:.1%})")
                pave.AddText(
                    f"N (peak < {value_cut:.2g}): {nele:.0f} ({nele/ntot:.1%})")
                line = ROOT.TLine(value_cut, 0, value_cut,
                                  hist_pv.GetMaximum())
                line.SetLineColor(ROOT.kRed)
                line.SetLineWidth(2)
                line.SetLineStyle(ROOT.kDashed)
                pm.plot_1d(hist_pv, f"{det}_peak_value", "Peak Value",
                           (pv_xmin, pv_xmax), yrange=(1, None),
                           style=_STYLE_SVC_1D_LOG, extraToDraw=[pave, line])

                hist_pv.GetXaxis().SetRange(0, hist_pv.GetNbinsX() + 1)
                hist_cdf = hist_pv.GetCumulative()
                hist_cdf.Scale(1.0 / hist_pv.Integral())
                line_cdf = ROOT.TLine(value_cut, 0, value_cut, 1)
                line_cdf.SetLineColor(ROOT.kRed)
                line_cdf.SetLineWidth(2)
                line_cdf.SetLineStyle(ROOT.kDashed)
                pm.plot_1d(hist_cdf, f"{det}_peak_value_cdf", "Peak Value",
                           (pv_xmin, pv_xmax), ylabel="Cumulative Fraction",
                           yrange=(0, 1.3), style=_STYLE_SVC_1D,
                           addOverflow=False, addUnderflow=False,
                           legendPos=[0.3, 0.80, 0.5, 0.85],
                           extraToDraw=line_cdf)
            else:
                pm.plot_1d(hist_pv, f"{det}_peak_value", "Peak Value",
                           (pv_xmin, pv_xmax), style=_STYLE_SVC_1D, leftlegend=True)

        hist_en = infile.Get(f"{det}_energy")
        if hist_en:
            if cut_method != "PeakValue":
                nhad = hist_en.Integral(hist_en.FindBin(value_cut), 10000)
                nele = hist_en.Integral(0, hist_en.FindBin(value_cut))
                ntot = nhad + nele + 1e-6
                pave = create_pave_text(0.23, 0.75, 0.55, 0.85)
                pave.SetFillColor(0)
                pave.AddText(
                    f"N (energy > {value_cut:.2g}): {nhad:.0f} ({nhad/ntot:.1%})")
                pave.AddText(
                    f"N (energy < {value_cut:.2g}): {nele:.0f} ({nele/ntot:.1%})")
                line = ROOT.TLine(value_cut, 0, value_cut,
                                  hist_en.GetMaximum())
                line.SetLineColor(ROOT.kRed)
                line.SetLineWidth(2)
                line.SetLineStyle(ROOT.kDashed)
                pm.plot_1d(hist_en, f"{det}_energy", "Energy (ADC)",
                           (en_xmin, en_xmax), yrange=(1, None),
                           style=_STYLE_SVC_1D_LOG, extraToDraw=[pave, line])

                hist_en.GetXaxis().SetRange(0, hist_en.GetNbinsX() + 1)
                hist_cdf = hist_en.GetCumulative()
                hist_cdf.Scale(1.0 / hist_en.Integral())
                line_cdf = ROOT.TLine(value_cut, 0, value_cut, 1)
                line_cdf.SetLineColor(ROOT.kRed)
                line_cdf.SetLineWidth(2)
                line_cdf.SetLineStyle(ROOT.kDashed)
                pm.plot_1d(hist_cdf, f"{det}_energy_cdf", "Energy (ADC)",
                           (en_xmin, en_xmax), ylabel="Cumulative Fraction",
                           yrange=(0, 1.3), style=_STYLE_SVC_1D,
                           addOverflow=False, addUnderflow=False,
                           legendPos=[0.3, 0.80, 0.5, 0.85],
                           extraToDraw=line_cdf)
            else:
                pm.plot_1d(hist_en, f"{det}_energy", "Energy (ADC)",
                           (en_xmin, en_xmax), yrange=(1, None),
                           style=_STYLE_SVC_1D_LOG)

        hist = infile.Get(f"{det}_integral_to_peak")
        if hist:
            pm.plot_1d(hist, f"{det}_integral_to_peak",
                       "Integral/Peak", (-5, 20), style=_STYLE_SVC_1D)

    intro_text = ("Pulse shape analysis for service DRS channels. "
                  "No selection applied unless specified.")
    return pm.generate_html(
        f"ServiceDRS/{suffix}.html", plots_per_row=8,
        intro_text=intro_text, title=f"{suffix.upper()} Analysis")


def _plot_pulse_correlations(channels, infile, pm, suffix):
    """2D correlation plots between detector pairs."""
    det_list = list(channels.keys())
    trigger_dets = ["KT1", "KT2", "T3", "T4"]
    det_list_notrigger = [d for d in det_list if d not in trigger_dets]

    output_htmls = []
    corr_categories = [("PID", det_list_notrigger), ("Trigger", trigger_dets)]

    for cat, tmp_list in corr_categories:
        pm.reset_plots()
        for idx1, det1 in enumerate(tmp_list):
            for idx2, det2 in enumerate(tmp_list):
                if idx2 <= idx1:
                    continue
                _, _, _, _, value_cut1, method1 = get_service_drs_cut(det1)
                _, _, _, _, value_cut2, method2 = get_service_drs_cut(det2)
                var1 = "peak_value" if method1 == "PeakValue" else "energy"
                var2 = "peak_value" if method2 == "PeakValue" else "energy"
                xmin, xmax = get_service_drs_processed_info_ranges(
                    det1, "peak_value" if method1 == "PeakValue" else "sum")
                ymin, ymax = get_service_drs_processed_info_ranges(
                    det2, "peak_value" if method2 == "PeakValue" else "sum")

                hist2d = infile.Get(f"{det1}_{var1}_vs_{det2}_{var2}")
                if not hist2d:
                    continue

                xPass = hist2d.GetXaxis().FindBin(value_cut1)
                yPass = hist2d.GetYaxis().FindBin(value_cut2)
                nPP = hist2d.Integral(xPass, 1000, yPass, 1000)
                nPF = hist2d.Integral(xPass, 1000, 0, yPass - 1)
                nFP = hist2d.Integral(0, xPass - 1, yPass, 1000)
                nFF = hist2d.Integral(0, xPass - 1, 0, yPass - 1)

                pave = create_pave_text(0.23, 0.20, 0.5, 0.40)
                pave.SetTextSize(0.03)
                name1 = det1.replace("Cerenkov", "Cer")
                name2 = det2.replace("Cerenkov", "Cer")
                pave.AddText(
                    f"N ({name1} > {value_cut1:.2g}, {name2} > {value_cut2:.2g}): {nPP:.0f}")
                pave.AddText(
                    f"N ({name1} > {value_cut1:.2g}, {name2} < {value_cut2:.2g}): {nPF:.0f}")
                pave.AddText(
                    f"N ({name1} < {value_cut1:.2g}, {name2} > {value_cut2:.2g}): {nFP:.0f}")
                pave.AddText(
                    f"N ({name1} < {value_cut1:.2g}, {name2} < {value_cut2:.2g}): {nFF:.0f}")

                line1 = ROOT.TLine(value_cut1, ymin, value_cut1, ymax)
                line2 = ROOT.TLine(xmin, value_cut2, xmax, value_cut2)
                for ln in [line1, line2]:
                    ln.SetLineWidth(2)
                    ln.SetLineStyle(ROOT.kDashed)
                    ln.SetLineColor(ROOT.kRed)

                pm.plot_2d(hist2d, f"{det1}_vs_{det2}_corr2D",
                           f"{det1} {var1}", (xmin, xmax),
                           f"{det2} {var2}", (ymin, ymax),
                           style=_STYLE_SVC_2D_LOG,
                           extraToDraw=[pave, line1, line2])

            pm.add_newline()

        intro_text = (f"Correlation plots of service DRS channels ({cat}). "
                      "No selection applied unless specified.")
        output_htmls.append(pm.generate_html(
            f"ServiceDRS/{cat}_correlation_{suffix}.html",
            plots_per_row=5, intro_text=intro_text))

    return output_htmls


def _plot_pulse(ctx, channels, suffix, include_correlations=True):
    """Combined plot runner: distributions + (optionally) correlations."""
    infile_name = f"{ctx.paths['root']}/drs_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open {infile_name}")

    with _pm(ctx) as pm:
        pm.set_output_dir(f"drs_{suffix}")
        output_htmls = [_plot_pulse_distributions(
            channels, infile, pm, suffix)]
        if include_correlations:
            output_htmls += _plot_pulse_correlations(
                channels, infile, pm, suffix)

    infile.Close()
    return output_htmls


def _plot_mcp_timing_diff(ctx, channels_mcp):
    """Gaussian fits to MCP CFD timing differences for all unique pairs."""
    infile_name = f"{ctx.paths['root']}/drs_mcp_timing_diff.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open {infile_name}")

    dets = list(channels_mcp.keys())

    with _pm(ctx) as pm:
        pm.set_output_dir("drs_mcp_timing_diff")
        for i, det1 in enumerate(dets):
            for det2 in dets[i + 1:]:
                hist = infile.Get(f"{det1}_cfd_diff_vs_{det2}")
                if not hist:
                    continue
                fit = ROOT.TF1("gaus_constrained", "gaus", -10, 10)
                mean_guess = hist.GetMean()
                sigma_guess = min(hist.GetRMS(), 1.9)
                fit.SetParameter(0, hist.GetMaximum())
                fit.SetParameter(1, mean_guess)
                fit.SetParameter(2, sigma_guess)
                fit.SetParLimits(1, mean_guess - 1, mean_guess + 1)
                fit.SetParLimits(2, 0.0, 2.0)
                hist.Fit(fit, "QB")
                fit.SetLineColor(ROOT.kRed)
                fit.SetLineWidth(2)
                sigma = abs(fit.GetParameter(2))
                n_gauss = (fit.Integral(hist.GetXaxis().GetXmin(),
                                        hist.GetXaxis().GetXmax())
                           / hist.GetBinWidth(1))
                n_total = hist.Integral(0, hist.GetNbinsX() + 1)
                pave = create_pave_text(0.55, 0.68, 0.90, 0.88)
                pave.AddText(f"Mean: {fit.GetParameter(1):.2f} TS")
                pave.AddText(f"Sigma: {sigma:.2f} TS")
                pave.AddText(f"N (Gaussian): {n_gauss:.0f}")
                pave.AddText(f"N (total): {n_total:.0f}")
                fitted_mean = fit.GetParameter(1)
                ymax = hist.GetMaximum() * 1.4
                pm.plot_1d(
                    hist, f"{det1}_cfd_diff_vs_{det2}",
                    f"#Delta t_{{CFD,ref}} ({det1} - {det2}) [TS]",
                    (fitted_mean - 2, fitted_mean + 2),
                    yrange=(0.9, ymax), ylabel="Counts", style=_STYLE_SVC_1D,
                    extraToDraw=[fit, pave])
            pm.add_newline()

        intro_text = (
            "Event-by-event MCP CFD time differences for all unique pairs. "
            "Reflects DRS+Method timing resolution.")
        output_html = pm.generate_html(
            "ServiceDRS/mcp_timing_diff.html",
            plots_per_row=len(channels_mcp) - 1,
            title="MCP CFD Timing Differences",
            intro_text=intro_text)

    infile.Close()
    return output_html


def _plot_hodo_peak(ctx):
    """Hodoscope peak position and timing distributions."""
    from channels.channel_map import build_hodo_pos_channels
    infile = ROOT.TFile(f"{ctx.paths['root']}/hodoscope_peaks.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError("Failed to open hodoscope_peaks.root")

    hodo_pos_channels = build_hodo_pos_channels(run=ctx.run_number)
    style_compare = PlotStyle(
        dology=False, drawoptions="HIST", mycolors=[1, 1])

    with _pm(ctx) as pm:
        pm.set_output_dir("DWC")
        for group, channels in hodo_pos_channels.items():
            labels = ["No Selection", "Pass HoleVeto"]
            linestyles = [1, 2]
            histos = {key: [] for key in [
                "diff", "diff_rel", "sum", "left", "right",
                "left_vs_right", "left_peak", "right_peak"
            ]}
            for sel in ["passNone", "is_HoleVeto_vetoed"]:
                cat = f"{group}_{sel}"
                histos["diff"].append(infile.Get(f"{cat}_delta_peak"))
                histos["diff_rel"].append(
                    infile.Get(f"{cat}_delta_peak_relative"))
                histos["sum"].append(infile.Get(f"{cat}_sum_peak"))
                histos["left"].append(infile.Get(f"{cat}_left_peak"))
                histos["right"].append(infile.Get(f"{cat}_right_peak"))
                histos["left_vs_right"].append(
                    infile.Get(f"{cat}_left_peak_vs_right_peak"))
                histos["left_peak"].append(
                    infile.Get(f"{cat}_left_peak_value"))
                histos["right_peak"].append(
                    infile.Get(f"{cat}_right_peak_value"))

            ls_style = PlotStyle(dology=False, drawoptions="HIST",
                                 mycolors=[1, 1], linestyles=linestyles)
            pm.plot_1d(histos["diff"], f"{group}_diff_peak",
                       "Peak Position Difference", (-250, 250),
                       legends=labels, style=ls_style)
            pm.plot_1d(histos["diff_rel"], f"{group}_diff_peak_relative",
                       "Peak Position Difference Relative", (-0.5, 0.5),
                       legends=labels, style=ls_style)
            pm.plot_1d(histos["sum"], f"{group}_sum_peak",
                       "Peak Position Sum", (600, 1200),
                       legends=labels, style=ls_style)
            pm.plot_1d(
                histos["left"] + histos["right"],
                f"{group}_peaks", "Peak Position", (300, 800),
                legends=["Left No Sel", "Left Pass Veto",
                         "Right No Sel", "Right Pass Veto"],
                style=PlotStyle(dology=False, drawoptions="HIST",
                                mycolors=[1, 1, 2, 2], linestyles=[1, 2, 1, 2]))
            for idx, cat_label in enumerate(["pass_NoSel", "pass_upstream_veto"]):
                if histos["left_vs_right"][idx]:
                    pm.plot_2d(histos["left_vs_right"][idx],
                               f"{group}_left_vs_right_peak_{cat_label}",
                               "Left Peak Position", (0, 1024),
                               "Right Peak Position", (0, 1024),
                               style=_STYLE_2D_LOG)
            pm.plot_1d(
                histos["left_peak"] + histos["right_peak"],
                f"{group}_peak_values", "Peak Value", (-1500, 100),
                legends=["Left No Sel", "Left Pass Veto",
                         "Right No Sel", "Right Pass Veto"],
                style=PlotStyle(dology=False, drawoptions="HIST",
                                mycolors=[1, 1, 2, 2], linestyles=[1, 2, 1, 2]))

        for sel in ["pass_NoSel", "pass_upstream_veto"]:
            for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
                cat = f"{dwc}_{sel}"
                hist2d = infile.Get(cat)
                if hist2d:
                    pm.plot_2d(hist2d, f"DWC_{dwc}_{sel}",
                               "LR Delta Peak Relative", (-0.4, 0.4),
                               "UD Delta Peak Relative", (-0.4, 0.4),
                               style=_STYLE_2D_LOG, prepend=True)
        pm.add_newline()
        output_html = pm.generate_html(
            "ServiceDRS/DWC.html", plots_per_row=7, title="DWC Positions")

    infile.Close()
    return output_html


# ---------------------------------------------------------------------------
# Public plot functions (ctx: CaloXAnalysisManager)
# ---------------------------------------------------------------------------

def plot_service_drs_pid(ctx):
    return _plot_pulse(ctx, get_pid_channels(ctx.run_number), suffix="services")


def plot_service_drs_mcp(ctx):
    return _plot_pulse(ctx, get_mcp_channels(ctx.run_number), suffix="mcp",
                       include_correlations=False)


def plot_service_drs_mcp_timing(ctx):
    return _plot_mcp_timing_diff(ctx, get_mcp_channels(ctx.run_number))


def plot_service_drs_hodo(ctx):
    return _plot_hodo_peak(ctx)
