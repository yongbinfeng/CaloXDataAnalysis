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
from configs.plot_config import (get_drs_plot_ranges, get_drs_cfd_finebins_range,
                                 get_drs_energy_range, get_drs_peak_value_range,
                                 get_drs_sum_vs_fers_ranges, get_fers_vs_drs_xmax,
                                 get_drs_prof_plot_ranges, get_drs_time_arr_ns_range,
                                 get_drs_time_ns_finebins_range,
                                 get_fers_1d_range, get_fers_max_range,
                                 get_fers_2d_hg_lg_range, get_fers_2d_mix_lg_range,
                                 getRangesForFERSEnergySums, get_ttu_hodo_ranges)
from utils.colors import colors
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle, STYLE_CER, STYLE_SCI, STYLE_CER_SCI
from plotting.calox_plot_helper import BoardPlotHelper, create_pave_text, create_board_info_pave
from utils.utils import number_to_string, round_up_to_1eN, get_channel_var, get_hist_mpv
from variables.drs import get_arr_name
from utils.visualization import (visualizeFERSBoards, visualizeDRSBoards,
                                 FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY,
                                 FERS_YMIN_DISPLAY, FERS_YMAX_DISPLAY,
                                 FERS_W_REF, FERS_H_REF)


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
    legendoptions="L",
)
_STYLE_1D_LOG = PlotStyle(
    dology=True,
    drawoptions="HIST",
    mycolors=[2, 6, 4, 8],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90],
    legendoptions="L",
)
_STYLE_1D_SMOOTH = PlotStyle(
    dology=False,
    drawoptions="C",
    mycolors=[2, 6, 4, 8],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90],
    legendoptions="L",
)


def _pm(ctx):
    """Convenience: create a PlotManager bound to ctx paths."""
    return PlotManager(
        ctx.paths["root"], ctx.paths["plots"], ctx.paths["html"], ctx.run_number,
        selection_text=getattr(ctx, "selection_summary", ""),
        use_jsroot=getattr(ctx.args, 'jsroot', False),
    )


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

def plot_fers_esum_vs_event(ctx):
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
    return DrawFERSBoards(run_number=ctx.run_number)


def plot_drs_mapping(ctx):
    return DrawDRSBoards(run_number=ctx.run_number)


# ---------------------------------------------------------------------------
# FERS: per-channel 1D distributions
# ---------------------------------------------------------------------------

def plot_fers_channels(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_1D")
        infile = pm._get_file("fers_all_channels_1d.root")

        _xrange = {g: r[1:] for g, r in get_fers_1d_range().items()}
        for gain in ["HG", "LG"]:
            for fersboard in ctx.fersboards.values():
                board_no = fersboard.board_no
                for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                    s_x = number_to_string(i_tower_x)
                    s_y = number_to_string(i_tower_y)
                    hist_c = infile.Get(
                        f"hist_FERS_Board{board_no}_Cer_{s_x}_{s_y}_{gain}")
                    hist_s = infile.Get(
                        f"hist_FERS_Board{board_no}_Sci_{s_x}_{s_y}_{gain}")
                    if not hist_c or not hist_s:
                        print(
                            f"Warning: Hists not found Board{board_no} Tower({i_tower_x},{i_tower_y}) {gain}")
                        continue
                    cer_ch = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=True).channel_no
                    sci_ch = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=False).channel_no
                    pave = create_pave_text(0.20, 0.65, 0.60, 0.90)
                    pave.AddText(
                        f"Board: {board_no}, Tower: ({i_tower_x}, {i_tower_y})")
                    pave.AddText(f"Cer Ch: {cer_ch}")
                    pave.AddText(f"Sci Ch: {sci_ch}")
                    pm.plot_1d(
                        [hist_c, hist_s],
                        f"Energy_{gain}_Board{board_no}_iTowerX{s_x}_iTowerY{s_y}",
                        f"Energy {gain} [ADC]", _xrange[gain],
                        ylabel="Counts", yrange=(1, 1e5),
                        legends=["Cer", "Sci"],
                        style=PlotStyle(
                            dology=True, drawoptions="HIST", mycolors=[2, 4]),
                        extraToDraw=pave)

        return pm.generate_html("FERS/ChannelADC.html")


# ---------------------------------------------------------------------------
# FERS: board-level statistics (mean, max, saturation, pedestal)
# ---------------------------------------------------------------------------

def plot_fers_stats(ctx, *, do_mean=True, do_max=True, do_satfreq=True, do_pedestal=True):
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

        helper = BoardPlotHelper(pm,
                                 xrange=(FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY),
                                 yrange=(FERS_YMIN_DISPLAY, FERS_YMAX_DISPLAY),
                                 W_ref=FERS_W_REF, H_ref=FERS_H_REF)
        plot_configs = [
            ("mean",     "FERS/Stat/Channel_Mean.html",
             0,   8000, 0, "Mean ADC",          do_mean),
            ("max",      "FERS/Stat/Channel_Max.html",
             0,   8000, 0, "Max ADC",            do_max),
            ("satfreq",  "FERS/Stat/Channel_SatFreq.html",
             0,   1,    2, "Saturation Rate",    do_satfreq),
            ("pedestal", "FERS/Stat/Channel_Pedestal.html",
             100, 300,  0, "Pedestal [ADC]",     do_pedestal),
        ]
        for stat, html_path, zmin, zmax, digits, zlabel, do_run in plot_configs:
            if not do_run:
                continue
            pm.reset_plots()
            for gain in ["HG", "LG"]:
                helper.plot_cer_sci_pair(
                    board_hists[f"{gain}_{stat}_Cer"],
                    board_hists[f"{gain}_{stat}_Sci"],
                    f"FERS_Boards_Run{ctx.run_number}_Stats_{gain}_{stat}",
                    zmin=zmin, zmax=zmax, nTextDigits=digits, zlabel=zlabel)
            output_htmls.append(pm.generate_html(html_path, plots_per_row=2))

    return output_htmls


# ---------------------------------------------------------------------------
# FERS: max-value distributions
# ---------------------------------------------------------------------------

def plot_fers_max(ctx):
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_MaxValues")
        infile = pm._get_file("fers_max_values.root")
        _, xmin, xmax = get_fers_max_range()
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
                    pave = create_pave_text(0.20, 0.70, 0.60, 0.90)
                    pave.AddText(f"Board: {board_no}, Ch: {chan.channel_no}")
                    pave.AddText(f"Tower: ({i_tower_x}, {i_tower_y})")
                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_x}_{s_y}_hg_VS_lg",
                        "HG [ADC]", get_fers_2d_hg_lg_range(
                        )[1:3], "LG [ADC]", get_fers_2d_hg_lg_range()[4:6],
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
# FERS: advanced physics plots (require define_physics_variables)
# ---------------------------------------------------------------------------

_FERS_GAIN_CALIBS = [("HG", False), ("LG", False), ("Mix", True)]


def plot_fers_energy_sum(ctx):
    HE = ctx.beam_energy >= 50
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_EnergySum")
        infile = pm._get_file("fers_energy_sum.root")
        for gain, calib in _FERS_GAIN_CALIBS:
            cfg = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE,
                run_number=ctx.run_number, beam_energy=ctx.beam_energy)
            gain_unit = cfg[f"title_{gain}"]
            gain_title = f"[{gain_unit}]" if gain == "Mix" else f"{gain} {gain_unit}"
            for cat in ["cer", "sci"]:
                is_cer = (cat == "cer")
                vname = ctx.fersboards.get_energy_sum_name(
                    gain=gain, isCer=is_cer, pdsub=True, calib=calib)
                h = infile.Get(f"hist_{vname}")
                if h:
                    style = STYLE_CER if is_cer else STYLE_SCI
                    pm.plot_1d(h, f"FERS_ESum_{gain}_{cat}",
                               f"{cat.capitalize()} {gain_title}",
                               (cfg["xmin_total"][f"{gain}_{cat}"],
                                cfg["xmax_total"][f"{gain}_{cat}"]),
                               style=style, prepend=(cat == "cer"))
        return pm.generate_html("FERS/ESum.html", plots_per_row=6,
                                title="FERS Energy Sums")


def plot_fers_energy_sum_board(ctx):
    """Plot per-board FERS energy sum 1D distributions (separate from the total)."""
    HE = ctx.beam_energy >= 50
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_EnergySum_Board")
        infile = pm._get_file("fers_energy_sum_board.root")
        for gain, calib in _FERS_GAIN_CALIBS:
            cfg = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE,
                run_number=ctx.run_number, beam_energy=ctx.beam_energy)
            gain_unit = cfg[f"title_{gain}"]
            gain_title = f"[{gain_unit}]" if gain == "Mix" else f"{gain} {gain_unit}"
            for cat in ["cer", "sci"]:
                is_cer = (cat == "cer")
                style = STYLE_CER if is_cer else STYLE_SCI
                for fb in ctx.fersboards.values():
                    vname = fb.get_energy_sum_name(
                        gain=gain, isCer=is_cer, pdsub=True, calib=calib)
                    h = infile.Get(f"hist_{vname}")
                    if h:
                        pm.plot_1d(h, f"FERS_ESum_Board{fb.board_no}_{gain}_{cat}",
                                   f"Board {fb.board_no} {cat.capitalize()} {gain_title}",
                                   (cfg["xmin_board"][f"{gain}_{cat}"],
                                    cfg["xmax_board"][f"{gain}_{cat}"]),
                                   style=style, prepend=(cat == "cer"))
        return pm.generate_html("FERS/ESum_Board.html", plots_per_row=6,
                                title="FERS Energy Sums (per board)")


def plot_fers_lg_vs_mix(ctx):
    """Plot FERS LG vs Mix 2D histograms per channel (requires Mix variables defined)."""
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_LG_vs_Mix")
        infile = pm._get_file("fers_lg_vs_mix_2d.root")
        for fersboard in ctx.fersboards.values():
            board_no = fersboard.board_no
            for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    hist = infile.Get(
                        f"hist_FERS_Board{board_no}_{var}_{s_x}_{s_y}_mix_VS_lg")
                    if not hist:
                        continue
                    pave = create_pave_text(0.20, 0.70, 0.60, 0.90)
                    pave.AddText(f"Board: {board_no}, Ch: {chan.channel_no}")
                    pave.AddText(f"Tower: ({i_tower_x}, {i_tower_y})")
                    pm.plot_2d(
                        hist,
                        f"FERS_Board{board_no}_{var}_{s_x}_{s_y}_mix_VS_lg",
                        "Mix [ADC]", get_fers_2d_mix_lg_range(
                        )[1:3], "LG [ADC]", get_fers_2d_mix_lg_range()[4:6],
                        style=_STYLE_2D_LOG, extraToDraw=pave)
        return pm.generate_html("FERS/LG_vs_Mix.html", plots_per_row=4)


def plot_fers_cer_vs_sci(ctx):
    HE = ctx.beam_energy >= 50
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_Cer_vs_Sci")
        infile = pm._get_file("fers_cer_vs_sci.root")
        for gain, calib in _FERS_GAIN_CALIBS:
            cfg = getRangesForFERSEnergySums(
                pdsub=True, calib=calib, clip=False, HE=HE,
                run_number=ctx.run_number, beam_energy=ctx.beam_energy)
            vc = ctx.fersboards.get_energy_sum_name(
                gain=gain, isCer=True,  pdsub=True, calib=calib)
            vs = ctx.fersboards.get_energy_sum_name(
                gain=gain, isCer=False, pdsub=True, calib=calib)
            hist = infile.Get(f"hist_{vc}_VS_{vs}")
            if not hist:
                continue
            extra = []
            gain_unit = cfg[f'title_{gain}']
            gain_title = f"[{gain_unit}]" if gain == "Mix" else f"{gain} {gain_unit}"
            if gain == "Mix":
                for fn, col in [("x", ROOT.kRed), ("0.5*x", ROOT.kBlue)]:
                    f = ROOT.TF1(f"f_{fn}_{gain}", fn, 0, 200)
                    f.SetLineColor(col)
                    f.SetLineStyle(2)
                    f.SetLineWidth(2)
                    extra.append(f)
            pm.plot_2d(hist, f"FERS_Total_Cer_VS_Sci_{gain}",
                       f"Sci {gain_title}",
                       (cfg["xmin_total"][f"{gain}_sci"],
                        cfg["xmax_total"][f"{gain}_sci"]),
                       f"Cer {gain_title}",
                       (cfg["xmin_total"][f"{gain}_cer"],
                        cfg["xmax_total"][f"{gain}_cer"]),
                       style=PlotStyle(dology=False, dologz=True,
                                       drawoptions="COLZ", zmin=1),
                       extraToDraw=extra or None, prepend=True, usePDF=True)
            if gain == "Mix":
                pm.add_newline()
                for cat, vx in [("Cer", vc), ("Sci", vs)]:
                    h = infile.Get(f"hist_{vx}_VS_PSD_Sum")
                    if h:
                        pm.plot_2d(h, f"FERS_Total_{cat}_VS_PSD_{gain}",
                                   "PSD Sum [ADC]", get_service_drs_processed_info_ranges(
                                       "PSD", "sum"),
                                   f"{cat} {gain_title}",
                                   (cfg[f"xmin_total"][f"{gain}_{cat.lower()}"],
                                    cfg[f"xmax_total"][f"{gain}_{cat.lower()}"]),
                                   style=PlotStyle(dology=False, dologz=True, drawoptions="COLZ", zmin=1))
                        prof = h.ProfileX(f"prof_{h.GetName()}", 1, -1)
                        pm.plot_1d(prof, f"FERS_Total_{cat}_VS_PSD_{gain}_Prof",
                                   "PSD Sum [ADC]", get_service_drs_processed_info_ranges(
                                       "PSD", "sum"),
                                   f"{cat} {gain_title}",
                                   style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[2 if cat == "Cer" else 4]))
        return pm.generate_html("FERS/Cer_VS_Sci.html", plots_per_row=3,
                                title="FERS Cer vs Sci correlations")


def plot_fers_dr(ctx):
    HE = ctx.beam_energy >= 50
    gain, calib = "Mix", True
    cfg = getRangesForFERSEnergySums(
        pdsub=True, calib=calib, clip=False, HE=HE,
        run_number=ctx.run_number, beam_energy=ctx.beam_energy)
    vc = ctx.fersboards.get_energy_sum_name(
        gain=gain, isCer=True,  pdsub=True, calib=calib)
    vs = ctx.fersboards.get_energy_sum_name(
        gain=gain, isCer=False, pdsub=True, calib=calib)
    xlo, xhi = cfg["xmin_total"][f"{gain}_sci"], cfg["xmax_total"][f"{gain}_sci"]

    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_DR")
        infile = pm._get_file("fers_dr.root")

        for var_sfx, xlabel, xrange in [("OuterRing", "Leakage Energy", (0, 20)),
                                        ("LeakCorr",  f"Leak-corrected [{cfg[f'title_{gain}']}]", (xlo, xhi))]:
            hl = [infile.Get(f"hist_{v}_{var_sfx}") for v in [vc, vs]]
            hl = [h for h in hl if h]
            if hl:
                pm.plot_1d(hl, f"FERS_DR_{var_sfx}", xlabel, xrange,
                           legends=["Cer", "Sci"], style=STYLE_CER_SCI)
        pm.add_newline()

        for vn, xlabel, xrange, style in [
            ("COverS", "C/S", (0, 2), STYLE_CER),
            ("fEM",    "f_{EM}", (0, 1.5), STYLE_SCI),
        ]:
            h = infile.Get(f"hist_{vn}")
            if h:
                pm.plot_1d(h, f"FERS_DR_{vn}", xlabel, xrange, style=style)

        vcersc = vc.replace("Cer", "CerSci")
        h = infile.Get(f"hist_{vcersc}")
        if h:
            pm.plot_1d(h, "FERS_DR_CerSci", f"Cer+Sci [{cfg[f'title_{gain}']}]", (xlo, xhi * 2.5),
                       style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[6]))
        pm.add_newline()

        hists_combined, legends_combined, tf1s, pave = [], [
        ], [], create_pave_text(0.20, 0.65, 0.60, 0.90)
        for label, vn in [("Cer", vc), ("Sci", vs),
                          ("DR", vc.replace("Cer", "DR")),
                          ("DR m2", vc.replace("Cer", "DR") + "_method2"),
                          ("DR m3", vc.replace("Cer", "DR") + "_method3")]:
            h = infile.Get(f"hist_{vn}" if label in ("DR", "DR m2", "DR m3") else
                           f"hist_{vn}")
            if h and h.Integral() > 0:
                fr = h.Fit("gaus", "SQ")
                if fr and fr.IsValid():
                    ff = h.GetFunction("gaus")
                    pave.AddText(f"{label}: #mu={ff.GetParameter(1):.1f}, "
                                 f"#sigma={ff.GetParameter(2):.1f}")
                    tf1s.append(ff)
                hists_combined.append(h)
                legends_combined.append(label)
        if hists_combined:
            pm.plot_1d(hists_combined, "FERS_DR_Energy_Combined",
                       "Energy [GeV]", (xlo, xhi),
                       legends=legends_combined,
                       style=PlotStyle(dology=False, drawoptions="HIST",
                                       mycolors=[2, 4, 1, 7, 8]),
                       extraToDraw=[pave] + tf1s, prepend=True)
        pm.add_newline()

        for sfx, cat, vy in [("Sci", "sci", vs), ("Cer", "cer", vc)]:
            for method in ["", "_method2", "_method3"]:
                vdr = vc.replace("Cer", "DR") + method
                h = infile.Get(f"hist_{vy}_VS_{vdr}")
                if h:
                    pm.plot_2d(h, f"FERS_DR_{sfx}_VS_DR{method}",
                               f"DR{method}", (xlo, xhi),
                               f"[{cfg[f'title_{gain}']}]", (cfg[f"xmin_total"][f"{gain}_{cat}"],
                                                             cfg[f"xmax_total"][f"{gain}_{cat}"]),
                               style=PlotStyle(
                                   dology=False, dologz=True, drawoptions="COLZ", zmin=1),
                               extra_text=sfx)

        return pm.generate_html("FERS/DR.html", plots_per_row=4, title="FERS Dual Readout")


def plot_fers_ewc(ctx):
    HE = ctx.beam_energy >= 50
    gain, calib = "Mix", True
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_EWC")
        infile = pm._get_file("fers_ewc.root")
        for cat in ["cer", "sci"]:
            is_cer = (cat == "cer")
            style = STYLE_CER if is_cer else STYLE_SCI
            vx = ctx.fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=is_cer, pdsub=True, calib=calib, isX=True)
            vy = ctx.fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=is_cer, pdsub=True, calib=calib, isX=False)
            for axis, vn in [("X", vx), ("Y", vy)]:
                h = infile.Get(f"hist_{vn}")
                if h:
                    pave = create_pave_text(0.20, 0.85, 0.60, 0.90)
                    pave.AddText(
                        f"Mean = {h.GetMean():.2f}, RMS = {h.GetRMS():.2f}")
                    pm.plot_1d(h, f"FERS_{cat.capitalize()}_EWC_{axis}",
                               f"EWC {axis} [cm]", (-15, 15),
                               style=style, extraToDraw=pave, extra_text=cat.capitalize())
            h2 = infile.Get(f"hist_{vy}_VS_{vx}")
            if h2:
                pm.plot_2d(h2, f"FERS_{cat.capitalize()}_EWC_Y_vs_X",
                           "EWC X [cm]", (-15, 15), "EWC Y [cm]", (-15, 15),
                           style=PlotStyle(drawoptions=["colz"], addOverflow=False,
                                           addUnderflow=False, zmin=1),
                           extra_text=cat.capitalize())
            hp = infile.Get(f"hprof_{vy}_VS_{vx}_WithEnergy")
            if hp:
                zmin = 0.7 * ctx.beam_energy if cat == "sci" else 0.5 * ctx.beam_energy
                zmax = 1.2 * ctx.beam_energy if cat == "sci" else 1.1 * ctx.beam_energy
                pm.plot_2d(hp, f"FERS_{cat.capitalize()}_EWC_Y_vs_X_WithEnergy",
                           "EWC X [cm]", (-15, 15), "EWC Y [cm]", (-15, 15),
                           style=PlotStyle(drawoptions=["colz"], addOverflow=False,
                                           addUnderflow=False, zmin=zmin, zmax=zmax,
                                           zlabel="Avg Energy [GeV]"),
                           extra_text=cat.capitalize())
        return pm.generate_html("FERS/EWC.html", plots_per_row=4,
                                title="FERS Energy Weighted Centre")


def plot_fers_ewc_vs_hodo(ctx):
    gain, calib = "Mix", True
    hodo_min, hodo_max, _ = get_ttu_hodo_ranges()
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_EWC_vs_Hodo")
        infile = pm._get_file("fers_ewc_vs_hodo.root")
        for cat in ["cer", "sci"]:
            is_cer = (cat == "cer")
            vx = ctx.fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=is_cer, pdsub=True, calib=calib, isX=True)
            vy = ctx.fersboards.get_energy_weighted_center_name(
                gain=gain, isCer=is_cer, pdsub=True, calib=calib, isX=False)
            zmin = 0.7 * ctx.beam_energy if cat == "sci" else 0.5 * ctx.beam_energy
            zmax = 1.2 * ctx.beam_energy if cat == "sci" else 1.1 * ctx.beam_energy
            for axis, ewc_v, hodo_ax in [("X", vx, "HodoX"), ("Y", vy, "HodoY")]:
                h = infile.Get(f"hist_{ewc_v}_VS_{hodo_ax}")
                if h:
                    pm.plot_2d(h, f"FERS_{cat.capitalize()}_EWC_{axis}_vs_{hodo_ax}",
                               f"Hodo {axis} [cm]", (hodo_min, hodo_max),
                               f"EWC {axis} [cm]", (-15, 15),
                               style=PlotStyle(drawoptions=["colz"], addOverflow=False,
                                               addUnderflow=False, zmin=1),
                               extra_text=cat.capitalize())
                hp = infile.Get(f"hprof_{ewc_v}_VS_{hodo_ax}_WithEnergy")
                if hp:
                    pm.plot_2d(hp, f"FERS_{cat.capitalize()}_EWC_{axis}_vs_{hodo_ax}_WithEnergy",
                               f"Hodo {axis} [cm]", (hodo_min, hodo_max),
                               f"EWC {axis} [cm]", (-15, 15),
                               style=PlotStyle(drawoptions=["colz"], addOverflow=False,
                                               addUnderflow=False, zmin=zmin, zmax=zmax,
                                               zlabel="Avg Energy [GeV]"),
                               extra_text=cat.capitalize(), usePDF=True)
            hp_hodo = infile.Get(f"hprof_HodoY_VS_HodoX_WithEnergy_{cat}")
            if hp_hodo:
                pm.plot_2d(hp_hodo, f"FERS_{cat.capitalize()}_HodoY_vs_HodoX_WithEnergy",
                           "Hodo X [cm]", (hodo_min, hodo_max),
                           "Hodo Y [cm]", (hodo_min, hodo_max),
                           style=PlotStyle(drawoptions=["colz"], addOverflow=False,
                                           addUnderflow=False, zmin=zmin, zmax=zmax,
                                           zlabel="Avg Energy [GeV]"),
                           extra_text=cat.capitalize(), usePDF=True)
        return pm.generate_html("FERS/EWC_vs_Hodo.html", plots_per_row=5,
                                title="EWC vs TTU Hodoscope")


def plot_fers_shower_shape(ctx):
    gain = "Mix"
    with _pm(ctx) as pm:
        pm.set_output_dir("FERS_ShowerShape")
        infile = pm._get_file("fers_shower_shape.root")
        hists_R = []
        for cat in ["cer", "sci"]:
            style = PlotStyle(dology=True, drawoptions="HIST",
                              mycolors=[2] if cat == "cer" else [4], donormalize=True)
            for axis in ["X", "Y"]:
                h = infile.Get(f"hist_Real{axis}_{gain}_{cat}")
                if h:
                    pm.plot_1d(h, f"FERS_ShowerShape_Real{axis}_{gain}_{cat}",
                               f"{axis} [cm]", (-20, 20),
                               ylabel="Frac. of Energy", yrange=(1e-4, 1), style=style,
                               extra_text=cat.capitalize())
            hR = infile.Get(f"hist_RealR_{gain}_{cat}")
            if hR:
                hists_R.append(hR)
        if hists_R:
            pm.plot_1d(hists_R, f"FERS_ShowerShape_RealR_{gain}",
                       "R [cm]", (0, 25), ylabel="Frac. of Energy",
                       yrange=(1e-4, 1), legends=["Cer", "Sci"],
                       style=PlotStyle(dology=True, drawoptions="HIST",
                                       mycolors=[2, 4], donormalize=True))
        return pm.generate_html("FERS/ShowerShape.html", plots_per_row=6,
                                title="FERS Shower Shape")


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
    ts_col = "ts" if mode == "raw" else get_arr_name(
        channel_name, use_mcp=(mode == "mcp"))
    pm.plot_2d(
        hist, plot_name, ts_col, (0, 1024), "DRS Output", (ymin, ymax),
        style=_STYLE_2D_LOG, extraToDraw=pave, extra_text=extra_text)


def plot_drs_waveforms(ctx, *, do_drs_vs_ts=True, do_mcp_vs_ts=False):
    """Plot DRS waveform 2D hists (TS correction modes). Returns list of HTML paths."""
    output_htmls = []

    if do_drs_vs_ts:
        # for mode in ["mcp"]:
        for mode in ["raw"]:
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
                            is_reference=chan.is_reference,
                            run_number=ctx.run_number)
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

    if do_mcp_vs_ts:
        list_mcp_channels = list(get_mcp_channels(ctx.run_number).values())
        with _pm(ctx) as pm:
            pm.set_output_dir("Service_DRS_VS_TS")
            infile = pm._get_file("drs_vs_ts.root")
            for channel_name in list_mcp_channels:
                ymin, ymax = get_drs_plot_ranges(
                    subtractMedian=True, is_mcp=True)
                pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                pave.AddText(channel_name)
                _plot_drs_channel_vs_ts(
                    pm, infile, channel_name, "raw", ymin, ymax, pave,
                    f"Service_DRS_ADC_VS_ts_raw_{channel_name}_blsub",
                    extra_text="MCP")
            output_htmls.append(
                pm.generate_html("DRS/Service_DRS_vs_ts_raw.html", plots_per_row=4))

    return output_htmls


def plot_drs_profiles(ctx, *, do_ts=True, do_time=True, do_mcp_only=False):
    """Plot DRS waveform profiles vs ts and vs time [ns]. Returns list of HTML paths."""
    output_htmls = []

    _VAR_ORDER = {"CerQuartz": 0, "CerPlastic": 1, "Sci": 2}
    _SIZES = ("3mm", "6mm")
    _VARS = ("CerQuartz", "CerPlastic", "Sci")
    _COLORS_3 = [ROOT.kRed+1, ROOT.kMagenta+1, ROOT.kBlue+1]
    _COLORS_6 = [ROOT.kRed+1, ROOT.kMagenta+1, ROOT.kBlue+1,
                 ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kCyan+2]

    def _data_yrange(hists, pad=0.25):
        lo, hi = float('inf'), float('-inf')
        for h in hists:
            for i in range(1, h.GetNbinsX() + 1):
                if h.GetBinError(i) > 0:
                    v = h.GetBinContent(i)
                    lo = min(lo, v)
                    hi = max(hi, v)
        if lo == float('inf'):
            return -1.0, 1.0
        span = hi - lo or 1.0
        return lo - pad * span, hi + pad * span

    def _align_peak(hists):
        result = []
        for h in hists:
            c = h.Clone()
            c.SetDirectory(0)
            peak = c.GetMaximum()
            if abs(peak) > 1e-10:
                c.Scale(1.0 / peak)
            result.append(c)
        return result

    def _combo_group_plots(pm, entries, group_name, xlabel, xrange_fn, prefix):
        if not entries:
            return
        hists_c = [e[0] for e in entries]
        labels_c = ([f"{e[2]} {e[1]}" for e in entries]
                    if group_name == "all" else [e[1] for e in entries])
        palette = _COLORS_6 if group_name == "all" else _COLORS_3
        xmin, xmax = xrange_fn(hists_c)
        ymin_c, ymax_c = _data_yrange(hists_c)
        pm.plot_1d(
            hists_c,
            f"{prefix}_{group_name}",
            xlabel, (xmin, xmax),
            "Mean DRS Output", (ymin_c, ymax_c),
            legends=labels_c,
            legendPos=[0.55, 0.70, 0.90, 0.90],
            style=_STYLE_1D_SMOOTH, mycolors=palette, extra_text=group_name,
            prepend=True)
        aligned = _align_peak(hists_c)
        ymin_a, ymax_a = _data_yrange(aligned)
        pm.plot_1d(
            aligned,
            f"{prefix}_{group_name}_aligned",
            xlabel, (xmin, xmax),
            "Mean DRS Output (a.u.)", (ymin_a, ymax_a),
            legends=labels_c,
            legendPos=[0.55, 0.70, 0.90, 0.90],
            style=_STYLE_1D_SMOOTH, mycolors=palette, extra_text=group_name,
            prepend=True)

    if do_ts:
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Prof_VS_TS")
            infile = pm._get_file("drs_profiles.root")
            pm.add_newline()

            for _, board in ctx.drsboards.items():
                for chan in board:
                    ch_blsub = chan.get_channel_name(blsub=True)

                    def _get_proj(name):
                        h = infile.Get(name)
                        return h.ProjectionX() if h else None

                    hist = _get_proj(f"prof_{ch_blsub}_VS_ts")
                    hist_ref = _get_proj(f"prof_{ch_blsub}_VS_ts_ref")
                    hist_mcp = _get_proj(f"prof_{ch_blsub}_VS_ts_mcp")
                    if not hist:
                        print(f"Warning: prof_{ch_blsub}_VS_ts not found")
                        continue
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

                    if do_mcp_only:
                        hists_to_plot = [h for h in [hist_mcp] if h]
                        legends_to_use = ["ts mcp"] if hist_mcp else []
                    else:
                        pairs = [(hist, "ts raw"), (hist_ref,
                                                    "ts ref"), (hist_mcp, "ts mcp")]
                        hists_to_plot = [h for h, _ in pairs if h]
                        legends_to_use = [lbl for h, lbl in pairs if h]

                    if hists_to_plot:
                        pm.plot_1d(
                            hists_to_plot,
                            f"DRS_ADC_Prof_VS_ts_{ch_blsub}_{var}",
                            "ts", (0, 1024),
                            "Mean DRS Output", (ymin_tmp, ymax_tmp),
                            legends=legends_to_use,
                            legendPos=[0.55, 0.70, 0.90, 0.90],
                            style=_STYLE_1D_SMOOTH, extraToDraw=pave, extra_text=var)

            combo_file = pm._get_file("drs_prof_combined.root")
            all_entries = []
            for size_tag in _SIZES:
                for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
                    h = combo_file.Get(f"combo_prof_{size_tag}_{var}_mcp")
                    if h:
                        all_entries.append((h, var, size_tag))
            grouped = {
                "3mm": [(h, v, s) for h, v, s in all_entries if s == "3mm"],
                "6mm": [(h, v, s) for h, v, s in all_entries if s == "6mm"],
                "all": all_entries,
            }
            for group_name, entries in [("3mm", grouped["3mm"]),
                                        ("6mm", grouped["6mm"]),
                                        ("all", grouped["all"])]:
                _combo_group_plots(pm, entries, group_name, "ts",
                                   lambda hs: (0, 1024),
                                   "DRS_ADC_Prof_VS_ts_Combined")

            output_htmls.append(
                pm.generate_html("DRS/DRS_Prof_vs_ts.html", plots_per_row=9))

    if do_time:
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Prof_VS_Time")
            infile = pm._get_file("drs_profiles.root")
            pm.add_newline()

            t_lo, t_hi = get_drs_time_arr_ns_range()
            for _, board in ctx.drsboards.items():
                for chan in board:
                    if chan.is_reference:
                        continue
                    ch_blsub = chan.get_channel_name(blsub=True)
                    hist_time_mcp = infile.Get(f"prof_{ch_blsub}_VS_time_mcp")
                    if not hist_time_mcp:
                        continue
                    ymin_tmp, ymax_tmp = get_drs_prof_plot_ranges(
                        subtractMedian=True,
                        is_amplified=chan.is_amplified,
                        is6mm=chan.is6mm,
                        is_reference=False,
                        is_cer=chan.isCer)
                    pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                    pave.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                    var = get_channel_var(chan)
                    pm.plot_1d(
                        [hist_time_mcp],
                        f"DRS_ADC_Prof_VS_time_mcp_{ch_blsub}_{var}",
                        "time [ns]", (t_lo, t_hi),
                        "Mean DRS Output", (ymin_tmp, ymax_tmp),
                        legends=["time"],
                        legendPos=[0.55, 0.70, 0.90, 0.90],
                        style=_STYLE_1D_SMOOTH, extraToDraw=pave, extra_text=var)

            try:
                corr_file = pm._get_file("drs_prof_corr_combined.root")
                corr_entries = []
                for size_tag in _SIZES:
                    for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
                        h = corr_file.Get(
                            f"combo_prof_corr_{size_tag}_{var}_mcp")
                        if h:
                            corr_entries.append((h, var, size_tag))
                corr_grouped = {
                    "3mm": [(h, v, s) for h, v, s in corr_entries if s == "3mm"],
                    "6mm": [(h, v, s) for h, v, s in corr_entries if s == "6mm"],
                    "all": corr_entries,
                }
                for group_name, entries in [("3mm", corr_grouped["3mm"]),
                                            ("6mm", corr_grouped["6mm"]),
                                            ("all", corr_grouped["all"])]:
                    _combo_group_plots(
                        pm, entries, group_name, "time [ns]",
                        lambda hs: (hs[0].GetXaxis().GetXmin(),
                                    hs[0].GetXaxis().GetXmax()),
                        "DRS_ADC_Prof_VS_corr_Combined")
            except FileNotFoundError:
                pass

            output_htmls.append(
                pm.generate_html("DRS/DRS_Prof_vs_time.html", plots_per_row=9))

    return output_htmls


# ---------------------------------------------------------------------------
# DRS: pulse statistics (peak, energy, timing)
# ---------------------------------------------------------------------------

def plot_drs_stats(ctx, *, do_peak=False, do_energy=True, do_timing=False,
                   do_timing_finebins_ts=True, do_timing_finebins_time=True):
    """Plot DRS peak, energy, and timing distributions. Returns list of HTML paths."""
    output_htmls = []

    if do_peak:
        # -- peak values --
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Stats")
            infile = pm._get_file("drs_stats.root")
            for _, board in ctx.drsboards.items():
                for chan in board:
                    if chan.is_reference:
                        continue
                    ch = chan.get_channel_name(blsub=False)
                    hist = infile.Get(f"hist_{ch}_peak_value")
                    if not hist:
                        print(f"Warning: hist_{ch}_peak_value not found")
                        continue
                    pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                    pave.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                    var = get_channel_var(chan)
                    pm.plot_1d(
                        hist, f"DRS_PeakValue_{ch}_{var}",
                        "DRS peak value", get_drs_peak_value_range()[
                            1:], "Counts", (0.9, None),
                        style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)
            output_htmls.append(pm.generate_html(
                "DRS/DRS_Peak.html", plots_per_row=9))

    if do_energy:
        # -- energy (CFD integral) --
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Stats")
            infile = pm._get_file("drs_stats.root")
            for _, board in ctx.drsboards.items():
                for chan in board:
                    if chan.is_reference:
                        continue
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
                        "CFD Energy", get_drs_energy_range(
                            chan.is6mm)[1:], "Counts", (0.9, None),
                        style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)
            output_htmls.append(pm.generate_html(
                "DRS/DRS_Sum.html", plots_per_row=9))

    if do_timing:
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
                    if not hp_ref or not hc_ref:
                        print(f"Warning: TS histograms for {ch} not found")
                        continue
                    var = get_channel_var(chan)
                    pave = create_pave_text(0.20, 0.80, 0.60, 0.90)
                    pave.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
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

    # shared constants for finebins pages
    _VAR_ORDER_FB = {"CerQuartz": 0, "CerPlastic": 1, "Sci": 2}
    _SIZE_ORDER_FB = {"3mm": 0, "6mm": 1}
    _VAR_IS_CER = {"CerQuartz": True, "CerPlastic": True, "Sci": False}
    _SIZES_FB = ("3mm", "6mm")
    _VARS_FB = ("CerQuartz", "CerPlastic", "Sci")
    _COLORS_3_FB = [ROOT.kRed+1, ROOT.kMagenta+1, ROOT.kBlue+1]
    _COLORS_6_FB = [ROOT.kRed+1, ROOT.kMagenta+1, ROOT.kBlue+1,
                    ROOT.kGreen+2, ROOT.kOrange+1, ROOT.kCyan+2]
    _NS_WINDOW = 1.2  # ±1.2 ns window for MPV fit

    if do_timing_finebins_ts:
        # -- fine-binned CFD TS (MCP-corrected) --
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Stats")
            infile = pm._get_file("drs_stats.root")

            try:
                combo_file = pm._get_file("drs_finebins_combined.root")
            except FileNotFoundError:
                combo_file = None

            all_entries = []
            if combo_file:
                for size_tag in _SIZES_FB:
                    for var in sorted(_VARS_FB, key=lambda v: _VAR_ORDER_FB[v]):
                        h = combo_file.Get(f"combo_finebins_{size_tag}_{var}")
                        if h:
                            all_entries.append((h, var, size_tag))

            grouped = {
                "3mm": [(h, v, s) for h, v, s in all_entries if s == "3mm"],
                "6mm": [(h, v, s) for h, v, s in all_entries if s == "6mm"],
                "all": all_entries,
            }
            for group_name in ("3mm", "6mm", "all"):
                entries = grouped[group_name]
                if not entries:
                    continue
                hists_combo = [e[0] for e in entries]
                labels_combo = ([f"{e[2]} {e[1]}" for e in entries]
                                if group_name == "all" else [e[1] for e in entries])
                palette = _COLORS_6_FB if group_name == "all" else _COLORS_3_FB
                pm.plot_1d(
                    hists_combo,
                    f"DRS_Time_FineBins_Combined_{group_name}",
                    "CFD TS (MCP-corrected)", (410, 490),
                    "Counts", (0.9, None),
                    legends=labels_combo,
                    legendPos=[0.55, 0.70, 0.90, 0.90],
                    style=_STYLE_1D_LOG, mycolors=palette, extra_text=group_name)

            pm.add_newline()

            for h, var, size_tag in sorted(
                    all_entries, key=lambda e: (_SIZE_ORDER_FB[e[2]], _VAR_ORDER_FB[e[1]])):
                xrange_fb = get_drs_cfd_finebins_range(_VAR_IS_CER[var])
                pave = create_pave_text(0.15, 0.75, 0.70, 0.90)
                pave.AddText(f"{size_tag} {var}")
                if h.GetEntries() >= 5:
                    mpv, mpv_err = get_hist_mpv(h)
                    pave.AddText(f"MPV: {mpv:.1f} #pm {mpv_err:.1f} TS")
                pm.plot_1d(
                    h, f"DRS_Time_FineBins_Summary_{size_tag}_{var}",
                    "CFD TS (MCP-corrected)", xrange_fb,
                    "Counts", (0.9, None),
                    style=_STYLE_1D_LOG, extraToDraw=pave,
                    extra_text=f"{size_tag} {var}")

            pm.add_newline()

            for _, board in ctx.drsboards.items():
                for chan in board:
                    if chan.is_reference:
                        continue
                    ch = chan.get_channel_name(blsub=False)
                    var = get_channel_var(chan)
                    hist = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
                    if not hist:
                        continue
                    pave = create_pave_text(0.15, 0.75, 0.70, 0.90)
                    pave.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    pave.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                    if hist.GetEntries() >= 5:
                        mpv, mpv_err = get_hist_mpv(hist)
                        pave.AddText(f"MPV: {mpv:.1f} #pm {mpv_err:.1f} TS")
                    pm.plot_1d(
                        hist, f"DRS_Time_FineBins_{ch}_{var}",
                        "CFD TS (MCP-corrected)", get_drs_cfd_finebins_range(chan.isCer),
                        "Counts", (0.9, None),
                        style=_STYLE_1D_LOG, extraToDraw=pave, extra_text=var)

            output_htmls.append(pm.generate_html(
                "DRS/DRS_Time_FineBins_TS.html", plots_per_row=9))

    if do_timing_finebins_time:
        # -- fine-binned CFD time [ns] (MCP + MPV corrected) --
        with _pm(ctx) as pm:
            pm.set_output_dir("DRS_Stats")
            infile = pm._get_file("drs_stats.root")

            try:
                corr_combo_file = pm._get_file(
                    "drs_finebins_corr_combined.root")
            except FileNotFoundError:
                corr_combo_file = None

            corr_entries = []
            if corr_combo_file:
                for size_tag in _SIZES_FB:
                    for var in sorted(_VARS_FB, key=lambda v: _VAR_ORDER_FB[v]):
                        h = corr_combo_file.Get(
                            f"combo_finebins_corr_{size_tag}_{var}")
                        if h:
                            corr_entries.append((h, var, size_tag))

            corr_grouped = {
                "3mm": [(h, v, s) for h, v, s in corr_entries if s == "3mm"],
                "6mm": [(h, v, s) for h, v, s in corr_entries if s == "6mm"],
                "all": corr_entries,
            }
            for group_name in ("3mm", "6mm", "all"):
                entries_c = corr_grouped[group_name]
                if not entries_c:
                    continue
                hists_combo = [e[0] for e in entries_c]
                labels_combo = ([f"{e[2]} {e[1]}" for e in entries_c]
                                if group_name == "all" else [e[1] for e in entries_c])
                palette = _COLORS_6_FB if group_name == "all" else _COLORS_3_FB
                pm.plot_1d(
                    hists_combo,
                    f"DRS_Time_FineBins_Corr_Combined_{group_name}",
                    "Time [ns]", (5, 15),
                    "Counts", (0.9, None),
                    legends=labels_combo,
                    legendPos=[0.55, 0.70, 0.90, 0.90],
                    style=_STYLE_1D_LOG, mycolors=palette, extra_text=group_name)

            pm.add_newline()

            for h, var, size_tag in sorted(
                    corr_entries, key=lambda e: (_SIZE_ORDER_FB[e[2]], _VAR_ORDER_FB[e[1]])):
                pave = create_pave_text(0.15, 0.75, 0.70, 0.90)
                pave.AddText(f"{size_tag} {var}")
                if h.GetEntries() >= 5:
                    mpv, mpv_err = get_hist_mpv(h, window_ts=_NS_WINDOW)
                    pave.AddText(f"MPV: {mpv:.2f} #pm {mpv_err:.2f} ns")
                pm.plot_1d(
                    h, f"DRS_Time_FineBins_Corr_Summary_{size_tag}_{var}",
                    "Time [ns]", (5, 15),
                    "Counts", (0.9, None),
                    style=_STYLE_1D_LOG, extraToDraw=pave,
                    extra_text=f"{size_tag} {var}")

            pm.add_newline()

            for _, board in ctx.drsboards.items():
                for chan in board:
                    if chan.is_reference:
                        continue
                    ch = chan.get_channel_name(blsub=False)
                    var = get_channel_var(chan)
                    hist_ns = infile.Get(f"hist_{ch}_Time_cfd_mcp_finebins")
                    if not hist_ns:
                        continue
                    tfb_lo, tfb_hi = get_drs_time_ns_finebins_range(chan.isCer)
                    pave_ns = create_pave_text(0.15, 0.75, 0.70, 0.90)
                    pave_ns.AddText(
                        f"B: {board.board_no}, G: {chan.group_no}, C: {chan.channel_no}")
                    pave_ns.AddText(
                        f"Tower: ({chan.i_tower_x}, {chan.i_tower_y})")
                    if hist_ns.GetEntries() >= 5:
                        mpv_ns, mpv_err_ns = get_hist_mpv(
                            hist_ns, window_ts=_NS_WINDOW)
                        pave_ns.AddText(
                            f"MPV: {mpv_ns:.2f} #pm {mpv_err_ns:.2f} ns")
                    pm.plot_1d(
                        hist_ns, f"DRS_Time_FineBins_ns_{ch}_{var}",
                        "Time [ns]", (tfb_lo, tfb_hi),
                        "Counts", (0.9, None),
                        style=_STYLE_1D_LOG, extraToDraw=pave_ns, extra_text=var)

            output_htmls.append(pm.generate_html(
                "DRS/DRS_Time_FineBins_Time.html", plots_per_row=9))

    return output_htmls


# ---------------------------------------------------------------------------
# DRS: CFD MPV board map (one value per channel)
# ---------------------------------------------------------------------------

def plot_drs_cfd_mpv(ctx):
    """Read the MPV of the fine-binned timing histograms and display as 2D board maps."""
    infile = ctx.hbook.open_file("drs_stats.root")

    mpv_map = {}       # all channels: raw TS MPV − 400
    mpv_map_cq = {}    # CerQuartz only
    mpv_map_cp = {}    # CerPlastic only
    ns_mpv_map = {}    # all channels: corrected Time [ns] MPV
    ns_mpv_map_cq = {}  # CerQuartz only
    ns_mpv_map_cp = {}  # CerPlastic only

    for _, board in ctx.drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            ch = chan.get_channel_name(blsub=False)
            var = get_channel_var(chan)
            hist_ts = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
            if hist_ts and hist_ts.GetEntries() >= 5:
                mpv, _ = get_hist_mpv(hist_ts)
                mpv_map[ch] = mpv - 400
                if var == "CerQuartz":
                    mpv_map_cq[ch] = mpv - 400
                elif var == "CerPlastic":
                    mpv_map_cp[ch] = mpv - 400
            else:
                mpv_map[ch] = 0.0
            hist_ns = infile.Get(f"hist_{ch}_Time_cfd_mcp_finebins")
            if hist_ns and hist_ns.GetEntries() >= 5:
                mpv_ns, _ = get_hist_mpv(hist_ns, window_ts=1.2)
                ns_mpv_map[ch] = round(mpv_ns, 2)
                if var == "CerQuartz":
                    ns_mpv_map_cq[ch] = round(mpv_ns, 2)
                elif var == "CerPlastic":
                    ns_mpv_map_cp[ch] = round(mpv_ns, 2)

    cer_hists, sci_hists = visualizeDRSBoards(
        ctx.drsboards, valuemaps=mpv_map,
        suffix=f"MPV_cfd_Run{ctx.run_number}")
    cq_hists, _ = visualizeDRSBoards(
        ctx.drsboards, valuemaps=mpv_map_cq,
        suffix=f"MPV_cfd_CerQuartz_Run{ctx.run_number}")
    cp_hists, _ = visualizeDRSBoards(
        ctx.drsboards, valuemaps=mpv_map_cp,
        suffix=f"MPV_cfd_CerPlastic_Run{ctx.run_number}")

    with _pm(ctx) as pm:
        pm.set_output_dir("DRS_CFD_MPV")
        helper = BoardPlotHelper(pm)

        helper.plot_cer_sci_pair(
            cer_hists, sci_hists,
            f"DRS_CFD_MPV_Run{ctx.run_number}",
            nTextDigits=1, zlabel="MPV [TS]")
        helper.plot_board_map(
            cq_hists, f"DRS_CFD_MPV_Run{ctx.run_number}_CerQuartz",
            extra_text="CerQuartz", nTextDigits=1, zlabel="MPV [TS]")
        helper.plot_board_map(
            cp_hists, f"DRS_CFD_MPV_Run{ctx.run_number}_CerPlastic",
            extra_text="CerPlastic", nTextDigits=1, zlabel="MPV [TS]")

        if ns_mpv_map:
            cer_hists_c, sci_hists_c = visualizeDRSBoards(
                ctx.drsboards, valuemaps=ns_mpv_map,
                suffix=f"MPV_cfd_corr_Run{ctx.run_number}")
            cq_hists_c, _ = visualizeDRSBoards(
                ctx.drsboards, valuemaps=ns_mpv_map_cq,
                suffix=f"MPV_cfd_corr_CerQuartz_Run{ctx.run_number}")
            cp_hists_c, _ = visualizeDRSBoards(
                ctx.drsboards, valuemaps=ns_mpv_map_cp,
                suffix=f"MPV_cfd_corr_CerPlastic_Run{ctx.run_number}")
            helper.plot_cer_sci_pair(
                cer_hists_c, sci_hists_c,
                f"DRS_CFD_MPV_corr_Run{ctx.run_number}",
                nTextDigits=2, zlabel="MPV [ns]")
            helper.plot_board_map(
                cq_hists_c, f"DRS_CFD_MPV_corr_Run{ctx.run_number}_CerQuartz",
                extra_text="CerQuartz", nTextDigits=2, zlabel="MPV [ns]")
            helper.plot_board_map(
                cp_hists_c, f"DRS_CFD_MPV_corr_Run{ctx.run_number}_CerPlastic",
                extra_text="CerPlastic", nTextDigits=2, zlabel="MPV [ns]")

        return pm.generate_html(
            "DRS/DRS_CFD_MPV.html", plots_per_row=4,
            title="DRS CFD MPV Board Map",
            intro_text=(
                "Per-channel MPV of the fine-binned MCP-corrected CFD time slice (TS board map, "
                "shifted by −400 TS) and the corrected time [ns] (Time_cfd_mcp_finebins). "
                "Channels with fewer than 5 entries are shown as 0."
            ))


# ---------------------------------------------------------------------------
# DRS: peak time-slice
# ---------------------------------------------------------------------------

def plot_drs_peak_ts(ctx, *, do_peak_ts=True, do_cer_vs_sci=True):
    output_htmls = []

    if do_peak_ts:
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

    if do_cer_vs_sci:
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
                            "Sci Peak TS", (400,
                                            600), "Cer Peak TS", (400, 600),
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
    _sum_ranges = get_drs_sum_vs_fers_ranges()
    _ranges = {"FERS": _sum_ranges["HG"],
               "FERSLG": _sum_ranges["LG"], "FERSMix": _sum_ranges["Mix"]}

    with _pm(ctx) as pm:
        pm.set_output_dir("DRSSum_VS_FERS")
        infile = pm._get_file("drssum_vs_fers.root")
        hists = []
        for _, board in ctx.drsboards.items():
            for i_tower_x, i_tower_y in board.get_list_of_towers():
                s_x = number_to_string(i_tower_x)
                s_y = number_to_string(i_tower_y)
                for var in ["Cer", "Sci"]:
                    for gain_tag in ["FERS", "FERSLG", "FERSMix"]:
                        hname = f"hist_DRSSum_VS_{gain_tag}_{var}_{s_x}_{s_y}"
                        hist = infile.Get(hname)
                        if not hist:
                            print(f"Warning: {hname} not found")
                            continue
                        hists.append(hist)
                        if ctx.do_detailed_plots:
                            zmax = round_up_to_1eN(
                                hist.Integral(0, 10000, 0, 10000))
                            tmp = _ranges[gain_tag][var]
                            pm.plot_2d(
                                hist, hname.replace("hist_", ""),
                                f"{gain_tag} [ADC]", (0, tmp[1]
                                                      ), "DRS Sum [ADC]", (0, tmp[0]),
                                style=PlotStyle(dologz=True, drawoptions="COLZ",
                                                zmin=1, zmax=zmax),
                                extra_text=var)

        for var in ["Cer", "Sci"]:
            for gain_tag in ["FERS", "FERSLG", "FERSMix"]:
                subset = [
                    h for h in hists if f"_{gain_tag}_{var}_" in h.GetName()]
                if subset:
                    combined = LHistos2Hist(
                        subset, f"hist_DRSSum_VS_{gain_tag}_{var}_Combined")
                    zmax = round_up_to_1eN(
                        combined.Integral(0, 10000, 0, 10000))
                    tmp = _ranges[gain_tag][var]
                    pm.plot_2d(
                        combined, f"DRSSum_VS_{gain_tag}_{var}_Combined",
                        f"{gain_tag} [ADC]", (0, tmp[1]
                                              ), "DRS Sum [ADC]", (0, tmp[0]),
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
                    for gain_tag, xmax in get_fers_vs_drs_xmax().items():
                        hname = f"hist_DRSPeak_VS_{gain_tag}_{var}_{s_x}_{s_y}"
                        hist = infile.Get(hname)
                        if not hist:
                            print(f"Warning: {hname} not found")
                            continue
                        pm.plot_2d(
                            hist, hname.replace("hist_", ""),
                            f"FERS {gain_tag} [ADC]", (0,
                                                       xmax), "DRS Peak [ADC]", (0, ymax),
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
        _, _, _, _, value_cut, cut_method = get_service_drs_cut(det, pm.run_number)

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
                # fail = below the cut bin; the cut bin itself goes to pass
                # (matches the 2D correlation convention; avoids double-count)
                nele = hist_pv.Integral(0, hist_pv.FindBin(value_cut) - 1)
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
                # fail = below the cut bin; the cut bin itself goes to pass
                # (matches the 2D correlation convention; avoids double-count)
                nele = hist_en.Integral(0, hist_en.FindBin(value_cut) - 1)
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


def _plot_pulse_correlations(channels, infile, pm, suffix, *, do_pid=True, do_trigger=True):
    """2D correlation plots between detector pairs."""
    det_list = list(channels.keys())
    trigger_dets = ["KT1", "KT2", "T3", "T4"]
    det_list_notrigger = [d for d in det_list if d not in trigger_dets]

    output_htmls = []
    corr_categories = [("PID", det_list_notrigger, do_pid),
                       ("Trigger", trigger_dets, do_trigger)]

    for cat, tmp_list, do_run in corr_categories:
        if not do_run:
            continue
        pm.reset_plots()
        for idx1, det1 in enumerate(tmp_list):
            for idx2, det2 in enumerate(tmp_list):
                if idx2 <= idx1:
                    continue
                _, _, _, _, value_cut1, method1 = get_service_drs_cut(det1, pm.run_number)
                _, _, _, _, value_cut2, method2 = get_service_drs_cut(det2, pm.run_number)
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


def _plot_pulse(ctx, channels, suffix, include_correlations=True, *,
                do_distributions=True, do_pid_correlations=True, do_trigger_correlations=True):
    """Combined plot runner: distributions + (optionally) correlations."""
    if not channels:
        print(f"No drs_{suffix} channels for run {ctx.run_number}; "
              f"nothing to plot.")
        return []

    infile_name = f"{ctx.paths['root']}/drs_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open {infile_name}")

    with _pm(ctx) as pm:
        pm.set_output_dir(f"drs_{suffix}")
        output_htmls = []
        if do_distributions:
            output_htmls.append(_plot_pulse_distributions(
                channels, infile, pm, suffix))
        if include_correlations:
            output_htmls += _plot_pulse_correlations(
                channels, infile, pm, suffix,
                do_pid=do_pid_correlations, do_trigger=do_trigger_correlations)

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

        # absolute CFD (ref-corrected) time per MCP (auto-zoom around the mean)
        for det in dets:
            hist = infile.Get(f"{det}_cfd_ref")
            if not hist:
                continue
            peak_pos = hist.GetBinCenter(hist.GetMaximumBin())
            pave = create_pave_text(0.55, 0.78, 0.90, 0.88)
            pave.AddText(f"Peak: {peak_pos:.2f} TS")
            pm.plot_1d(
                hist, f"{det}_cfd_ref",
                f"t_{{CFD,ref}} ({det}) [TS]", (420, 480),
                yrange=(0.9, hist.GetMaximum() * 1.4), ylabel="Counts",
                style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[1],
                                addOverflow=True, addUnderflow=True),
                extraToDraw=[pave])
        pm.add_newline()

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
                ymax = hist.GetMaximum() * 1.4
                pm.plot_1d(
                    hist, f"{det1}_cfd_diff_vs_{det2}",
                    f"#Delta t_{{CFD,ref}} ({det1} - {det2}) [TS]",
                    (-10, 10),
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

    hodo_pos_channels = build_hodo_pos_channels(run_number=ctx.run_number)
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

def plot_service_drs_pid(ctx, *, do_distributions=True, do_pid_correlations=True, do_trigger_correlations=True):
    return _plot_pulse(ctx, get_pid_channels(ctx.run_number), suffix="services",
                       do_distributions=do_distributions,
                       do_pid_correlations=do_pid_correlations,
                       do_trigger_correlations=do_trigger_correlations)


def plot_service_drs_mcp(ctx):
    return _plot_pulse(ctx, get_mcp_channels(ctx.run_number), suffix="mcp",
                       include_correlations=False)


def plot_service_drs_mcp_timing(ctx):
    return _plot_mcp_timing_diff(ctx, get_mcp_channels(ctx.run_number))


def plot_service_drs_hodo(ctx):
    return _plot_hodo_peak(ctx)


def _print_ttu_hodo_table(ctx, counts):
    """Print a table of event counts after each TTU-hodo requirement combination."""
    total = counts.get("all", 0)
    width = max((len(l) for l in counts), default=12)
    print(f"\n{'='*(width+30)}")
    print(f"{'TTU Hodo selection counts (Run ' + str(ctx.run_number) + ')':^{width+30}}")
    print(f"{'='*(width+30)}")
    print(f"{'Requirement':<{width}} | {'Events':>10} | {'Frac.':>8}")
    print(f"{'-'*(width+30)}")
    for label, n in counts.items():
        frac = (100.0 * n / total) if total else 0.0
        print(f"{label.replace('_', ' '):<{width}} | {n:>10} | {frac:>7.2f}%")
    print(f"{'='*(width+30)}\n")


def _ttu_hodo_table_html(counts):
    """Build a single-line HTML table of TTU-hodo counts/fractions for the page."""
    total = counts.get("all", 0)
    th = ("padding:3px 14px;border:1px solid #ccc;text-align:right")
    th_l = ("padding:3px 14px;border:1px solid #ccc;text-align:left")
    rows = (f"<tr><th style='{th_l}'>Requirement</th>"
            f"<th style='{th}'>Events</th><th style='{th}'>Fraction</th></tr>")
    for label, n in counts.items():
        frac = (100.0 * n / total) if total else 0.0
        rows += (f"<tr><td style='{th_l}'>{label.replace('_', ' ')}</td>"
                 f"<td style='{th}'>{n}</td><td style='{th}'>{frac:.2f}%</td></tr>")
    return f"<table style='border-collapse:collapse;font-size:13px'>{rows}</table>"


def plot_service_ttu_hodo(ctx):
    """Plot TTU hodoscope XY over all detector requirement combinations + table."""
    from analysis.hist_functions import _ttu_hodo_labels
    hodo_min, hodo_max, _ = get_ttu_hodo_ranges()
    labels = _ttu_hodo_labels()
    counts = {}
    with _pm(ctx) as pm:
        pm.set_output_dir("TTU_Hodo")
        infile = pm._get_file("ttu_hodo.root")
        for label in labels:
            h = infile.Get(f"TTU_Hodo_XY_{label}")
            if not h:
                continue
            # weighted integral (incl. overflow) = #events passing the requirement
            counts[label] = int(round(
                h.Integral(0, h.GetNbinsX() + 1, 0, h.GetNbinsY() + 1)))
            pave = create_pave_text(0.18, 0.82, 0.70, 0.90)
            pave.AddText(label.replace("_", " "))
            pm.plot_2d(h, f"TTU_Hodo_XY_{label}",
                       "Hodo X [cm]", (hodo_min, hodo_max),
                       "Hodo Y [cm]", (hodo_min, hodo_max),
                       style=_STYLE_2D_LOG, extraToDraw=pave)
        _print_ttu_hodo_table(ctx, counts)
        return pm.generate_html("ServiceDRS/TTU_Hodo.html", plots_per_row=4,
                                title="TTU Hodoscope XY by requirement combination",
                                intro_text=_ttu_hodo_table_html(counts))
