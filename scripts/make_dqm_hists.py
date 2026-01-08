import json
import ROOT
from channels.channel_map import (
    build_hodo_pos_channels,
    build_time_reference_channels,
    get_mcp_channels,
    get_service_drs_channels,
)
from configs.plot_ranges import get_drs_plot_ranges, get_service_drs_plot_ranges
from core.analysis_manager import CaloXAnalysisManager
from utils.root_setup import setup_root
from utils.parser import get_args
from utils.plot_helper import get_run_paths, save_hists_to_file
from utils.timing import auto_timer
from utils.utils import number_to_string
from variables.drs import get_drs_stats
from variables.fers import get_fers_energy_max, get_fers_energy_sum

auto_timer("Total Execution Time")

setup_root()

args = get_args()

analysis = CaloXAnalysisManager(args).prepare()
rdf = analysis.get_rdf()
run_number = analysis.run_number
DRSBoards = analysis.drsboards
fersboards = analysis.fersboards

rdf = get_fers_energy_max(rdf, fersboards, gain="HG")
rdf = get_fers_energy_max(rdf, fersboards, gain="LG")
rdf = get_drs_stats(rdf, run_number, DRSBoards, 0, 1000, 9)

rdf = get_fers_energy_sum(rdf, fersboards, gain="HG")
rdf = get_fers_energy_sum(rdf, fersboards, gain="LG")


debug_drs = False
n_events = 60000
n_bins_event = 500


def monitor_conditions():
    # monitor the V, I, T conditions
    hprofs_Condition_VS_Event = []
    for fersboard in fersboards.values():
        hprof_SipmHV = rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event",
            f"FERS Board - SipmHV VS Event;Event;SipmHV (V)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_sipm_hv_name()
        )
        hprof_SipmI = rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_i_name()}_VS_Event",
            f"FERS Board - SipmI VS Event;Event;SipmI (mA)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_sipm_i_name()
        )
        hprof_TempDET = rdf.Profile1D((
            f"hprof_{fersboard.get_temp_det_name()}_VS_Event",
            f"FERS Board - TempDET VS Event;Event;TempDET (C)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_temp_det_name()
        )
        hprof_TempFPGA = rdf.Profile1D((
            f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event",
            f"FERS Board - TempFPGA VS Event;Event;TempFPGA (C)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_temp_fpga_name()
        )
        hprofs_Condition_VS_Event.append(hprof_SipmHV)
        hprofs_Condition_VS_Event.append(hprof_SipmI)
        hprofs_Condition_VS_Event.append(hprof_TempDET)
        hprofs_Condition_VS_Event.append(hprof_TempFPGA)

    return hprofs_Condition_VS_Event


def monitor_fers_energy_sum():
    hprofs_EnergySum_VS_Event = []
    for fersboard in fersboards.values():
        # monitor FERS energy sum
        for var in ["Cer", "Sci"]:
            for gain in ["HG", "LG"]:
                hprof_EnergySum = rdf.Profile1D((
                    f'hprof_{fersboard.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))}_VS_Event',
                    f"FERS Board - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                    n_bins_event, 0, n_events),
                    "event_n", fersboard.get_energy_sum_name(
                        gain=gain, isCer=(var == "Cer"))
                )
                hprofs_EnergySum_VS_Event.append(hprof_EnergySum)

    for var in ["Cer", "Sci"]:
        for gain in ["HG", "LG"]:
            hprof_EnergySum = rdf.Profile1D((
                f'hprof_{fersboards.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))}_VS_Event',
                f"FERS - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                n_bins_event, 0, n_events),
                "event_n", fersboards.get_energy_sum_name(
                    gain=gain, isCer=(var == "Cer"))
            )
            hprofs_EnergySum_VS_Event.append(hprof_EnergySum)

    return hprofs_EnergySum_VS_Event


def makeFERS1DHists(gain="HG"):
    hists1d_FERS = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist = rdf.Histo1D((
                    f"hist_FERS_Board{board_no}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS Board {board_no} - {var} i_tower_x {sTowerX} i_tower_y {sTowerY};{var} Energy;Counts",
                    3000, 0, 9000),
                    chan.get_channel_name(gain=gain)
                )
                hists1d_FERS.append(hist)

    return hists1d_FERS


def collectFERSPedestals(hists1d_FERS, gain="HG"):
    pedestals = {}
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                channelName = chan.get_channel_name(gain=gain)

                hname = f"hist_FERS_Board{board_no}_{var}_{sTowerX}_{sTowerY}"
                hist = None
                for h in hists1d_FERS:
                    if h.GetName() == hname:
                        hist = h
                        break
                if hist is None:
                    print(
                        f"Warning: Histogram {hname} not found in hists1d_FERS")
                    pedestals[channelName] = None
                    continue
                ibinmin = hist.GetXaxis().FindBin(100.0)
                ibinmax = hist.GetXaxis().FindBin(500.0)
                maxVal = -1
                ibinPed = -1
                for ibin in range(ibinmin, ibinmax):
                    if hist.GetBinContent(ibin) > maxVal:
                        maxVal = hist.GetBinContent(ibin)
                        ibinPed = ibin
                # pedestal = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())
                pedestal = hist.GetXaxis().GetBinCenter(ibinPed)
                pedestals[channelName] = pedestal
    return pedestals


def collectFERSStats():
    from configs.plot_ranges import getFERSSaturationValue
    saturation_value = getFERSSaturationValue()
    stats = {}
    # mean, max,
    # and how frequent the saturation value is reached
    for fersboard in fersboards.values():
        for chan in fersboard:
            channelName_HG = chan.get_channel_name(gain="HG")
            stats[channelName_HG] = (
                rdf.Mean(channelName_HG),
                rdf.Max(channelName_HG),
                rdf.Filter(f"{channelName_HG} >= {saturation_value}").Count()
            )
            channelName_LG = chan.get_channel_name(gain="LG")
            stats[channelName_LG] = (
                rdf.Mean(channelName_LG),
                rdf.Max(channelName_LG),
                rdf.Filter(f"{channelName_LG} >= {saturation_value}").Count()
            )

    return stats


def makeFERSMaxValueHists():
    """
    get the max FERS readout per board and per event 
    """
    hists_board_cer_HG_max = []
    hists_board_cer_LG_max = []
    hists_board_sci_HG_max = []
    hists_board_sci_LG_max = []
    nbins = 100
    xmin = 7500
    xmax = 8500
    for fersboard in fersboards.values():
        hist_board_cer_HG_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='HG', isCer=True)}",
            f"FERS Board - CER Energy HG max;CER Energy HG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='HG', isCer=True)
        )
        hist_board_cer_LG_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='LG', isCer=True)}",
            f"FERS Board - CER Energy LG max;CER Energy LG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='LG', isCer=True)
        )
        hist_board_sci_HG_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='HG', isCer=False)}",
            f"FERS Board - SCI Energy HG max;SCI Energy HG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='HG', isCer=False)
        )
        hist_board_sci_LG_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='LG', isCer=False)}",
            f"FERS Board - SCI Energy LG max;SCI Energy LG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='LG', isCer=False)
        )
        hists_board_cer_HG_max.append(hist_board_cer_HG_max)
        hists_board_cer_LG_max.append(hist_board_cer_LG_max)
        hists_board_sci_HG_max.append(hist_board_sci_HG_max)
        hists_board_sci_LG_max.append(hist_board_sci_LG_max)

    hist_cer_HG_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='HG', isCer=True)}",
        "FERS - CER Energy HG max;CER Energy HG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='HG', isCer=True)
    )
    hist_cer_LG_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='LG', isCer=True)}",
        "FERS - CER Energy LG max;CER Energy LG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='LG', isCer=True)
    )
    hist_sci_HG_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='HG', isCer=False)}",
        "FERS - SCI Energy HG max;SCI Energy HG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='HG', isCer=False)
    )
    hist_sci_LG_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='LG', isCer=False)}",
        "FERS - SCI Energy LG max;SCI Energy LG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='LG', isCer=False)
    )
    return hists_board_cer_HG_max + hists_board_cer_LG_max + hists_board_sci_HG_max + hists_board_sci_LG_max + [hist_cer_HG_max, hist_cer_LG_max, hist_sci_HG_max, hist_sci_LG_max]


def trackFERSHists():
    hists2d_FERS_VS_Event = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist = rdf.Histo2D((
                    f"hist_FERS_Board{board_no}_{var}_VS_Event_{sTowerX}_{sTowerY}",
                    f"FERS Board {board_no} - Event VS {var} {chan.channel_no} in i_tower_x {sTowerX} i_tower_y {sTowerY};Event;{var} Energy HG",
                    n_bins_event, 0, n_events, 1000, 0, 9000),
                    "event_n", chan.get_channel_name(gain="HG")
                )
                hists2d_FERS_VS_Event.append(hist)
    return hists2d_FERS_VS_Event


def makeFERS2DHists():
    hists2d_FERS = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)
            chan_Cer = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            chan_Sci = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=False)

            iCer = chan_Cer.channel_no
            iSci = chan_Sci.channel_no
            # high gain
            hist = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{sTowerX}_{sTowerY}",
                f"CER {iCer} VS SCI {iSci} in i_tower_x {sTowerX} i_tower_y {sTowerY};Sci Energy HG;Cer Energy HG",
                300, 0, 9000, 300, 0, 9000),
                chan_Sci.get_channel_name(gain="HG"),
                chan_Cer.get_channel_name(gain="HG")
            )
            hist_zoomed = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{sTowerX}_{sTowerY}_zoom",
                f"CER {iCer} VS SCI {iSci} in i_tower_x {sTowerX} i_tower_y {sTowerY} (zoomed);SCI Energy HG;CER Energy HG",
                300, 0, 1000, 200, 0, 2000),
                chan_Sci.get_channel_name(gain="HG"),
                chan_Cer.get_channel_name(gain="HG")
            )
            hists2d_FERS.append(hist)
            hists2d_FERS.append(hist_zoomed)

            # low gain (Y) VS high gain (X) for Sci
            hist_sci_lg_vs_hg = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Sci_{sTowerX}_{sTowerY}_hg_VS_lg",
                f"SCI {iSci} LG VS HG;SCI Energy HG;SCI Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_Sci.get_channel_name(gain="HG"),
                chan_Sci.get_channel_name(gain="LG")
            )
            hists2d_FERS.append(hist_sci_lg_vs_hg)
            hist_cer_lg_vs_hg = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_{sTowerX}_{sTowerY}_hg_VS_lg",
                f"CER {iCer} LG VS HG;CER Energy HG;CER Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_Cer.get_channel_name(gain="HG"),
                chan_Cer.get_channel_name(gain="LG")
            )
            hists2d_FERS.append(hist_cer_lg_vs_hg)
    return hists2d_FERS


def makeDRS2DHists(debug=False):
    hists2d_DRS_VS_TS = []
    if debug:
        hists2d_DRS_VS_RTSpos = []
        hists2d_DRS_VS_RTSneg = []
    for _, DRSBoard in DRSBoards.items():
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan is None:
                    continue
                channelName = chan.get_channel_name(blsub=True)
                ymin, ymax = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm)
                hist_subtractMedian = rdf.Histo2D((
                    f"hist_{channelName}_VS_TS",
                    "DRS values (subtract baseline);TS;DRS values",
                    1024, 0, 1024, 50, ymin, ymax),
                    "TS", channelName
                )
                hists2d_DRS_VS_TS.append(hist_subtractMedian)

                if debug:
                    hist_subtractMedian_RTSpos = rdf.Histo2D((
                        f"hist_{channelName}_VS_TS_RTSpos",
                        "DRS values (subtract baseline) VS RTS pos;RTS pos;DRS values",
                        1024, 0, 1024, 50, ymin, ymax),
                        f"RTS_pos_{channelName}", channelName
                    )
                    hist_subtractMedian_RTSneg = rdf.Histo2D((
                        f"hist_{channelName}_VS_TS_RTSneg",
                        "DRS values (subtract baseline) VS RTS neg;RTS neg;DRS values",
                        1024, 0, 1024, 50, ymin, ymax),
                        f"RTS_neg_{channelName}", channelName
                    )
                    hists2d_DRS_VS_RTSpos.append(hist_subtractMedian_RTSpos)
                    hists2d_DRS_VS_RTSneg.append(hist_subtractMedian_RTSneg)
    if debug:
        return hists2d_DRS_VS_TS, hists2d_DRS_VS_RTSpos, hists2d_DRS_VS_RTSneg
    return hists2d_DRS_VS_TS


def compareDRSChannels(channels_to_compare, isServiceDRS=False):
    if isinstance(channels_to_compare, dict):
        channels_to_compare = list(channels_to_compare.values())
    hists_trigger = []
    for chan_name in channels_to_compare:
        ymin = -2500
        ymax = 500
        if isServiceDRS:
            ymin, ymax = get_service_drs_plot_ranges(
                chan_name, subtractMedian=True)
        hist_subtractMedian = rdf.Histo2D((
            f"hist_{chan_name}_blsub",
            f"{chan_name} (subtract baseline);TS;DRS values",
            1024, 0, 1024,
            300, ymin, ymax),
            "TS", chan_name + "_blsub"
        )
        hists_trigger.append(hist_subtractMedian)
    return hists_trigger


def checkDRSSumVSFERS():
    xymax = {
        "Cer": (20000, 8500),
        "Sci": (30000, 8500)
    }
    xymax_LG = {
        "Cer": (20000, 2000),
        "Sci": (30000, 4000)
    }

    # correlate DRS sum and FERS outputs
    h2s_DRSSum_VS_FERS = []
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue
                chan_FERS = None
                for fersboard in fersboards.values():
                    chan_FERS = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_FERS is not None:
                        break
                if chan_FERS is None:
                    print(
                        f"Warning: FERS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                h2_DRSSum_VS_FERS = rdf.Histo2D((
                    f"hist_DRSSum_VS_FERS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS sum VS FERS energy correlation for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}; FERS Energy HG; DRS Sum",
                    100, 0, xymax[var][1], 100, 0, xymax[var][0]
                ),
                    chan_FERS.get_channel_name(gain="HG"),
                    chan_DRS.get_channel_sum_name(),
                )
                h2s_DRSSum_VS_FERS.append(h2_DRSSum_VS_FERS)

                h2_DRSSum_VS_FERSLG = rdf.Histo2D((
                    f"hist_DRSSum_VS_FERSLG_{var}_{sTowerX}_{sTowerY}",
                    f"DRS sum VS FERS LG energy correlation for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, 0, xymax_LG[var][1], 100, 0, xymax_LG[var][0]
                ),
                    chan_FERS.get_channel_name(gain="LG"),
                    chan_DRS.get_channel_sum_name(),
                )
                h2s_DRSSum_VS_FERS.append(h2_DRSSum_VS_FERSLG)

    return h2s_DRSSum_VS_FERS


def checkDRSPeakVSFERS():
    h2s_DRSPeak_VS_FERS = []
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue
                chan_FERS = None
                for fersboard in fersboards.values():
                    chan_FERS = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_FERS is not None:
                        break
                if chan_FERS is None:
                    print(
                        f"Warning: FERS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                _, ymax = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan_DRS.is_amplified, is6mm=chan_DRS.is6mm)
                h2_DRSPeak_VS_FERS = rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak VS FERS energy correlation for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}; FERS Energy HG; DRS Peak",
                    100, 0, 9000, 100, 0, ymax
                ),
                    chan_FERS.get_channel_name(gain="HG"),
                    chan_DRS.get_channel_peak_name(),
                )
                h2s_DRSPeak_VS_FERS.append(h2_DRSPeak_VS_FERS)

                h2_DRSPeak_VS_FERSLG = rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERSLG_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak VS FERS LG energy correlation for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}; FERS Energy LG; DRS Peak",
                    100, 0, 3000, 100, 0, ymax
                ),
                    chan_FERS.get_channel_name(gain="LG"),
                    chan_DRS.get_channel_peak_name(),
                )
                h2s_DRSPeak_VS_FERS.append(h2_DRSPeak_VS_FERSLG)

    return h2s_DRSPeak_VS_FERS


def checkDRSPeakTS():
    h1s_DRSPeakTS = {}
    h1s_DRSPeakTS["Cer"] = []
    h1s_DRSPeakTS["Sci"] = []
    h2s_DRSPeakTS_Cer_VS_Sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_DRS is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.get_channel_peak_ts_name()
                channelNames[var] = channelName

                h1_DRSPeakTS = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak TS for Board{board_no}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
                    1000, 0, 1000),
                    f"{channelName}_good"
                )
                h1s_DRSPeakTS[var].append(h1_DRSPeakTS)

            if len(channelNames) < 2:
                # print(
                #    f"Warning: Not enough channels found for Board{boardNo}, Tower({sTowerX}, {sTowerY})")
                continue

            h2_DRSPeak_Cer_VS_Sci = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}",
                f"DRS Peak TS - CER VS SCI for Board{board_no}, Tower({sTowerX}, {sTowerY});SCI Peak TS;CER Peak TS",
                1000, 0, 1000, 1000, 0, 1000),
                f'{channelNames["Sci"]}_good',
                f'{channelNames["Cer"]}_good',
            )
            h2s_DRSPeakTS_Cer_VS_Sci.append(h2_DRSPeak_Cer_VS_Sci)
    return h1s_DRSPeakTS["Cer"], h1s_DRSPeakTS["Sci"], h2s_DRSPeakTS_Cer_VS_Sci


def main():
    paths = get_run_paths(run_number)

    hprofs_conditions = monitor_conditions()
    hprofs_energysum = monitor_fers_energy_sum()

    hists1d_FERS = makeFERS1DHists()
    hists1d_FRRS_LG = makeFERS1DHists(gain="LG")

    # hists2d_FERS = makeFERS2DHists()
    # hists2d_FERS_VS_Event = trackFERSHists()

    hists2d_DRS_VS_RTSpos = None
    hists2d_DRS_VS_RTSneg = None
    if debug_drs:
        hists2d_DRS_VS_TS, hists2d_DRS_VS_RTSpos, hists2d_DRS_VS_RTSneg = makeDRS2DHists(
            debug=True)
    else:
        hists2d_DRS_VS_TS = makeDRS2DHists(debug=False)

    hists2d_DRSPeak_VS_FERS = checkDRSPeakVSFERS()

    hists1d_DRSPeakTS_Cer, hists1d_DRSPeakTS_Sci, hists2d_DRSPeakTS_Cer_VS_Sci = checkDRSPeakTS()

    time_reference_channels = build_time_reference_channels(run=run_number)
    hists2d_time_reference = compareDRSChannels(time_reference_channels)

    hists2d_hodo_trigger = None
    # hodo_trigger_channels = build_hodo_trigger_channels(run=runNumber)
    # hists2d_hodo_trigger = compareDRSChannels(hodo_trigger_channels)

    service_drs_channels = get_service_drs_channels(run=run_number)
    hists2d_service_drs = compareDRSChannels(
        service_drs_channels, isServiceDRS=True)

    hodo_pos_channels = build_hodo_pos_channels(run=run_number)
    channels = [channel for channels in hodo_pos_channels.values()
                for channel in channels]
    hists2d_hodo_pos = compareDRSChannels(channels)

    mcp_channels_map = get_mcp_channels(run=run_number)
    mcp_channels = [channel for channels in mcp_channels_map.values()
                    for channel in channels]
    hists2d_mcp = compareDRSChannels(mcp_channels)

    hists2d_DRSSum_VS_FERS = checkDRSSumVSFERS()

    hists_FERS_max = makeFERSMaxValueHists()

    stats = collectFERSStats()

    # pedestals
    # event loop should be triggered here
    pedestals_HG = collectFERSPedestals(hists1d_FERS, gain="HG")
    pedestals_LG = collectFERSPedestals(hists1d_FRRS_LG, gain="LG")

    # print(f"Loops executed after pedestals: {rdf.GetNRuns()}")

    print("\033[94mSave results\033[0m")

    # dump stats into a json file
    stats_results = {}
    for channelName, (mean, max_value, sat_frequency) in stats.items():
        stats_results[channelName] = (mean.GetValue(), max_value.GetValue(), float(
            sat_frequency.GetValue()) / n_events)
    with open(f"{paths['root']}/fers_stats.json", "w") as f:
        json.dump(stats_results, f, indent=4)

    with open(f"{paths['root']}/fers_pedestals_hg.json", "w") as f:
        json.dump(pedestals_HG, f, indent=4)
    with open(f"{paths['root']}/fers_pedestals_lg.json", "w") as f:
        json.dump(pedestals_LG, f, indent=4)

    # Save histograms to an output ROOT file
    save_hists_to_file(hprofs_conditions,
                       f"{paths['root']}/conditions_vs_event.root")

    save_hists_to_file(hprofs_energysum,
                       f"{paths['root']}/fers_energysum_vs_event.root")

    save_hists_to_file(hists1d_FERS,
                       f"{paths['root']}/fers_all_channels_1d.root")

    save_hists_to_file(hists2d_DRS_VS_TS, f"{paths['root']}/drs_vs_ts.root")
    if debug_drs:
        save_hists_to_file(hists2d_DRS_VS_RTSpos,
                           f"{paths['root']}/drs_VS_RTSpos.root")
        save_hists_to_file(hists2d_DRS_VS_RTSneg,
                           f"{paths['root']}/drs_VS_RTSneg.root")

    save_hists_to_file(hists2d_DRSPeak_VS_FERS,
                       f"{paths['root']}/drspeak_vs_fers.root")

    save_hists_to_file(hists1d_DRSPeakTS_Cer + hists1d_DRSPeakTS_Sci +
                       hists2d_DRSPeakTS_Cer_VS_Sci, f"{paths['root']}/drspeakts.root")

    save_hists_to_file(hists2d_time_reference,
                       f"{paths['root']}/time_reference_channels.root")

    save_hists_to_file(hists2d_hodo_trigger,
                       f"{paths['root']}/hodo_trigger_channels.root")

    save_hists_to_file(hists2d_service_drs,
                       f"{paths['root']}/service_drs_channels.root")

    save_hists_to_file(hists2d_mcp, f"{paths['root']}/mcp_channels.root")

    save_hists_to_file(
        hists2d_hodo_pos, f"{paths['root']}/hodo_pos_channels.root")

    save_hists_to_file(hists2d_DRSSum_VS_FERS,
                       f"{paths['root']}/drssum_vs_fers.root")

    save_hists_to_file(hists_FERS_max, f"{paths['root']}/fers_max_values.root")


if __name__ == "__main__":
    main()
