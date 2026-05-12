import json
from channels.channel_map import get_mcp_channels
from configs.plot_config import get_drs_plot_ranges, get_fers_saturation_value
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
    hprofs_condition_vs_event = []
    for fersboard in fersboards.values():
        hprof_sipm_hv = rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event",
            f"FERS Board - SipmHV VS Event;Event;SipmHV (V)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_sipm_hv_name()
        )
        hprof_sipm_i = rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_i_name()}_VS_Event",
            f"FERS Board - SipmI VS Event;Event;SipmI (mA)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_sipm_i_name()
        )
        hprof_temp_det = rdf.Profile1D((
            f"hprof_{fersboard.get_temp_det_name()}_VS_Event",
            f"FERS Board - TempDET VS Event;Event;TempDET (C)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_temp_det_name()
        )
        hprof_temp_fpga = rdf.Profile1D((
            f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event",
            f"FERS Board - TempFPGA VS Event;Event;TempFPGA (C)",
            n_bins_event, 0, n_events),
            "event_n", fersboard.get_temp_fpga_name()
        )
        hprofs_condition_vs_event.append(hprof_sipm_hv)
        hprofs_condition_vs_event.append(hprof_sipm_i)
        hprofs_condition_vs_event.append(hprof_temp_det)
        hprofs_condition_vs_event.append(hprof_temp_fpga)

    return hprofs_condition_vs_event


def monitor_fers_energy_sum():
    hprofs_energy_sum_vs_event = []
    for fersboard in fersboards.values():
        # monitor FERS energy sum
        for var in ["Cer", "Sci"]:
            for gain in ["HG", "LG"]:
                hprof_energy_sum = rdf.Profile1D((
                    f'hprof_{fersboard.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))}_VS_Event',
                    f"FERS Board - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                    n_bins_event, 0, n_events),
                    "event_n", fersboard.get_energy_sum_name(
                        gain=gain, isCer=(var == "Cer"))
                )
                hprofs_energy_sum_vs_event.append(hprof_energy_sum)

    for var in ["Cer", "Sci"]:
        for gain in ["HG", "LG"]:
            hprof_energy_sum = rdf.Profile1D((
                f'hprof_{fersboards.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))}_VS_Event',
                f"FERS - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                n_bins_event, 0, n_events),
                "event_n", fersboards.get_energy_sum_name(
                    gain=gain, isCer=(var == "Cer"))
            )
            hprofs_energy_sum_vs_event.append(hprof_energy_sum)

    return hprofs_energy_sum_vs_event


def make_fers_1d_hists(gain="HG"):
    hists1d_fers = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist = rdf.Histo1D((
                    f"hist_FERS_Board{board_no}_{var}_{s_tower_x}_{s_tower_y}",
                    f"FERS Board {board_no} - {var} i_tower_x {s_tower_x} i_tower_y {s_tower_y};{var} Energy;Counts",
                    3000, 0, 9000),
                    chan.get_channel_name(gain=gain)
                )
                hists1d_fers.append(hist)

    return hists1d_fers


def collect_fers_pedestals(hists1d_fers, gain="HG"):
    pedestals = {}
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                channel_name = chan.get_channel_name(gain=gain)

                hname = f"hist_FERS_Board{board_no}_{var}_{s_tower_x}_{s_tower_y}"
                hist = None
                for h in hists1d_fers:
                    if h.GetName() == hname:
                        hist = h
                        break
                if hist is None:
                    print(
                        f"Warning: Histogram {hname} not found in hists1d_fers")
                    pedestals[channel_name] = None
                    continue
                ibinmin = hist.GetXaxis().FindBin(100.0)
                ibinmax = hist.GetXaxis().FindBin(500.0)
                max_val = -1
                ibin_ped = -1
                for ibin in range(ibinmin, ibinmax):
                    if hist.GetBinContent(ibin) > max_val:
                        max_val = hist.GetBinContent(ibin)
                        ibin_ped = ibin
                # pedestal = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())
                pedestal = hist.GetXaxis().GetBinCenter(ibin_ped)
                pedestals[channel_name] = pedestal
    return pedestals


def collect_fers_stats():
    saturation_value = get_fers_saturation_value()
    stats = {}
    # mean, max,
    # and how frequent the saturation value is reached
    for fersboard in fersboards.values():
        for chan in fersboard:
            channel_name_hg = chan.get_channel_name(gain="HG")
            stats[channel_name_hg] = (
                rdf.Mean(channel_name_hg),
                rdf.Max(channel_name_hg),
                rdf.Filter(f"{channel_name_hg} >= {saturation_value}").Count()
            )
            channel_name_lg = chan.get_channel_name(gain="LG")
            stats[channel_name_lg] = (
                rdf.Mean(channel_name_lg),
                rdf.Max(channel_name_lg),
                rdf.Filter(f"{channel_name_lg} >= {saturation_value}").Count()
            )

    return stats


def make_fers_max_value_hists():
    """
    get the max FERS readout per board and per event 
    """
    hists_board_cer_hg_max = []
    hists_board_cer_lg_max = []
    hists_board_sci_hg_max = []
    hists_board_sci_lg_max = []
    nbins = 100
    xmin = 7500
    xmax = 8500
    for fersboard in fersboards.values():
        hist_board_cer_hg_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='HG', isCer=True)}",
            f"FERS Board - CER Energy HG max;CER Energy HG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='HG', isCer=True)
        )
        hist_board_cer_lg_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='LG', isCer=True)}",
            f"FERS Board - CER Energy LG max;CER Energy LG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='LG', isCer=True)
        )
        hist_board_sci_hg_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='HG', isCer=False)}",
            f"FERS Board - SCI Energy HG max;SCI Energy HG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='HG', isCer=False)
        )
        hist_board_sci_lg_max = rdf.Histo1D((
            f"hist_{fersboard.get_energy_max_name(gain='LG', isCer=False)}",
            f"FERS Board - SCI Energy LG max;SCI Energy LG max per board;Counts",
            nbins, xmin, xmax),
            fersboard.get_energy_max_name(gain='LG', isCer=False)
        )
        hists_board_cer_hg_max.append(hist_board_cer_hg_max)
        hists_board_cer_lg_max.append(hist_board_cer_lg_max)
        hists_board_sci_hg_max.append(hist_board_sci_hg_max)
        hists_board_sci_lg_max.append(hist_board_sci_lg_max)

    hist_cer_hg_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='HG', isCer=True)}",
        "FERS - CER Energy HG max;CER Energy HG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='HG', isCer=True)
    )
    hist_cer_lg_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='LG', isCer=True)}",
        "FERS - CER Energy LG max;CER Energy LG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='LG', isCer=True)
    )
    hist_sci_hg_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='HG', isCer=False)}",
        "FERS - SCI Energy HG max;SCI Energy HG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='HG', isCer=False)
    )
    hist_sci_lg_max = rdf.Histo1D((
        f"hist_{fersboards.get_energy_max_name(gain='LG', isCer=False)}",
        "FERS - SCI Energy LG max;SCI Energy LG max;Counts",
        nbins, xmin, xmax),
        fersboards.get_energy_max_name(gain='LG', isCer=False)
    )
    return hists_board_cer_hg_max + hists_board_cer_lg_max + hists_board_sci_hg_max + hists_board_sci_lg_max + [hist_cer_hg_max, hist_cer_lg_max, hist_sci_hg_max, hist_sci_lg_max]


def track_fers_hists():
    hists2d_fers_vs_event = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hist = rdf.Histo2D((
                    f"hist_FERS_Board{board_no}_{var}_VS_Event_{s_tower_x}_{s_tower_y}",
                    f"FERS Board {board_no} - Event VS {var} {chan.channel_no} in i_tower_x {s_tower_x} i_tower_y {s_tower_y};Event;{var} Energy HG",
                    n_bins_event, 0, n_events, 1000, 0, 9000),
                    "event_n", chan.get_channel_name(gain="HG")
                )
                hists2d_fers_vs_event.append(hist)
    return hists2d_fers_vs_event


def make_fers_2d_hists():
    hists2d_fers = []
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)
            chan_cer = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            chan_sci = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=False)

            i_cer = chan_cer.channel_no
            i_sci = chan_sci.channel_no
            # high gain
            hist = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{s_tower_x}_{s_tower_y}",
                f"CER {i_cer} VS SCI {i_sci} in i_tower_x {s_tower_x} i_tower_y {s_tower_y};Sci Energy HG;Cer Energy HG",
                300, 0, 9000, 300, 0, 9000),
                chan_sci.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="HG")
            )
            hist_zoomed = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{s_tower_x}_{s_tower_y}_zoom",
                f"CER {i_cer} VS SCI {i_sci} in i_tower_x {s_tower_x} i_tower_y {s_tower_y} (zoomed);SCI Energy HG;CER Energy HG",
                300, 0, 1000, 200, 0, 2000),
                chan_sci.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="HG")
            )
            hists2d_fers.append(hist)
            hists2d_fers.append(hist_zoomed)

            # low gain (Y) VS high gain (X) for Sci
            hist_sci_lg_vs_hg = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Sci_{s_tower_x}_{s_tower_y}_hg_VS_lg",
                f"SCI {i_sci} LG VS HG;SCI Energy HG;SCI Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_sci.get_channel_name(gain="HG"),
                chan_sci.get_channel_name(gain="LG")
            )
            hists2d_fers.append(hist_sci_lg_vs_hg)
            hist_cer_lg_vs_hg = rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_{s_tower_x}_{s_tower_y}_hg_VS_lg",
                f"CER {i_cer} LG VS HG;CER Energy HG;CER Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_cer.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="LG")
            )
            hists2d_fers.append(hist_cer_lg_vs_hg)
    return hists2d_fers


def make_drs_2d_hists():
    hists2d_DRS_VS_TS = []
    for _, DRSBoard in DRSBoards.items():
        for chan in DRSBoard:
            if chan is None:
                continue
            channelName = chan.get_channel_name(blsub=True)
            ymin, ymax = get_drs_plot_ranges(
                subtractMedian=True, is_amplified=chan.is_amplified, is6mm=chan.is6mm, is_reference=chan.is_reference)
            hist_subtractMedian = rdf.Histo2D((
                f"hist_{channelName}_VS_TS",
                "DRS values (subtract baseline);TS;DRS values",
                1024, 0, 1024, 50, ymin, ymax),
                "TS", channelName
            )
            hists2d_DRS_VS_TS.append(hist_subtractMedian)
            hprof_subtractMedian = rdf.Profile1D((
                f"prof_{channelName}_VS_TS",
                "DRS values (subtract baseline) VS TS;TS;Mean DRS values",
                1024, 0, 1024),
                "TS", channelName
            )
            hists2d_DRS_VS_TS.append(hprof_subtractMedian)
            hist_subtractMedian_TSCalib = rdf.Histo2D((
                f"hist_{channelName}_VS_AlignedTS",
                "DRS values (subtract baseline) VS AlignedTS;AlignedTS;DRS values",
                1024, 0, 1024, 50, ymin, ymax),
                f"{channelName}_AlignedTS", channelName
            )
            hists2d_DRS_VS_TS.append(hist_subtractMedian_TSCalib)
            hprof_subtractMedian_TSCalib = rdf.Profile1D((
                f"prof_{channelName}_VS_AlignedTS",
                "DRS values (subtract baseline) VS AlignedTS;AlignedTS;Mean DRS values",
                1024, 0, 1024),
                f"{channelName}_AlignedTS", channelName
            )
            hists2d_DRS_VS_TS.append(hprof_subtractMedian_TSCalib)
            hist_subtractMedian_CFDTSCalib = rdf.Histo2D((
                f"hist_{channelName}_VS_CFDAlignedTS",
                "DRS values (subtract baseline) VS CFD AlignedTS;CFD AlignedTS;DRS values",
                1024, 0, 1024, 50, ymin, ymax),
                f"{channelName}_cfdalignedts", channelName
            )
            hists2d_DRS_VS_TS.append(hist_subtractMedian_CFDTSCalib)
            hprof_subtractMedian_CFDTSCalib = rdf.Profile1D((
                f"prof_{channelName}_VS_CFDAlignedTS",
                "DRS values (subtract baseline) VS CFD AlignedTS;CFD AlignedTS;Mean DRS values",
                1024, 0, 1024),
                f"{channelName}_cfdalignedts", channelName
            )
            hists2d_DRS_VS_TS.append(hprof_subtractMedian_CFDTSCalib)
            hist_subtractMedian_CFDTSCalib_MCP = rdf.Histo2D((
                f"hist_{channelName}_VS_CFDAlignedTS_MCP",
                "DRS values (subtract baseline) VS CFD AlignedTS (calibrated to MCP);CFD AlignedTS calibrated to MCP;DRS values",
                1024, 0, 1024, 50, ymin, ymax),
                f"{channelName}_cfdalignedts_mcp", channelName
            )
            hists2d_DRS_VS_TS.append(hist_subtractMedian_CFDTSCalib_MCP)
            hprof_subtractMedian_CFDTSCalib_MCP = rdf.Profile1D((
                f"prof_{channelName}_VS_CFDAlignedTS_MCP",
                "DRS values (subtract baseline) VS CFD AlignedTS (calibrated to MCP);CFD AlignedTS calibrated to MCP;Mean DRS values",
                1024, 0, 1024),
                f"{channelName}_cfdalignedts_mcp", channelName
            )
            hists2d_DRS_VS_TS.append(hprof_subtractMedian_CFDTSCalib_MCP)

    return hists2d_DRS_VS_TS


def make_drs_service_2d_hists(channel_list):
    hists2d_DRS_service = []
    for channelName in channel_list:
        hist_subtractMedian_TSCalib = rdf.Histo2D((
            f"hist_{channelName}_VS_TS",
            "DRS values (subtract baseline) VS AlignedTS;AlignedTS;DRS values",
            1024, 0, 1024, 50, -1500, 500),
            f"TS", channelName
        )
        hists2d_DRS_service.append(hist_subtractMedian_TSCalib)
    return hists2d_DRS_service


def make_drs_stats_hists():
    # plot DRS sum and peak values
    hists_DRSStats = []
    for _, DRSBoard in DRSBoards.items():
        for chan in DRSBoard:
            channelSumName = chan.get_channel_sum_name()
            channelPeakName = chan.get_channel_peak_name()
            channelName_blsub = chan.get_channel_name(blsub=True)

            hist_DRS_Peak = rdf.Histo1D((
                f"hist_{channelPeakName}",
                "DRS Peak;Peak;Counts",
                200, 0, 800),
                channelPeakName
            )
            hists_DRSStats.append(hist_DRS_Peak)

            # two ways of getting the DRS sum: bare integral and cfd integral
            hist_DRS_Sum = rdf.Histo1D((
                f"hist_{channelSumName}",
                "DRS Sum; DRS Sum;Counts",
                650, -500, 6000),
                channelSumName
            )
            hists_DRSStats.append(hist_DRS_Sum)
            hist_DRS_Sum_cfd = rdf.Histo1D((
                f"hist_{channelName_blsub}_cfd",
                "DRS CFD Integral;DRS CFD Integral;Counts",
                650, -500, 6000),
                f"{channelName_blsub}_cfd"
            )
            hists_DRSStats.append(hist_DRS_Sum_cfd)

            # two ways of getting the DRS peak time: peak and cfd
            hists_DRS_RelPeakTS = rdf.Histo1D((
                f"hist_{channelName_blsub}_RelPeakTS",
                "DRS Relative Peak TS;Relative Peak TS (TS);Counts",
                1024, 0, 1024),
                f"{channelName_blsub}_RelPeakTS"
            )
            hists_DRSStats.append(hists_DRS_RelPeakTS)
            hists_DRS_RelCFDTS = rdf.Histo1D((
                f"hist_{channelName_blsub}_RelCFDTS",
                "DRS Relative CFD TS;Relative CFD TS (TS);Counts",
                1024, 0, 1024),
                f"{channelName_blsub}_cfdrelts"
            )
            hists_DRSStats.append(hists_DRS_RelCFDTS)

            # peak ts with respect to the MCP signal
            hists_DRS_RelPeakTS_MCP = rdf.Histo1D((
                f"hist_{channelName_blsub}_RelPeakTS_MCP",
                "DRS Relative Peak TS (calibrated to MCP);Relative Peak TS (TS);Counts",
                1024, 0, 1024),
                f"{channelName_blsub}_RelPeakTS_MCP"
            )
            hists_DRSStats.append(hists_DRS_RelPeakTS_MCP)
            hists_DRS_RelCFDTS_MCP = rdf.Histo1D((
                f"hist_{channelName_blsub}_RelCFDTS_MCP",
                "DRS Relative CFD TS (calibrated to MCP);Relative CFD TS (TS);Counts",
                1024, 0, 1024),
                f"{channelName_blsub}_cfdrelts_mcp"
            )
            hists_DRSStats.append(hists_DRS_RelCFDTS_MCP)

    return hists_DRSStats


def check_drs_sum_vs_fers():
    xymax = {
        "Cer": (20000, 8500),
        "Sci": (30000, 8500)
    }
    xymax_LG = {
        "Cer": (20000, 2000),
        "Sci": (30000, 4000)
    }

    # correlate DRS sum and FERS outputs
    h2s_drs_sum_vs_fers = []
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan_drs = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_drs is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}")
                    continue
                chan_fers = None
                for fersboard in fersboards.values():
                    chan_fers = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_fers is not None:
                        break
                if chan_fers is None:
                    print(
                        f"Warning: FERS Channel not found for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}")
                    continue

                h2_drs_sum_vs_fers = rdf.Histo2D((
                    f"hist_DRSSum_VS_FERS_{var}_{s_tower_x}_{s_tower_y}",
                    f"DRS sum VS FERS energy correlation for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}; FERS Energy HG; DRS Sum",
                    100, 0, xymax[var][1], 100, 0, xymax[var][0]
                ),
                    chan_fers.get_channel_name(gain="HG"),
                    chan_drs.get_channel_sum_name(),
                )
                h2s_drs_sum_vs_fers.append(h2_drs_sum_vs_fers)

                h2_drs_sum_vs_fers_lg = rdf.Histo2D((
                    f"hist_DRSSum_VS_FERSLG_{var}_{s_tower_x}_{s_tower_y}",
                    f"DRS sum VS FERS LG energy correlation for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}",
                    100, 0, xymax_LG[var][1], 100, 0, xymax_LG[var][0]
                ),
                    chan_fers.get_channel_name(gain="LG"),
                    chan_drs.get_channel_sum_name(),
                )
                h2s_drs_sum_vs_fers.append(h2_drs_sum_vs_fers_lg)

    return h2s_drs_sum_vs_fers


def check_drs_peak_vs_fers():
    h2s_drs_peak_vs_fers = []
    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                chan_drs = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_drs is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}")
                    continue
                chan_fers = None
                for fersboard in fersboards.values():
                    chan_fers = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_fers is not None:
                        break
                if chan_fers is None:
                    print(
                        f"Warning: FERS Channel not found for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}")
                    continue

                _, ymax = get_drs_plot_ranges(
                    subtractMedian=True, is_amplified=chan_drs.is_amplified, is6mm=chan_drs.is6mm)
                h2_drs_peak_vs_fers = rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERS_{var}_{s_tower_x}_{s_tower_y}",
                    f"DRS Peak VS FERS energy correlation for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}; FERS Energy HG; DRS Peak",
                    100, 0, 9000, 100, 0, ymax
                ),
                    chan_fers.get_channel_name(gain="HG"),
                    chan_drs.get_channel_peak_name(),
                )
                h2s_drs_peak_vs_fers.append(h2_drs_peak_vs_fers)

                h2_drs_peak_vs_fers_lg = rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERSLG_{var}_{s_tower_x}_{s_tower_y}",
                    f"DRS Peak VS FERS LG energy correlation for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}; FERS Energy LG; DRS Peak",
                    100, 0, 3000, 100, 0, ymax
                ),
                    chan_fers.get_channel_name(gain="LG"),
                    chan_drs.get_channel_peak_name(),
                )
                h2s_drs_peak_vs_fers.append(h2_drs_peak_vs_fers_lg)

    return h2s_drs_peak_vs_fers


def check_drs_peak_ts():
    h1s_drs_peak_ts = {}
    h1s_drs_peak_ts["Cer"] = []
    h1s_drs_peak_ts["Sci"] = []
    h2s_drs_peak_ts_cer_vs_sci = []

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            s_tower_x = number_to_string(i_tower_x)
            s_tower_y = number_to_string(i_tower_y)

            channel_names = {}
            for var in ["Cer", "Sci"]:
                chan_drs = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_drs is None:
                    if var == "Cer":
                        print(
                            f"Warning: DRS Channel not found for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var}")
                    continue

                channel_name = chan_drs.get_channel_peak_ts_name()
                channel_names[var] = channel_name

                h1_drs_peak_ts = rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{s_tower_x}_{s_tower_y}",
                    f"DRS Peak TS for Board{board_no}, Tower({s_tower_x}, {s_tower_y}), {var};Peak TS;Counts",
                    1000, 0, 1000),
                    f"{channel_name}_good"
                )
                h1s_drs_peak_ts[var].append(h1_drs_peak_ts)

            if len(channel_names) < 2:
                # print(
                #    f"Warning: Not enough channels found for Board{boardNo}, Tower({s_tower_x}, {s_tower_y})")
                continue

            h2_drs_peak_cer_vs_sci = rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{s_tower_x}_{s_tower_y}",
                f"DRS Peak TS - CER VS SCI for Board{board_no}, Tower({s_tower_x}, {s_tower_y});SCI Peak TS;CER Peak TS",
                1000, 0, 1000, 1000, 0, 1000),
                f'{channel_names["Sci"]}_good',
                f'{channel_names["Cer"]}_good',
            )
            h2s_drs_peak_ts_cer_vs_sci.append(h2_drs_peak_cer_vs_sci)
    return h1s_drs_peak_ts["Cer"], h1s_drs_peak_ts["Sci"], h2s_drs_peak_ts_cer_vs_sci


def main():
    paths = get_run_paths(run_number)

    hprofs_conditions = monitor_conditions()
    hprofs_energysum = monitor_fers_energy_sum()

    hists1d_fers = make_fers_1d_hists()
    hists1d_fers_lg = make_fers_1d_hists(gain="LG")

    # hists2d_fers = make_fers_2d_hists()
    # hists2d_fers_vs_event = track_fers_hists()

    hists2d_drs_vs_ts = make_drs_2d_hists()
    hists_drs_stats = make_drs_stats_hists()

    map_mcp_channels = get_mcp_channels(run_number)
    list_mcp_channels = [
        channel + "_blsub" for channels in map_mcp_channels.values() for channel in channels]
    hists2d_drs_service = make_drs_service_2d_hists(list_mcp_channels)

    # hists2d_DRSPeak_VS_FERS = check_drs_peak_vs_fers()

    # hists1d_DRSPeakTS_Cer, hists1d_DRSPeakTS_Sci, hists2d_DRSPeakTS_Cer_VS_Sci = check_drs_peak_ts()

    # hists2d_DRSSum_VS_FERS = check_drs_sum_vs_fers()

    hists_fers_max = make_fers_max_value_hists()

    stats = collect_fers_stats()

    # pedestals
    # event loop should be triggered here
    pedestals_hg = collect_fers_pedestals(hists1d_fers, gain="HG")
    pedestals_lg = collect_fers_pedestals(hists1d_fers_lg, gain="LG")

    # print(f"Loops executed after pedestals: {rdf.GetNRuns()}")

    print("\033[94mSave results\033[0m")

    # dump stats into a json file
    stats_results = {}
    for channel_name, (mean, max_value, sat_frequency) in stats.items():
        stats_results[channel_name] = (mean.GetValue(), max_value.GetValue(), float(
            sat_frequency.GetValue()) / n_events)
    with open(f"{paths['root']}/fers_stats.json", "w") as f:
        json.dump(stats_results, f, indent=4)

    with open(f"{paths['root']}/fers_pedestals_hg.json", "w") as f:
        json.dump(pedestals_hg, f, indent=4)
    with open(f"{paths['root']}/fers_pedestals_lg.json", "w") as f:
        json.dump(pedestals_lg, f, indent=4)

    # Save histograms to an output ROOT file
    save_hists_to_file(hprofs_conditions,
                       f"{paths['root']}/conditions_vs_event.root")

    save_hists_to_file(hprofs_energysum,
                       f"{paths['root']}/fers_energysum_vs_event.root")

    save_hists_to_file(hists1d_fers,
                       f"{paths['root']}/fers_all_channels_1d.root")

    save_hists_to_file(hists2d_drs_vs_ts, f"{paths['root']}/drs_vs_ts.root")

    save_hists_to_file(hists_drs_stats, f"{paths['root']}/drs_stats.root")

    save_hists_to_file(hists2d_drs_service,
                       f"{paths['root']}/drs_service_channels.root")

    # save_hists_to_file(hists2d_DRSPeak_VS_FERS,
    #                   f"{paths['root']}/drspeak_vs_fers.root")

    # save_hists_to_file(hists1d_DRSPeakTS_Cer + hists1d_DRSPeakTS_Sci +
    #                   hists2d_DRSPeakTS_Cer_VS_Sci, f"{paths['root']}/drspeakts.root")

    # save_hists_to_file(hists2d_DRSSum_VS_FERS,
    #                   f"{paths['root']}/drssum_vs_fers.root")

    save_hists_to_file(hists_fers_max, f"{paths['root']}/fers_max_values.root")


if __name__ == "__main__":
    main()
