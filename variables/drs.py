# collect all functions related to DRS here
import re
from channels.channel_map import findFanoutTimeReferenceDelay, findDRSTriggerMap, get_mcp_channels, get_service_drs_channels

TS_END = 1024


def preProcessDRSBoards(rdf, drsboards=None):
    import re
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"^DRS_Board\d+_Group\d+_Channel\d+$")
    drs_branches = [b for b in branches if pattern.search(b)]

    # Create an array of indices for DRS outputs
    rdf = rdf.Define("TS", "FillIndices(1024)")

    drs_amplified_channels = []
    if drsboards is not None:
        for _, DRSBoard in drsboards.items():
            for channel in DRSBoard:
                if channel.is6mm and channel.is_amplified:
                    drs_amplified_channels.append(
                        channel.get_channel_name(blsub=False))
        # print("6mm amplified channels:", drs_amplified_channels)

    # find the baseline of each DRS channel (median here)
    # and subtract it from the DRS outputs
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_bl",
            f"compute_baseline_median({varname}, 0, 200)"
        )
        if varname in drs_amplified_channels:
            # for 6mm amplified channels, flip the signal
            rdf = rdf.Define(
                f"{varname}_blsub",
                f"-({varname} - {varname}_bl)"
            )
        else:
            rdf = rdf.Define(
                f"{varname}_blsub",
                f"{varname} - {varname}_bl"
            )

    return rdf


def processDRSRefChannels(rdf, DRSBoards):
    # find the reference time of each DRS board group:
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            if not channel.is_reference:
                continue
            channelName_blsub = channel.get_channel_name(blsub=True)
            rdf = rdf.Define(
                f"{channelName_blsub}_REFHit",
                f"process_dynamic_led({channelName_blsub}, 0)"
            )
            rdf = rdf.Define(
                f"{channelName_blsub}_refts",
                f"{channelName_blsub}_REFHit.time_slice")

            # use the peak
            rdf = rdf.Define(
                f"{channelName_blsub}_PeakTS",
                f"ArgMinRange({channelName_blsub}, 0, 1024, -500.0)"
            )

    return rdf


def processDRSChannelsCFD(rdf, DRSBoards, TS_start=0, TS_end=TS_END):
    # get the mean of DRS outputs per channel
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            channelName_blsub = channel.get_channel_name(blsub=True)
            rdf = rdf.Define(
                f"{channelName_blsub}_CFD",
                f"compute_cfd_integral({channelName_blsub})"
            )
            rdf = rdf.Define(f"{channelName_blsub}_cfd",
                             f"{channelName_blsub}_CFD.energy")
            rdf = rdf.Define(f"{channelName_blsub}_cfdts",
                             f"{channelName_blsub}_CFD.time_slice")

            channel_TS = re.sub(
                r"_Channel[0-7]", "_Channel8", channelName_blsub)
            rdf = rdf.Define(
                f"{channelName_blsub}_cfdrelts", f"790 + (int){channelName_blsub}_cfdts - (int){channel_TS}_refts"
            )
            rdf = rdf.Define(
                f"{channelName_blsub}_cfdalignedts", f"790 + TS - (int){channel_TS}_refts - 0"
            )
    return rdf


def processDRSChannelsPeak(rdf, DRSBoards, TSminDRS=0, TSmaxDRS=TS_END, threshold=1.0):
    # jitter offset
    value_diffcorrs = [0, 0, 0, 0, 0, 0, 0]

    # calibration all DRS channels to MCP US channel 0
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            channelName_blsub = channel.get_channel_name(blsub=True)

            # get peak value
            channelPeakName = channel.get_channel_peak_name()

            rdf = rdf.Define(
                channelPeakName,
                f"MaxRange({channelName_blsub}, {TSminDRS}, {TSmaxDRS})"
            )

            # bare integral
            channelSumName = channel.get_channel_sum_name()
            rdf = rdf.Define(
                channelSumName,
                f"SumRange({channelName_blsub}, 0, -1)"
            )

            # define the relative peak TS with respect to the reference channel
            channel_TS = re.sub(
                r"_Channel[0-7]", "_Channel8", channelName_blsub)
            if not channel.is_reference:
                rdf = rdf.Define(
                    f"{channelName_blsub}_PeakTS",
                    f"ArgMaxRange({channelName_blsub}, {TSminDRS}, {TSmaxDRS}, {threshold})"
                )
            rdf = rdf.Define(
                f"{channelName_blsub}_RelPeakTS", f"790 + (int){channelName_blsub}_PeakTS - (int){channel_TS}_PeakTS")

            # align TS with respect to the trigger
            rdf = rdf.Define(
                f"{channelName_blsub}_AlignedTS", f"790 + TS - (int){channel_TS}_PeakTS - 0 - {value_diffcorrs[DRSBoard.board_no]}"
            )

    return rdf


def get_psd_energy_deposit(rdf, run_number, TS_start=100, TS_end=400):
    channel_name = get_service_drs_channels(run_number).get("PSD")
    rdf = rdf.Define(
        "PSD_Sum", f"SumRange({channel_name}_blsub, {TS_start}, {TS_end})")
    rdf = rdf.Define(
        "PSD_Energy_Cer", "PSD_Sum * 5.15852e-05 * (-1)").Define(
        "PSD_Energy_Sci", "PSD_Sum * 6.53095e-05 * (-1)")
    return rdf


def calibrateDRSPeakTS(rdf, run, DRSBoards, TSminMCP=500, TSmaxMCP=600, TSminDRS=0, TSmaxDRS=400, threshold=1.0):
    map_mcp_channels = get_mcp_channels(run)

    # get time reference channel TS
    channels_TF = [
        "DRS_Board0_Group0_Channel8",
        "DRS_Board1_Group0_Channel8",
        "DRS_Board2_Group0_Channel8",
        "DRS_Board3_Group0_Channel8",
        "DRS_Board4_Group0_Channel8",
        "DRS_Board5_Group0_Channel8",
        "DRS_Board6_Group0_Channel8",
        "DRS_Board0_Group1_Channel8",
        "DRS_Board1_Group1_Channel8",
        "DRS_Board2_Group1_Channel8",
        "DRS_Board3_Group1_Channel8",
        "DRS_Board4_Group1_Channel8",
        "DRS_Board5_Group1_Channel8",
        "DRS_Board6_Group1_Channel8",
        "DRS_Board0_Group2_Channel8",
        "DRS_Board1_Group2_Channel8",
        "DRS_Board2_Group2_Channel8",
        "DRS_Board3_Group2_Channel8",
        "DRS_Board4_Group2_Channel8",
        "DRS_Board5_Group2_Channel8",
        "DRS_Board6_Group2_Channel8",
        "DRS_Board0_Group3_Channel8",
        "DRS_Board1_Group3_Channel8",
        "DRS_Board2_Group3_Channel8",
        "DRS_Board3_Group3_Channel8",
        "DRS_Board4_Group3_Channel8",
        "DRS_Board5_Group3_Channel8",
        "DRS_Board6_Group3_Channel8",
    ]

    # jitter offset
    value_diffcorrs = [0, 0, 0, 0, 0, 0, 0]

    for channel in channels_TF:
        rdf = rdf.Define(f"{channel}_PeakTS",
                         f"ArgMinRange({channel}_blsub, 0, 1024, -70.0)")

    for det, channels in map_mcp_channels.items():
        for idx, channel in enumerate(channels):
            rdf = rdf.Define(f"{channel}_PeakTS",
                             f"ArgMinRange({channel}_blsub, {TSminMCP}, {TSmaxMCP}, -300.0)")
            rdf = rdf.Define(
                f"{channel}_Peak", f"MinRange({channel}_blsub, {TSminMCP}, {TSmaxMCP})")
            # define the relative peak TS with respect to the reference channel
            channel_TS = re.sub(r"_Channel[0-7]", "_Channel8", channel)
            rdf = rdf.Define(
                f"{channel}_RelPeakTS", f"(int){channel}_PeakTS - (int){channel_TS}_PeakTS")

    # calibration all DRS channels to MCP US channel 0
    for _, DRSBoard in DRSBoards.items():
        # if DRSBoard.board_no > 3:
        #    continue
        for channel in DRSBoard:
            channelName = channel.get_channel_name(blsub=False)

            # define the relative peak TS with respect to the reference channel
            channel_TS = re.sub(r"_Channel[0-7]", "_Channel8", channelName)
            rdf = rdf.Define(
                f"{channelName}_PeakTS",
                f"ArgMaxRange({channelName}_blsub, {TSminDRS}, {TSmaxDRS}, {threshold})"
            )
            rdf = rdf.Define(
                f"{channelName}_RelPeakTS", f"(int){channelName}_PeakTS - (int){channel_TS}_PeakTS")

            # define the difference of relative peak TS with respect to the first channel of US
            rdf = rdf.Define(
                f"{channelName}_DiffRelPeakTS_US", f"(int){channelName}_RelPeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[DRSBoard.board_no]}")

            # align TS
            rdf = rdf.Define(
                f"{channelName}_AlignedTS", f"TS - (int){channel_TS}_PeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[DRSBoard.board_no]}"
            )

            # map aligned TS to real distance in z (cm)
            # t0TS = 116.5 + 2 + 31
            # if channel.isQuartz:
            #    refrac = 1.468 * 1.468
            #    t0TS = 145
            # else:
            #    refrac = 1.621
            #    refrac = 1.58 * 1.58
            #    t0TS = 154
            # rdf = rdf.Define(
            #    f"{channelName}_MeasuredZ", f"{refrac} / ({refrac} - 1) * 250 - 6.0 / ({refrac} - 1) * ({channelName}_AlignedTS + {t0TS})")
            if channel.isQuartz:
                intercetp = 11.2
                slope = 0.016
            else:
                intercetp = 10.3
                slope = 0.018
            rdf = rdf.Define(f"{channelName}_MeasuredZ",
                             f"-1.0/{slope} * ({channelName}_AlignedTS * 0.2 + {intercetp})")

            # channel sum
            rdf = rdf.Define(
                f"{channelName}_Sum", f"SumRange({channelName}_blsub, {TSminDRS}, {TSmaxDRS})")

    return rdf


def get_drs_stats(rdf, run, DRSBoards, TS_start=0, TS_end=400, threshold=1.0):
    """
    Get the statistics of DRS outputs per channel.
    """
    rdf = processDRSRefChannels(rdf, DRSBoards)
    rdf = processDRSChannelsCFD(rdf, DRSBoards, TS_start, TS_end)
    rdf = processDRSChannelsPeak(rdf, DRSBoards, TS_start, TS_end)
    # rdf = getDRSSum(rdf, DRSBoards, TS_start, TS_end)
    # rdf = getDRSPeakTS(rdf, run, DRSBoards, TS_start, TS_end, threshold)
    # rdf = getDRSPeak(rdf, DRSBoards, TS_start, TS_end)
    return rdf
