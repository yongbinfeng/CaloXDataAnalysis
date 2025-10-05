# collect all functions related to DRS here
import os
import re
import json
from utils.utils import getBranchStats
from channels.channel_map import findFanoutTimeReferenceDelay, findDRSTriggerMap, getMCPChannels


def preProcessDRSBoards(rdf, debug=False, drsboards=None, runNumber=None):
    import re
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"^DRS_Board\d+_Group\d+_Channel\d+$")
    drs_branches = [b for b in branches if pattern.search(b)]
    # check if the drs stats file already exists
    stats_file = f"results/root/Run{runNumber}/drs_stat_preprocessed.json"
    if not os.path.exists(stats_file):
        # get the statistics of DRS branches
        print(f"Calculating DRS stats and saving to {stats_file}")
        stats = getBranchStats(rdf, drs_branches)

        print("DRS branches statistics:")
        for br, res in stats.items():
            print(f"{br}: mean = {res['mean'].GetValue():.4f}, "
                  f"min = {res['min'].GetValue():.4f}, "
                  f"max = {res['max'].GetValue():.4f}")
            stats[br] = {
                "mean": res['mean'].GetValue(),
                "min": res['min'].GetValue(),
                "max": res['max'].GetValue()
            }
        # save the stats to a json file
        os.makedirs(os.path.dirname(stats_file), exist_ok=True)
        with open(stats_file, "w") as f:
            json.dump(stats, f, indent=4)

    else:
        print(f"Loading DRS stats from {stats_file}")
        with open(stats_file, "r") as f:
            stats = json.load(f)

    # Create an array of indices for DRS outputs
    rdf = rdf.Define("TS", "FillIndices(1024)")

    # get the list of 6mm amplified channels
    drs_amplified_channels = []
    if drsboards is not None:
        for _, DRSBoard in drsboards.items():
            for channel in DRSBoard:
                if channel.is6mm and channel.isAmplified:
                    drs_amplified_channels.append(
                        channel.GetChannelName(blsub=False))
        # print("6mm amplified channels:", drs_amplified_channels)

    # find the baseline of each DRS channel (median here)
    # and subtract it from the DRS outputs
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_bl",
            f"compute_median({varname})"
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
    if debug:
        # define relative TS with respect to the StartIndexCell
        for varrname in drs_branches:
            # replace the string "Channel[0-9]+" with "StartIndexCell"
            var_StartIndexCell = re.sub(
                r"Channel[0-9]+", "StartIndexCell", varrname)
            rdf = rdf.Define(
                f"RTS_pos_{varrname}",
                f"(TS + {var_StartIndexCell}) % 1024"
            )
            rdf = rdf.Define(
                f"RTS_neg_{varrname}",
                f"((TS - {var_StartIndexCell}) % 1024 + 1024) % 1024"
            )

    return rdf


def getDRSSum(rdf, DRSBoards, TS_start=0, TS_end=400):
    # get the mean of DRS outputs per channel
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            channelName_blsub = channel.GetChannelName(blsub=True)
            channelSumName = channel.GetChannelSumName()
            rdf = rdf.Define(
                channelSumName,
                f"SumRange({channelName_blsub}, {TS_start}, {TS_end})"
            )
            # rdf = rdf.Define(
            #    f"{varname}_subtractMedian_positive",
            #    f"clipToZero({varname}_subtractMedian)"
            # )
            # rdf = rdf.Define(
            #    f"{varname}_sum",
            #    f"SumRange({varname}_subtractMedian_positive, {TS_start}, {TS_end})"
            # )
    return rdf


def getDRSPeakTS(rdf, run, DRSBoards, TS_start=0, TS_end=400, threshold=1.0):
    # get the peak TS of DRS outputs per channel
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            channelName_sub = channel.GetChannelName(blsub=True)
            channelPeakTSName = channel.GetChannelPeakTSName()

            # find the trigger channel for this DRS channel
            channelName = channel.GetChannelName(blsub=False)
            triggerName = findDRSTriggerMap(channelName, run=run)
            triggerPeakTSName = triggerName + "_peakTS"
            triggerDelay = findFanoutTimeReferenceDelay(triggerName, run=run)
            # convert to time samples
            triggerDelayTS = int((triggerDelay / 0.2) + 800.0)

            rdf = rdf.Define(
                channelPeakTSName,
                f"ArgMaxRange({channelName_sub}, {TS_start}, {TS_end}, {threshold})"
            )
            if (not rdf.HasColumn(triggerPeakTSName)):
                rdf = rdf.Define(
                    triggerPeakTSName,
                    f"ArgMinRange({triggerName}_blsub, 500, 1024)"
                )
            rdf = rdf.Define(
                f"{channelPeakTSName}_good",
                f"({channelPeakTSName} - {triggerPeakTSName}) + {triggerDelayTS}"
            )
    return rdf


def calibrateDRSPeakTS(rdf, run, DRSBoards, TSminMCP=500, TSmaxMCP=600, TSminDRS=0, TSmaxDRS=400, threshold=1.0):
    map_mcp_channels = getMCPChannels(run)

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
        if DRSBoard.boardNo > 3:
            continue
        for channel in DRSBoard:
            channelName = channel.GetChannelName(blsub=False)

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
                f"{channelName}_DiffRelPeakTS_US", f"(int){channelName}_RelPeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[DRSBoard.boardNo]}")

            # align TS
            rdf = rdf.Define(
                f"{channelName}_AlignedTS", f"TS - (int){channel_TS}_PeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[DRSBoard.boardNo]}"
            )

            # map aligned TS to real distance in z (cm)
            t0TS = 116.5 + 2
            if channel.isQuartz:
                refrac = 1.468
            else:
                refrac = 1.621
                refrac = 1.58
            rdf = rdf.Define(
                f"{channelName}_MeasuredZ", f"{refrac} / ({refrac} - 1) * 250 - 6.0 / ({refrac} - 1) * ({channelName}_AlignedTS + {t0TS})")

            # channel sum
            rdf = rdf.Define(
                f"{channelName}_Sum", f"SumRange({channelName}_blsub, {TSminDRS}, {TSmaxDRS})")

    return rdf


def getDRSPeak(rdf, DRSBoards, TS_start=0, TS_end=400):
    """
    Get the peak value of DRS outputs per channel.
    """
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            channelName_sub = channel.GetChannelName(blsub=True)
            channelPeakName = channel.GetChannelPeakName()
            rdf = rdf.Define(
                channelPeakName,
                f"MaxRange({channelName_sub}, {TS_start}, {TS_end})"
            )
    return rdf


def getDRSStats(rdf, run, DRSBoards, TS_start=0, TS_end=400, threshold=1.0):
    """
    Get the statistics of DRS outputs per channel.
    """
    rdf = getDRSSum(rdf, DRSBoards, TS_start, TS_end)
    rdf = getDRSPeakTS(rdf, run, DRSBoards, TS_start, TS_end, threshold)
    rdf = getDRSPeak(rdf, DRSBoards, TS_start, TS_end)
    return rdf
