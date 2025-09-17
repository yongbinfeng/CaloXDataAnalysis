# collect all functions related to DRS here
from utils.utils import getBranchStats
from channels.channel_map import findFanoutTimeReferenceDelay, findDRSTriggerMap

def preProcessDRSBoards(rdf, debug=False):
    import re
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"^DRS_Board\d+_Group\d+_Channel\d+$")
    drs_branches = [b for b in branches if pattern.search(b)]
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

    # Create an array of indices for DRS outputs
    rdf = rdf.Define("TS", "FillIndices(1024)")

    # find the baseline of each DRS channel (median here)
    # and subtract it from the DRS outputs
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_bl",
            f"compute_median({varname})"
        )
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
            triggerDelayTS = int((triggerDelay / 0.2) + 800.0) # convert to time samples

            rdf = rdf.Define(
                channelPeakTSName,
                f"ArgMaxRange({channelName_sub}, {TS_start}, {TS_end}, {threshold})"
            )
            if(not rdf.HasColumn(triggerPeakTSName)):
                rdf = rdf.Define(
                    triggerPeakTSName,
                    f"ArgMinRange({triggerName}_blsub, 500, 1024)"
                )
            rdf = rdf.Define(
                f"{channelPeakTSName}_good",
                f"({channelPeakTSName} - {triggerPeakTSName}) + {triggerDelayTS}"
            )
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
