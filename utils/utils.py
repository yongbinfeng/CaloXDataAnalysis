from utils.channel_map import findDRSTriggerMap, findTimeReferenceDelay

def number2string(n):
    s = str(n)
    return s.replace('-', 'm').replace('.', 'p')


def string2number(s):
    return float(s.replace('m', '-').replace('p', '.'))


def round_up_to_1eN(x):
    import math
    """
    Round a number up to the nearest 10^N.
    """
    if x <= 0:
        return 0
    return 10 ** math.ceil(math.log10(x))


def getBranchStats(rdf, branches):
    stats = {
        br: {
            "mean": rdf.Mean(br),
            "min": rdf.Min(br),
            "max": rdf.Max(br)
        } for br in branches
    }
    return stats

def preProcessTimeCorrections(rdf, DRSBoards, runNumber):
    for _, DRSBoard in DRSBoards.items():
        for chan in DRSBoard:
            channelName = chan.GetChannelName()
            triggerName = findDRSTriggerMap(channelName, run=runNumber)
            triggerDelay = findTimeReferenceDelay(triggerName, run=runNumber)
            channelTimingName = chan.GetChannelTimeName()
            triggerTimingName = f"{triggerName}_LP2_50"
            # rdf = rdf.Define(f"{channelTimingName}_goodTime", f"{channelTimingName} - {triggerTimingName} + {triggerDelay}")
            rdf = rdf.Define(f"{channelTimingName}_goodTime", f"{channelTimingName} - {triggerTimingName} + {triggerDelay} - (DRS_Board4_Group3_Channel6_LP2_50 - DRS_Board4_Group3_Channel8_LP2_50 - 1.735) + 35.0")

    rdf = rdf.Filter(f"DRS_Board4_Group3_Channel6_LP2_50 > 1")

    return rdf
