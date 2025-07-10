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


def getDataFile(runNumber):
    runNum = str(runNumber)
    import json
    jsonFile = "data/datafiles.json"
    with open(jsonFile, 'r') as f:
        data = json.load(f)
    if runNum in data:
        return data[runNum]
    else:
        raise ValueError(f"Run number {runNum} not found in datafiles.json")


def getBranchStats(rdf, branches):
    stats = {
        br: {
            "mean": rdf.Mean(br),
            "min": rdf.Min(br),
            "max": rdf.Max(br)
        } for br in branches
    }
    return stats


def vectorizeFERS(rdf, FERSBoards):
    # FRES board outputs
    # define variables as RDF does not support reading vectors
    # with indices directly
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for channel in FERSBoard:
            rdf = rdf.Define(
                f"FERS_Board{boardNo}_energyHG_{channel.channelNo}",
                f"FERS_Board{boardNo}_energyHG[{channel.channelNo}]")
            rdf = rdf.Define(
                f"FERS_Board{boardNo}_energyLG_{channel.channelNo}",
                f"FERS_Board{boardNo}_energyLG[{channel.channelNo}]"
            )
    return rdf


def processDRSBoards(rdf, debug=False):
    import re
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"DRS.*Group.*Channel.*")
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

    # get the mean of DRS outputs per channel
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_median",
            f"compute_median({varname})"
        )
        rdf = rdf.Define(
            f"{varname}_subtractMedian",
            f"{varname} - {varname}_median"
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


def calculateEnergySumFERS(rdf, FERSBoards):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    boardNos = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        boardNos.append(boardNo)

        channels_Cer = FERSBoard.GetCerChannels()
        channels_Sci = FERSBoard.GetSciChannels()

        string_CerEnergyHG = "+".join(
            chan.GetHGChannelName() for chan in channels_Cer
        )
        string_CerEnergyLG = "+".join(
            chan.GetLGChannelName() for chan in channels_Cer
        )
        string_SciEnergyHG = "+".join(
            chan.GetHGChannelName() for chan in channels_Sci
        )
        string_SciEnergyLG = "+".join(
            chan.GetLGChannelName() for chan in channels_Sci
        )
        # per-board energy sum
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_CerEnergyHG",
            f"({string_CerEnergyHG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_CerEnergyLG",
            f"({string_CerEnergyLG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_SciEnergyHG",
            f"({string_SciEnergyHG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_SciEnergyLG",
            f"({string_SciEnergyLG})")
    # per-event energy sum
    string_CerEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_CerEnergyHG" for boardNo in boardNos
    )
    string_CerEnergyLG_Total = "+".join(
        f"FERS_Board{boardNo}_CerEnergyLG" for boardNo in boardNos
    )
    string_SciEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_SciEnergyHG" for boardNo in boardNos
    )
    string_SciEnergyLG_Total = "+".join(
        f"FERS_Board{boardNo}_SciEnergyLG" for boardNo in boardNos
    )
    rdf = rdf.Define("FERS_CerEnergyHG", f"({string_CerEnergyHG_Total})")
    rdf = rdf.Define("FERS_CerEnergyLG", f"({string_CerEnergyLG_Total})")
    rdf = rdf.Define("FERS_SciEnergyHG", f"({string_SciEnergyHG_Total})")
    rdf = rdf.Define("FERS_SciEnergyLG", f"({string_SciEnergyLG_Total})")

    return rdf


def getDRSSum(rdf, DRSBoards, TS_start=0, TS_end=400):
    # get the mean of DRS outputs per channel
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            varname = channel.GetChannelName()
            rdf = rdf.Define(
                f"{varname}_subtractMedian_positive",
                f"clipToZero({varname}_subtractMedian)"
            )
            rdf = rdf.Define(
                f"{varname}_sum",
                f"SumRange({varname}_subtractMedian_positive, {TS_start}, {TS_end})"
            )
    return rdf


def getDRSPeakTS(rdf, DRSBoards, TS_start=0, TS_end=400, threshold=1.0):
    # get the peak TS of DRS outputs per channel
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            varname = channel.GetChannelName()
            rdf = rdf.Define(
                f"{varname}_peakTS",
                f"ArgMaxRange({varname}_subtractMedian_positive, {TS_start}, {TS_end}, {threshold})"
            )
    return rdf


def loadRDF(runNumber, firstEvent=0, lastEvent=-1):
    import ROOT
    # Open the input ROOT file
    ifile = getDataFile(runNumber)
    infile = ROOT.TFile(ifile, "READ")
    rdf_org = ROOT.RDataFrame("EventTree", infile)
    nevents = rdf_org.Count().GetValue()
    if lastEvent < 0 or lastEvent > nevents:
        lastEvent = nevents
    print(
        f"\033[94mfiltering events from {firstEvent} to {lastEvent} in run {runNumber}, total {nevents} events in the file.\033[0m")
    # Apply the event range filter
    rdf = rdf_org.Filter(f"event_n >= {firstEvent} && event_n < {lastEvent}")

    return rdf, rdf_org


def filterPrefireEvents(rdf, TS=350):
    # use the hodo trigger to filter prefire events
    from utils.channel_map import buildHodoTriggerChannels
    trigger_name_top, trigger_name_bottom = buildHodoTriggerChannels()
    # index of the minimum value in the trigger channels
    rdf = rdf.Define(
        "TS_fired_up", f"ROOT::VecOps::ArgMin({trigger_name_top})")
    rdf = rdf.Define(
        "TS_fired_down", f"ROOT::VecOps::ArgMin({trigger_name_bottom})")

    rdf = rdf.Define(
        "NormalFired", f"(TS_fired_up >= {TS}) && (TS_fired_down >= {TS})")

    rdf_prefilter = rdf
    rdf = rdf.Filter("NormalFired == 1")

    return rdf, rdf_prefilter
