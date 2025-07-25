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


def IsScanRun(runNumber):
    import json
    f_scanruns = "data/scanruns.json"
    with open(f_scanruns, 'r') as f:
        temp = json.load(f)
        scanruns = temp["scanruns"]
    return runNumber in scanruns


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


def subtractFERSPedestal(rdf, FERSBoards, pedestals):
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for channel in FERSBoard:
            channelNo = channel.channelNo
            channelNameHG = channel.GetHGChannelName()
            pedestal = pedestals[channelNameHG]
            # subtract pedestal from HG and LG energies
            rdf = rdf.Define(
                f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted",
                f"FERS_Board{boardNo}_energyHG_{channelNo} - {pedestal}"
            )
            # not on LG yet
            # rdf = rdf.Define(
            #    f"FERS_Board{boardNo}_energyLG_{channelNo}_subtracted",
            #    f"FERS_Board{boardNo}_energyLG_{channelNo} - 0."
            # )
    return rdf


def calibrateFERSChannels(rdf, FERSBoards, file_gains, file_pedestals):
    """
    Calibrate FERS channels using gains and pedestals from the provided files.
    """
    import json
    with open(file_gains, 'r') as f:
        gains = json.load(f)
    with open(file_pedestals, 'r') as f:
        pedestals = json.load(f)

    # Subtract pedestal and apply gain calibration
    rdf = subtractFERSPedestal(rdf, FERSBoards, pedestals)

    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for channel in FERSBoard:
            channelNo = channel.channelNo
            channelNameHG = channel.GetHGChannelName()
            gain = gains[channelNameHG]
            rdf = rdf.Define(
                f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calibrated",
                f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted / {gain}"
            )
            rdf = rdf.Define(
                f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calibrated_clipped",
                f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calibrated > 0.8 ? FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calibrated : 0"
            )
            # Not calibrating LG yet
            # rdf = rdf.Define(
            #    f"FERS_Board{boardNo}_energyLG_{channelNo}_calibrated",
            #    f"FERS_Board{boardNo}_energyLG_{channelNo}_subtracted"
            # )

    return rdf


def preProcessDRSBoards(rdf, debug=False):
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

    # find the baseline of each DRS channel
    # and subtract it from the DRS outputs
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


def calculateEnergySumFERS(rdf, FERSBoards, subtractPedestal=False, calibrate=False, clip=False):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    suffix = ""
    if subtractPedestal:
        suffix = "_subtracted"
    if calibrate:
        suffix += "_calibrated"
    if clip:
        suffix += "_clipped"
    boardNos = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        boardNos.append(boardNo)

        channels_Cer = FERSBoard.GetCerChannels()
        channels_Sci = FERSBoard.GetSciChannels()

        string_CerEnergyHG = "+".join(
            chan.GetHGChannelName() + suffix for chan in channels_Cer
        )
        # string_CerEnergyLG = "+".join(
        #    chan.GetLGChannelName() for chan in channels_Cer
        # )
        string_SciEnergyHG = "+".join(
            chan.GetHGChannelName() + suffix for chan in channels_Sci
        )
        # string_SciEnergyLG = "+".join(
        #    chan.GetLGChannelName() for chan in channels_Sci
        # )
        # per-board energy sum
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_CerEnergyHG" + suffix,
            f"({string_CerEnergyHG})")
        # rdf = rdf.Define(
        #    f"FERS_Board{boardNo}_CerEnergyLG" + suffix,
        #    f"({string_CerEnergyLG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_SciEnergyHG" + suffix,
            f"({string_SciEnergyHG})")
        # rdf = rdf.Define(
        #    f"FERS_Board{boardNo}_SciEnergyLG" + suffix,
        #    f"({string_SciEnergyLG})")
    # per-event energy sum
    string_CerEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_CerEnergyHG" + suffix for boardNo in boardNos
    )
    # string_CerEnergyLG_Total = "+".join(
    #    f"FERS_Board{boardNo}_CerEnergyLG" + suffix for boardNo in boardNos
    # )
    string_SciEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_SciEnergyHG" + suffix for boardNo in boardNos
    )
    # string_SciEnergyLG_Total = "+".join(
    #    f"FERS_Board{boardNo}_SciEnergyLG" + suffix for boardNo in boardNos
    # )
    rdf = rdf.Define("FERS_CerEnergyHG" + suffix,
                     f"({string_CerEnergyHG_Total})")
    # rdf = rdf.Define("FERS_CerEnergyLG" + suffix,
    #                 f"({string_CerEnergyLG_Total})")
    rdf = rdf.Define("FERS_SciEnergyHG" + suffix,
                     f"({string_SciEnergyHG_Total})")
    # rdf = rdf.Define("FERS_SciEnergyLG" + suffix,
    #                 f"({string_SciEnergyLG_Total})")

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


def getDRSPeak(rdf, DRSBoards, TS_start=0, TS_end=400):
    """
    Get the peak value of DRS outputs per channel.
    """
    TS_start = int(TS_start)
    TS_end = int(TS_end)
    for _, DRSBoard in DRSBoards.items():
        for channel in DRSBoard:
            varname = channel.GetChannelName()
            rdf = rdf.Define(
                f"{varname}_peak",
                f"MaxRange({varname}_subtractMedian_positive, {TS_start}, {TS_end})"
            )
    return rdf


def prepareDRSStats(rdf, DRSBoards, TS_start=0, TS_end=400, threshold=1.0):
    """
    Get the statistics of DRS outputs per channel.
    """
    rdf = getDRSSum(rdf, DRSBoards, TS_start, TS_end)
    rdf = getDRSPeakTS(rdf, DRSBoards, TS_start, TS_end, threshold)
    rdf = getDRSPeak(rdf, DRSBoards, TS_start, TS_end)
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


def filterPrefireEvents(rdf, runNumber, TS=350):
    # use the hodo trigger to filter prefire events
    from utils.channel_map import buildHodoTriggerChannels
    trigger_names = buildHodoTriggerChannels(runNumber)
    if not trigger_names:
        return rdf, rdf  # No hodo trigger channels available for this run

    trigger_name_top, trigger_name_bottom = trigger_names[0], trigger_names[1]
    print(
        f"Filtering prefire events with TS >= {TS} using triggers: {trigger_name_top}, {trigger_name_bottom}")
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
