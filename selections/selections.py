def filterPrefireEvents(rdf, runNumber, TS=350):
    # use the hodo trigger to filter prefire events
    # in cosmic runs
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


def vetoMuonCounter(rdf, TSmin=400, TSmax=600, cut=-80):
    from utils.channel_map import getDownStreamMuonChannel
    muon_channel = getDownStreamMuonChannel()
    if muon_channel is None:
        print("Muon counter channel not found, skipping veto.")
        return rdf

    print(
        f"Vetoing events with muon counter signal between {TSmin} and {TSmax} with cut {cut}")
    rdf = rdf.Define("MuonCounterMin",
                     f"MinRange({muon_channel}_subtractMedian, {TSmin}, {TSmax})")
    rdf = rdf.Define("HasMuonCounterSignal",
                     f"MuonCounterMin < {cut}")

    rdf_prefilter = rdf
    rdf = rdf.Filter("HasMuonCounterSignal == 0")
    print(
        f"Events before and after muon counter veto: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter


def PSDSelection(rdf, runNumber, isHadron=False):
    from utils.channel_map import getPreShowerChannel
    preshower_channel = getPreShowerChannel(runNumber)
    if preshower_channel is None:
        print("Pre-shower channel not found, skipping PSD selection.")
        return rdf

    print("Applying PSD selection based on pre-shower channel.")
    rdf = rdf.Define(f"{preshower_channel}_peak_value",
                     f"MinRange({preshower_channel}_subtractMedian, 100, 400)")
    if not isHadron:
        rdf = rdf.Define("pass_psd_selection",
                         f"({preshower_channel}_peak_value < -200.0)")
    else:
        rdf = rdf.Define("pass_psd_selection",
                         f"({preshower_channel}_peak_value > -200.0)")

    rdf_prefilter = rdf
    rdf = rdf.Filter("pass_psd_selection == 1")
    print(
        f"Events before and after PSD selection: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter


def checkUpstreamVeto(rdf, runNumber):
    from utils.channel_map import getUpstreamVetoChannel
    chan_upveto = getUpstreamVetoChannel(runNumber)

    rdf = rdf.Define(f"{chan_upveto}_peak_position",
                     f"ROOT::VecOps::ArgMin({chan_upveto}_subtractMedian)")
    rdf = rdf.Define(f"{chan_upveto}_peak_value",
                     f"ROOT::VecOps::Min({chan_upveto}_subtractMedian)")

    rdf = rdf.Define(f"pass_upstream_veto",
                     f"({chan_upveto}_peak_value > -1000.0)") \
        .Define("pass_NoSel", "1.0")

    return rdf


def applyUpstreamVeto(rdf, runNumber):
    rdf = checkUpstreamVeto(rdf, runNumber)
    rdf_prefilter = rdf
    rdf = rdf.Filter("pass_upstream_veto == 1")
    print(
        f"Events before and after upstream veto: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter
