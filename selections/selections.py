def filterPrefireEvents(rdf, runNumber, TS=350):
    # use the hodo trigger to filter prefire events
    # in cosmic runs
    from channels.channel_map import buildHodoTriggerChannels
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


def getPSDSumCutValue():
    return -1e3


def getCC1SumCutValue():
    return -1e3


def vetoMuonCounter(rdf, TSmin=400, TSmax=600, cut=-80):
    from channels.channel_map import getDownStreamMuonChannel
    muon_channel = getDownStreamMuonChannel()
    if muon_channel is None:
        print("Muon counter channel not found, skipping veto.")
        return rdf

    print(
        f"Vetoing events with muon counter signal between {TSmin} and {TSmax} with cut {cut}")
    rdf = rdf.Define("MuonCounterMin",
                     f"MinRange({muon_channel}_blsub, {TSmin}, {TSmax})")
    rdf = rdf.Define("HasMuonCounterSignal",
                     f"MuonCounterMin < {cut}")

    rdf_prefilter = rdf
    rdf = rdf.Filter("HasMuonCounterSignal == 0")
    print(
        f"Events before and after muon counter veto: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter


def applyPSDSelection(rdf, runNumber, isHadron=False, applyCut=True):
    from channels.channel_map import getPreShowerChannel
    preshower_channel = getPreShowerChannel(runNumber)
    if preshower_channel is None:
        print("Pre-shower channel not found, skipping PSD selection.")
        return rdf

    print("Applying PSD selection based on pre-shower channel.")
    rdf = rdf.Define(f"{preshower_channel}_peak_value",
                     f"MinRange({preshower_channel}_blsub, 100, 400)")
    rdf = rdf.Define(f"{preshower_channel}_sum",
                     f"SumRange({preshower_channel}_blsub, 100, 400)")

    valCut = getPSDSumCutValue()
    rdf = rdf.Define("pass_PSDEle_selection",
                     f"({preshower_channel}_sum < {valCut})")

    # rdf = rdf.Define("pass_PSDEle_selection", f"({preshower_channel}_peak_value < -200.0)")

    if not applyCut:
        return rdf

    rdf_prefilter = rdf
    if not isHadron:
        rdf = rdf.Filter("pass_PSDEle_selection == 1")
    else:
        rdf = rdf.Filter("pass_PSDEle_selection == 0")
    print(
        f"Events before and after PSD selection: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter


def applyCC1Selection(rdf, runNumber, isHadron=False, applyCut=True):
    from channels.channel_map import getCerenkovCounters
    cerenkov_channels = getCerenkovCounters(runNumber)
    if cerenkov_channels is None or len(cerenkov_channels) == 0:
        print("Cerenkov channels not found, skipping CC1 selection.")
        return rdf

    # Use the first Cerenkov channel for CC1 selection
    cc1_channel = cerenkov_channels[0]

    print("Applying CC1 selection based on Cerenkov1 channel.")
    rdf = rdf.Define(f"{cc1_channel}_peak_value",
                     f"MinRange({cc1_channel}_blsub, 600, 800)")
    rdf = rdf.Define(f"{cc1_channel}_sum",
                     f"SumRange({cc1_channel}_blsub, 600, 800)")

    valCut = getCC1SumCutValue()
    rdf = rdf.Define("pass_CC1Ele_selection",
                     f"({cc1_channel}_sum < {valCut})")

    #    rdf = rdf.Define("pass_cc1_selection",
    #                     f"({cerenkov1_channel}_peak_value > -200.0)")
    if not applyCut:
        return rdf

    rdf_prefilter = rdf
    if not isHadron:
        rdf = rdf.Filter("pass_CC1Ele_selection == 1")
    else:
        rdf = rdf.Filter("pass_CC1Ele_selection == 0")
    print(
        f"Events before and after CC1 selection: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter


def applyUpstreamVeto(rdf, runNumber, applyCut=True):
    from channels.channel_map import getUpstreamVetoChannel
    chan_upveto = getUpstreamVetoChannel(runNumber)

    rdf = rdf.Define(f"{chan_upveto}_peak_position",
                     f"ROOT::VecOps::ArgMin({chan_upveto}_blsub)")
    rdf = rdf.Define(f"{chan_upveto}_peak_value",
                     f"ROOT::VecOps::Min({chan_upveto}_blsub)")

    rdf = rdf.Define(f"pass_upstream_veto",
                     f"({chan_upveto}_peak_value > -1000.0)")
    rdf = rdf.Define("pass_NoSel", "1.0")
    if not applyCut:
        return rdf

    rdf_prefilter = rdf
    rdf = rdf.Filter("pass_upstream_veto == 1")
    print(
        f"Events before and after upstream veto: {rdf_prefilter.Count().GetValue()}, {rdf.Count().GetValue()}")
    return rdf, rdf_prefilter
