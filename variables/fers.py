# collect all the FERS variables on the fly

def vectorizeFERS(rdf, fersboards):
    """
    FRES board outputs
    define variables as RDF does not support reading vectors
    with indices directly
    """
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        for channel in fersboard:
            rdf = rdf.Define(
                channel.GetChannelName(useHG=True),
                f"FERS_Board{boardNo}_energyHG[{channel.channelNo}]")
            rdf = rdf.Define(
                channel.GetChannelName(useHG=False),
                f"FERS_Board{boardNo}_energyLG[{channel.channelNo}]"
            )
    return rdf


def subtractFERSPedestal(rdf, fersboards, pedestalsHG: dict = None, pedestalsLG: dict = None):
    for fersboard in fersboards.values():
        for channel in fersboard:
            if pedestalsHG is not None:
                channelName = channel.GetChannelName(useHG=True)
                channelNamePDSub = channel.GetChannelName(
                    useHG=True, pdsub=True)
                if f"{channelNamePDSub}" in rdf.GetColumnNames():
                    # already defined
                    continue
                if channelName not in pedestalsHG:
                    raise ValueError(
                        f"Pedestal for channel {channelName} not found in pedestalsHG.")
                pedestal = pedestalsHG[channelName]
                rdf = rdf.Define(channelNamePDSub,
                                 f"{channelName} - {pedestal}")

            if pedestalsLG is not None:
                channelName = channel.GetChannelName(useHG=False)
                channelNamePDSub = channel.GetChannelName(
                    useHG=False, pdsub=True)
                if f"{channelNamePDSub}" in rdf.GetColumnNames():
                    # already defined
                    continue
                if channelName not in pedestalsLG:
                    raise ValueError(
                        f"Pedestal for channel {channelName} not found in pedestalsLG.")
                pedestal_LG = pedestalsLG[channelName]
                rdf = rdf.Define(channelNamePDSub,
                                 f"{channelName} - {pedestal_LG}")
    return rdf


def calibrateFERSChannels(rdf, fersboards, file_gains: str, file_pedestals: str):
    """
    Calibrate FERS channels using gains and pedestals from the provided files.
    """
    import json
    with open(file_gains, 'r') as f:
        gains = json.load(f)
    with open(file_pedestals, 'r') as f:
        pedestals = json.load(f)

    # Subtract pedestal and apply gain calibration
    rdf = subtractFERSPedestal(rdf, fersboards, pedestals)

    for fersboard in fersboards.values():
        for channel in fersboard:
            channelName = channel.GetChannelName(useHG=True)
            channelNamePDSub = channel.GetChannelName(
                useHG=True, pdsub=True)
            channelNamePDSubCal = channel.GetChannelName(
                useHG=True, pdsub=True, calib=True)
            if f"{channelNamePDSubCal}" in rdf.GetColumnNames():
                # already defined
                continue
            if channelName not in gains:
                raise ValueError(
                    f"Gain for channel {channelName} not found in gains.")
            gain = gains[channelName]
            rdf = rdf.Define(channelNamePDSubCal,
                             f"{channelNamePDSub} / {gain}")
            # rdf = rdf.Define(
            #    f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calib_clipped",
            #    f"FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calib > 0.8 ? FERS_Board{boardNo}_energyHG_{channelNo}_subtracted_calib : 0"
            # )
            # Not calibrating LG yet

    return rdf


def getFERSEnergyMax(rdf, fersboards):
    """
    Define per board max energy
    and per event max energy
    for Cer and Sci, HG and LG
    """
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.GetEnergyMaxName(useHG=True, isCer=True),
                         f"std::max({{{', '.join(chan.GetChannelName(useHG=True) for chan in fersboard.GetCerChannels())}}})")
        rdf = rdf.Define(fersboard.GetEnergyMaxName(useHG=False, isCer=True),
                         f"std::max({{{', '.join(chan.GetChannelName(useHG=False) for chan in fersboard.GetCerChannels())}}})")
        rdf = rdf.Define(fersboard.GetEnergyMaxName(useHG=True, isCer=False),
                         f"std::max({{{', '.join(chan.GetChannelName(useHG=True) for chan in fersboard.GetSciChannels())}}})")
        rdf = rdf.Define(fersboard.GetEnergyMaxName(useHG=False, isCer=False),
                         f"std::max({{{', '.join(chan.GetChannelName(useHG=False) for chan in fersboard.GetSciChannels())}}})")

    rdf = rdf.Define(fersboards.GetEnergyMaxName(useHG=True, isCer=True),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(useHG=True, isCer=True) for fersboard in fersboards.values())}}})")
    rdf = rdf.Define(fersboards.GetEnergyMaxName(useHG=False, isCer=True),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(useHG=False, isCer=True) for fersboard in fersboards.values())}}})")
    rdf = rdf.Define(fersboards.GetEnergyMaxName(useHG=True, isCer=False),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(useHG=True, isCer=False) for fersboard in fersboards.values())}}})")
    rdf = rdf.Define(fersboards.GetEnergyMaxName(useHG=False, isCer=False),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(useHG=False, isCer=False) for fersboard in fersboards.values())}}})")
    return rdf


def getFERSEnergySum(rdf, fersboards, pdsub=False, calib=False):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    # per-board energy sum
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=True, isCer=True, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=True, pdsub=pdsub, calib=calib) for chan in fersboard.GetCerChannels())})")
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=False, isCer=True, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=False, pdsub=pdsub, calib=calib) for chan in fersboard.GetCerChannels())})")
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=True, isCer=False, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=True, pdsub=pdsub, calib=calib) for chan in fersboard.GetSciChannels())})")
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=False, isCer=False, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=False, pdsub=pdsub, calib=calib) for chan in fersboard.GetSciChannels())})")

    # per-event energy sum
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=True, isCer=True, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=True, isCer=True, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=False, isCer=True, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=False, isCer=True, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=True, isCer=False, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=True, isCer=False, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=False, isCer=False, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=False, isCer=False, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    return rdf


def getFERSEnergyWeightedCenter(rdf, fersboards, pdsub=False, calib=False):
    """
    Calculate the weighted center of the energy distribution for FERS boards.
    """
    for gain in ["HG", "LG"]:
        for var in ["Cer", "Sci"]:
            x_center = "0."
            y_center = "0."
            for fersboard in fersboards.values():
                for channel in fersboard.GetListOfChannels(isCer=(var == "Cer")):
                    channelName = channel.GetChannelName(
                        useHG=(gain == "HG"), pdsub=pdsub, calib=calib)
                    x_center += f"+ {channel.iTowerX} * {channelName}"
                    y_center += f"+ {channel.iTowerY} * {channelName}"

            rdf = rdf.Define(fersboards.GetEnergyWeightedCenterName(useHG=(gain == "HG"), isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=True),
                             f"({x_center}) / ({fersboards.GetEnergySumName(useHG=(gain == 'HG'), isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")
            rdf = rdf.Define(fersboards.GetEnergyWeightedCenterName(useHG=(gain == "HG"), isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=False),
                             f"({y_center}) / ({fersboards.GetEnergySumName(useHG=(gain == 'HG'), isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")

    return rdf
