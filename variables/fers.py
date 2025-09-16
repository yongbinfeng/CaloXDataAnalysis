# collect all the FERS variables on the fly
from configs.plotranges import getFERSSaturationValue


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


def addFERSPosXY(rdf, fersboards):
    """
    Add iTowerX and iTowerY for each FERS channel
    """
    for fersboard in fersboards.values():
        for channel in fersboard:
            rdf = rdf.Define(
                channel.GetPosName(isX=True),
                f"{channel.iTowerX}"
            )
            rdf = rdf.Define(
                channel.GetPosName(isX=False),
                f"{channel.iTowerY}"
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


def mixFERSHGLG(rdf, fersboards, file_HG2LG: str):
    """
    Mix FERS HG and LG channels using HG2LG ratios from the provided file.
    """
    import json
    with open(file_HG2LG, 'r') as f:
        HG2LG_ratios = json.load(f)

    for fersboard in fersboards.values():
        for channel in fersboard:
            channelName_HG = channel.GetChannelName(useHG=True, pdsub=True)
            channelName_LG = channel.GetChannelName(useHG=False, pdsub=True)
            channelName_Mixed = channel.GetChannelName(
                useHG=None)  # None for mixed
            channelName_key = f"FERS_Board{fersboard.boardNo}_{channel.channelNo}"
            if f"{channelName_Mixed}" in rdf.GetColumnNames():
                # already defined
                continue
            if channelName_key not in HG2LG_ratios:
                raise ValueError(
                    f"HG2LG ratio for channel {channelName_HG} not found in {file_HG2LG}.")
            incep, ratio = HG2LG_ratios[channelName_key]
            rdf = rdf.Define(channelName_Mixed,
                             f"({channelName_HG} < {getFERSSaturationValue() - 150.0}) ? {channelName_HG} : ({channelName_LG} - {incep})/{ratio}")

    return rdf


def calibrateFERSChannels(rdf, fersboards, file_gains: str, file_pedestals: str, useMix=False, file_HG2LG: str = None, file_pedestals_LG: str = None):
    """
    Calibrate FERS channels using gains and pedestals from the provided files.
    """
    import json
    with open(file_gains, 'r') as f:
        gains = json.load(f)
    with open(file_pedestals, 'r') as f:
        pedestals = json.load(f)
    pedestals_LG = None
    if file_pedestals_LG is not None:
        with open(file_pedestals_LG, 'r') as f:
            pedestals_LG = json.load(f)

    # Subtract pedestal
    rdf = subtractFERSPedestal(rdf, fersboards, pedestals, pedestals_LG)

    if useMix:
        if file_HG2LG is None:
            raise ValueError(
                "file_HG2LG must be provided when useMix is True.")
        rdf = mixFERSHGLG(rdf, fersboards, file_HG2LG)

    for fersboard in fersboards.values():
        for channel in fersboard:
            if not useMix:
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
            else:
                channelName = channel.GetChannelName(useHG=None)
                channelNamePDSubCal = channel.GetChannelName(
                    useHG=None, pdsub=True, calib=True)
                if f"{channelNamePDSubCal}" in rdf.GetColumnNames():
                    # already defined
                    continue
                if channelName not in gains:
                    # raise ValueError(
                    #    f"Gain for channel {channelName} not found in gains.")
                    # masked
                    print(
                        f"\033[93mGain for channel {channelName} not found in gains. Channel Masked.\033[0m")
                    rdf = rdf.Define(channelNamePDSubCal, f"0.")
                else:
                    gain = gains[channelName]["response"]
                    rdf = rdf.Define(channelNamePDSubCal,
                                     f"{channelName} * {gain}")
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


def getFERSEnergySum(rdf, fersboards, useHG=True, pdsub=False, calib=False):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    # per-board energy sum
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=useHG, isCer=True, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=useHG, pdsub=pdsub, calib=calib) for chan in fersboard.GetCerChannels())})")
        rdf = rdf.Define(fersboard.GetEnergySumName(useHG=useHG, isCer=False, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(useHG=useHG, pdsub=pdsub, calib=calib) for chan in fersboard.GetSciChannels())})")

    # per-event energy sum
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=useHG, isCer=True, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=useHG, isCer=True, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.GetEnergySumName(useHG=useHG, isCer=False, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(useHG=useHG, isCer=False, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
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
