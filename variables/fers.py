# collect all the FERS variables on the fly
import json
from configs.plotranges import getFERSSaturationValue
from channels.gainvalidator import enforce_gain


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
                channel.GetChannelName(gain="HG"),
                f"FERS_Board{boardNo}_energyHG[{channel.channelNo}]")
            rdf = rdf.Define(
                channel.GetChannelName(gain="LG"),
                f"FERS_Board{boardNo}_energyLG[{channel.channelNo}]"
            )
    return rdf


def addFERSPosXY(rdf, fersboards):
    """
    Add iTowerX, iTowerY, RealX, and RealY for each FERS channel
    """
    for fersboard in fersboards.values():
        for channel in fersboard:
            rdf = rdf.Define(
                channel.GetTowerPosName(isX=True),
                f"{channel.iTowerX}"
            )
            # unit cm
            rdf = rdf.Define(
                channel.GetRealPosName(isX=True),
                f"{channel.iTowerX * 1.2}"
            )
            rdf = rdf.Define(
                channel.GetTowerPosName(isX=False),
                f"{channel.iTowerY}"
            )
            # unit cm
            rdf = rdf.Define(
                channel.GetRealPosName(isX=False),
                f"{channel.iTowerY * 1.6}"
            )
    return rdf


@enforce_gain
def subtractFERSPedestal(rdf, fersboards, gain="HG", file_pedestals: dict = None):
    if file_pedestals is None:
        print("\033[93mNo pedestal provided, skip pedestal subtraction.\033[0m")
        return rdf
    with open(file_pedestals, 'r') as f:
        pedestals = json.load(f)

    for fersboard in fersboards.values():
        for channel in fersboard:
            channelName = channel.GetChannelName(gain=gain)
            channelNamePDSub = channel.GetChannelName(
                gain=gain, pdsub=True)
            if f"{channelNamePDSub}" in rdf.GetColumnNames():
                # already defined
                continue
            if channelName not in pedestals:
                raise ValueError(
                    f"Pedestal for channel {channelName} not found in pedestals.")
            pedestal = pedestals[channelName]
            rdf = rdf.Define(channelNamePDSub,
                             f"{channelName} - {pedestal}")

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
            channelName_HG = channel.GetChannelName(gain="HG", pdsub=True)
            channelName_LG = channel.GetChannelName(gain="LG", pdsub=True)
            channelName_Mix = channel.GetChannelName(gain="Mix")
            channelName_key = f"FERS_Board{fersboard.boardNo}_{channel.channelNo}"
            if f"{channelName_Mix}" in rdf.GetColumnNames():
                # already defined
                continue
            if channelName_key not in HG2LG_ratios:
                raise ValueError(
                    f"HG2LG ratio for channel {channelName_HG} not found in {file_HG2LG}.")
            incep, ratio = HG2LG_ratios[channelName_key]
            rdf = rdf.Define(channelName_Mix,
                             f"({channelName_HG} < {getFERSSaturationValue()}) ? {channelName_HG} : ({channelName_LG} - {incep})/{ratio}")

    return rdf


def calibrateFERSChannels(rdf, fersboards, file_calibrations: str, gain="HG", file_deadchannels: str = None, toyCalib: bool = False):
    """
    Calibrate FERS channels using gains and pedestals from the provided files.
    """
    import json
    with open(file_calibrations, 'r') as f:
        calibrations = json.load(f)

    deadchannels = {}
    if file_deadchannels is not None:
        with open(file_deadchannels, 'r') as f:
            deadchannels = json.load(f)
        print(
            f"\033[93mThe following channels are marked as dead and will be set to 0:\033[0m")

    for fersboard in fersboards.values():
        for channel in fersboard:
            channelNamePDSub = channel.GetChannelName(
                gain=gain, pdsub=True)
            channelNamePDSubCal = channel.GetChannelName(
                gain=gain, pdsub=True, calib=True)
            if f"{channelNamePDSubCal}" in rdf.GetColumnNames():
                # already defined
                continue
            if str(fersboard.boardNo) in deadchannels and str(channel.channelNo) in deadchannels[str(fersboard.boardNo)]:
                print(
                    f"\033[93mChannel {channelNamePDSub} is marked as dead. Channel Masked.\033[0m")
                rdf = rdf.Define(channelNamePDSubCal, "0.")
            elif channelNamePDSub not in calibrations:
                print(
                    f"\033[93mGain for channel {channelNamePDSub} not found in gains. Channel Masked.\033[0m")
                rdf = rdf.Define(channelNamePDSubCal, "0.")
            else:
                if not toyCalib:
                    calibration = calibrations[channelNamePDSub]["response"]
                else:
                    if channel.isCer and channel.isQuartz:
                        calibration = 80.0 / 62000.0
                    elif channel.isCer and not channel.isQuartz:
                        calibration = 80.0 / 210000.0
                    elif not channel.isCer:
                        calibration = 80.0 / 300000.0
                    else:
                        print(
                            "\033[93mUnknown channel type for toy calibration. Channel Masked.\033[0m")
                rdf = rdf.Define(channelNamePDSubCal,
                                 f"({channelNamePDSub} * {calibration})")
    return rdf


def getFERSEnergyMax(rdf, fersboards, gain="HG"):
    """
    Define per board max energy
    and per event max energy
    for Cer and Sci
    """
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.GetEnergyMaxName(gain=gain, isCer=True),
                         f"std::max({{{', '.join(chan.GetChannelName(gain=gain) for chan in fersboard.GetListOfChannels(isCer=True))}}})")
        rdf = rdf.Define(fersboard.GetEnergyMaxName(gain=gain, isCer=False),
                         f"std::max({{{', '.join(chan.GetChannelName(gain=gain) for chan in fersboard.GetListOfChannels(isCer=False))}}})")

    rdf = rdf.Define(fersboards.GetEnergyMaxName(gain=gain, isCer=True),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(gain=gain, isCer=True) for fersboard in fersboards.values())}}})")
    rdf = rdf.Define(fersboards.GetEnergyMaxName(gain=gain, isCer=False),
                     f"std::max({{{', '.join(fersboard.GetEnergyMaxName(gain=gain, isCer=False) for fersboard in fersboards.values())}}})")
    return rdf


def getFERSEnergySum(rdf, fersboards, gain="HG", pdsub=False, calib=False):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    # per-board energy sum
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.GetEnergySumName(gain=gain, isCer=True, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(gain=gain, pdsub=pdsub, calib=calib) for chan in fersboard.GetListOfChannels(isCer=True))})")
        rdf = rdf.Define(fersboard.GetEnergySumName(gain=gain, isCer=False, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.GetChannelName(gain=gain, pdsub=pdsub, calib=calib) for chan in fersboard.GetListOfChannels(isCer=False))})")

    # per-event energy sum
    rdf = rdf.Define(fersboards.GetEnergySumName(gain=gain, isCer=True, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(gain=gain, isCer=True, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.GetEnergySumName(gain=gain, isCer=False, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.GetEnergySumName(gain=gain, isCer=False, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    return rdf


def getFERSEnergyWeightedCenter(rdf, fersboards, gain="HG", pdsub=False, calib=False, useRealPos=True):
    """
    Calculate the weighted center of the energy distribution for FERS boards.
    """
    scaleX = 1.2 if useRealPos else 1.0
    scaleY = 1.6 if useRealPos else 1.0
    for var in ["Cer", "Sci"]:
        x_center = "0."
        y_center = "0."
        for fersboard in fersboards.values():
            for channel in fersboard.GetListOfChannels(isCer=(var == "Cer")):
                channelName = channel.GetChannelName(
                    gain=gain, pdsub=pdsub, calib=calib)
                x_center += f"+ {channel.iTowerX * scaleX} * {channelName}"
                y_center += f"+ {channel.iTowerY * scaleY} * {channelName}"

        rdf = rdf.Define(fersboards.GetEnergyWeightedCenterName(gain=gain, isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=True),
                         f"({x_center}) / ({fersboards.GetEnergySumName(gain=gain, isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")
        rdf = rdf.Define(fersboards.GetEnergyWeightedCenterName(gain=gain, isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=False),
                         f"({y_center}) / ({fersboards.GetEnergySumName(gain=gain, isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")

    return rdf
