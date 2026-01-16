# collect all the FERS variables on the fly
import json
from configs.plot_config import get_fers_saturation_value
from channels.gain_validator import enforce_gain


def vectorizeFERS(rdf, fersboards):
    """
    FRES board outputs
    define variables as RDF does not support reading vectors
    with indices directly
    """
    for fersboard in fersboards.values():
        board_no = fersboard.board_no
        for channel in fersboard:
            rdf = rdf.Define(
                channel.get_channel_name(gain="HG"),
                f"FERS_Board{board_no}_energyHG[{channel.channel_no}]")
            rdf = rdf.Define(
                channel.get_channel_name(gain="LG"),
                f"FERS_Board{board_no}_energyLG[{channel.channel_no}]"
            )
    return rdf


def addFERSPosXY(rdf, fersboards):
    """
    Add iTowerX, iTowerY, RealX, and RealY for each FERS channel
    """
    for fersboard in fersboards.values():
        for channel in fersboard:
            rdf = rdf.Define(
                channel.get_tower_pos_name(isX=True),
                f"{channel.i_tower_x}"
            )
            # unit cm
            rdf = rdf.Define(
                channel.get_real_pos_name(isX=True),
                f"{channel.i_tower_x * 1.2}"
            )
            rdf = rdf.Define(
                channel.get_tower_pos_name(isX=False),
                f"{channel.i_tower_y}"
            )
            # unit cm
            rdf = rdf.Define(
                channel.get_real_pos_name(isX=False),
                f"{channel.i_tower_y * 1.6}"
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
            channelName = channel.get_channel_name(gain=gain)
            channelNamePDSub = channel.get_channel_name(
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
            channelName_HG = channel.get_channel_name(gain="HG", pdsub=True)
            channelName_LG = channel.get_channel_name(gain="LG", pdsub=True)
            channelName_Mix = channel.get_channel_name(gain="Mix")
            channelName_key = f"FERS_Board{fersboard.board_no}_{channel.channel_no}"
            if f"{channelName_Mix}" in rdf.GetColumnNames():
                # already defined
                continue
            if channelName_key not in HG2LG_ratios:
                raise ValueError(
                    f"HG2LG ratio for channel {channelName_HG} not found in {file_HG2LG}.")
            incep, ratio = HG2LG_ratios[channelName_key]
            rdf = rdf.Define(channelName_Mix,
                             f"({channelName_HG} < {get_fers_saturation_value()}) ? {channelName_HG} : ({channelName_LG} - {incep})/{ratio}")

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
            channelNamePDSub = channel.get_channel_name(
                gain=gain, pdsub=True)
            channelNamePDSubCal = channel.get_channel_name(
                gain=gain, pdsub=True, calib=True)
            if f"{channelNamePDSubCal}" in rdf.GetColumnNames():
                # already defined
                continue
            if str(fersboard.board_no) in deadchannels and str(channel.channel_no) in deadchannels[str(fersboard.board_no)]:
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


def get_fers_energy_max(rdf, fersboards, gain="HG"):
    """
    Define per board max energy
    and per event max energy
    for Cer and Sci
    """
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.get_energy_max_name(gain=gain, isCer=True),
                         f"std::max({{{', '.join(chan.get_channel_name(gain=gain) for chan in fersboard.get_list_of_channels(isCer=True))}}})")
        rdf = rdf.Define(fersboard.get_energy_max_name(gain=gain, isCer=False),
                         f"std::max({{{', '.join(chan.get_channel_name(gain=gain) for chan in fersboard.get_list_of_channels(isCer=False))}}})")

    rdf = rdf.Define(fersboards.get_energy_max_name(gain=gain, isCer=True),
                     f"std::max({{{', '.join(fersboard.get_energy_max_name(gain=gain, isCer=True) for fersboard in fersboards.values())}}})")
    rdf = rdf.Define(fersboards.get_energy_max_name(gain=gain, isCer=False),
                     f"std::max({{{', '.join(fersboard.get_energy_max_name(gain=gain, isCer=False) for fersboard in fersboards.values())}}})")
    return rdf


def get_fers_energy_sum(rdf, fersboards, gain="HG", pdsub=False, calib=False):
    """
    Calculate the Sci and Cer energy sum for FERS boards, per board and per event.
    """
    # per-board energy sum
    for fersboard in fersboards.values():
        rdf = rdf.Define(fersboard.get_energy_sum_name(gain=gain, isCer=True, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.get_channel_name(gain=gain, pdsub=pdsub, calib=calib) for chan in fersboard.get_list_of_channels(isCer=True))})")
        rdf = rdf.Define(fersboard.get_energy_sum_name(gain=gain, isCer=False, pdsub=pdsub, calib=calib),
                         f"({' + '.join(chan.get_channel_name(gain=gain, pdsub=pdsub, calib=calib) for chan in fersboard.get_list_of_channels(isCer=False))})")

    # per-event energy sum
    rdf = rdf.Define(fersboards.get_energy_sum_name(gain=gain, isCer=True, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.get_energy_sum_name(gain=gain, isCer=True, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")
    rdf = rdf.Define(fersboards.get_energy_sum_name(gain=gain, isCer=False, pdsub=pdsub, calib=calib),
                     f"({' + '.join(fersboard.get_energy_sum_name(gain=gain, isCer=False, pdsub=pdsub, calib=calib) for fersboard in fersboards.values())})")

    # calculate adhoc lateral leakage correction
    # E_leak = E_{out ring} * f, f=1.0
    for cat in ["Cer", "Sci"]:
        channels_outerring = []
        for fersboard in fersboards.values():
            for channel in fersboard.get_list_of_channels(isCer=(cat == "Cer")):
                if abs(channel.i_tower_x) >= 8 or abs(channel.i_tower_y) >= 6:
                    channels_outerring.append(
                        channel.get_channel_name(gain=gain, pdsub=pdsub, calib=calib))

        name_total = fersboards.get_energy_sum_name(
            gain=gain, isCer=(cat == "Cer"), pdsub=pdsub, calib=calib)
        name_outer = name_total + "_OuterRing"
        rdf = rdf.Define(f"{name_outer}",
                         f"({' + '.join(channels_outerring)})")
        name_corr = name_total + "_LeakCorr"
        rdf = rdf.Define(f"{name_corr}",
                         f"{name_total} + {name_outer} * 1.0")

    return rdf


def getFERSEnergyDR(rdf, fersboards, energy=0.0):
    #
    # for dual-readout calculation
    #
    name_cer = fersboards.get_energy_sum_name(
        gain="Mix", isCer=True, pdsub=True, calib=True)
    name_sci = fersboards.get_energy_sum_name(
        gain="Mix", isCer=False, pdsub=True, calib=True)
    name_cer_corr = name_cer + "_LeakCorr"
    name_sci_corr = name_sci + "_LeakCorr"

    # per-event energy ratio
    rdf = rdf.Define("COverS", f"{name_cer} / ({name_sci} + 1e-9)")

    eOverh_C = 5.0
    eOverh_S = 1.3
    hOvere_S = 1.0 / eOverh_S
    hOvere_C = 1.0 / eOverh_C

    # fEM = (S * hOvere_C - C * hOvere_S) / ( (1 - hOvere_S) * C - (1 - hOvere_C) * S)
    # rdf = rdf.Define(
    #    "fEM", f"({name_sci} / {eOverh_C} - {name_cer} / {eOverh_S}) / ( (1.0-1.0/{eOverh_S})*{name_cer} - (1.0-1.0/{eOverh_C})*{name_sci})")
    rdf = rdf.Define(
        "fEM", f"({name_sci_corr} * {hOvere_C} - {name_cer} * {hOvere_S}) / ((1.0 - {hOvere_S}) * {name_cer} - (1.0 - {hOvere_C}) * {name_sci_corr})"
    )

    chi = 0.29
    name_dr = name_cer.replace("Cer", "DR")
    name_sum = name_cer.replace("Cer", "CerSci")
    rdf = rdf.Define(name_sum,
                     f"({name_cer} + {name_sci})")
    # DR: E = (E_sci - E_cer * chi) / (1 - chi)
    rdf = rdf.Define(name_dr,
                     f"({name_sci} - {name_cer} * {chi}) / (1 - {chi})")
    # DR method 2 (energy dependent correction):
    # E = S + k * (b * E - (S+C))
    slope = 0.482
    rdf = rdf.Define(name_dr + "_method2",
                     f"{name_sci} + {slope} * (1.766 * {energy}  - {name_sum})")
    slope = 0.789
    rdf = rdf.Define(name_dr + "_method3",
                     f"{name_sci} + 1.0 / {slope} * (60.0 / 80.0 * {energy}  - {name_cer})")
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
            for channel in fersboard.get_list_of_channels(isCer=(var == "Cer")):
                channelName = channel.get_channel_name(
                    gain=gain, pdsub=pdsub, calib=calib)
                x_center += f"+ {channel.i_tower_x * scaleX} * {channelName}"
                y_center += f"+ {channel.i_tower_y * scaleY} * {channelName}"

        rdf = rdf.Define(fersboards.get_energy_weighted_center_name(gain=gain, isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=True),
                         f"({x_center}) / ({fersboards.get_energy_sum_name(gain=gain, isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")
        rdf = rdf.Define(fersboards.get_energy_weighted_center_name(gain=gain, isCer=(var == "Cer"), pdsub=pdsub, calib=calib, isX=False),
                         f"({y_center}) / ({fersboards.get_energy_sum_name(gain=gain, isCer=(var == 'Cer'), pdsub=pdsub, calib=calib)} + 1e-9)")

    return rdf


def buildTTUHodo(rdf, val_cut=4000):
    """
    Build TTU Hodoscope hit information in RDF.
    """
    ttuhodo_x_name = "FERS_Board1_energyHG"
    ttuhodo_y_name = "FERS_Board0_energyHG"
    # find number of hits
    # by summing boolean array where energy deposit > threshold
    hit_x_expr = f"ROOT::VecOps::Sum( ( {ttuhodo_x_name} > {val_cut} ) )"
    hit_y_expr = f"ROOT::VecOps::Sum( ( {ttuhodo_y_name} > {val_cut} ) )"
    rdf = rdf.Define("TTU_Hodo_nHitX", hit_x_expr)
    rdf = rdf.Define("TTU_Hodo_nHitY", hit_y_expr)
    # get the channel index with max energy deposit
    rdf = rdf.Define("TTU_Hodo_MaxChannelX",
                     f"ROOT::VecOps::ArgMax( {ttuhodo_x_name} )")
    rdf = rdf.Define("TTU_Hodo_MaxChannelY",
                     f"ROOT::VecOps::ArgMax( {ttuhodo_y_name} )")
    # map to the real hodo channel
    rdf = rdf.Define("TTU_Hodo_X",
                     f"get_x_index( TTU_Hodo_MaxChannelX )")
    rdf = rdf.Define("TTU_Hodo_Y",
                     f"get_y_index( TTU_Hodo_MaxChannelY )")
    return rdf
