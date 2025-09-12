# collect all the FERS variables on the fly

def vectorizeFERS(rdf, FERSBoards):
    """
    FRES board outputs
    define variables as RDF does not support reading vectors
    with indices directly
    """
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


def subtractFERSPedestal(rdf, FERSBoards, pedestalsHG, pedestalsLG=None):
    for _, FERSBoard in FERSBoards.items():
        for channel in FERSBoard:
            channelNameHG = channel.GetChannelName(useHG=True)
            pedestal = pedestalsHG[channelNameHG]

            # subtract pedestal from HG and LG energies
            rdf = rdf.Define(
                f"{channelNameHG}_subtracted",
                f"{channelNameHG} - {pedestal}"
            )

            if pedestalsLG is not None:
                channelNameLG = channel.GetChannelName(useHG=False)
                pedestal_LG = pedestalsLG[channelNameLG]
                rdf = rdf.Define(
                    f"{channelNameLG}_subtracted",
                    f"{channelNameLG} - {pedestal_LG}"
                )
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


def getFERSBoardMax(rdf, FERSBoards):
    """
    Define per board max energy
    and per event max energy
    for cer and sci, HG and LG
    """
    branches_cer_HG = []
    branches_cer_LG = []
    branches_sci_HG = []
    branches_sci_LG = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        branches_board_cer_HG = []
        branches_board_cer_LG = []
        branches_board_sci_HG = []
        branches_board_sci_LG = []
        for chan in FERSBoard:
            if chan.isCer:
                branches_board_cer_HG.append(chan.GetHGChannelName())
                branches_board_cer_LG.append(chan.GetLGChannelName())
            else:
                branches_board_sci_HG.append(chan.GetHGChannelName())
                branches_board_sci_LG.append(chan.GetLGChannelName())
        rdf = rdf.Define(f"FERS_Board{boardNo}_energy_cer_HG_max",
                         f"std::max({{{', '.join(branches_board_cer_HG)}}})")
        rdf = rdf.Define(f"FERS_Board{boardNo}_energy_cer_LG_max",
                         f"std::max({{{', '.join(branches_board_cer_LG)}}})")
        rdf = rdf.Define(f"FERS_Board{boardNo}_energy_sci_HG_max",
                         f"std::max({{{', '.join(branches_board_sci_HG)}}})")
        rdf = rdf.Define(f"FERS_Board{boardNo}_energy_sci_LG_max",
                         f"std::max({{{', '.join(branches_board_sci_LG)}}})")
        branches_cer_HG.append(f"FERS_Board{boardNo}_energy_cer_HG_max")
        branches_cer_LG.append(f"FERS_Board{boardNo}_energy_cer_LG_max")
        branches_sci_HG.append(f"FERS_Board{boardNo}_energy_sci_HG_max")
        branches_sci_LG.append(f"FERS_Board{boardNo}_energy_sci_LG_max")

    rdf = rdf.Define("FERS_energy_cer_HG_max",
                     f"std::max({{{', '.join(branches_cer_HG)}}})")
    rdf = rdf.Define("FERS_energy_cer_LG_max",
                     f"std::max({{{', '.join(branches_cer_LG)}}})")
    rdf = rdf.Define("FERS_energy_sci_HG_max",
                     f"std::max({{{', '.join(branches_sci_HG)}}})")
    rdf = rdf.Define("FERS_energy_sci_LG_max",
                     f"std::max({{{', '.join(branches_sci_LG)}}})")
    return rdf


def getFERSEnergySum(rdf, FERSBoards, subtractPedestal=False, calibrate=False, clip=False):
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
        string_CerEnergyLG = "+".join(
            chan.GetLGChannelName() + suffix for chan in channels_Cer
        )
        string_SciEnergyHG = "+".join(
            chan.GetHGChannelName() + suffix for chan in channels_Sci
        )
        string_SciEnergyLG = "+".join(
            chan.GetLGChannelName() + suffix for chan in channels_Sci
        )
        # per-board energy sum
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_CerEnergyHG" + suffix,
            f"({string_CerEnergyHG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_CerEnergyLG" + suffix,
            f"({string_CerEnergyLG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_SciEnergyHG" + suffix,
            f"({string_SciEnergyHG})")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_SciEnergyLG" + suffix,
            f"({string_SciEnergyLG})")
    # per-event energy sum
    string_CerEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_CerEnergyHG" + suffix for boardNo in boardNos
    )
    string_CerEnergyLG_Total = "+".join(
        f"FERS_Board{boardNo}_CerEnergyLG" + suffix for boardNo in boardNos
    )
    string_SciEnergyHG_Total = "+".join(
        f"FERS_Board{boardNo}_SciEnergyHG" + suffix for boardNo in boardNos
    )
    string_SciEnergyLG_Total = "+".join(
        f"FERS_Board{boardNo}_SciEnergyLG" + suffix for boardNo in boardNos
    )
    rdf = rdf.Define("FERS_CerEnergyHG" + suffix,
                     f"({string_CerEnergyHG_Total})")
    rdf = rdf.Define("FERS_CerEnergyLG" + suffix,
                     f"({string_CerEnergyLG_Total})")
    rdf = rdf.Define("FERS_SciEnergyHG" + suffix,
                     f"({string_SciEnergyHG_Total})")
    rdf = rdf.Define("FERS_SciEnergyLG" + suffix,
                     f"({string_SciEnergyLG_Total})")

    return rdf
