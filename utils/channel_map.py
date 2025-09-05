import numpy as np
from utils.CaloXChannel import FERSBoard, DRSBoard
import json

f_scanruns = "data/scanruns.json"
with open(f_scanruns, 'r') as f:
    temp = json.load(f)
    scanruns = temp["scanruns"]
    print("Loaded scan runs from", f_scanruns)
print("Scan runs:", scanruns)

f_triggerdelay = "data/triggerdelay.json"
with open(f_triggerdelay, 'r') as f:
    triggerdelay = json.load(f)


def buildFERSBoards(run=316):
    """
    Build a map for ixy and FERS channels for both boards.
    """
    base_FERSBoard_6mm = FERSBoard(boardNo=-1, is6mm=True)
    base_FERSBoard_3mm = FERSBoard(boardNo=-1, is6mm=False)
    FERSBoards = {}
    if run == 316:
        # 5 FERS board in 316
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board3"] = base_FERSBoard_6mm.copy(boardNo=3)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board5"] = base_FERSBoard_6mm.copy(boardNo=5)

        # positions
        FERSBoards["Board1"].MoveTo(6.5, -0.5)
        FERSBoards["Board2"].MoveTo(2.5, -0.5)
        FERSBoards["Board3"].MoveTo(-1.5, -2.5)
        FERSBoards["Board4"].MoveTo(-5.5, -0.5)
        FERSBoards["Board5"].MoveTo(-9.5, -0.5)

    elif run == 571:
        # 10 FERS boards in 571
        FERSBoards["Board0"] = base_FERSBoard_6mm.copy(boardNo=0)
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board3"] = base_FERSBoard_6mm.copy(boardNo=3)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board5"] = base_FERSBoard_6mm.copy(boardNo=5)
        FERSBoards["Board8"] = base_FERSBoard_6mm.copy(boardNo=8)
        FERSBoards["Board9"] = base_FERSBoard_6mm.copy(boardNo=9)
        FERSBoards["Board10"] = base_FERSBoard_6mm.copy(boardNo=10)
        FERSBoards["Board11"] = base_FERSBoard_6mm.copy(boardNo=11)
        # FERSBoards["Board12"] = base_FERSBoard_6mm.copy(boardNo=12)

        FERSBoards["Board0"].MoveTo(-13.5, 3.5)
        FERSBoards["Board1"].MoveTo(-9.5, 7.5)
        FERSBoards["Board2"].MoveTo(-5.5, 7.5)
        FERSBoards["Board3"].MoveTo(-1.5, 9.5)
        FERSBoards["Board4"].MoveTo(2.5, 7.5)
        FERSBoards["Board5"].MoveTo(6.5, 7.5)
        FERSBoards["Board8"].MoveTo(-9.5, -0.5)
        FERSBoards["Board9"].MoveTo(-5.5, -0.5)
        FERSBoards["Board10"].MoveTo(-1.5, -2.5)
        FERSBoards["Board11"].MoveTo(2.5, -0.5)
        # FERSBoards["Board12"].MoveTo(6.5, -0.5)

    elif run >= 583 and run < 685:
        # include 3mm FERS board in 583
        FERSBoards["Board0"] = base_FERSBoard_6mm.copy(boardNo=0)
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board3"] = base_FERSBoard_6mm.copy(boardNo=3)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board8"] = base_FERSBoard_6mm.copy(boardNo=8)
        FERSBoards["Board9"] = base_FERSBoard_6mm.copy(boardNo=9)
        FERSBoards["Board10"] = base_FERSBoard_6mm.copy(boardNo=10)
        FERSBoards["Board11"] = base_FERSBoard_6mm.copy(boardNo=11)
        # FERSBoards["Board12"] = base_FERSBoard_3mm.copy(boardNo=12)

        FERSBoards["Board5"] = base_FERSBoard_3mm.copy(boardNo=5)

        FERSBoards["Board0"].MoveTo(-13.5, 3.5)
        FERSBoards["Board1"].MoveTo(-9.5, 7.5)
        FERSBoards["Board2"].MoveTo(-5.5, 7.5)
        FERSBoards["Board3"].MoveTo(-1.5, 9.5)
        FERSBoards["Board4"].MoveTo(2.5, 7.5)
        FERSBoards["Board8"].MoveTo(-9.5, -0.5)
        FERSBoards["Board9"].MoveTo(-5.5, -0.5)
        FERSBoards["Board10"].MoveTo(-1.5, -2.5)
        FERSBoards["Board11"].MoveTo(2.5, -0.5)
        # FERSBoards["Board12"].MoveTo(6.5, -0.5)
        FERSBoards["Board5"].MoveTo(-1.5, 1.875)

    elif run >= 685 and run < 895:
        FERSBoards["Board0"] = base_FERSBoard_6mm.copy(boardNo=0)
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board5"] = base_FERSBoard_6mm.copy(boardNo=5)
        FERSBoards["Board6"] = base_FERSBoard_6mm.copy(boardNo=6)

        FERSBoards["Board3"] = base_FERSBoard_3mm.copy(boardNo=3)
        FERSBoards["Board7"] = base_FERSBoard_3mm.copy(boardNo=7)

        FERSBoards["Board0"].MoveTo(-5.5, 7.5)
        FERSBoards["Board1"].MoveTo(-1.5, 9.5)
        FERSBoards["Board2"].MoveTo(2.5, 7.5)
        FERSBoards["Board4"].MoveTo(-5.5, -0.5)
        FERSBoards["Board5"].MoveTo(-1.5, -2.5)
        FERSBoards["Board6"].MoveTo(2.5, -0.5)

        FERSBoards["Board3"].MoveTo(-1.5, 1.875)
        FERSBoards["Board7"].MoveTo(0.5, 1.875)

    elif run >= 895 and run < 1100:
        FERSBoards["Board0"] = base_FERSBoard_6mm.copy(boardNo=0)
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board5"] = base_FERSBoard_6mm.copy(boardNo=5)
        FERSBoards["Board6"] = base_FERSBoard_6mm.copy(boardNo=6)
        FERSBoards["Board7"] = base_FERSBoard_6mm.copy(boardNo=7)
        FERSBoards["Board8"] = base_FERSBoard_6mm.copy(boardNo=8)
        FERSBoards["Board9"] = base_FERSBoard_6mm.copy(boardNo=9)
        FERSBoards["Board10"] = base_FERSBoard_6mm.copy(boardNo=10)
        FERSBoards["Board12"] = base_FERSBoard_6mm.copy(boardNo=12)
        FERSBoards["Board13"] = base_FERSBoard_6mm.copy(boardNo=13)

        FERSBoards["Board3"] = base_FERSBoard_3mm.copy(boardNo=3)
        FERSBoards["Board11"] = base_FERSBoard_3mm.copy(boardNo=11)

        FERSBoards["Board0"].MoveTo(-13.5, 3.5)
        FERSBoards["Board1"].MoveTo(-9.5, 7.5)
        FERSBoards["Board2"].MoveTo(-5.5, 7.5)
        FERSBoards["Board4"].MoveTo(-1.5, 9.5)
        FERSBoards["Board5"].MoveTo(2.5, 7.5)
        FERSBoards["Board6"].MoveTo(6.5, 7.5)
        FERSBoards["Board7"].MoveTo(10.5, 3.5)
        FERSBoards["Board8"].MoveTo(-9.5, -0.5)
        FERSBoards["Board9"].MoveTo(-5.5, -0.5)
        FERSBoards["Board10"].MoveTo(-1.5, -2.5)
        FERSBoards["Board12"].MoveTo(2.5, -0.5)
        FERSBoards["Board13"].MoveTo(6.5, -0.5)

        FERSBoards["Board3"].MoveTo(-1.5, 1.875)
        FERSBoards["Board11"].MoveTo(0.5, 1.875)
    elif run >= 1173:
        # test beam
        FERSBoards["Board0"] = base_FERSBoard_6mm.copy(boardNo=0)
        FERSBoards["Board1"] = base_FERSBoard_6mm.copy(boardNo=1)
        FERSBoards["Board2"] = base_FERSBoard_6mm.copy(boardNo=2)
        FERSBoards["Board3"] = base_FERSBoard_6mm.copy(boardNo=3)
        FERSBoards["Board4"] = base_FERSBoard_6mm.copy(boardNo=4)
        FERSBoards["Board5"] = base_FERSBoard_6mm.copy(boardNo=5)
        FERSBoards["Board6"] = base_FERSBoard_6mm.copy(boardNo=6)
        FERSBoards["Board9"] = base_FERSBoard_6mm.copy(boardNo=9)
        FERSBoards["Board10"] = base_FERSBoard_6mm.copy(boardNo=10)
        FERSBoards["Board11"] = base_FERSBoard_6mm.copy(boardNo=11)
        FERSBoards["Board12"] = base_FERSBoard_6mm.copy(boardNo=12)
        FERSBoards["Board13"] = base_FERSBoard_6mm.copy(boardNo=13)

        FERSBoards["Board7"] = base_FERSBoard_3mm.copy(boardNo=7)
        FERSBoards["Board8"] = base_FERSBoard_3mm.copy(boardNo=8)

        FERSBoards["Board0"].MoveTo(-13.5, 3.5)
        FERSBoards["Board1"].MoveTo(-9.5, 7.5)
        FERSBoards["Board2"].MoveTo(-5.5, 7.5)
        FERSBoards["Board3"].MoveTo(-1.5, 9.5)
        FERSBoards["Board4"].MoveTo(2.5, 7.5)
        FERSBoards["Board5"].MoveTo(6.5, 7.5)
        FERSBoards["Board6"].MoveTo(10.5, 3.5)
        FERSBoards["Board9"].MoveTo(-9.5, -0.5)
        FERSBoards["Board10"].MoveTo(-5.5, -0.5)
        FERSBoards["Board11"].MoveTo(-1.5, -2.5)
        FERSBoards["Board12"].MoveTo(2.5, -0.5)
        FERSBoards["Board13"].MoveTo(6.5, -0.5)

        FERSBoards["Board7"].MoveTo(-1.5, 1.875)
        FERSBoards["Board8"].MoveTo(0.5, 1.875)

    else:
        raise ValueError(f"Unsupported run number {run} for FERS boards.")
    return FERSBoards


def buildDRSBoards(run=316):
    """
    Build a map for ixy and DRS channels.
    DRS boards are the same for all runs for now.
    """
    base_DRSBoard_6mm = DRSBoard(boardNo=-1, is6mm=True)
    base_DRSBoard_3mm = DRSBoard(boardNo=-1, is6mm=False)
    DRSBoards = {}
    if run in scanruns:
        # no DRS boards in scan runs
        # only FERS
        return DRSBoards

    if run < 685:
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(boardNo=2)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 5)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 6)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 7)
        DRSBoards["Board2"].MoveTo(-1.5, -2.5)
        channel = DRSBoards["Board2"].GetChannelByGroupChannel(3, 4)
        channel.iTowerX = 1.5
        channel.iTowerY = -9.5
        channel.isCer = False

        DRSBoards["Board0"] = base_DRSBoard_6mm.copy(boardNo=0)
        DRSBoards["Board0"].MoveTo(-1.5, -6.5)
        channels = DRSBoards["Board0"].GetListOfChannels()
        for channel in channels:
            if channel.isCer:
                channel.iTowerY += 1
                channel.isCer = False
            else:
                channel.isCer = True
        prechannel_towerX = -999
        prechannel_towerY = -999
        prechannel_isCer = 0
        for group in [3, 2]:
            for channel in [7, 6, 5, 4, 3, 2, 1, 0]:
                channel = DRSBoards["Board0"].GetChannelByGroupChannel(
                    group, channel)
                temp_channel = channel.__copy__()
                channel.iTowerX = prechannel_towerX
                channel.iTowerY = prechannel_towerY
                channel.isCer = prechannel_isCer

                prechannel_towerX = temp_channel.iTowerX
                prechannel_towerY = temp_channel.iTowerY
                prechannel_isCer = temp_channel.isCer
        channel = DRSBoards["Board0"].GetChannelByGroupChannel(1, 7)
        channel.iTowerX = prechannel_towerX
        channel.iTowerY = prechannel_towerY
        channel.isCer = prechannel_isCer
        DRSBoards["Board0"].RemoveChannelByGroupChannel(3, 7)

    elif run >= 685 and run < 1003:
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(boardNo=1)
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(boardNo=2)

        # remove channels that are not used
        DRSBoards["Board1"].MoveTo(-1.5, 9.5)
        DRSBoards["Board2"].MoveTo(-1.5, -6.5)
        # no more cables..
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 4)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 5)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 6)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 7)
    elif run >= 1003 and run < 1100:
        # 3 DRS boards in 1003
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(boardNo=1)
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(boardNo=2)
        DRSBoards["Board3"] = base_DRSBoard_6mm.copy(boardNo=3)

        DRSBoards["Board1"].MoveTo(-1.5, 9.5)
        DRSBoards["Board2"].MoveTo(-1.5, 5.5)
        DRSBoards["Board3"].MoveTo(-1.5, -6.5)
        # all channels are connected.
        # no need to remove
        if run >= 1033:
            # two amplified channels in board 1
            channel = DRSBoards["Board1"].GetChannelByGroupChannel(0, 0)
            channel.isAmplified = True
            channel = DRSBoards["Board1"].GetChannelByGroupChannel(0, 1)
            channel.isAmplified = True
    elif run >= 1173:
        # test beam
        DRSBoards["Board0"] = base_DRSBoard_3mm.copy(boardNo=0)
        DRSBoards["Board1"] = base_DRSBoard_3mm.copy(boardNo=1)
        DRSBoards["Board2"] = base_DRSBoard_3mm.copy(boardNo=2)
        DRSBoards["Board3"] = base_DRSBoard_3mm.copy(boardNo=3)

        DRSBoards["Board0"].MoveTo(-1.5, 1.875)
        DRSBoards["Board1"].MoveTo(-1.5, -0.125)
        DRSBoards["Board2"].MoveTo(0.5, 1.875)
        DRSBoards["Board3"].MoveTo(0.5, -0.125)

        DRSBoards["Board4"] = buildDRSBoardTestBeam(boardNo=4)
        DRSBoards["Board5"] = buildDRSBoardTestBeam(boardNo=5)
        DRSBoards["Board6"] = buildDRSBoardTestBeam(boardNo=6)
    else:
        raise ValueError(f"Unsupported run number {run} for DRS boards.")
    return DRSBoards


def buildDRSBoardTestBeam(boardNo=4):
    # hacky way to build the confusing DRS mapping for the test beam
    from utils.CaloXChannel import DRSChannel, drs_map
    if boardNo == 4:
        channels_DRS = []
        for ix in range(0, 4):
            channels_DRS_one_row = []
            for iy in range(0, 8):
                if iy < 4:
                    base_iX = -5.5
                    base_iY = 3.5
                else:
                    base_iX = -1.5
                    base_iY = 5.5 + 4
                channelNo = drs_map[ix, iy]
                groupNo = (channelNo // 8)
                chanNo = (channelNo % 8)
                isCer = True
                channel = DRSChannel(base_iX + ix, base_iY - iy, ix, iy, isCer,
                                     chanNo, groupNo, boardNo, is6mm=True, isAmplified=False)
                channels_DRS_one_row.append(channel)
            channels_DRS.append(channels_DRS_one_row)
            drsboard = DRSBoard(boardNo=boardNo, channels=channels_DRS)
            drsboard.RemoveChannelByGroupChannel(3, 7)
            drsboard.RemoveChannelByGroupChannel(3, 6)
        return drsboard
    elif boardNo == 5:
        channels_DRS = []
        channels_DRS_one_row = []
        # two remnant channels
        channel = DRSChannel(-0.5, 2.5, 0, 0, True, 0, 0,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channel = DRSChannel(-1.5, 2.5, 0, 1, True, 1, 0,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channels_DRS.append(channels_DRS_one_row)

        for ix in range(0, 4):
            channels_DRS_one_row = []
            for iy in range(0, 4):
                base_iX = 2.5
                base_iY = 3.5
                channelNo = drs_map[ix, iy] + 2
                groupNo = (channelNo // 8)
                chanNo = (channelNo % 8)
                isCer = True
                channel = DRSChannel(base_iX + ix, base_iY - iy, ix, iy, isCer,
                                     chanNo, groupNo, boardNo, is6mm=True, isAmplified=False)
                channels_DRS_one_row.append(channel)
            channels_DRS.append(channels_DRS_one_row)

        for ix in range(0, 4):
            channels_DRS_one_row = []
            for iy in range(4, 6):
                base_iX = -3.5
                if iy == 4:
                    base_iY = -1.5
                else:
                    base_iY = -3.5
                channelNo = drs_map[ix, iy] + 2
                groupNo = (channelNo // 8)
                chanNo = (channelNo % 8)
                isCer = True
                channel = DRSChannel(base_iX + (ix % 2), base_iY + ix//2, ix, iy, isCer,
                                     chanNo, groupNo, boardNo, is6mm=True, isAmplified=False)
                channels_DRS_one_row.append(channel)
            channels_DRS.append(channels_DRS_one_row)

        channels_DRS_one_row = []
        channel = DRSChannel(1.5, -2.5, 0, 0, True, 2, 3,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channel = DRSChannel(0.5, -2.5, 0, 1, True, 3, 3,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channel = DRSChannel(-0.5, -2.5, 0, 2, True, 4, 3,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channel = DRSChannel(-1.5, -2.5, 0, 3, True, 5, 3,
                             boardNo, is6mm=True, isAmplified=False)
        channels_DRS_one_row.append(channel)
        channels_DRS.append(channels_DRS_one_row)

        drsboard = DRSBoard(boardNo=boardNo, channels=channels_DRS)

        return drsboard

    elif boardNo == 6:
        base_iX_R = -1.5
        base_iY_R = -2.5
        base_iX_T = 2.5
        base_iY_T = -0.5
        maps_board6 = {
            (base_iX_R + 3, base_iY_R - 1): (0, 2),
            (base_iX_R + 2, base_iY_R - 1): (0, 3),
            (base_iX_R + 1, base_iY_R - 1): (0, 4),
            (base_iX_R + 0, base_iY_R - 1): (0, 5),
            (base_iX_R + 3, base_iY_R - 2): (0, 6),
            (base_iX_R + 2, base_iY_R - 2): (0, 7),
            (base_iX_R + 1, base_iY_R - 2): (2, 2),
            (base_iX_R + 0, base_iY_R - 2): (2, 3),
            (base_iX_R + 3, base_iY_R - 3): (1, 0),
            (base_iX_R + 2, base_iY_R - 3): (1, 1),
            (base_iX_R + 1, base_iY_R - 3): (2, 4),
            (base_iX_R + 0, base_iY_R - 3): (2, 5),

            (base_iX_T + 1, base_iY_T): (1, 2),
            (base_iX_T + 0, base_iY_T): (1, 3),
            (base_iX_T + 1, base_iY_T - 1): (1, 4),
            (base_iX_T + 0, base_iY_T - 1): (1, 5),
            (base_iX_T + 1, base_iY_T - 2): (1, 6),
            (base_iX_T + 0, base_iY_T - 2): (1, 7),
            (base_iX_T + 1, base_iY_T - 3): (2, 0),
            (base_iX_T + 0, base_iY_T - 3): (2, 1),
        }
        channels_DRS = []
        for iX in range(0, 4):
            channels_DRS_one_row = []
            channels_DRS.append(channels_DRS_one_row)

        for (iTowerX, iTowerY), (groupNo, chanNo) in maps_board6.items():
            isCer = True
            channel = DRSChannel(iTowerX, iTowerY, groupNo, chanNo, isCer,
                                 chanNo, groupNo, boardNo, is6mm=True, isAmplified=False)
            channels_DRS[groupNo].append(channel)

        drsboard = DRSBoard(boardNo=boardNo, channels=channels_DRS)
        return drsboard


def buildTimeReferenceChannels(run=316):
    """
    Returns a list of time reference channels.
    """
    time_reference_channels = []
    if run in scanruns:
        # no time reference channels in scan runs
        # since no drs boards
        return time_reference_channels

    if run < 685:
        time_reference_channels.append("DRS_Board0_Group3_Channel7")
        time_reference_channels.append("DRS_Board2_Group3_Channel7")
        time_reference_channels.append("DRS_Board1_Group0_Channel0")
    elif run >= 685:
        time_reference_channels.append("DRS_Board1_Group0_Channel8")
        time_reference_channels.append("DRS_Board1_Group1_Channel8")
        time_reference_channels.append("DRS_Board1_Group2_Channel8")
        time_reference_channels.append("DRS_Board1_Group3_Channel8")
        time_reference_channels.append("DRS_Board2_Group0_Channel8")
        time_reference_channels.append("DRS_Board2_Group1_Channel8")
        time_reference_channels.append("DRS_Board2_Group2_Channel8")
        time_reference_channels.append("DRS_Board2_Group3_Channel8")
        time_reference_channels.append("DRS_Board0_Group0_Channel8")
    else:
        raise ValueError(
            f"Unsupported run number {run} for time reference channels.")

    return time_reference_channels


def buildHodoTriggerChannels(run=316):
    """
    Returns a list of hodoscope trigger channels.
    """
    hodo_trigger_channels = []
    if run in scanruns:
        # no hodoscope trigger channels in scan runs
        # since no drs boards
        return hodo_trigger_channels

    if run < 685:
        hodo_trigger_channels.append("DRS_Board1_Group2_Channel0")
        hodo_trigger_channels.append("DRS_Board1_Group2_Channel1")
    elif run >= 685:
        hodo_trigger_channels.append("DRS_Board0_Group2_Channel0")
        hodo_trigger_channels.append("DRS_Board0_Group2_Channel1")
    else:
        raise ValueError(
            f"Unsupported run number {run} for hodoscope trigger channels.")

    return hodo_trigger_channels


def buildHodoPosChannels(run=316):
    """
    Returns a dictionary containing the hodoscope channels for the position measurements
    """
    hodoscope_channels = {}
    if run in scanruns:
        # no hodoscope position channels in scan runs
        # since no drs boards
        return hodoscope_channels

    hodoscope_channels["TopX"] = [
        "DRS_Board1_Group0_Channel1",
        "DRS_Board1_Group0_Channel2",
    ]
    hodoscope_channels["TopZ"] = [
        "DRS_Board1_Group0_Channel3",
        "DRS_Board1_Group0_Channel4",
    ]
    hodoscope_channels["BottomX"] = [
        "DRS_Board1_Group0_Channel5",
        "DRS_Board1_Group0_Channel6",
    ]
    if run < 583:
        # bottom z seems to have flipped left and right
        hodoscope_channels["BottomZ"] = [
            "DRS_Board1_Group1_Channel0",
            "DRS_Board1_Group0_Channel7",
        ]
    elif run >= 583:
        hodoscope_channels["BottomZ"] = [
            "DRS_Board1_Group0_Channel7",
            "DRS_Board1_Group1_Channel0",
        ]

    if run >= 685 and run < 1170:
        # For runs >= 685, board 0 is used for hodoscope position channels
        hodoscope_channels["TopX"] = [
            "DRS_Board0_Group0_Channel0",
            "DRS_Board0_Group0_Channel1",
        ]
        hodoscope_channels["TopZ"] = [
            "DRS_Board0_Group0_Channel2",
            "DRS_Board0_Group0_Channel3",
        ]
        hodoscope_channels["BottomX"] = [
            "DRS_Board0_Group0_Channel4",
            "DRS_Board0_Group0_Channel5",
        ]
        hodoscope_channels["BottomZ"] = [
            "DRS_Board0_Group0_Channel6",
            "DRS_Board0_Group0_Channel7",
        ]

    if run >= 1170:
        # test beam
        # Left and right; top and bottom for Hodoscope2
        # 1170 is not exactly the first run number
        hodoscope_channels = {}
        hodoscope_channels["LR"] = [
            "DRS_Board7_Group0_Channel4",
            "DRS_Board7_Group0_Channel5",
        ]
        hodoscope_channels["UD"] = [
            "DRS_Board7_Group0_Channel6",
            "DRS_Board7_Group0_Channel7",
        ]

    return hodoscope_channels


def findTimeReferenceDelay(channel, run=1040):
    if str(run) not in triggerdelay.keys:
        return triggerdelay["default"][channel]
    else:
        return triggerdelay[str(run)][channel]


def getUpstreamVetoChannel(run=1184):
    if run < 1183:
        return None
    else:
        return "DRS_Board7_Group1_Channel6"


def getDownStreamMuonChannel(run=1184):
    if run < 1183:
        return None
    else:
        return "DRS_Board7_Group1_Channel0"


def getPreShowerChannel(run=1184):
    """
    Returns the pre-shower channel.
    """
    if run < 1183:
        return None
    else:
        return "DRS_Board7_Group1_Channel1"


def getCerenkovCounters(run=1184):
    """
    Returns a list of Cerenkov counter channels.
    """
    if run < 1183:
        return []
    else:
        return [
            "DRS_Board7_Group1_Channel2",
            "DRS_Board7_Group1_Channel3",
            "DRS_Board7_Group1_Channel4",
            "DRS_Board7_Group1_Channel5",
        ]


def getServiceDRSChannels(run=1184):
    """
    Returns a list of service DRS channels.
    """
    if run < 1183:
        return []
    elif run < 1260:
        return [
            "DRS_Board7_Group0_Channel0",
            "DRS_Board7_Group0_Channel1",
            "DRS_Board7_Group0_Channel2",
            "DRS_Board7_Group0_Channel3",
            "DRS_Board7_Group0_Channel4",
            "DRS_Board7_Group0_Channel5",
            "DRS_Board7_Group0_Channel6",
            "DRS_Board7_Group0_Channel7",
            "DRS_Board7_Group1_Channel0",
            "DRS_Board7_Group1_Channel1",
            "DRS_Board7_Group1_Channel2",
            "DRS_Board7_Group1_Channel3",
            "DRS_Board7_Group1_Channel4",
            "DRS_Board7_Group1_Channel5",
            "DRS_Board7_Group1_Channel6",
            "DRS_Board7_Group1_Channel7",
        ]
    else:
        return [
            "DRS_Board7_Group0_Channel0",
            "DRS_Board7_Group0_Channel1",
            "DRS_Board7_Group0_Channel2",
            "DRS_Board7_Group0_Channel3",
            "DRS_Board7_Group0_Channel4",
            "DRS_Board7_Group0_Channel5",
            "DRS_Board7_Group0_Channel6",
            "DRS_Board7_Group0_Channel7",
            "DRS_Board7_Group1_Channel0",
            "DRS_Board7_Group1_Channel1",
            "DRS_Board7_Group1_Channel2",
            "DRS_Board7_Group1_Channel3",
            "DRS_Board7_Group1_Channel4",
            "DRS_Board7_Group1_Channel5",
            "DRS_Board7_Group1_Channel6",
            "DRS_Board7_Group1_Channel7",
            "DRS_Board7_Group2_Channel0",
            "DRS_Board7_Group2_Channel1",
        ]


if __name__ == "__main__":
    # Example usage
    run_number = 583
    fers_boards = buildFERSBoards(run=run_number)
    drs_boards = buildDRSBoards(run=run_number)

    print("FERS Boards:")
    for board_name, board in fers_boards.items():
        print(f"{board_name}: {board}")

    print("\nDRS Boards:")
    for board_name, board in drs_boards.items():
        print(f"{board_name}: {board}")

    print("\nHodoscope Position Channels:")
    hodo_channels = buildHodoPosChannels(run=run_number)
    for hodo_type, channels in hodo_channels.items():
        print(f"{hodo_type}: {channels}")
