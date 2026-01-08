from channels.calox_channel import FERSBoard, DRSBoard, DRSChannel, drs_map, FERSBoards
from utils.data_loader import is_scan_run
import json
from collections import OrderedDict

f_triggerdelay = "data/triggerdelay.json"
with open(f_triggerdelay, 'r') as f:
    triggerdelay = json.load(f)

f_drstriggermap = "data/drstriggermap.json"
with open(f_drstriggermap, 'r') as f:
    triggermap = json.load(f)


def build_fers_boards(run=316):
    """
    Build a map for ixy and FERS channels for both boards.
    """
    base_FERSBoard_6mm = FERSBoard(board_no=-1, is6mm=True)
    base_FERSBoard_3mm = FERSBoard(board_no=-1, is6mm=False)
    fersboards = FERSBoards()
    if run == 316:
        # 5 FERS board in 316
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)

        # positions
        fersboards["Board1"].move_to(6.5, -0.5)
        fersboards["Board2"].move_to(2.5, -0.5)
        fersboards["Board3"].move_to(-1.5, -2.5)
        fersboards["Board4"].move_to(-5.5, -0.5)
        fersboards["Board5"].move_to(-9.5, -0.5)

    elif run == 571:
        # 10 FERS boards in 571
        fersboards["Board0"] = base_FERSBoard_6mm.copy(board_no=0)
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board8"] = base_FERSBoard_6mm.copy(board_no=8)
        fersboards["Board9"] = base_FERSBoard_6mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_6mm.copy(board_no=10)
        fersboards["Board11"] = base_FERSBoard_6mm.copy(board_no=11)
        # fersboards["Board12"] = base_FERSBoard_6mm.copy(boardNo=12)

        fersboards["Board0"].move_to(-13.5, 3.5)
        fersboards["Board1"].move_to(-9.5, 7.5)
        fersboards["Board2"].move_to(-5.5, 7.5)
        fersboards["Board3"].move_to(-1.5, 9.5)
        fersboards["Board4"].move_to(2.5, 7.5)
        fersboards["Board5"].move_to(6.5, 7.5)
        fersboards["Board8"].move_to(-9.5, -0.5)
        fersboards["Board9"].move_to(-5.5, -0.5)
        fersboards["Board10"].move_to(-1.5, -2.5)
        fersboards["Board11"].move_to(2.5, -0.5)
        # fersboards["Board12"].move_to(6.5, -0.5)

    elif run >= 583 and run < 685:
        # include 3mm FERS board in 583
        fersboards["Board0"] = base_FERSBoard_6mm.copy(board_no=0)
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board8"] = base_FERSBoard_6mm.copy(board_no=8)
        fersboards["Board9"] = base_FERSBoard_6mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_6mm.copy(board_no=10)
        fersboards["Board11"] = base_FERSBoard_6mm.copy(board_no=11)
        # fersboards["Board12"] = base_FERSBoard_3mm.copy(boardNo=12)

        fersboards["Board5"] = base_FERSBoard_3mm.copy(board_no=5)

        fersboards["Board0"].move_to(-13.5, 3.5)
        fersboards["Board1"].move_to(-9.5, 7.5)
        fersboards["Board2"].move_to(-5.5, 7.5)
        fersboards["Board3"].move_to(-1.5, 9.5)
        fersboards["Board4"].move_to(2.5, 7.5)
        fersboards["Board8"].move_to(-9.5, -0.5)
        fersboards["Board9"].move_to(-5.5, -0.5)
        fersboards["Board10"].move_to(-1.5, -2.5)
        fersboards["Board11"].move_to(2.5, -0.5)
        # fersboards["Board12"].move_to(6.5, -0.5)
        fersboards["Board5"].move_to(-1.5, 1.875)

    elif run >= 685 and run < 895:
        fersboards["Board0"] = base_FERSBoard_6mm.copy(board_no=0)
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board6"] = base_FERSBoard_6mm.copy(board_no=6)

        fersboards["Board3"] = base_FERSBoard_3mm.copy(board_no=3)
        fersboards["Board7"] = base_FERSBoard_3mm.copy(board_no=7)

        fersboards["Board0"].move_to(-5.5, 7.5)
        fersboards["Board1"].move_to(-1.5, 9.5)
        fersboards["Board2"].move_to(2.5, 7.5)
        fersboards["Board4"].move_to(-5.5, -0.5)
        fersboards["Board5"].move_to(-1.5, -2.5)
        fersboards["Board6"].move_to(2.5, -0.5)

        fersboards["Board3"].move_to(-1.5, 1.875)
        fersboards["Board7"].move_to(0.5, 1.875)

    elif run >= 895 and run < 1100:
        fersboards["Board0"] = base_FERSBoard_6mm.copy(board_no=0)
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board6"] = base_FERSBoard_6mm.copy(board_no=6)
        fersboards["Board7"] = base_FERSBoard_6mm.copy(board_no=7)
        fersboards["Board8"] = base_FERSBoard_6mm.copy(board_no=8)
        fersboards["Board9"] = base_FERSBoard_6mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_6mm.copy(board_no=10)
        fersboards["Board12"] = base_FERSBoard_6mm.copy(board_no=12)
        fersboards["Board13"] = base_FERSBoard_6mm.copy(board_no=13)

        fersboards["Board3"] = base_FERSBoard_3mm.copy(board_no=3)
        fersboards["Board11"] = base_FERSBoard_3mm.copy(board_no=11)

        fersboards["Board0"].move_to(-13.5, 3.5)
        fersboards["Board1"].move_to(-9.5, 7.5)
        fersboards["Board2"].move_to(-5.5, 7.5)
        fersboards["Board4"].move_to(-1.5, 9.5)
        fersboards["Board5"].move_to(2.5, 7.5)
        fersboards["Board6"].move_to(6.5, 7.5)
        fersboards["Board7"].move_to(10.5, 3.5)
        fersboards["Board8"].move_to(-9.5, -0.5)
        fersboards["Board9"].move_to(-5.5, -0.5)
        fersboards["Board10"].move_to(-1.5, -2.5)
        fersboards["Board12"].move_to(2.5, -0.5)
        fersboards["Board13"].move_to(6.5, -0.5)

        fersboards["Board3"].move_to(-1.5, 1.875)
        fersboards["Board11"].move_to(0.5, 1.875)
    elif run >= 1173 and run < 1327:
        # test beam
        fersboards["Board0"] = base_FERSBoard_6mm.copy(board_no=0)
        fersboards["Board1"] = base_FERSBoard_6mm.copy(board_no=1)
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board6"] = base_FERSBoard_6mm.copy(board_no=6)
        fersboards["Board9"] = base_FERSBoard_6mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_6mm.copy(board_no=10)
        fersboards["Board11"] = base_FERSBoard_6mm.copy(board_no=11)
        fersboards["Board12"] = base_FERSBoard_6mm.copy(board_no=12)
        fersboards["Board13"] = base_FERSBoard_6mm.copy(board_no=13)

        fersboards["Board7"] = base_FERSBoard_3mm.copy(board_no=7)
        fersboards["Board8"] = base_FERSBoard_3mm.copy(board_no=8)

        fersboards["Board0"].move_to(-13.5, 3.5)
        fersboards["Board1"].move_to(-9.5, 7.5)
        fersboards["Board2"].move_to(-5.5, 7.5)
        fersboards["Board3"].move_to(-1.5, 9.5)
        fersboards["Board4"].move_to(2.5, 7.5)
        fersboards["Board5"].move_to(6.5, 7.5)
        fersboards["Board6"].move_to(10.5, 3.5)
        fersboards["Board9"].move_to(-9.5, -0.5)
        fersboards["Board10"].move_to(-5.5, -0.5)
        fersboards["Board11"].move_to(-1.5, -2.5)
        fersboards["Board12"].move_to(2.5, -0.5)
        fersboards["Board13"].move_to(6.5, -0.5)

        fersboards["Board7"].move_to(-1.5, 1.875)
        fersboards["Board8"].move_to(0.5, 1.875)

    elif run >= 1342:
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board6"] = base_FERSBoard_6mm.copy(board_no=6)
        fersboards["Board7"] = base_FERSBoard_6mm.copy(board_no=7)
        fersboards["Board8"] = base_FERSBoard_6mm.copy(board_no=8)
        fersboards["Board11"] = base_FERSBoard_6mm.copy(board_no=11)
        fersboards["Board12"] = base_FERSBoard_6mm.copy(board_no=12)
        fersboards["Board13"] = base_FERSBoard_6mm.copy(board_no=13)
        fersboards["Board14"] = base_FERSBoard_6mm.copy(board_no=14)
        fersboards["Board15"] = base_FERSBoard_6mm.copy(board_no=15)

        fersboards["Board9"] = base_FERSBoard_3mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_3mm.copy(board_no=10)

        fersboards["Board2"].move_to(-13.5, 3.5)
        fersboards["Board3"].move_to(-9.5, 7.5)
        fersboards["Board4"].move_to(-5.5, 7.5)
        fersboards["Board5"].move_to(-1.5, 9.5)
        fersboards["Board6"].move_to(2.5, 7.5)
        fersboards["Board7"].move_to(6.5, 7.5)
        fersboards["Board8"].move_to(10.5, 3.5)
        fersboards["Board11"].move_to(-9.5, -0.5)
        fersboards["Board12"].move_to(-5.5, -0.5)
        fersboards["Board13"].move_to(-1.5, -2.5)
        fersboards["Board14"].move_to(2.5, -0.5)
        fersboards["Board15"].move_to(6.5, -0.5)

        # 3mm
        fersboards["Board9"].move_to(-1.5, 1.875)
        fersboards["Board10"].move_to(0.5, 1.875)

    else:
        raise ValueError(f"Unsupported run number {run} for FERS boards.")
    update_quartz_channels(fersboards)
    return fersboards


def build_drs_boards(run=316):
    """
    Build a map for ixy and DRS channels.
    DRS boards are the same for all runs for now.
    """
    base_DRSBoard_6mm = DRSBoard(board_no=-1, is6mm=True)
    base_DRSBoard_3mm = DRSBoard(board_no=-1, is6mm=False)
    DRSBoards = {}
    if is_scan_run(run):
        # no DRS boards in scan runs
        # only FERS
        return DRSBoards

    if run < 685:
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(board_no=2)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 5)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 7)
        DRSBoards["Board2"].move_to(-1.5, -2.5)
        channel = DRSBoards["Board2"].get_channel_by_group_channel(3, 4)
        channel.i_tower_x = 1.5
        channel.i_tower_y = -9.5
        channel.isCer = False

        DRSBoards["Board0"] = base_DRSBoard_6mm.copy(board_no=0)
        DRSBoards["Board0"].move_to(-1.5, -6.5)
        channels = DRSBoards["Board0"].get_list_of_channels()
        for channel in channels:
            if channel.isCer:
                channel.i_tower_y += 1
                channel.isCer = False
            else:
                channel.isCer = True
        prechannel_towerX = -999
        prechannel_towerY = -999
        prechannel_isCer = 0
        for group in [3, 2]:
            for channel in [7, 6, 5, 4, 3, 2, 1, 0]:
                channel = DRSBoards["Board0"].get_channel_by_group_channel(
                    group, channel)
                temp_channel = channel.__copy__()
                channel.i_tower_x = prechannel_towerX
                channel.i_tower_y = prechannel_towerY
                channel.isCer = prechannel_isCer

                prechannel_towerX = temp_channel.i_tower_x
                prechannel_towerY = temp_channel.i_tower_y
                prechannel_isCer = temp_channel.isCer
        channel = DRSBoards["Board0"].get_channel_by_group_channel(1, 7)
        channel.i_tower_x = prechannel_towerX
        channel.i_tower_y = prechannel_towerY
        channel.isCer = prechannel_isCer
        DRSBoards["Board0"].remove_channel_by_group_channel(3, 7)

    elif run >= 685 and run < 1003:
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(board_no=2)

        # remove channels that are not used
        DRSBoards["Board1"].move_to(-1.5, 9.5)
        DRSBoards["Board2"].move_to(-1.5, -6.5)
        # no more cables..
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 4)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 5)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 7)
    elif run >= 1003 and run < 1100:
        # 3 DRS boards in 1003
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(board_no=2)
        DRSBoards["Board3"] = base_DRSBoard_6mm.copy(board_no=3)

        DRSBoards["Board1"].move_to(-1.5, 9.5)
        DRSBoards["Board2"].move_to(-1.5, 5.5)
        DRSBoards["Board3"].move_to(-1.5, -6.5)
        # all channels are connected.
        # no need to remove
        if run >= 1033:
            # two amplified channels in board 1
            channel = DRSBoards["Board1"].get_channel_by_group_channel(0, 0)
            channel.is_amplified = True
            channel = DRSBoards["Board1"].get_channel_by_group_channel(0, 1)
            channel.is_amplified = True
    elif run >= 1173 and run < 1327:
        # test beam
        DRSBoards["Board0"] = base_DRSBoard_3mm.copy(board_no=0)
        DRSBoards["Board1"] = base_DRSBoard_3mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_3mm.copy(board_no=2)
        DRSBoards["Board3"] = base_DRSBoard_3mm.copy(board_no=3)

        DRSBoards["Board0"].move_to(-1.5, 1.875)
        DRSBoards["Board1"].move_to(-1.5, -0.125)
        DRSBoards["Board2"].move_to(0.5, 1.875)
        DRSBoards["Board3"].move_to(0.5, -0.125)

        DRSBoards["Board4"] = buildDRSBoardTestBeam(board_no=4)
        DRSBoards["Board5"] = buildDRSBoardTestBeam(board_no=5)
        DRSBoards["Board6"] = buildDRSBoardTestBeam(board_no=6)
    elif run >= 1342:
        # september 2nd test beam channel maps
        DRSBoards["Board0"] = base_DRSBoard_3mm.copy(board_no=0)
        DRSBoards["Board1"] = base_DRSBoard_3mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_3mm.copy(board_no=2)
        DRSBoards["Board3"] = base_DRSBoard_3mm.copy(board_no=3)

        DRSBoards["Board0"].move_to(-1.5, 1.875)
        DRSBoards["Board1"].move_to(-1.5, -0.125)
        DRSBoards["Board2"].move_to(0.5, 1.875)
        DRSBoards["Board3"].move_to(0.5, -0.125)
        # channels for MCP
        DRSBoards["Board0"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board0"].remove_channel_by_group_channel(3, 7)
        DRSBoards["Board1"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board1"].remove_channel_by_group_channel(3, 7)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board2"].remove_channel_by_group_channel(3, 7)
        DRSBoards["Board3"].remove_channel_by_group_channel(3, 6)
        DRSBoards["Board3"].remove_channel_by_group_channel(3, 7)

        DRSBoards["Board4"] = buildDRSBoardTestBeamSep(board_no=4)
        DRSBoards["Board5"] = buildDRSBoardTestBeamSep(board_no=5)
        DRSBoards["Board6"] = buildDRSBoardTestBeamSep(board_no=6)

    else:
        raise ValueError(f"Unsupported run number {run} for DRS boards.")
    update_quartz_channels(DRSBoards)
    return DRSBoards


def buildDRSBoardTestBeam(board_no=4):
    # hacky way to build the confusing DRS mapping for the test beam
    if board_no == 4:
        channels_DRS = []
        for ix in range(0, 4):
            for iy in range(0, 8):
                if iy < 4:
                    base_iX = -5.5
                    base_iY = 3.5
                else:
                    base_iX = -1.5
                    base_iY = 5.5 + 4
                channel_no = drs_map[ix, iy]
                group_no = (channel_no // 8)
                chan_no = (channel_no % 8)
                isCer = True
                channel = DRSChannel(base_iX + ix, base_iY - iy, isCer,
                                     chan_no, group_no, board_no, is6mm=True, is_amplified=False)
                channels_DRS.append(channel)
            drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)
            drsboard.remove_channel_by_group_channel(3, 7)
            drsboard.remove_channel_by_group_channel(3, 6)
        return drsboard
    elif board_no == 5:
        channels_DRS = []
        # two remnant channels
        channel = DRSChannel(-0.5, 2.5, True, 0, 0,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)
        channel = DRSChannel(-1.5, 2.5, True, 1, 0,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)

        for ix in range(0, 4):
            for iy in range(0, 4):
                base_iX = 2.5
                base_iY = 3.5
                channel_no = drs_map[ix, iy] + 2
                group_no = (channel_no // 8)
                chan_no = (channel_no % 8)
                isCer = True
                channel = DRSChannel(base_iX + ix, base_iY - iy, isCer,
                                     chan_no, group_no, board_no, is6mm=True, is_amplified=False)
                channels_DRS.append(channel)

        for ix in range(0, 4):
            for iy in range(4, 6):
                base_iX = -3.5
                if iy == 4:
                    base_iY = -1.5
                else:
                    base_iY = -3.5
                channel_no = drs_map[ix, iy] + 2
                group_no = (channel_no // 8)
                chan_no = (channel_no % 8)
                isCer = True
                channel = DRSChannel(base_iX + (ix % 2), base_iY + ix//2, isCer,
                                     chan_no, group_no, board_no, is6mm=True, is_amplified=False)
                channels_DRS.append(channel)

        channel = DRSChannel(1.5, -2.5, True, 2, 3,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)
        channel = DRSChannel(0.5, -2.5, True, 3, 3,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)
        channel = DRSChannel(-0.5, -2.5, True, 4, 3,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)
        channel = DRSChannel(-1.5, -2.5, True, 5, 3,
                             board_no, is6mm=True, is_amplified=False)
        channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)

        return drsboard

    elif board_no == 6:
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

        for (i_tower_x, i_tower_y), (group_no, chan_no) in maps_board6.items():
            isCer = True
            channel = DRSChannel(i_tower_x, i_tower_y, isCer, chan_no,
                                 group_no, board_no, is6mm=True, is_amplified=False)
            channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)
        return drsboard


def buildDRSBoardTestBeamSep(board_no=4):
    if board_no == 4:
        # top left
        base_iX = -1.5
        base_iY = -2.5
        maps_board4 = {
            (base_iX + 3, base_iY): (0, 0),
            (base_iX + 2, base_iY): (0, 1),
            (base_iX + 1, base_iY): (0, 2),
            (base_iX + 0, base_iY): (0, 3),
            (base_iX + 3, base_iY - 1): (0, 4),
            (base_iX + 2, base_iY - 1): (0, 5),
            (base_iX + 1, base_iY - 1): (0, 6),
            (base_iX + 0, base_iY - 1): (0, 7),
            (base_iX + 3, base_iY - 2): (1, 0),
            (base_iX + 2, base_iY - 2): (1, 1),
            (base_iX + 1, base_iY - 2): (1, 2),
            (base_iX + 0, base_iY - 2): (1, 3),
            (base_iX + 3, base_iY - 3): (1, 4),
            (base_iX + 2, base_iY - 3): (1, 5),
            (base_iX + 1, base_iY - 3): (1, 6),
            (base_iX + 0, base_iY - 3): (1, 7),
            # two channels for Sci
            (base_iX + 3, base_iY - 4): (2, 2),
            (base_iX + 2, base_iY - 4): (2, 3),
            (base_iX + 1, base_iY - 4): (2, 4),
            (base_iX + 0, base_iY - 4): (2, 5),
            # only two Cer channels for this row
            (base_iX + 2, base_iY - 5): (2, 6),
            (base_iX + 1, base_iY - 5): (2, 7),
            (base_iX + 3, base_iY - 6): (3, 0),
            (base_iX + 2, base_iY - 6): (3, 1),
            (base_iX + 1, base_iY - 6): (3, 2),
            (base_iX + 0, base_iY - 6): (3, 3),
            # two remnant 6mm cer channels
            (2.5, -3.5): (3, 6),
            (-2.5, -3.5): (3, 7),
        }

        channels_DRS = []

        for (i_tower_x, i_tower_y), (group_no, chan_no) in maps_board4.items():
            isCer = True
            channel = DRSChannel(i_tower_x, i_tower_y, isCer, chan_no,
                                 group_no, board_no, is6mm=True, is_amplified=True)
            channels_DRS.append(channel)

        # add four channels for Sci
        channel = DRSChannel(0.5, -5.5, False, 0, 2,
                             board_no, is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(-0.5, -5.5, False, 1, 2, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(0.5, -8.5, False, 4, 3, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(-0.5, -8.5, False, 5, 3, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)

    elif board_no == 6:
        # top middle
        base_iX = -1.5
        base_iY = 9.5
        maps_board6 = {
            (base_iX + 3, base_iY - 1): (0, 0),
            (base_iX + 2, base_iY - 1): (0, 1),
            (base_iX + 1, base_iY - 1): (0, 2),
            (base_iX + 0, base_iY - 1): (0, 3),
            # leave two channels for Sci
            # only two Cer channels for this row
            (base_iX + 2, base_iY - 2): (0, 6),
            (base_iX + 1, base_iY - 2): (0, 7),
            (base_iX + 3, base_iY - 3): (1, 0),
            (base_iX + 2, base_iY - 3): (1, 1),
            (base_iX + 1, base_iY - 3): (1, 2),
            (base_iX + 0, base_iY - 3): (1, 3),
            (base_iX + 3, base_iY - 4): (1, 4),
            (base_iX + 2, base_iY - 4): (1, 5),
            (base_iX + 1, base_iY - 4): (1, 6),
            (base_iX + 0, base_iY - 4): (1, 7),
            # leave two channels for Sci
            (base_iX + 3, base_iY - 5): (2, 2),
            (base_iX + 2, base_iY - 5): (2, 3),
            (base_iX + 1, base_iY - 5): (2, 4),
            (base_iX + 0, base_iY - 5): (2, 5),
            (base_iX + 3, base_iY - 6): (2, 6),
            (base_iX + 2, base_iY - 6): (2, 7),
            (base_iX + 1, base_iY - 6): (3, 0),
            (base_iX + 0, base_iY - 6): (3, 1),
            (base_iX + 3, base_iY - 7): (3, 2),
            (base_iX + 2, base_iY - 7): (3, 3),
            (base_iX + 1, base_iY - 7): (3, 4),
            (base_iX + 0, base_iY - 7): (3, 5),
            # two remnant 6mm cer channels
            (2.5, 3.5): (3, 6),
            (-2.5, 3.5): (3, 7),
        }

        channels_DRS = []

        for (i_tower_x, i_tower_y), (group_no, chan_no) in maps_board6.items():
            isCer = True
            channel = DRSChannel(i_tower_x, i_tower_y, isCer, chan_no,
                                 group_no, board_no, is6mm=True, is_amplified=True)
            channels_DRS.append(channel)

        # add four channels for Sci
        channel = DRSChannel(0.5, 8.5, False, 4, 0, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(-0.5, 8.5, False, 5, 0, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(0.5, 5.5, False, 0, 2, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)
        channel = DRSChannel(-0.5, 5.5, False, 1, 2, board_no,
                             is6mm=True, is_amplified=True)
        channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)

    elif board_no == 5:
        maps_board5 = {
            # iTowerX, iTowerY: (groupNo, chanNo)
            # right
            (3.5, 2.5): (0, 0),
            (2.5, 2.5): (0, 1),
            (3.5, 1.5): (0, 2),
            (2.5, 1.5): (0, 3),
            (3.5, 0.5): (0, 4),
            (2.5, 0.5): (0, 5),
            (3.5, -0.5): (0, 6),
            (2.5, -0.5): (0, 7),
            (3.5, -1.5): (1, 0),
            (2.5, -1.5): (1, 1),
            (3.5, -2.5): (1, 2),
            (2.5, -2.5): (1, 3),
            # left
            (-2.5, 2.5): (1, 4),
            (-3.5, 2.5): (1, 5),
            (-2.5, 1.5): (1, 6),
            (-3.5, 1.5): (1, 7),
            (-2.5, 0.5): (2, 0),
            (-3.5, 0.5): (2, 1),
            (-2.5, -0.5): (2, 2),
            (-3.5, -0.5): (2, 3),
            (-2.5, -1.5): (2, 4),
            (-3.5, -1.5): (2, 5),
            (-2.5, -2.5): (2, 6),
            (-3.5, -2.5): (2, 7),
        }
        channels_DRS = []

        for (i_tower_x, i_tower_y), (group_no, chan_no) in maps_board5.items():
            isCer = True
            channel = DRSChannel(i_tower_x, i_tower_y, isCer, chan_no,
                                 group_no, board_no, is6mm=True, is_amplified=True)
            channels_DRS.append(channel)

        # 8 channels for 3mm Cer and Sci
        # NOTE
        # not sure if isCer starts with True or False
        # !!!!!
        isCer = True
        maps_board5_3mm = [
            ((-1.5, 0.125), (3, 0, isCer)),
            ((-1.5, 0.125), (3, 1, not isCer)),
            ((-1.5, -1.875), (3, 2, isCer)),
            ((-1.5, -1.875), (3, 3, not isCer)),
            ((0.5, 0.125), (3, 4, isCer)),
            ((0.5, 0.125), (3, 5, not isCer)),
            ((0.5, -1.875), (3, 6, isCer)),
            ((0.5, -1.875), (3, 7, not isCer)),
        ]
        for (i_tower_x, i_tower_y), (group_no, chan_no, isCer) in maps_board5_3mm:
            channel = DRSChannel(i_tower_x, i_tower_y, isCer, chan_no,
                                 group_no, board_no, is6mm=False, is_amplified=True, isQuartz=isCer)
            channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)

    else:
        raise ValueError(f"Unsupported board number {board_no} for test beam.")

    return drsboard


def build_time_reference_channels(run=316):
    """
    Returns a list of time reference channels.
    """
    time_reference_channels = []
    if is_scan_run(run):
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


def build_hodo_trigger_channels(run=316):
    """
    Returns a list of hodoscope trigger channels.
    """
    hodo_trigger_channels = []
    if is_scan_run(run):
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


def build_hodo_pos_channels(run=316):
    """
    Returns a dictionary containing the hodoscope channels for the position measurements
    """
    hodoscope_channels = {}
    if is_scan_run(run):
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

    if run >= 1170 and run < 1327:
        # test beam
        # Left and right; top and bottom for Hodoscope2
        # 1170 is not exactly the first run number
        hodoscope_channels = {}
        hodoscope_channels["LR1"] = [
            "DRS_Board7_Group0_Channel4",
            "DRS_Board7_Group0_Channel5",
        ]
        hodoscope_channels["UD1"] = [
            "DRS_Board7_Group0_Channel6",
            "DRS_Board7_Group0_Channel7",
        ]

    if run >= 1342:
        # Left and right; top and bottom for Hodoscope2
        hodoscope_channels = {}
        hodoscope_channels["LR1"] = [
            "DRS_Board7_Group0_Channel0",
            "DRS_Board7_Group0_Channel1",
        ]
        hodoscope_channels["UD1"] = [
            "DRS_Board7_Group0_Channel2",
            "DRS_Board7_Group0_Channel3",
        ]
        hodoscope_channels["LR2"] = [
            "DRS_Board7_Group0_Channel4",
            "DRS_Board7_Group0_Channel5",
        ]
        hodoscope_channels["UD2"] = [
            "DRS_Board7_Group0_Channel6",
            "DRS_Board7_Group0_Channel7",
        ]

    return hodoscope_channels


def findFanoutTimeReferenceDelay(channel, run=1040):
    if str(run) not in triggerdelay.keys():
        return triggerdelay["default"][channel]
    else:
        return triggerdelay[str(run)][channel]


def findDRSTriggerMap(channel, run=1040):
    result = "_".join(channel.split("_")[:3])
    # print(triggermap.keys())
    if str(run) not in triggermap.keys():
        return triggermap["default"][result]
    else:
        return triggermap[str(run)][result]


def get_hole_veto_channel(run=1184):
    if run < 1183:
        return None
    elif run < 1327:
        return "DRS_Board7_Group1_Channel6"
    else:
        return get_service_drs_channels(run)["HoleVeto"]


def get_downstream_muon_channel(run=1184):
    if run < 1183:
        return None
    elif run < 1327:
        return "DRS_Board7_Group1_Channel0"
    else:
        return "DRS_Board7_Group1_Channel0"


def get_downstream_ttu_muon_channel(run=1184):
    if run < 1183:
        return None
    elif run < 1327:
        return None
    else:
        return get_service_drs_channels(run)["TTUMuonVeto"]


def get_pre_shower_channel(run=1184):
    """
    Returns the pre-shower channel.
    """
    if run < 1183:
        return None
    elif run < 1327:
        return "DRS_Board7_Group1_Channel1"
    else:
        return get_service_drs_channels(run)["PSD"]


def get_cerenkov_counters(run=1184):
    """
    Returns a list of Cerenkov counter channels.
    """
    if run < 1183:
        return []
    elif run < 1327:
        return [
            "DRS_Board7_Group1_Channel2",
            "DRS_Board7_Group1_Channel3",
            "DRS_Board7_Group1_Channel4",
            "DRS_Board7_Group1_Channel5",
        ]
    else:
        return [
            "DRS_Board7_Group2_Channel5",
            "DRS_Board7_Group2_Channel6",
            "DRS_Board7_Group2_Channel7",
        ]


def get_mcp_channels(run=1184):
    """
    Returns a list of MCP channels.
    """
    # todo: add MCP channels for August test beam
    if run < 1342:
        return []
    else:
        return {
            "DS": [
                "DRS_Board0_Group3_Channel6",
                "DRS_Board1_Group3_Channel6",
                "DRS_Board2_Group3_Channel6",
                "DRS_Board3_Group3_Channel6",
            ],
            "US": [
                "DRS_Board0_Group3_Channel7",
                "DRS_Board1_Group3_Channel7",
                "DRS_Board2_Group3_Channel7",
                "DRS_Board3_Group3_Channel7",
            ]
        }


def get_service_drs_channels(run=1184):
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
        return {
            "DWC1Left": "DRS_Board7_Group0_Channel0",
            "DWC1Right": "DRS_Board7_Group0_Channel1",
            "DWC1Up": "DRS_Board7_Group0_Channel2",
            "DWC1Down": "DRS_Board7_Group0_Channel3",
            "DWC2Left": "DRS_Board7_Group0_Channel4",
            "DWC2Right": "DRS_Board7_Group0_Channel5",
            "DWC2Up": "DRS_Board7_Group0_Channel6",
            "DWC2Down": "DRS_Board7_Group0_Channel7",
            "MuonVeto": "DRS_Board7_Group1_Channel0",
            "PSD": "DRS_Board7_Group1_Channel1",
            "HoleVeto": "DRS_Board7_Group1_Channel6",
            "NC": "DRS_Board7_Group1_Channel7",
            "T3": "DRS_Board7_Group2_Channel0",
            "T4": "DRS_Board7_Group2_Channel1",
            "KT1": "DRS_Board7_Group2_Channel2",
            "KT2": "DRS_Board7_Group2_Channel3",
            "TTUMuonVeto": "DRS_Board7_Group2_Channel4",
            "Cer474": "DRS_Board7_Group2_Channel5",
            "Cer519": "DRS_Board7_Group2_Channel6",
            "Cer537": "DRS_Board7_Group2_Channel7",
        }


def get_quartz_channel_list():
    channels_quartz = []
    # fmt off
    # 3mm region
    channels_quartz += [
        (-1.5, 1.875),
        (-1.5, 1.625), (-0.5, 1.625), (0.5, 1.625), (1.5, 1.625),
        (-1.5, 1.125), (-0.5, 1.125), (0.5, 1.125), (1.5, 1.125),
        (-1.5, 0.625), (-0.5, 0.625), (0.5, 0.625), (1.5, 0.625),
        (-0.5, 0.375),
        (-1.5, 0.125), (-0.5, 0.125), (0.5, 0.125), (1.5, 0.125),
        (-1.5, -0.375), (-0.5, -0.375), (0.5, -0.375), (1.5, -0.375),
        (-1.5, -0.875), (-0.5, -0.875), (0.5, -0.875), (1.5, -0.875),
        (-1.5, -1.375), (-0.5, -1.375), (0.5, -1.375), (1.5, -1.375),
        (-1.5, -1.875), (-0.5, -1.875), (0.5, -1.875), (1.5, -1.875),
    ]
    # 6mm region
    channels_quartz += [
        (-3.5, 5.5), (-2.5, 5.5), (-1.5, 5.5), (-0.5, 5.5), (0.5, 5.5), (1.5, 5.5), (2.5, 5.5), (3.5, 5.5),  # noqa
        (-5.5, 4.5), (-4.5, 4.5), (-3.5, 4.5), (-2.5, 4.5), (-1.5, 4.5), (-0.5, 4.5), (0.5, 4.5), (1.5, 4.5), (2.5, 4.5), (3.5, 4.5), (4.5, 4.5), (5.5, 4.5),  # noqa
        (-6.5, 3.5), (-5.5, 3.5), (-4.5, 3.5), (-3.5, 3.5), (-2.5, 3.5), (-1.5, 3.5), (-0.5, 3.5), (0.5, 3.5), (1.5, 3.5), (2.5, 3.5), (3.5, 3.5), (4.5, 3.5), (5.5, 3.5), (6.5, 3.5),  # noqa
        (-7.5, 2.5), (-6.5, 2.5), (-5.5, 2.5), (-4.5, 2.5), (-3.5, 2.5), (-2.5, 2.5), (-1.5, 2.5), (-0.5, 2.5), (0.5, 2.5), (1.5, 2.5), (2.5, 2.5), (3.5, 2.5), (4.5, 2.5), (5.5, 2.5), (6.5, 2.5), (7.5, 2.5),  # noqa
        (-7.5, 1.5), (-6.5, 1.5), (-5.5, 1.5), (-4.5, 1.5), (-3.5, 1.5), (-2.5, 1.5), (-1.5, 1.5), (-0.5, 1.5), (0.5, 1.5), (1.5, 1.5), (2.5, 1.5), (3.5, 1.5), (4.5, 1.5), (5.5, 1.5), (6.5, 1.5), (7.5, 1.5),  # noqa
        (-7.5, 0.5), (-6.5, 0.5), (-5.5, 0.5), (-4.5, 0.5), (-3.5, 0.5), (-2.5, 0.5), (-1.5, 0.5), (-0.5, 0.5), (0.5, 0.5), (1.5, 0.5), (2.5, 0.5), (3.5, 0.5), (4.5, 0.5), (5.5, 0.5), (6.5, 0.5), (7.5, 0.5),  # noqa
        (-7.5, -0.5), (-6.5, -0.5), (-5.5, -0.5), (-4.5, -0.5), (-3.5, -0.5), (-2.5, -0.5), (-1.5, -0.5), (-0.5, -0.5), (0.5, -0.5), (1.5, -0.5), (2.5, -0.5), (3.5, -0.5), (4.5, -0.5), (5.5, -0.5), (6.5, -0.5), (7.5, -0.5),  # noqa
        (-7.5, -1.5), (-6.5, -1.5), (-5.5, -1.5), (-4.5, -1.5), (-3.5, -1.5), (-2.5, -1.5), (-1.5, -1.5), (-0.5, -1.5), (0.5, -1.5), (1.5, -1.5), (2.5, -1.5), (3.5, -1.5), (4.5, -1.5), (5.5, -1.5), (6.5, -1.5), (7.5, -1.5),  # noqa
        (-7.5, -2.5), (-6.5, -2.5), (-5.5, -2.5), (-4.5, -2.5), (-3.5, -2.5), (-2.5, -2.5), (-1.5, -2.5), (-0.5, -2.5), (0.5, -2.5), (1.5, -2.5), (2.5, -2.5), (3.5, -2.5), (4.5, -2.5), (5.5, -2.5), (6.5, -2.5), (7.5, -2.5),  # noqa
        (-6.5, -3.5), (-5.5, -3.5), (-4.5, -3.5), (-3.5, -3.5), (-2.5, -3.5), (-1.5, -3.5), (-0.5, -3.5), (0.5, -3.5), (1.5, -3.5), (2.5, -3.5), (3.5, -3.5), (4.5, -3.5), (5.5, -3.5), (6.5, -3.5),  # noqa
        (-5.5, -4.5), (-4.5, -4.5), (-3.5, -4.5), (-2.5, -4.5), (-1.5, -4.5), (-0.5, -4.5), (0.5, -4.5), (1.5, -4.5), (2.5, -4.5), (3.5, -4.5), (4.5, -4.5), (5.5, -4.5),  # noqa
        (-3.5, -5.5), (-2.5, -5.5), (-1.5, -5.5), (-0.5, -5.5), (0.5, -5.5), (1.5, -5.5), (2.5, -5.5), (3.5, -5.5),  # noqa

    ]
    # fmt on

    return channels_quartz


def update_quartz_channels(boards):
    """
    Update the quartz channels to match the FERS/DRS board positions.
    """
    channels_quartz = get_quartz_channel_list()
    for board in boards.values():
        for channel in board.get_list_of_channels(isCer=True):
            if (channel.i_tower_x, channel.i_tower_y) in channels_quartz:
                channel.isQuartz = True


if __name__ == "__main__":
    # Example usage
    run_number = 583
    fers_boards = build_fers_boards(run=run_number)
    drs_boards = build_drs_boards(run=run_number)

    print("FERS Boards:")
    for board_name, board in fers_boards.items():
        print(f"{board_name}: {board}")

    print("\nDRS Boards:")
    for board_name, board in drs_boards.items():
        print(f"{board_name}: {board}")

    print("\nHodoscope Position Channels:")
    hodo_channels = build_hodo_pos_channels(run=run_number)
    for hodo_type, channels in hodo_channels.items():
        print(f"{hodo_type}: {channels}")
