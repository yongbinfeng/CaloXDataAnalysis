import numpy as np
from utils.CaloXChannel import FERSBoard, DRSBoard


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

    elif run >= 895:
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

    else:
        raise ValueError(f"Unsupported run number {run} for FERS boards.")
    return FERSBoards


def buildDRSBoards(run=316):
    """
    Build a map for ixy and DRS channels.
    DRS boards are the same for all runs for now.
    """
    base_DRSBoard = DRSBoard(boardNo=-1)
    DRSBoards = {}
    if run < 685:
        DRSBoards["Board2"] = base_DRSBoard.copy(boardNo=2)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 5)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 6)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 7)
        DRSBoards["Board2"].MoveTo(-1.5, -2.5)
        channel = DRSBoards["Board2"].GetChannelByGroupChannel(3, 4)
        channel.iTowerX = 1.5
        channel.iTowerY = -9.5
        channel.isCer = False

        DRSBoards["Board0"] = base_DRSBoard.copy(boardNo=0)
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

    elif run >= 685:
        DRSBoards["Board1"] = base_DRSBoard.copy(boardNo=1)
        DRSBoards["Board2"] = base_DRSBoard.copy(boardNo=2)

        # remove channels that are not used
        DRSBoards["Board1"].MoveTo(-1.5, 9.5)
        DRSBoards["Board2"].MoveTo(-1.5, -6.5)
        # no more cables..
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 4)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 5)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 6)
        DRSBoards["Board2"].RemoveChannelByGroupChannel(3, 7)
    else:
        raise ValueError(f"Unsupported run number {run} for DRS boards.")
    return DRSBoards


def buildTimeReferenceChannels(run=316):
    """
    Returns a list of time reference channels.
    """
    time_reference_channels = []
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

    if run >= 685:
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

    return hodoscope_channels


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
