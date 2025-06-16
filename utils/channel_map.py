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
        for tower_ix in range(0, 4):
            for tower_iy in range(0, 16):
                FERSBoards["Board2"][tower_ix, tower_iy].iTowerX -= 4
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerX -= 8
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerY -= 4
                FERSBoards["Board4"][tower_ix, tower_iy].iTowerX -= 12
                FERSBoards["Board5"][tower_ix, tower_iy].iTowerX -= 16

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
        # FERSBoards["Board12"] = base_FERSBoard_3mm

        for tower_ix in range(0, 4):
            for tower_iy in range(0, 16):
                FERSBoards["Board0"][tower_ix, tower_iy].iTowerX -= 12
                FERSBoards["Board0"][tower_ix, tower_iy].iTowerY += 6
                FERSBoards["Board1"][tower_ix, tower_iy].iTowerX -= 8
                FERSBoards["Board1"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board2"][tower_ix, tower_iy].iTowerX -= 4
                FERSBoards["Board2"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerX += 0
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerY += 12
                FERSBoards["Board4"][tower_ix, tower_iy].iTowerX += 4
                FERSBoards["Board4"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board5"][tower_ix, tower_iy].iTowerX += 8
                FERSBoards["Board5"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board8"][tower_ix, tower_iy].iTowerX -= 8
                FERSBoards["Board8"][tower_ix, tower_iy].iTowerY += 2
                FERSBoards["Board9"][tower_ix, tower_iy].iTowerX -= 4
                FERSBoards["Board9"][tower_ix, tower_iy].iTowerY += 2
                FERSBoards["Board10"][tower_ix, tower_iy].iTowerX += 0
                FERSBoards["Board10"][tower_ix, tower_iy].iTowerY += 0
                FERSBoards["Board11"][tower_ix, tower_iy].iTowerX += 4
                FERSBoards["Board11"][tower_ix, tower_iy].iTowerY += 2
                # FERSBoards["Board12"][tower_ix, tower_iy].iTowerX += 8
                # FERSBoards["Board12"][tower_ix, tower_iy].iTowerY += 2

    elif run == 583:
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

        for tower_ix in range(0, 4):
            for tower_iy in range(0, 16):
                FERSBoards["Board0"][tower_ix, tower_iy].iTowerX -= 12
                FERSBoards["Board0"][tower_ix, tower_iy].iTowerY += 6
                FERSBoards["Board1"][tower_ix, tower_iy].iTowerX -= 8
                FERSBoards["Board1"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board2"][tower_ix, tower_iy].iTowerX -= 4
                FERSBoards["Board2"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerX += 0
                FERSBoards["Board3"][tower_ix, tower_iy].iTowerY += 12
                FERSBoards["Board4"][tower_ix, tower_iy].iTowerX += 4
                FERSBoards["Board4"][tower_ix, tower_iy].iTowerY += 10
                FERSBoards["Board8"][tower_ix, tower_iy].iTowerX -= 8
                FERSBoards["Board8"][tower_ix, tower_iy].iTowerY += 2
                FERSBoards["Board9"][tower_ix, tower_iy].iTowerX -= 4
                FERSBoards["Board9"][tower_ix, tower_iy].iTowerY += 2
                FERSBoards["Board10"][tower_ix, tower_iy].iTowerX += 0
                FERSBoards["Board10"][tower_ix, tower_iy].iTowerY += 0
                FERSBoards["Board11"][tower_ix, tower_iy].iTowerX += 4
                FERSBoards["Board11"][tower_ix, tower_iy].iTowerY += 2
                # FERSBoards["Board12"][tower_ix, tower_iy].iTowerX += 8
                # FERSBoards["Board12"][tower_ix, tower_iy].iTowerY += 2

                FERSBoards["Board5"][tower_ix, tower_iy].iTowerX += 0
                FERSBoards["Board5"][tower_ix, tower_iy].iTowerY += 4.375
    else:
        raise ValueError(f"Unsupported run number {run} for FERS boards.")
    return FERSBoards


def buildDRSBoards(run=316):
    """
    Build a map for ixy and DRS channels.
    """
    base_DRSBoard = DRSBoard(boardNo=-1)
    DRSBoards = {}
    # 2 DRS boards in 316
    DRSBoards["Board0"] = base_DRSBoard.copy(boardNo=0)
    DRSBoards["Board2"] = base_DRSBoard.copy(boardNo=2)
    for tower_ix in range(0, 4):
        for tower_iy in range(0, 8):
            DRSBoards["Board2"][tower_ix, tower_iy].iTowerX -= 0
            DRSBoards["Board2"][tower_ix, tower_iy].iTowerY -= 4
    return DRSBoards


def buildHodoChannels(run=316):
    """
    Returns a dictionary containing the hodoscope channels.
    """
    hodoscope_channels = {}
    hodoscope_channels["trigger"] = ["DRS_Board1_Group0_Channel0"]
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

    return hodoscope_channels


def DrawFERSBoards(run=316):
    """
    Draws the FERS boards for a given run in a TH2D
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)  # Disable interactive mode
    xmax = 11.5
    xmin = -12.5
    ymax = 15.5
    ymin = -8.5
    fers_boards = buildFERSBoards(run)
    h2_Cer = ROOT.TH2D(f"h_FERSBoards_Run{run}_Cer", f"FERS Board Cer Channels for Run {run}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Sci = ROOT.TH2D(f"h_FERSBoards_Run{run}_Sci", f"FERS Board Sci Channels for Run {run}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Cer_3mm = ROOT.TH2D(f"h_FERSBoards_Run{run}_Cer_3mm", f"FERS Board Cer Channels for Run {run} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Sci_3mm = ROOT.TH2D(f"h_FERSBoards_Run{run}_Sci_3mm", f"FERS Board Sci Channels for Run {run} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Cer_3mm.SetMarkerSize(0.70)
    h2_Sci_3mm.SetMarkerSize(0.70)
    bboxes = []
    colors = [
        ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1,
        ROOT.kCyan+2, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kViolet+1, ROOT.kAzure+1,
        ROOT.kSpring+3, ROOT.kTeal+2, ROOT.kPink+2, ROOT.kGray+1, ROOT.kOrange+7
    ]

    for _, FERSBoard in fers_boards.items():
        towerx_min = 99
        towerx_max = -99
        towery_min = 99
        towery_max = -99
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            channel_Cer = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            channel_Sci = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)
            if FERSBoard.Is3mm():
                h2_Cer_3mm.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                                channel_Cer.channelNo)
                h2_Sci_3mm.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                                channel_Sci.channelNo)
                extra_x_side = 0.5
                extra_y_side = 0.5 / 4
            else:
                h2_Cer.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                            channel_Cer.channelNo)
                h2_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                            channel_Sci.channelNo)
                extra_x_side = 0.5
                extra_y_side = 0.5
            towerx_min = min(towerx_min, channel_Cer.iTowerX)
            towerx_max = max(towerx_max, channel_Cer.iTowerX)
            towery_min = min(towery_min, channel_Cer.iTowerY)
            towery_max = max(towery_max, channel_Cer.iTowerY)

        bbox = ROOT.TBox(towerx_min - extra_x_side, towery_min - extra_y_side,
                         towerx_max + extra_x_side, towery_max + extra_y_side)
        bbox.SetFillColorAlpha(colors[len(bboxes)], 0.1)
        bboxes.append(bbox)

    extraToDraw = ROOT.TPaveText(0.01, 0.73, 0.12, 0.88, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Run: {run}")
    extraToDraw_Cer = extraToDraw.Clone("extraToDraw_2")
    extraToDraw_Cer.AddText("Cerenkov")
    extraToDraw_Sci = extraToDraw.Clone("extraToDraw_3")
    extraToDraw_Sci.AddText("Scintillator")

    import sys
    sys.path.append("CMSPLOTS")  # noqa
    from myFunction import DrawHistos
    output_dir = f"plots/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=1000, extraToDraw=bboxes + [extraToDraw_Cer])
    DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=1000, extraToDraw=bboxes + [extraToDraw_Sci])


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

    print("\nHodoscope Channels:")
    hodo_channels = buildHodoChannels(run=run_number)
    for hodo_type, channels in hodo_channels.items():
        print(f"{hodo_type}: {channels}")

    DrawFERSBoards(run=run_number)
