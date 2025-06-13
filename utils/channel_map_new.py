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
                FERSBoards["Board5"][tower_ix, tower_iy].iTowerY += 11.75
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


if __name__ == "__main__":
    # Example usage
    run_number = 316
    fers_boards = buildFERSBoards(run=run_number)
    drs_boards = buildDRSBoards(run=run_number)

    print("FERS Boards:")
    for board_name, board in fers_boards.items():
        print(f"{board_name}: {board}")

    print("\nDRS Boards:")
    for board_name, board in drs_boards.items():
        print(f"{board_name}: {board}")
