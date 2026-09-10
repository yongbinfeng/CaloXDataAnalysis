"""FERS board layouts and the geometry of the calorimeter face."""
from channels.calox_channel import (FERSBoard, FERSBoards, A5202_map,
                                    A5205_map_3mm)
from channels.maps.eras import _for_run
from channels.maps.quartz import update_quartz_channels


# The calorimeter face is chamfered at the four corners: a few towers along the
# outer edge of the side arms and of the outermost top/bottom boards are not
# where a plain 4x8 board rectangle would put them, but follow the octagonal
# detector edge. The corrections are expressed as absolute tower positions, so
# they apply to whichever board happens to cover them in a given run.
def _build_corner_chamfer():
    """Map (from) -> (to) tower position for the four chamfered corners.

    For the top-left corner (the others are the mirror images in x and/or y):
      side arm: (-13.5, -12.5) x (+3.5, +2.5) -> (-11.5, -10.5) x (+5.5, +4.5)
      shoulder: (-9.5, +7.5) -> (-8.5, +8.5),  (-9.5, +6.5) -> (-7.5, +8.5)
    """
    chamfer = {}
    for sx in (-1, 1):
        for sy in (-1, 1):
            # the 2x2 block at the outer end of the side arm, shifted two
            # columns inwards and two rows further out
            for dx in (0, 1):
                for dy in (0, 1):
                    chamfer[(sx * (13.5 - dx), sy * (3.5 - dy))] = \
                        (sx * (11.5 - dx), sy * (5.5 - dy))
            # the two towers at the outer top (bottom) corner of the shoulder
            # board, folded around the corner onto the row above (below)
            chamfer[(sx * 9.5, sy * 7.5)] = (sx * 8.5, sy * 8.5)
            chamfer[(sx * 9.5, sy * 6.5)] = (sx * 7.5, sy * 8.5)
    return chamfer


_CORNER_CHAMFER = _build_corner_chamfer()


def apply_corner_chamfer(boards):
    """Move the corner towers of the outermost boards (see above).

    Source and target positions are disjoint, so this is safe to apply once to
    a freshly positioned set of boards. Boards that cover none of the corner
    positions (the 3mm boards, and all boards in the smaller run configs) are
    left untouched.
    """
    for board in boards.values():
        for channel in board.channels:
            new_pos = _CORNER_CHAMFER.get(
                (channel.i_tower_x, channel.i_tower_y))
            if new_pos is not None:
                channel.i_tower_x, channel.i_tower_y = new_pos


# FERS board layout per run era: {board_no: (granularity, iTowerX, iTowerY)}.
# The position is where the board's first channel is moved to; a 6mm board
# spans 4x8 towers from there, a 3mm board 2x4 (at 4x finer granularity in y).
# Ranges are [run_min, run_max); run_max None means "and everything after".
_FERS_LAYOUTS = (
    ((316, 317), {  # 5 FERS board in 316
        1: ("6mm", 6.5, -0.5),
        2: ("6mm", 2.5, -0.5),
        3: ("6mm", -1.5, -2.5),
        4: ("6mm", -5.5, -0.5),
        5: ("6mm", -9.5, -0.5),
    }),
    ((571, 572), {  # 10 FERS boards in 571
        0: ("6mm", -13.5, 3.5),
        1: ("6mm", -9.5, 7.5),
        2: ("6mm", -5.5, 7.5),
        3: ("6mm", -1.5, 9.5),
        4: ("6mm", 2.5, 7.5),
        5: ("6mm", 6.5, 7.5),
        8: ("6mm", -9.5, -0.5),
        9: ("6mm", -5.5, -0.5),
        10: ("6mm", -1.5, -2.5),
        11: ("6mm", 2.5, -0.5),
    }),
    ((583, 685), {  # include 3mm FERS board in 583
        0: ("6mm", -13.5, 3.5),
        1: ("6mm", -9.5, 7.5),
        2: ("6mm", -5.5, 7.5),
        3: ("6mm", -1.5, 9.5),
        4: ("6mm", 2.5, 7.5),
        8: ("6mm", -9.5, -0.5),
        9: ("6mm", -5.5, -0.5),
        10: ("6mm", -1.5, -2.5),
        11: ("6mm", 2.5, -0.5),
        5: ("3mm", -1.5, 1.875),
    }),
    ((685, 895), {
        0: ("6mm", -5.5, 7.5),
        1: ("6mm", -1.5, 9.5),
        2: ("6mm", 2.5, 7.5),
        4: ("6mm", -5.5, -0.5),
        5: ("6mm", -1.5, -2.5),
        6: ("6mm", 2.5, -0.5),
        3: ("3mm", -1.5, 1.875),
        7: ("3mm", 0.5, 1.875),
    }),
    ((895, 1100), {
        0: ("6mm", -13.5, 3.5),
        1: ("6mm", -9.5, 7.5),
        2: ("6mm", -5.5, 7.5),
        4: ("6mm", -1.5, 9.5),
        5: ("6mm", 2.5, 7.5),
        6: ("6mm", 6.5, 7.5),
        7: ("6mm", 10.5, 3.5),
        8: ("6mm", -9.5, -0.5),
        9: ("6mm", -5.5, -0.5),
        10: ("6mm", -1.5, -2.5),
        12: ("6mm", 2.5, -0.5),
        13: ("6mm", 6.5, -0.5),
        3: ("3mm", -1.5, 1.875),
        11: ("3mm", 0.5, 1.875),
    }),
    ((1173, 1327), {  # test beam
        0: ("6mm", -13.5, 3.5),
        1: ("6mm", -9.5, 7.5),
        2: ("6mm", -5.5, 7.5),
        3: ("6mm", -1.5, 9.5),
        4: ("6mm", 2.5, 7.5),
        5: ("6mm", 6.5, 7.5),
        6: ("6mm", 10.5, 3.5),
        9: ("6mm", -9.5, -0.5),
        10: ("6mm", -5.5, -0.5),
        11: ("6mm", -1.5, -2.5),
        12: ("6mm", 2.5, -0.5),
        13: ("6mm", 6.5, -0.5),
        7: ("3mm", -1.5, 1.875),
        8: ("3mm", 0.5, 1.875),
    }),
    ((1342, 1600), {  # September 2024 test beam
        2: ("6mm", -13.5, 3.5),
        3: ("6mm", -9.5, 7.5),
        4: ("6mm", -5.5, 7.5),
        5: ("6mm", -1.5, 9.5),
        6: ("6mm", 2.5, 7.5),
        7: ("6mm", 6.5, 7.5),
        8: ("6mm", 10.5, 3.5),
        11: ("6mm", -9.5, -0.5),
        12: ("6mm", -5.5, -0.5),
        13: ("6mm", -1.5, -2.5),
        14: ("6mm", 2.5, -0.5),
        15: ("6mm", 6.5, -0.5),
        9: ("3mm", -1.5, 1.875),
        10: ("3mm", 0.5, 1.875),
    }),
    ((1720, None), {  # 2026 test beam
        2: ("6mm", -9.5, -0.5),
        3: ("6mm", -5.5, -0.5),
        4: ("6mm", -1.5, -2.5),
        5: ("6mm", 2.5, -0.5),
        8: ("6mm", 6.5, -0.5),
        9: ("6mm", 10.5, 3.5),
        10: ("6mm", 2.5, 7.5),
        11: ("6mm", 6.5, 7.5),
        12: ("6mm", -13.5, 3.5),
        13: ("6mm", -9.5, 7.5),
        14: ("6mm", -5.5, 7.5),
        15: ("6mm", -1.5, 9.5),
        6: ("3mm", -1.5, 1.875),
        7: ("3mm", 0.5, 1.875),
    }),
)


def _fers_layout(run_number):
    """The FERS board layout for a run (see _FERS_LAYOUTS)."""
    layout = _for_run(_FERS_LAYOUTS, run_number)
    if layout is None:
        raise ValueError(f"Unsupported run_number {run_number} for FERS boards.")
    return layout


def build_fers_boards(run_number=316):
    """
    Build a map for ixy and FERS channels for both boards.
    """
    base = {"6mm": FERSBoard(board_no=-1, is6mm=True),
            "3mm": FERSBoard(board_no=-1, is6mm=False)}
    fersboards = FERSBoards()
    for board_no, (size, i_tower_x, i_tower_y) in _fers_layout(run_number).items():
        board = base[size].copy(board_no=board_no)
        board.move_to(i_tower_x, i_tower_y)
        fersboards[f"Board{board_no}"] = board
    apply_corner_chamfer(fersboards)
    update_quartz_channels(fersboards)
    return fersboards


def physical_to_fers_channel(physical, is6mm):
    """Map a FERS 'Physical' position (1-64) to the FERS readout channel number.

    The A5202_map (6mm) / A5205_map_3mm (3mm) arrays store, at [ix, iy], the
    readout channel number. Within each group of 4, the physical positions fill
    ix right-to-left: physical 1,2,3,4 -> ix 2,3,0,1 (i.e. ix = ((p-1)+2) % 4),
    iy = (p-1)//4. With this ordering the CSV materials (Cer/Sci and
    quartz/plastic) agree exactly with the FERS map. Same convention for 6mm
    and 3mm.
    """
    fers_map = A5202_map if is6mm else A5205_map_3mm
    p = int(physical) - 1
    return int(fers_map[(p + 2) % 4, p // 4])
