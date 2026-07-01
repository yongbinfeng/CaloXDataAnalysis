from channels.calox_channel import FERSBoard, DRSBoard, DRSChannel, drs_map, FERSBoards, add_drs_reference_channel, A5202_map, A5205_map_3mm
from utils.data_loader import is_scan_run
import json
import os
import csv
import re
from collections import OrderedDict

# Run number at which DRS branches gained the "Brg{N}_" prefix.
_DRS_BRG_RUN = 1700
# Run number at which DRS boards switched from 3mm to 6mm fibers.
_DRS_6MM_RUN = 1748
_DRS_FULL_RUN_TB2026 = 1828
_DRS_TESTFIBER_RUN_TB2026 = 1994

# CSV mapping each DRS readout channel (DRS_ROOT_Bridge_Mapping) to its FERS
# channel (fers_root_board, Physical-1). Used only for the from-FERS era
# (run >= _DRS_FULL_RUN_TB2026); earlier runs use the function-based board maps.
# The mapping was rewired between phases, so the CSV is run-dependent:
#   1828 <= run < 1896  : DRS_PHASE_1.csv
#   run >= 1896         : DRS_PHASE_2.csv
_CHANNEL_MAP_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "channel_maps")
_DRS_FERS_CSV_PHASE1 = os.path.join(_CHANNEL_MAP_DIR, "DRS_PHASE_1.csv")
_DRS_FERS_CSV_PHASE2 = os.path.join(_CHANNEL_MAP_DIR, "DRS_PHASE_2.csv")
_DRS_PHASE2_RUN = 1896


def _drs_fers_csv(run_number):
    """DRS->FERS channel-map CSV for a run (see ranges above)."""
    if run_number is not None and run_number >= _DRS_PHASE2_RUN:
        return _DRS_FERS_CSV_PHASE2
    return _DRS_FERS_CSV_PHASE1


def _drs(board, group, ch, brg=None):
    """Build a DRS channel branch name, optionally with a bridge prefix."""
    if brg is not None:
        return f"DRS_Brg{brg}_Board{board}_Group{group}_Channel{ch}"
    return f"DRS_Board{board}_Group{group}_Channel{ch}"

f_triggerdelay = "data/triggerdelay.json"
with open(f_triggerdelay, 'r') as f:
    triggerdelay = json.load(f)

f_drstriggermap = "data/drstriggermap.json"
with open(f_drstriggermap, 'r') as f:
    triggermap = json.load(f)


def build_fers_boards(run_number=316):
    """
    Build a map for ixy and FERS channels for both boards.
    """
    base_FERSBoard_6mm = FERSBoard(board_no=-1, is6mm=True)
    base_FERSBoard_3mm = FERSBoard(board_no=-1, is6mm=False)
    fersboards = FERSBoards()
    if run_number == 316:
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

    elif run_number == 571:
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

    elif run_number >= 583 and run_number < 685:
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

    elif run_number >= 685 and run_number < 895:
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

    elif run_number >= 895 and run_number < 1100:
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
    elif run_number >= 1173 and run_number < 1327:
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

    elif run_number >= 1342 and run_number < 1600:
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
    elif run_number >= 1720:
        fersboards["Board2"] = base_FERSBoard_6mm.copy(board_no=2)
        fersboards["Board3"] = base_FERSBoard_6mm.copy(board_no=3)
        fersboards["Board4"] = base_FERSBoard_6mm.copy(board_no=4)
        fersboards["Board5"] = base_FERSBoard_6mm.copy(board_no=5)
        fersboards["Board8"] = base_FERSBoard_6mm.copy(board_no=8)
        fersboards["Board9"] = base_FERSBoard_6mm.copy(board_no=9)
        fersboards["Board10"] = base_FERSBoard_6mm.copy(board_no=10)
        fersboards["Board11"] = base_FERSBoard_6mm.copy(board_no=11)
        fersboards["Board12"] = base_FERSBoard_6mm.copy(board_no=12)
        fersboards["Board13"] = base_FERSBoard_6mm.copy(board_no=13)
        fersboards["Board14"] = base_FERSBoard_6mm.copy(board_no=14)
        fersboards["Board15"] = base_FERSBoard_6mm.copy(board_no=15)

        fersboards["Board6"] = base_FERSBoard_3mm.copy(board_no=6)
        fersboards["Board7"] = base_FERSBoard_3mm.copy(board_no=7)

        fersboards["Board12"].move_to(-13.5, 3.5)
        fersboards["Board13"].move_to(-9.5, 7.5)
        fersboards["Board14"].move_to(-5.5, 7.5)
        fersboards["Board15"].move_to(-1.5, 9.5)
        fersboards["Board10"].move_to(2.5, 7.5)
        fersboards["Board11"].move_to(6.5, 7.5)
        fersboards["Board9"].move_to(10.5, 3.5)
        fersboards["Board2"].move_to(-9.5, -0.5)
        fersboards["Board3"].move_to(-5.5, -0.5)
        fersboards["Board4"].move_to(-1.5, -2.5)
        fersboards["Board5"].move_to(2.5, -0.5)
        fersboards["Board8"].move_to(6.5, -0.5)

        # 3mm
        fersboards["Board6"].move_to(-1.5, 1.875)
        fersboards["Board7"].move_to(0.5, 1.875)

    else:
        raise ValueError(f"Unsupported run_number {run_number} for FERS boards.")
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


def load_drs_fers_lookup(csv_path=_DRS_FERS_CSV_PHASE1):
    """Parse the DRS->FERS mapping CSV.

    Returns a list of dicts, one per signal DRS readout channel:
        bridge, board, group, channel : from the DRS_ROOT_Bridge_Mapping column
                                        (DRS_Bridge{b}_Board{B}_Group{G}_Channel{C},
                                         which equals the data branch DRS_Brg{b}_...)
        fers_board : FERS ROOT board number  (fers_root_board column)
        physical   : FERS physical position  (Physical column, 1-64); converted to
                     a FERS channel number with physical_to_fers_channel()
        isCer      : True for Plastic/Quartz (Cherenkov), False for Scintillation
        isQuartz   : True for Quartz
    Service-DRS rows (channels not wired to the calorimeter readout) are simply
    absent from the CSV, so they are naturally excluded.
    """
    entries = []
    with open(csv_path) as f:
        # Parse by column name (layout differs between phases, e.g. PHASE_2 has
        # no MCP/blank columns) rather than by fixed position.
        reader = csv.DictReader(f)
        for row in reader:
            mapping = (row.get("DRS_ROOT_Bridge_Mapping") or "").strip()
            m = re.match(
                r"DRS_Bridge(\d+)_Board(\d+)_Group(\d+)_Channel(\d+)", mapping)
            if not m:
                continue
            brg, bd, gr, ch = (int(m.group(i)) for i in range(1, 5))
            typ = (row.get("type") or "").strip()
            entries.append({
                "bridge": brg, "board": bd, "group": gr, "channel": ch,
                "fers_board": int(row["fers_root_board"]),
                "physical": int(row["Physical"]),
                "isCer": typ in ("Plastic", "Quartz"),
                "isQuartz": typ == "Quartz",
            })
    return entries


def build_drs_boards_from_fers(fersboards, csv_path=_DRS_FERS_CSV_PHASE1):
    """Build DRS boards whose channel (x, y) positions come from the FERS map.

    Each DRS readout channel is matched, via the CSV, to a FERS channel:
    fers_root_board + the FERS 'Physical' position, mapped through the
    A5202_map/A5205_map_3mm arrays (physical_to_fers_channel) to a FERS channel
    number. The DRS channel inherits that FERS channel's tower position AND its
    Cer/Sci + quartz/plastic type, so the DRS map matches the FERS map exactly
    (the FERS map already encodes the special quartz/plastic pattern of the few
    core channels that differ from the regular 3mm layout). The CSV material
    column agrees with the FERS type and is kept only as a cross-check.

    All channels are amplified. The 6mm/3mm flag is inherited from the matched
    FERS channel (the core DRS boards read the 3mm region, the outer ones the
    6mm region), so the channels land in the correct granularity map. A reference
    channel (Channel8, position (-999,-999)) is added per group, as the old
    builder did. Boards are keyed "Brg{bridge}_Board{board}". FERS maps are not
    modified.
    """
    fers_by = {(board.board_no, ch.channel_no): ch
               for board in fersboards.values() for ch in board}
    # Granularity is a board-level property; the per-channel FERSChannel.is6mm
    # flag is unreliable (build_fers_base does not set it), so take it from the
    # FERS board.
    fers_board_is6mm = {board.board_no: board.is6mm
                        for board in fersboards.values()}

    grouped = OrderedDict()
    for e in load_drs_fers_lookup(csv_path):
        is6mm = fers_board_is6mm[e["fers_board"]]
        fers_channel = physical_to_fers_channel(e["physical"], is6mm)
        fers_ch = fers_by.get((e["fers_board"], fers_channel))
        if fers_ch is None:
            raise ValueError(
                f"DRS_Brg{e['bridge']}_Board{e['board']}_Group{e['group']}_"
                f"Channel{e['channel']} maps to FERS board {e['fers_board']} "
                f"physical {e['physical']} (channel {fers_channel}), "
                f"which does not exist.")
        chan = DRSChannel(
            fers_ch.i_tower_x, fers_ch.i_tower_y, fers_ch.isCer,
            e["channel"], e["group"], e["board"],
            is_amplified=True, is6mm=is6mm,
            isQuartz=fers_ch.isQuartz, bridge_no=e["bridge"])
        grouped.setdefault((e["bridge"], e["board"]), []).append(chan)

    DRSBoards = {}
    for (brg, bd), channels in sorted(grouped.items()):
        # All signal channels of a board share one granularity; the reference
        # channel inherits it.
        board_is6mm = channels[0].is6mm
        for group_no in sorted({c.group_no for c in channels}):
            channels.append(add_drs_reference_channel(
                group_no, bd, is6mm=board_is6mm, bridge_no=brg))
        DRSBoards[f"Brg{brg}_Board{bd}"] = DRSBoard(
            board_no=bd, channels=channels, bridge_no=brg)
    return DRSBoards


def build_drs_boards(run_number=316):
    """
    Build a map for ixy and DRS channels.
    DRS boards are the same for all runs for now.
    """
    base_DRSBoard_6mm = DRSBoard(board_no=-1, is6mm=True)
    base_DRSBoard_3mm = DRSBoard(board_no=-1, is6mm=False)
    DRSBoards = {}
    if is_scan_run(run_number):
        # no DRS boards in scan runs
        # only FERS
        return DRSBoards

    if run_number < 685:
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

    elif run_number >= 685 and run_number < 1003:
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
    elif run_number >= 1003 and run_number < 1100:
        # 3 DRS boards in 1003
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_6mm.copy(board_no=2)
        DRSBoards["Board3"] = base_DRSBoard_6mm.copy(board_no=3)

        DRSBoards["Board1"].move_to(-1.5, 9.5)
        DRSBoards["Board2"].move_to(-1.5, 5.5)
        DRSBoards["Board3"].move_to(-1.5, -6.5)
        # all channels are connected.
        # no need to remove
        if run_number >= 1033:
            # two amplified channels in board 1
            channel = DRSBoards["Board1"].get_channel_by_group_channel(0, 0)
            channel.is_amplified = True
            channel = DRSBoards["Board1"].get_channel_by_group_channel(0, 1)
            channel.is_amplified = True
    elif run_number >= 1173 and run_number < 1327:
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
    elif run_number >= _DRS_FULL_RUN_TB2026:
        # 6mm boards (run_number >= 1748): the DRS channel positions are taken
        # from the corresponding FERS channels via the DRS->FERS CSV mapping
        # (build_drs_boards_from_fers). FERS maps are unchanged. The per-channel
        # bridge numbering and quartz/plastic flags come from the CSV, so we skip
        # the set_bridge_no / update_quartz_channels post-processing below.
        fersboards = build_fers_boards(run_number=run_number)
        return build_drs_boards_from_fers(
            fersboards, csv_path=_drs_fers_csv(run_number))
    elif run_number >= _DRS_6MM_RUN:
        # 6mm boards (run_number >= 1748): 2 calo boards (0, 1), all 6mm, no MCP
        # All channels are amplified (inverted), so is_amplified=True on every channel.
        DRSBoards["Board0"] = base_DRSBoard_6mm.copy(board_no=0)
        DRSBoards["Board1"] = base_DRSBoard_6mm.copy(board_no=1)

        DRSBoards["Board0"].move_to(-1.5, 1.875)
        DRSBoards["Board1"].move_to(-1.5, -0.125)

        for board in DRSBoards.values():
            for channel in board:
                channel.is_amplified = True

    elif run_number >= _DRS_BRG_RUN:
        # 3mm boards with bridge numbering (run_number 1700–1747): 3 calo boards, no MCP
        DRSBoards["Board0"] = base_DRSBoard_3mm.copy(board_no=0)
        DRSBoards["Board1"] = base_DRSBoard_3mm.copy(board_no=1)
        DRSBoards["Board2"] = base_DRSBoard_3mm.copy(board_no=2)

        DRSBoards["Board0"].move_to(-1.5, 1.875)
        DRSBoards["Board1"].move_to(-1.5, -0.125)
        DRSBoards["Board2"].move_to(0.5, 1.875)

    elif run_number >= 1342:
        # September 2024 test beam: 4 calo boards + 3 test beam boards
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
        raise ValueError(f"Unsupported run_number {run_number} for DRS boards.")

    if run_number >= _DRS_6MM_RUN:
        for board in DRSBoards.values():
            board.set_bridge_no(0)
    elif run_number >= _DRS_BRG_RUN:
        for board in DRSBoards.values():
            board.set_bridge_no(1)

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

        # add the reference channels (Channel 8)
        for igroup in range(4):
            channel = add_drs_reference_channel(igroup, board_no, is6mm=True)
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

        # add the reference channels (Channel 8)
        for igroup in range(4):
            channel = add_drs_reference_channel(igroup, board_no, is6mm=True)
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

        # add the reference channels (Channel 8)
        for igroup in range(4):
            channel = add_drs_reference_channel(igroup, board_no, is6mm=True)
            channels_DRS.append(channel)

        drsboard = DRSBoard(board_no=board_no, channels=channels_DRS)

    else:
        raise ValueError(f"Unsupported board number {board_no} for test beam.")

    return drsboard


def build_time_reference_channels(run_number=316):
    """
    Deprecated: use the reference channels in the DRS boards instead.
    """
    time_reference_channels = []
    if is_scan_run(run_number):
        # no time reference channels in scan runs
        # since no drs boards
        return time_reference_channels

    brg = 1 if run_number >= _DRS_BRG_RUN else None
    if run_number < 685:
        time_reference_channels.append(_drs(0, 3, 7))
        time_reference_channels.append(_drs(2, 3, 7))
        time_reference_channels.append(_drs(1, 0, 0))
    elif run_number >= 685:
        time_reference_channels.append(_drs(1, 0, 8, brg))
        time_reference_channels.append(_drs(1, 1, 8, brg))
        time_reference_channels.append(_drs(1, 2, 8, brg))
        time_reference_channels.append(_drs(1, 3, 8, brg))
        time_reference_channels.append(_drs(2, 0, 8, brg))
        time_reference_channels.append(_drs(2, 1, 8, brg))
        time_reference_channels.append(_drs(2, 2, 8, brg))
        time_reference_channels.append(_drs(2, 3, 8, brg))
        time_reference_channels.append(_drs(0, 0, 8, brg))
    else:
        raise ValueError(
            f"Unsupported run_number {run_number} for time reference channels.")

    return time_reference_channels


def build_hodo_trigger_channels(run_number=316):
    """
    Returns a list of hodoscope trigger channels.
    """
    hodo_trigger_channels = []
    if is_scan_run(run_number):
        # no hodoscope trigger channels in scan runs
        # since no drs boards
        return hodo_trigger_channels

    brg = 1 if run_number >= _DRS_BRG_RUN else None
    if run_number < 685:
        hodo_trigger_channels.append(_drs(1, 2, 0))
        hodo_trigger_channels.append(_drs(1, 2, 1))
    elif run_number >= 685:
        hodo_trigger_channels.append(_drs(0, 2, 0, brg))
        hodo_trigger_channels.append(_drs(0, 2, 1, brg))
    else:
        raise ValueError(
            f"Unsupported run_number {run_number} for hodoscope trigger channels.")

    return hodo_trigger_channels


def build_hodo_pos_channels(run_number=316):
    """
    Returns a dictionary containing the hodoscope channels for the position measurements
    """
    hodoscope_channels = {}
    if is_scan_run(run_number):
        # no hodoscope position channels in scan runs
        # since no drs boards
        return hodoscope_channels

    brg = 1 if run_number >= _DRS_BRG_RUN else None
    hodoscope_channels["TopX"] = [_drs(1, 0, 1), _drs(1, 0, 2)]
    hodoscope_channels["TopZ"] = [_drs(1, 0, 3), _drs(1, 0, 4)]
    hodoscope_channels["BottomX"] = [_drs(1, 0, 5), _drs(1, 0, 6)]
    if run_number < 583:
        hodoscope_channels["BottomZ"] = [_drs(1, 1, 0), _drs(1, 0, 7)]
    elif run_number >= 583:
        hodoscope_channels["BottomZ"] = [_drs(1, 0, 7), _drs(1, 1, 0)]

    if run_number >= 685 and run_number < 1170:
        hodoscope_channels["TopX"] = [_drs(0, 0, 0, brg), _drs(0, 0, 1, brg)]
        hodoscope_channels["TopZ"] = [_drs(0, 0, 2, brg), _drs(0, 0, 3, brg)]
        hodoscope_channels["BottomX"] = [_drs(0, 0, 4, brg), _drs(0, 0, 5, brg)]
        hodoscope_channels["BottomZ"] = [_drs(0, 0, 6, brg), _drs(0, 0, 7, brg)]

    if run_number >= 1170 and run_number < 1327:
        hodoscope_channels = {}
        hodoscope_channels["LR1"] = [_drs(7, 0, 4), _drs(7, 0, 5)]
        hodoscope_channels["UD1"] = [_drs(7, 0, 6), _drs(7, 0, 7)]

    if run_number >= 1342:
        hodoscope_channels = {}
        hodoscope_channels["LR1"] = [_drs(7, 0, 0, brg), _drs(7, 0, 1, brg)]
        hodoscope_channels["UD1"] = [_drs(7, 0, 2, brg), _drs(7, 0, 3, brg)]
        hodoscope_channels["LR2"] = [_drs(7, 0, 4, brg), _drs(7, 0, 5, brg)]
        hodoscope_channels["UD2"] = [_drs(7, 0, 6, brg), _drs(7, 0, 7, brg)]

    return hodoscope_channels


def findFanoutTimeReferenceDelay(channel, run_number=1040):
    if str(run_number) not in triggerdelay.keys():
        return triggerdelay["default"][channel]
    else:
        return triggerdelay[str(run_number)][channel]


def findDRSTriggerMap(channel, run_number=1040):
    result = "_".join(channel.split("_")[:3])
    # print(triggermap.keys())
    if str(run_number) not in triggermap.keys():
        return triggermap["default"][result]
    else:
        return triggermap[str(run_number)][result]


def get_hole_veto_channel(run_number=1184):
    brg = 1 if run_number >= _DRS_BRG_RUN else None
    if run_number < 1183:
        return None
    elif run_number < 1327:
        return _drs(7, 1, 6)
    else:
        return get_service_drs_channels(run_number)["HoleVeto"]


def get_downstream_muon_channel(run_number=1184):
    brg = 1 if run_number >= _DRS_BRG_RUN else None
    if run_number < 1183:
        return None
    else:
        return _drs(7, 1, 0, brg)


def get_downstream_ttu_muon_channel(run_number=1184):
    if run_number < 1183:
        return None
    elif run_number < 1327:
        return None
    else:
        return get_service_drs_channels(run_number)["TTUMuonVeto"]


def get_pre_shower_channel(run_number=1184):
    """
    Returns the pre-shower channel.
    """
    if run_number < 1183:
        return None
    elif run_number < 1327:
        return _drs(7, 1, 1)
    else:
        return get_service_drs_channels(run_number)["PSD"]


def get_cerenkov_counters(run_number=1184):
    """
    Returns a list of Cerenkov counter channels.
    """
    brg = 1 if run_number >= _DRS_BRG_RUN else None
    if run_number < 1183:
        return []
    elif run_number < 1327:
        return [_drs(7, 1, 2), _drs(7, 1, 3), _drs(7, 1, 4), _drs(7, 1, 5)]
    else:
        return [_drs(7, 2, 5, brg), _drs(7, 2, 6, brg), _drs(7, 2, 7, brg)]


def get_mcp_channels(run_number=1184):
    """
    Returns a dict of MCP detector name -> channel branch name.
    """
    if run_number < 1342 or (run_number >= _DRS_BRG_RUN and run_number < 1764):
        # no MCP for runs before Sep 2024 test beam or 2025+ (all channels are signal)
        return {}
    elif run_number < 1600:
        # 4 calo boards, Sep 2025 test beam
        return {
            "MCP_DS_0": _drs(0, 3, 6, None),
            "MCP_DS_1": _drs(1, 3, 6, None),
            "MCP_DS_2": _drs(2, 3, 6, None),
            "MCP_DS_3": _drs(3, 3, 6, None),
            "MCP_US_0": _drs(0, 3, 7, None),
            "MCP_US_1": _drs(1, 3, 7, None),
            "MCP_US_2": _drs(2, 3, 7, None),
            "MCP_US_3": _drs(3, 3, 7, None),
        }
    elif run_number > 1700 and run_number < _DRS_FULL_RUN_TB2026:
        # only going to server drs
        return {
            "MCP_1": _drs(0, 0, 5, 1),
            "MCP_2": _drs(0, 0, 6, 1),
        }
    elif run_number < _DRS_TESTFIBER_RUN_TB2026:
        return {
            "MCP_DS_0": _drs(3, 3, 7, 1),
            "MCP_US_0": _drs(3, 3, 6, 1),
            "MCP_DS_1": _drs(3, 3, 7, 0),
            "MCP_US_1": _drs(3, 3, 6, 0),
        }
    else:
        return {
            "MCP_DS_0": _drs(2, 3, 7, 0),
            "MCP_US_0": _drs(2, 3, 6, 0),
            "MCP_DS_1": _drs(0, 3, 7, 0),
            "MCP_US_1": _drs(0, 3, 6, 0),
        }


def get_service_drs_channels(run_number=1184):
    """
    Returns a list of service DRS channels.
    """
    if run_number < 1183:
        return []
    elif run_number < 1260:
        return [_drs(7, g, c) for g in range(2) for c in range(8)]
    elif run_number < 1600:
        return {
            "DWC1Left":   _drs(7, 0, 0, None),
            "DWC1Right":  _drs(7, 0, 1, None),
            "DWC1Up":     _drs(7, 0, 2, None),
            "DWC1Down":   _drs(7, 0, 3, None),
            "DWC2Left":   _drs(7, 0, 4, None),
            "DWC2Right":  _drs(7, 0, 5, None),
            "DWC2Up":     _drs(7, 0, 6, None),
            "DWC2Down":   _drs(7, 0, 7, None),
            "MuonVeto":   _drs(7, 1, 0, None),
            "PSD":        _drs(7, 1, 1, None),
            "HoleVeto":   _drs(7, 1, 6, None),
            "NC":         _drs(7, 1, 7, None),
            "T3":         _drs(7, 2, 0, None),
            "T4":         _drs(7, 2, 1, None),
            "KT1":        _drs(7, 2, 2, None),
            "KT2":        _drs(7, 2, 3, None),
            "TTUMuonVeto":_drs(7, 2, 4, None),
            "Cer474":     _drs(7, 2, 5, None),
            "Cer519":     _drs(7, 2, 6, None),
            "Cer537":     _drs(7, 2, 7, None),
        }
    else:
        # tb 2026
        channels = {
            "TailCatcher": _drs(0, 0, 0, 1),
            "TTUMuonVeto": _drs(0, 0, 1, 1),
            "Cer474":      _drs(0, 0, 2, 1),
            "Cer519":      _drs(0, 0, 3, 1),
            "Cer537":      _drs(0, 0, 4, 1),
            "HoleVeto":    _drs(0, 0, 7, 1),
        }
        if run_number >= 1824:
            channels["ST1"] = _drs(0, 1, 0, 1)
            channels["ST3"] = _drs(0, 1, 1, 1)
        channels.update(get_mcp_channels(run_number))
        return channels


def get_pid_channels(run_number=1184):
    """Return an OrderedDict of PID detector name -> channel name for service DRS analysis."""
    from collections import OrderedDict
    service = get_service_drs_channels(run_number=run_number)
    pid_dets = [
        "HoleVeto", "PSD", "TTUMuonVeto",
        "Cer474", "Cer519", "Cer537",
        "KT1", "KT2", "T3", "T4",
        "MCP_1", "MCP_2",
        "MCP_DS_0", "MCP_US_0", "MCP_DS_1", "MCP_US_1",
        "ST1", "ST3",
        "TailCatcher"
    ]
    return OrderedDict((det, service.get(det)) for det in pid_dets)


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
    fers_boards = build_fers_boards(run_number=run_number)
    drs_boards = build_drs_boards(run_number=run_number)

    print("FERS Boards:")
    for board_name, board in fers_boards.items():
        print(f"{board_name}: {board}")

    print("\nDRS Boards:")
    for board_name, board in drs_boards.items():
        print(f"{board_name}: {board}")

    print("\nHodoscope Position Channels:")
    hodo_channels = build_hodo_pos_channels(run_number=run_number)
    for hodo_type, channels in hodo_channels.items():
        print(f"{hodo_type}: {channels}")
