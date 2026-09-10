"""DRS calorimeter maps.

Runs from _DRS_FULL_RUN_TB2026 on derive their channel positions from the FERS
map through the DRS->FERS CSVs; earlier eras read a frozen JSON map from
data/channel_maps/drs/.
"""
import csv
import json
import os
import re

from channels.calox_channel import (DRSBoard, DRSChannel,
                                    add_drs_reference_channel)
from channels.maps.eras import (_DRS_FULL_RUN_TB2026, _DRS_PHASE2_RUN,
                                _for_run)
from channels.maps.fers import build_fers_boards, physical_to_fers_channel
from channels.maps.quartz import update_quartz_channels
from utils.data_loader import is_scan_run

from collections import OrderedDict


# CSV mapping each DRS readout channel (DRS_ROOT_Bridge_Mapping) to its FERS
# channel (fers_root_board, Physical-1). Used only for the from-FERS era
# (run >= _DRS_FULL_RUN_TB2026); earlier runs use the function-based board maps.
# The mapping was rewired between phases, so the CSV is run-dependent:
#   1828 <= run < 1896  : DRS_PHASE_1.csv
#   run >= 1896         : DRS_PHASE_2.csv
_CHANNEL_MAP_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "channel_maps")


_DRS_FERS_CSV_PHASE1 = os.path.join(_CHANNEL_MAP_DIR, "DRS_PHASE_1.csv")


_DRS_FERS_CSV_PHASE2 = os.path.join(_CHANNEL_MAP_DIR, "DRS_PHASE_2.csv")


def _drs_fers_csv(run_number):
    """DRS->FERS channel-map CSV for a run (see ranges above)."""
    if run_number is not None and run_number >= _DRS_PHASE2_RUN:
        return _DRS_FERS_CSV_PHASE2
    return _DRS_FERS_CSV_PHASE1


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


# Per-era DRS maps live as JSON under data/channel_maps/drs/. Ranges are
# [run_min, run_max); run_min None means "everything before run_max". Runs from
# _DRS_FULL_RUN_TB2026 on are not listed here: those are derived from the FERS
# map via the DRS->FERS CSVs (build_drs_boards_from_fers).
_DRS_MAP_DIR = os.path.join(_CHANNEL_MAP_DIR, "drs")


_DRS_LAYOUTS = (
    ((None, 685), "run_0_685.json"),
    ((685, 1003), "run_685_1003.json"),
    ((1003, 1033), "run_1003_1033.json"),
    ((1033, 1100), "run_1033_1100.json"),
    ((1173, 1327), "run_1173_1327.json"),
    ((1342, 1700), "run_1342_1700.json"),   # September 2024 test beam
    ((1700, 1748), "run_1700_1748.json"),
    ((1748, 1828), "run_1748_1828.json"),
)


def _drs_layout_file(run_number):
    """Path of the DRS map JSON for a run (see _DRS_LAYOUTS)."""
    fname = _for_run(_DRS_LAYOUTS, run_number)
    if fname is None:
        raise ValueError(f"Unsupported run_number {run_number} for DRS boards.")
    return os.path.join(_DRS_MAP_DIR, fname)


def load_drs_boards(path):
    """Read a DRS map written by dump_drs_boards.

    Quartz flags are not stored: they are a function of the tower position and
    are applied by update_quartz_channels after loading, so that the quartz
    region stays defined in exactly one place (get_quartz_channel_list).
    """
    with open(path, 'r') as f:
        spec = json.load(f)
    boards = {}
    for board in spec["boards"]:
        channels = [
            DRSChannel(c["x"], c["y"], c["cer"], c["channel"], c["group"],
                       board["board_no"], is_amplified=c["amplified"],
                       is6mm=c["is6mm"], is_reference=c["reference"])
            for c in board["channels"]]
        boards[board["key"]] = DRSBoard(board_no=board["board_no"],
                                        channels=channels,
                                        bridge_no=board["bridge_no"])
    return boards


def dump_drs_boards(boards, path, run_range=None, note=None):
    """Write a DRS map in the format load_drs_boards expects.

    Use this to add a new era: build the boards once however the cabling
    dictates, dump them, and add the file to _DRS_LAYOUTS.
    """
    spec = {
        "run_range": list(run_range) if run_range else None,
        "note": note,
        "boards": [
            {"key": key,
             "board_no": int(board.board_no),
             "bridge_no": board.bridge_no,
             "channels": [
                 {"group": int(c.group_no), "channel": int(c.channel_no),
                  "x": c.i_tower_x, "y": c.i_tower_y, "cer": c.isCer,
                  "amplified": c.is_amplified, "is6mm": c.is6mm,
                  "reference": c.is_reference}
                 for c in board.channels]}
            for key, board in boards.items()],
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(spec, f, indent=1)
        f.write("\n")


def build_drs_boards(run_number=316):
    """
    Build a map for ixy and DRS channels.

    Runs from _DRS_FULL_RUN_TB2026 on take their channel positions from the
    FERS map through the DRS->FERS CSVs; earlier runs read the frozen per-era
    JSON map named by _DRS_LAYOUTS.
    """
    if is_scan_run(run_number):
        # no DRS boards in scan runs
        # only FERS
        return {}

    if run_number >= _DRS_FULL_RUN_TB2026:
        # The per-channel bridge numbering and quartz/plastic flags come from
        # the CSV, so the update_quartz_channels post-processing is skipped.
        fersboards = build_fers_boards(run_number=run_number)
        return build_drs_boards_from_fers(
            fersboards, csv_path=_drs_fers_csv(run_number))

    DRSBoards = load_drs_boards(_drs_layout_file(run_number))
    update_quartz_channels(DRSBoards)
    return DRSBoards
