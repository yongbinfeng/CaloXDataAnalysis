"""Auxiliary (non-calorimeter) channels.

Time reference, hodoscope, MCP and the service DRS board. Each table maps a
run range to named channels given as (board, group, channel, bridge) tuples;
_resolve turns those into DRS branch names for a particular run.
"""
import json

from collections import OrderedDict

from channels.maps.eras import (_DRS_BRG_RUN, _DRS_FULL_RUN_TB2026,
                                _DRS_TESTFIBER_RUN_TB2026, _for_run)
from utils.data_loader import is_scan_run


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


# Auxiliary (non-calorimeter) DRS channels: time reference, hodoscope, MCP and
# the service DRS board. Each block below is a table of
# (run_min, run_max) -> channels, with the channels given as
# (board, group, channel, bridge) and resolved to branch names by _resolve.
# Ranges are [run_min, run_max); None means open-ended on that side.
_AUTO_BRG = "auto"   # bridge follows the run: 1 from _DRS_BRG_RUN on, else none


def _resolve(spec, run_number):
    """Turn a (board, group, channel, bridge) tuple into a DRS branch name."""
    board, group, channel, brg = spec
    if brg == _AUTO_BRG:
        brg = 1 if run_number >= _DRS_BRG_RUN else None
    return _drs(board, group, channel, brg)


_TIME_REFERENCE_CHANNELS = (
    ((None, 685), [(0, 3, 7, None), (2, 3, 7, None), (1, 0, 0, None)]),
    ((685, None), [(1, 0, 8, _AUTO_BRG), (1, 1, 8, _AUTO_BRG),
                   (1, 2, 8, _AUTO_BRG), (1, 3, 8, _AUTO_BRG),
                   (2, 0, 8, _AUTO_BRG), (2, 1, 8, _AUTO_BRG),
                   (2, 2, 8, _AUTO_BRG), (2, 3, 8, _AUTO_BRG),
                   (0, 0, 8, _AUTO_BRG)]),
)


_HODO_TRIGGER_CHANNELS = (
    ((None, 685), [(1, 2, 0, None), (1, 2, 1, None)]),
    ((685, None), [(0, 2, 0, _AUTO_BRG), (0, 2, 1, _AUTO_BRG)]),
)


# Runs 1327-1342 fall back to the pre-685 hodoscope wiring: that is what the
# original if-chain did, and build_drs_boards has no map for those runs anyway.
_HODO_POS_EARLY = {
    "TopX": [(1, 0, 1, None), (1, 0, 2, None)],
    "TopZ": [(1, 0, 3, None), (1, 0, 4, None)],
    "BottomX": [(1, 0, 5, None), (1, 0, 6, None)],
    "BottomZ": [(1, 1, 0, None), (1, 0, 7, None)],
}


_HODO_POS_LEGACY = dict(_HODO_POS_EARLY,
                        BottomZ=[(1, 0, 7, None), (1, 1, 0, None)])


_HODO_POS_CHANNELS = (
    ((None, 583), _HODO_POS_EARLY),
    ((583, 685), _HODO_POS_LEGACY),
    ((685, 1170), {
        "TopX": [(0, 0, 0, _AUTO_BRG), (0, 0, 1, _AUTO_BRG)],
        "TopZ": [(0, 0, 2, _AUTO_BRG), (0, 0, 3, _AUTO_BRG)],
        "BottomX": [(0, 0, 4, _AUTO_BRG), (0, 0, 5, _AUTO_BRG)],
        "BottomZ": [(0, 0, 6, _AUTO_BRG), (0, 0, 7, _AUTO_BRG)],
    }),
    ((1170, 1327), {
        "LR1": [(7, 0, 4, None), (7, 0, 5, None)],
        "UD1": [(7, 0, 6, None), (7, 0, 7, None)],
    }),
    ((1327, 1342), _HODO_POS_LEGACY),
    ((1342, None), {
        "LR1": [(7, 0, 0, _AUTO_BRG), (7, 0, 1, _AUTO_BRG)],
        "UD1": [(7, 0, 2, _AUTO_BRG), (7, 0, 3, _AUTO_BRG)],
        "LR2": [(7, 0, 4, _AUTO_BRG), (7, 0, 5, _AUTO_BRG)],
        "UD2": [(7, 0, 6, _AUTO_BRG), (7, 0, 7, _AUTO_BRG)],
    }),
)


_MCP_TB2026 = {
    "MCP_DS_0": (3, 3, 7, 1), "MCP_US_0": (3, 3, 6, 1),
    "MCP_DS_1": (3, 3, 7, 0), "MCP_US_1": (3, 3, 6, 0),
}


_MCP_CHANNELS = (
    ((None, 1342), {}),                 # before the Sep 2024 test beam
    ((1342, 1600), {                    # 4 calo boards, Sep 2024 test beam
        "MCP_DS_0": (0, 3, 6, None), "MCP_DS_1": (1, 3, 6, None),
        "MCP_DS_2": (2, 3, 6, None), "MCP_DS_3": (3, 3, 6, None),
        "MCP_US_0": (0, 3, 7, None), "MCP_US_1": (1, 3, 7, None),
        "MCP_US_2": (2, 3, 7, None), "MCP_US_3": (3, 3, 7, None),
    }),
    ((1600, 1700), _MCP_TB2026),        # kept from the original if-chain; no
                                        # real runs here (see _DRS_LAYOUTS)
    ((1700, 1764), {}),                 # all channels are signal
    ((1764, 1828), {"MCP_1": (0, 0, 5, 1), "MCP_2": (0, 0, 6, 1)}),
    ((_DRS_FULL_RUN_TB2026, _DRS_TESTFIBER_RUN_TB2026), _MCP_TB2026),
    ((_DRS_TESTFIBER_RUN_TB2026, None), {
        "MCP_DS_0": (2, 3, 7, 0), "MCP_US_0": (2, 3, 6, 0),
        "MCP_DS_1": (0, 3, 7, 0), "MCP_US_1": (0, 3, 6, 0),
    }),
)


# Service DRS board. For runs 1183-1260 the original code listed these 16
# channels unnamed; the names below are carried over from the 1260-1600 table,
# which the accessors of that era corroborate for MuonVeto (7,1,0), PSD
# (7,1,1) and HoleVeto (7,1,6). The four Cerenkov counters of that era sit at
# (7,1,2)-(7,1,5) per the original get_cerenkov_counters, but which optical
# filter each one carried is not recorded, hence Cer1-Cer4.
_SERVICE_DRS_CHANNELS = (
    ((None, 1183), {}),
    ((1183, 1260), {
        "DWC1Left": (7, 0, 0, None), "DWC1Right": (7, 0, 1, None),
        "DWC1Up": (7, 0, 2, None), "DWC1Down": (7, 0, 3, None),
        "DWC2Left": (7, 0, 4, None), "DWC2Right": (7, 0, 5, None),
        "DWC2Up": (7, 0, 6, None), "DWC2Down": (7, 0, 7, None),
        "MuonVeto": (7, 1, 0, None), "PSD": (7, 1, 1, None),
        "Cer1": (7, 1, 2, None), "Cer2": (7, 1, 3, None),
        "Cer3": (7, 1, 4, None), "Cer4": (7, 1, 5, None),
        "HoleVeto": (7, 1, 6, None), "NC": (7, 1, 7, None),
    }),
    ((1260, 1600), {
        "DWC1Left": (7, 0, 0, None), "DWC1Right": (7, 0, 1, None),
        "DWC1Up": (7, 0, 2, None), "DWC1Down": (7, 0, 3, None),
        "DWC2Left": (7, 0, 4, None), "DWC2Right": (7, 0, 5, None),
        "DWC2Up": (7, 0, 6, None), "DWC2Down": (7, 0, 7, None),
        "MuonVeto": (7, 1, 0, None), "PSD": (7, 1, 1, None),
        "HoleVeto": (7, 1, 6, None), "NC": (7, 1, 7, None),
        "T3": (7, 2, 0, None), "T4": (7, 2, 1, None),
        "KT1": (7, 2, 2, None), "KT2": (7, 2, 3, None),
        "TTUMuonVeto": (7, 2, 4, None), "Cer474": (7, 2, 5, None),
        "Cer519": (7, 2, 6, None), "Cer537": (7, 2, 7, None),
    }),
    ((1600, 1824), {                    # tb 2026
        "TailCatcher": (0, 0, 0, 1), "TTUMuonVeto": (0, 0, 1, 1),
        "Cer474": (0, 0, 2, 1), "Cer519": (0, 0, 3, 1),
        "Cer537": (0, 0, 4, 1), "HoleVeto": (0, 0, 7, 1),
    }),
    ((1824, None), {                    # tb 2026, scintillator telescope added
        "TailCatcher": (0, 0, 0, 1), "TTUMuonVeto": (0, 0, 1, 1),
        "Cer474": (0, 0, 2, 1), "Cer519": (0, 0, 3, 1),
        "Cer537": (0, 0, 4, 1), "HoleVeto": (0, 0, 7, 1),
        "ST1": (0, 1, 0, 1), "ST3": (0, 1, 1, 1),
    }),
)


# From this run on the service DRS board also carries the MCP channels.
_SERVICE_DRS_WITH_MCP_RUN = 1600


def build_time_reference_channels(run_number=316):
    """
    Deprecated: use the reference channels in the DRS boards instead.
    """
    if is_scan_run(run_number):
        # no time reference channels in scan runs, since no drs boards
        return []
    specs = _for_run(_TIME_REFERENCE_CHANNELS, run_number)
    if specs is None:
        raise ValueError(
            f"Unsupported run_number {run_number} for time reference channels.")
    return [_resolve(s, run_number) for s in specs]


def build_hodo_trigger_channels(run_number=316):
    """
    Returns a list of hodoscope trigger channels.
    """
    if is_scan_run(run_number):
        # no hodoscope trigger channels in scan runs, since no drs boards
        return []
    specs = _for_run(_HODO_TRIGGER_CHANNELS, run_number)
    if specs is None:
        raise ValueError(
            f"Unsupported run_number {run_number} for hodoscope trigger channels.")
    return [_resolve(s, run_number) for s in specs]


def build_hodo_pos_channels(run_number=316):
    """
    Returns a dictionary containing the hodoscope channels for the position measurements
    """
    if is_scan_run(run_number):
        # no hodoscope position channels in scan runs, since no drs boards
        return {}
    layout = _for_run(_HODO_POS_CHANNELS, run_number, {})
    return {plane: [_resolve(s, run_number) for s in specs]
            for plane, specs in layout.items()}


def findFanoutTimeReferenceDelay(channel, run_number=1040):
    if str(run_number) not in triggerdelay.keys():
        return triggerdelay["default"][channel]
    else:
        return triggerdelay[str(run_number)][channel]


def findDRSTriggerMap(channel, run_number=1040):
    result = "_".join(channel.split("_")[:3])
    if str(run_number) not in triggermap.keys():
        return triggermap["default"][result]
    else:
        return triggermap[str(run_number)][result]


def get_mcp_channels(run_number=1184):
    """
    Returns a dict of MCP detector name -> channel branch name.
    """
    return {name: _resolve(spec, run_number)
            for name, spec in _for_run(_MCP_CHANNELS, run_number, {}).items()}


def get_service_drs_channels(run_number=1184):
    """
    Returns a dict of service detector name -> channel branch name.
    """
    channels = {name: _resolve(spec, run_number)
                for name, spec in
                _for_run(_SERVICE_DRS_CHANNELS, run_number, {}).items()}
    if run_number >= _SERVICE_DRS_WITH_MCP_RUN:
        channels.update(get_mcp_channels(run_number))
    return channels


def get_hole_veto_channel(run_number=1184):
    return get_service_drs_channels(run_number).get("HoleVeto")


def get_downstream_muon_channel(run_number=1184):
    return get_service_drs_channels(run_number).get("MuonVeto")


def get_downstream_ttu_muon_channel(run_number=1184):
    return get_service_drs_channels(run_number).get("TTUMuonVeto")


def get_pre_shower_channel(run_number=1184):
    """
    Returns the pre-shower channel.
    """
    return get_service_drs_channels(run_number).get("PSD")


def get_cerenkov_counters(run_number=1184):
    """
    Returns a list of Cerenkov counter channels.
    """
    return [channel for name, channel
            in get_service_drs_channels(run_number).items()
            if name.startswith("Cer")]


def get_pid_channels(run_number=1184):
    """Return an OrderedDict of PID detector name -> channel name for service DRS analysis."""
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
