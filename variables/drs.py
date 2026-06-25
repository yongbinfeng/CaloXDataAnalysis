# collect all functions related to DRS here
import re
import json
import os
from channels.channel_map import get_mcp_channels, get_pid_channels, get_service_drs_channels
from configs.selection_config import get_service_drs_cut
from utils.utils import get_channel_var

TS_END = 1024
MCP_REF = "MCP_DS_0"
MCP_REF = "MCP_1"

# Column naming convention — uppercase prefix = scalar, lowercase prefix = RVec array:
#   _ref_TS          scalar: LED discriminator crossing time slice of the reference channel (Channel8)
#   _ref_TS_peak     scalar: peak position time slice of the reference channel
#   _TS_cfd          scalar float: CFD crossing ti
# me slice from compute_cfd_integral
#   _TS_peak         scalar float: peak position time slice from compute_cfd_integral
#   _TS_cfd_ref      scalar: _TS_cfd corrected for the reference channel timing
#   _TS_peak_ref     scalar: _TS_peak corrected for the reference channel timing
#   _TS_cfd_mcp      scalar: _TS_cfd_ref further corrected for MCP timing  [TS]
#   _TS_peak_mcp     scalar: _TS_peak_ref further corrected for MCP timing [TS]
#   _ts_ref          RVec:   per-event TS array shifted by the reference channel time slice
#   _ts_mcp          RVec:   per-event TS array shifted by reference + MCP time slice
#   _Time_cfd_mcp    scalar: _TS_cfd_mcp group-MPV-corrected, converted to ns
#   _time_mcp        RVec:   _ts_mcp group-MPV-corrected, converted to ns


def get_drs_branches(rdf):
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"^DRS_(Brg\d+_)?Board\d+_Group\d+_Channel\d+$")
    drs_branches = [b for b in branches if pattern.search(b)]
    return drs_branches


def get_drs_branches_to_flip(run_number, drs_channels_ref=None, drsboards=None):
    channels_services = list(get_service_drs_channels(run_number).values())
    channels_mcp = list(get_mcp_channels(run_number).values())

    # Service and MCP channels on bridge 0 are not inverted (opposite readout
    # polarity) — EXCEPT for runs > 1839, where all service/MCP channels are
    # flipped. This exception applies ONLY to service/MCP channels; the calo
    # ("real") channels below still need to be double-checked.
    def _bridge_no(name):
        m = re.search(r"_Brg(\d+)_", name)
        return int(m.group(1)) if m else None
    if run_number is not None and run_number >= 1839:
        drs_channels_to_flip = list(channels_services + channels_mcp)
    else:
        drs_channels_to_flip = [c for c in channels_services + channels_mcp
                                if _bridge_no(c) != 0]

    if drs_channels_ref is not None:
        drs_channels_to_flip += drs_channels_ref

    if drsboards is not None:
        for _, DRSBoard in drsboards.items():
            for channel in DRSBoard:
                if (channel.is6mm and channel.is_amplified):
                    drs_channels_to_flip.append(
                        channel.get_channel_name(blsub=False))

    return drs_channels_to_flip


def get_ref_ts_name(drs_channel_name, use_peak=False):
    """Return the reference-channel column name (_ref_TS or _ref_TS_peak) for a given DRS channel."""
    board_group_name = drs_channel_name.rsplit("_Channel", 1)[0]
    return f"{board_group_name}_ref_TS_peak" if use_peak else f"{board_group_name}_ref_TS"


def subtract_baseline(rdf, drs_branches, drs_channels_to_flip=None,
                      linear=False, pre_window=(0, 120), post_window=(800, 1000)):
    """Define {ch}_bl and {ch}_blsub for each DRS branch.

    By default (linear=False) the baseline is a single constant — the median of
    time slices [0, 200) — subtracted from the whole record.

    With linear=True a two-anchor *linear* baseline is used instead: the medians
    of a clean pre-pulse window and a clean post-pulse window define a line that
    is subtracted sample-by-sample. This removes the post-pulse pedestal droop
    seen on large pulses (e.g. service-DRS HoleVeto), at the cost of one extra
    median per channel. pre_window/post_window must straddle the pulse region.
    """
    if drs_channels_to_flip is None:
        drs_channels_to_flip = []
    existing = [str(c) for c in rdf.GetColumnNames()]
    if "ts" not in existing:
        rdf = rdf.Define("ts", "FillIndices(1024)")

    pre_c = 0.5 * (pre_window[0] + pre_window[1])
    post_c = 0.5 * (post_window[0] + post_window[1])

    for channel_name in drs_branches:
        if linear:
            rdf = rdf.Define(
                f"{channel_name}_bl_pre",
                f"compute_baseline_median({channel_name}, {pre_window[0]}, {pre_window[1]})")
            rdf = rdf.Define(
                f"{channel_name}_bl_post",
                f"compute_baseline_median({channel_name}, {post_window[0]}, {post_window[1]})")
            # constant scalar kept for downstream code (pre-pulse level)
            rdf = rdf.Define(f"{channel_name}_bl", f"{channel_name}_bl_pre")
            bl_lo, bl_hi = pre_window
            # per-sample baseline line straddling the pulse
            rdf = rdf.Define(
                f"{channel_name}_bl_line",
                f"{channel_name}_bl_pre + (({channel_name}_bl_post - {channel_name}_bl_pre)"
                f" / double({post_c} - {pre_c})) * (ts - {pre_c})")
            baseline_expr = f"{channel_name}_bl_line"
        else:
            rdf = rdf.Define(
                f"{channel_name}_bl",
                f"compute_baseline_median({channel_name}, 0, 200)")
            bl_lo, bl_hi = 0, 200
            baseline_expr = f"{channel_name}_bl"

        # per-event baseline noise (RMS over the baseline window)
        rdf = rdf.Define(
            f"{channel_name}_blrms",
            f"compute_baseline_rms({channel_name}, {bl_lo}, {bl_hi})")

        if channel_name in drs_channels_to_flip:
            rdf = rdf.Define(
                f"{channel_name}_blsub",
                f"-({channel_name} - {baseline_expr})")
        else:
            rdf = rdf.Define(
                f"{channel_name}_blsub",
                f"{channel_name} - {baseline_expr}")

    return rdf


def process_reference_channels(rdf, drs_channels_ref):
    """Define _ref_TS and _ref_TS_peak for each reference channel (Channel8)."""
    for channel_name in drs_channels_ref:
        if not channel_name.endswith("Channel8"):
            raise ValueError(
                f"Reference channel {channel_name} does not end with 'Channel8'. Check the channel naming convention.")
        board_group_name = channel_name.replace("_Channel8", "")
        rdf = rdf.Define(
            f"{channel_name}_REFHit", f"process_dynamic_led({channel_name}_blsub)")
        rdf = rdf.Define(f"{board_group_name}_ref_TS",
                         f"{channel_name}_REFHit.time_slice")
        rdf = rdf.Define(f"{board_group_name}_ref_TS_peak",
                         f"{channel_name}_REFHit.peak_position")

    return rdf


def process_mcp_channels(rdf, run_number):
    """Define _TS_cfd, _TS_peak, and their ref-corrected variants for each MCP detector."""
    existing = [str(c) for c in rdf.GetColumnNames()]
    map_mcp_channels = get_mcp_channels(run_number)
    for det, channel_name in map_mcp_channels.items():
        channel_name_blsub = channel_name + "_blsub"
        if channel_name_blsub not in existing:
            continue
        ts_begin, ts_end, window_pre, window_post, _, _ = get_service_drs_cut(
            det, run_number)

        rdf = rdf.Define(f"{det}_CFD",
                         f"compute_cfd_integral({channel_name_blsub}, 1, 10, {ts_begin}, {ts_end}, {window_pre}, {window_post})")
        rdf = rdf.Define(f"{det}_energy", f"{det}_CFD.energy")
        rdf = rdf.Define(f"{det}_peak_value", f"{det}_CFD.peak_value")
        rdf = rdf.Define(f"{det}_TS_cfd", f"{det}_CFD.time_slice")
        rdf = rdf.Define(f"{det}_TS_peak", f"{det}_CFD.peak_position")
        rdf = rdf.Define(f"{det}_integral_to_peak",
                         f"{det}_CFD.integral_to_peak")

        ref_TS = get_ref_ts_name(channel_name, use_peak=False)
        rdf = rdf.Define(f"{det}_TS_cfd_ref",
                         f"790 + {det}_TS_cfd - {ref_TS}")

        ref_TS_peak = get_ref_ts_name(channel_name, use_peak=True)
        rdf = rdf.Define(f"{det}_TS_peak_ref",
                         f"{det}_TS_peak - {ref_TS_peak}")

    return rdf


def process_pid_channels(rdf, run_number):
    """Define _CFD-derived columns for each PID service-DRS detector."""
    existing = [str(c) for c in rdf.GetColumnNames()]
    for det, channel in get_pid_channels(run_number).items():
        if channel is None:
            continue
        channel_blsub = channel + "_blsub"
        if channel_blsub not in existing:
            continue
        # Skip detectors already processed elsewhere: when MCP is enabled,
        # MCP_1/MCP_2 are defined by process_mcp_channels (which runs first)
        # and also appear in the PID list, so their CFD columns already exist.
        if f"{det}_CFD" in existing:
            continue
        ts_begin, ts_end, window_pre, window_post, _, _ = get_service_drs_cut(
            det, run_number)

        rdf = rdf.Define(f"{det}_CFD",
                         f"compute_cfd_integral({channel_blsub}, 1, 10.0, {ts_begin}, {ts_end}, {window_pre}, {window_post})")
        rdf = rdf.Define(f"{det}_energy",          f"{det}_CFD.energy")
        rdf = rdf.Define(f"{det}_peak_value",
                         f"{det}_CFD.peak_value")
        rdf = rdf.Define(f"{det}_TS_cfd",
                         f"{det}_CFD.time_slice")
        rdf = rdf.Define(f"{det}_TS_peak",
                         f"{det}_CFD.peak_position")
        rdf = rdf.Define(f"{det}_integral_to_peak",
                         f"{det}_CFD.integral_to_peak")

    return rdf


def get_arr_name(drs_channel_name, use_mcp=False, in_ns=False):
    """Return the time array column name for a DRS channel.

    in_ns=False (default): board-group TS array (_ts_ref or _ts_mcp) [TS]
    in_ns=True:            per-channel corrected time array (_time_mcp) [ns]
    use_mcp: selects _ts_mcp vs _ts_ref (only relevant when in_ns=False)
    """
    if in_ns:
        return f"{drs_channel_name}_time_mcp"
    board_group = drs_channel_name.rsplit("_Channel", 1)[0]
    return f"{board_group}_ts_mcp" if use_mcp else f"{board_group}_ts_ref"


def update_ts(rdf, drs_channels_ref, mcp_det="MCP_US_0"):
    """Define per-event calibrated ts arrays (_ts_ref and _ts_mcp) for each board.

    When mcp_det is None (no MCP available), _ts_mcp is set equal to _ts_ref
    so that downstream code using _ts_mcp still works with zero MCP correction.
    """
    for channel_name in drs_channels_ref:
        board_group_name = channel_name.replace("_Channel8", "")
        ref_TS = f"{board_group_name}_ref_TS"

        ts_ref = get_arr_name(channel_name)
        rdf = rdf.Define(ts_ref, f"ts - {ref_TS} + 790")
        if mcp_det is not None:
            rdf = rdf.Define(
                get_arr_name(channel_name, use_mcp=True),
                f"{ts_ref} - {mcp_det}_TS_cfd_ref + 500")
        else:
            rdf = rdf.Define(get_arr_name(channel_name, use_mcp=True), ts_ref)
    return rdf


def process_drs_channels(rdf, drs_channel_names, mcp_det="MCP_US_0"):
    """Define _TS_cfd, _TS_peak, and their ref/mcp-corrected variants for each DRS channel."""
    for channel_name in drs_channel_names:
        channelName_blsub = channel_name + "_blsub"

        #rdf = rdf.Define(f"{channel_name}_CFD",
        #                 f"compute_cfd_integral({channelName_blsub}, 1, 3.0, 200, 220, 3, 15)")
        rdf = rdf.Define(f"{channel_name}_CFD",
                         f"compute_cfd_integral({channelName_blsub})")
        rdf = rdf.Define(f"{channel_name}_energy",
                         f"{channel_name}_CFD.energy")
        rdf = rdf.Define(f"{channel_name}_TS_cfd",
                         f"{channel_name}_CFD.time_slice")
        rdf = rdf.Define(f"{channel_name}_TS_peak",
                         f"{channel_name}_CFD.peak_position")
        rdf = rdf.Define(f"{channel_name}_peak_value",
                         f"{channel_name}_CFD.peak_value")
        rdf = rdf.Define(f"{channel_name}_integral_to_peak",
                         f"{channel_name}_CFD.integral_to_peak")

        ref_TS = get_ref_ts_name(channel_name, use_peak=False)
        rdf = rdf.Define(f"{channel_name}_TS_cfd_ref",
                         f"790 + {channel_name}_TS_cfd - {ref_TS}")
        ref_TS_peak = get_ref_ts_name(channel_name, use_peak=True)
        rdf = rdf.Define(f"{channel_name}_TS_peak_ref",
                         f"790 + {channel_name}_TS_peak - {ref_TS_peak}")

        if mcp_det is not None:
            rdf = rdf.Define(f"{channel_name}_TS_cfd_mcp",
                             f"{channel_name}_TS_cfd_ref - {mcp_det}_TS_cfd_ref + 500")
            rdf = rdf.Define(f"{channel_name}_TS_peak_mcp",
                             f"{channel_name}_TS_peak_ref - {mcp_det}_TS_cfd_ref + 500")
        else:
            rdf = rdf.Define(f"{channel_name}_TS_cfd_mcp",
                             f"{channel_name}_TS_cfd_ref")
            rdf = rdf.Define(f"{channel_name}_TS_peak_mcp",
                             f"{channel_name}_TS_peak_ref")
    return rdf


def get_psd_energy_deposit(rdf, run_number, TS_start=100, TS_end=400):
    channel_name = get_service_drs_channels(run_number).get("PSD")
    if channel_name is None:
        return rdf
    existing = [str(c) for c in rdf.GetColumnNames()]
    if f"{channel_name}_blsub" not in existing:
        return rdf
    rdf = rdf.Define(
        "PSD_Sum", f"SumRange({channel_name}_blsub, {TS_start}, {TS_end})")
    rdf = rdf.Define(
        "PSD_Energy_Cer", "PSD_Sum * 5.15852e-05 * (-1)").Define(
        "PSD_Energy_Sci", "PSD_Sum * 6.53095e-05 * (-1)")
    return rdf


def process_drs_data(rdf, run_number, drsboards, do_mcp=True,
                     linear_baseline=False,
                     baseline_pre_window=(0, 120),
                     baseline_post_window=(800, 1000)):
    drs_branches = get_drs_branches(rdf)
    drs_channels_ref = [ch for ch in drs_branches if ch.endswith("Channel8")]
    drs_branches_to_flip = get_drs_branches_to_flip(
        run_number, drs_channels_ref=drs_channels_ref, drsboards=drsboards)

    # Only correct for MCP timing if the reference MCP channel actually exists
    # for this run (e.g. runs >= _DRS_BRG_RUN have no MCP channels).
    mcp_available = MCP_REF in get_mcp_channels(run_number)
    mcp_det = MCP_REF if (do_mcp and mcp_available) else None

    rdf = subtract_baseline(
        rdf, drs_branches, drs_channels_to_flip=drs_branches_to_flip,
        linear=linear_baseline, pre_window=baseline_pre_window,
        post_window=baseline_post_window)
    rdf = process_reference_channels(rdf, drs_channels_ref)
    if do_mcp:
        rdf = process_mcp_channels(rdf, run_number)
        rdf = process_pid_channels(rdf, run_number)
    rdf = update_ts(rdf, drs_channels_ref, mcp_det=mcp_det)
    drs_channels_physics = [
        ch for ch in drs_branches if not ch.endswith("Channel8")]
    rdf = process_drs_channels(rdf, drs_channels_physics, mcp_det=mcp_det)
    rdf = define_time_ns(rdf, drsboards)
    return rdf


_MPV_GROUP_JSON_PATH = "data/drs/drs_cfd_mpv_at_1500mm_by_group.json"
_TS_TO_NS = 0.2


def _grp_key(chan, board):
    """Return the group JSON key for a channel (board needed for B5 special case)."""
    var = get_channel_var(chan)
    size = "6mm" if chan.is6mm else "3mm"
    cer_key = "Cer" if var in ("CerQuartz", "CerPlastic") else "Sci"
    if not chan.is6mm and var == "CerQuartz" and board.board_no == 5:
        return "3mm_CerQuartz_B5"
    return f"{size}_{cer_key}"


def subtract_type_mpv(drsboards, raw_mpv_map, json_path=_MPV_GROUP_JSON_PATH):
    """Subtract the group-average MPV at 1500 mm from per-channel raw MPV values.

    Uses the group JSON (drs_cfd_mpv_at_1500mm_by_group.json) with keys:
        3mm_Cer, 3mm_Sci, 6mm_Cer, 6mm_Sci, 3mm_CerQuartz_B5

    Both CerQuartz and CerPlastic channels use the inclusive Cer reference.
    Board-5 3mm CerQuartz uses the dedicated 3mm_CerQuartz_B5 reference.
    Returns a dict {ch: corrected_value} or None if the JSON file does not exist.
    """
    if not os.path.exists(json_path):
        return None

    with open(json_path) as _f:
        group_mpv = json.load(_f)

    corr_map = {}
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            ch = chan.get_channel_name(blsub=False)
            if ch not in raw_mpv_map:
                continue
            key = _grp_key(chan, board)
            if key not in group_mpv:
                continue
            corr_map[ch] = raw_mpv_map[ch] - group_mpv[key]

    return corr_map


def define_time_ns(rdf, drsboards, json_path=_MPV_GROUP_JSON_PATH):
    """Define per-channel time_cfd_mcp and time_mcp columns in physical units [ns].

    For each non-reference DRS channel adds:
        {ch}_Time_cfd_mcp  scalar: ({ch}_TS_cfd_mcp - group_mpv) * 0.2  [ns]
        {ch}_time_mcp      RVec:   ({board_group}_ts_mcp - group_mpv) * 0.2  [ns]

    Channels whose group key is absent from the JSON are skipped.
    """
    if not os.path.exists(json_path):
        return rdf

    with open(json_path) as _f:
        group_mpv = json.load(_f)

    existing = [str(c) for c in rdf.GetColumnNames()]
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            key = _grp_key(chan, board)
            ref_ts = group_mpv.get(key)
            if ref_ts is None:
                continue
            ch = chan.get_channel_name(blsub=False)
            if f"{ch}_TS_cfd_mcp" not in existing:
                continue
            rdf = rdf.Define(f"{ch}_Time_cfd_mcp",
                             f"({ch}_TS_cfd_mcp - {ref_ts}) * {_TS_TO_NS}")
            rdf = rdf.Define(get_arr_name(ch, in_ns=True),
                             f"({get_arr_name(ch, use_mcp=True)} - {ref_ts}) * {_TS_TO_NS}")

    return rdf
