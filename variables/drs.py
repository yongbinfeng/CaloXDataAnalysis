# collect all functions related to DRS here
import re
from channels.channel_map import get_mcp_channels, get_service_drs_channels

TS_END = 1024

# Column naming convention (all time quantities are time-slice indices, not physical time):
#   _ref_TS          scalar: LED discriminator crossing time slice of the reference channel (Channel8)
#   _ref_TS_peak     scalar: peak position time slice of the reference channel
#   _TS_cfd          scalar float: CFD crossing time slice from compute_cfd_integral
#   _TS_peak         scalar float: peak position time slice from compute_cfd_integral
#   _TS_cfd_ref      scalar: _TS_cfd corrected for the reference channel timing
#   _TS_peak_ref     scalar: _TS_peak corrected for the reference channel timing
#   _TS_cfd_mcp      scalar: _TS_cfd_ref further corrected for MCP timing
#   _TS_peak_mcp     scalar: _TS_peak_ref further corrected for MCP timing
#   _ts_ref          RVec:   per-event TS array shifted by the reference channel time slice
#   _ts_mcp          RVec:   per-event TS array shifted by reference + MCP time slice


def get_drs_branches(rdf):
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"^DRS_Board\d+_Group\d+_Channel\d+$")
    drs_branches = [b for b in branches if pattern.search(b)]
    return drs_branches


def get_drs_branches_to_flip(run_number, drs_channels_ref=None, drsboards=None):
    channels_services = get_service_drs_channels(run_number).values()
    channels_mcp = get_mcp_channels(run_number).values()

    drs_channels_to_flip = list(channels_services) + list(channels_mcp)
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


def subtract_baseline(rdf, drs_branches, drs_channels_to_flip=None):
    existing = [str(c) for c in rdf.GetColumnNames()]
    if "ts" not in existing:
        rdf = rdf.Define("ts", "FillIndices(1024)")

    for channel_name in drs_branches:
        rdf = rdf.Define(
            f"{channel_name}_bl",
            f"compute_baseline_median({channel_name}, 0, 200)"
        )
        if channel_name in drs_channels_to_flip:
            rdf = rdf.Define(
                f"{channel_name}_blsub",
                f"-({channel_name} - {channel_name}_bl)"
            )
        else:
            rdf = rdf.Define(
                f"{channel_name}_blsub",
                f"{channel_name} - {channel_name}_bl"
            )

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
    map_mcp_channels = get_mcp_channels(run_number)
    for det, channel_name in map_mcp_channels.items():
        channel_name_blsub = channel_name + "_blsub"

        rdf = rdf.Define(f"{det}_CFD",
                         f"compute_cfd_integral({channel_name_blsub}, 1, 200.0)")
        rdf = rdf.Define(f"{det}_TS_cfd", f"{det}_CFD.time_slice")
        rdf = rdf.Define(f"{det}_TS_peak", f"{det}_CFD.peak_position")

        ref_TS = get_ref_ts_name(channel_name, use_peak=False)
        rdf = rdf.Define(f"{det}_TS_cfd_ref",
                         f"790 + {det}_TS_cfd - {ref_TS}")

        ref_TS_peak = get_ref_ts_name(channel_name, use_peak=True)
        rdf = rdf.Define(f"{det}_TS_peak_ref",
                         f"{det}_TS_peak - {ref_TS_peak}")

    return rdf


def get_ts_arr_name(drs_channel_name, use_mcp=False):
    """Return the calibrated ts array column name (_ts_ref or _ts_mcp) for a given DRS channel."""
    board_group_name = drs_channel_name.rsplit("_Channel", 1)[0]
    return f"{board_group_name}_ts_mcp" if use_mcp else f"{board_group_name}_ts_ref"


def update_ts(rdf, drs_channels_ref, mcp_det="MCP_US_0"):
    """Define per-event calibrated ts arrays (_ts_ref and _ts_mcp) for each board."""
    for channel_name in drs_channels_ref:
        board_group_name = channel_name.replace("_Channel8", "")
        ref_TS = f"{board_group_name}_ref_TS"

        rdf = rdf.Define(get_ts_arr_name(channel_name, use_mcp=False),
                         f"ts - {ref_TS} + 790")
        if mcp_det is not None:
            rdf = rdf.Define(
                get_ts_arr_name(channel_name, use_mcp=True),
                f"ts - {ref_TS} - {mcp_det}_TS_cfd_ref + 1300")
    return rdf


def process_drs_channels(rdf, drs_channel_names, mcp_det="MCP_US_0"):
    """Define _TS_cfd, _TS_peak, and their ref/mcp-corrected variants for each DRS channel."""
    for channel_name in drs_channel_names:
        channelName_blsub = channel_name + "_blsub"

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
    return rdf


def get_psd_energy_deposit(rdf, run_number, TS_start=100, TS_end=400):
    channel_name = get_service_drs_channels(run_number).get("PSD")
    rdf = rdf.Define(
        "PSD_Sum", f"SumRange({channel_name}_blsub, {TS_start}, {TS_end})")
    rdf = rdf.Define(
        "PSD_Energy_Cer", "PSD_Sum * 5.15852e-05 * (-1)").Define(
        "PSD_Energy_Sci", "PSD_Sum * 6.53095e-05 * (-1)")
    return rdf


def process_drs_data(rdf, run_number, drsboards):
    drs_branches = get_drs_branches(rdf)
    drs_channels_ref = [ch for ch in drs_branches if ch.endswith("Channel8")]
    drs_branches_to_flip = get_drs_branches_to_flip(
        run_number, drs_channels_ref=drs_channels_ref, drsboards=drsboards)

    rdf = subtract_baseline(
        rdf, drs_branches, drs_channels_to_flip=drs_branches_to_flip)
    rdf = process_reference_channels(rdf, drs_channels_ref)
    rdf = process_mcp_channels(rdf, run_number)
    rdf = update_ts(rdf, drs_channels_ref, mcp_det="MCP_US_0")
    drs_channels_physics = [
        ch for ch in drs_branches if not ch.endswith("Channel8")]
    rdf = process_drs_channels(rdf, drs_channels_physics, mcp_det="MCP_US_0")
    return rdf
