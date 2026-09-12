"""Channel maps for CaloX.

The maps themselves live in channels.maps (see that package's docstring for the
split); this module re-exports them so the rest of the code base has one import
site. Prefer importing from here.
"""
from channels.maps.eras import (  # noqa: F401
    _DRS_BRG_RUN, _DRS_6MM_RUN, _DRS_FULL_RUN_TB2026,
    _DRS_TESTFIBER_RUN_TB2026, _DRS_PHASE2_RUN, _for_run,
)
from channels.maps.quartz import (  # noqa: F401
    get_quartz_channel_list, update_quartz_channels,
)
from channels.maps.fers import (  # noqa: F401
    apply_corner_chamfer, build_fers_boards, physical_to_fers_channel,
    _CORNER_CHAMFER, _FERS_LAYOUTS, _fers_layout,
)
from channels.maps.drs import (  # noqa: F401
    build_drs_boards, build_drs_boards_from_fers, dump_drs_boards,
    load_drs_boards, load_drs_fers_lookup,
    _DRS_FERS_CSV_PHASE1, _DRS_FERS_CSV_PHASE2, _DRS_LAYOUTS, _DRS_MAP_DIR,
    _drs_fers_csv, _drs_layout_file,
)
from channels.maps.services import (  # noqa: F401
    build_hodo_pos_channels, build_hodo_trigger_channels,
    build_time_reference_channels, findDRSTriggerMap,
    findFanoutTimeReferenceDelay, get_cerenkov_counters,
    get_downstream_muon_channel, get_downstream_ttu_muon_channel,
    get_hole_veto_channel, get_mcp_channels, get_mcp_reference, get_pid_channels,
    get_pre_shower_channel, get_service_drs_channels,
    _MCP_CHANNELS, _SERVICE_DRS_CHANNELS, _drs, triggerdelay, triggermap,
)


if __name__ == "__main__":
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
