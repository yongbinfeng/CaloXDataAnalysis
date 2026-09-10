"""Channel maps for CaloX, split by concern.

    eras      run-number boundaries shared by every map, and the range lookup
    quartz    the quartz-fibre region of the calorimeter face
    fers      FERS board layouts and face geometry
    drs       DRS calorimeter maps (frozen JSON, and the from-FERS CSV path)
    services  auxiliary channels: time reference, hodoscope, MCP, service DRS

channels.channel_map re-exports all of this and stays the import site for the
rest of the code base.
"""
