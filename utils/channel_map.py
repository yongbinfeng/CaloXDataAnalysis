# map of fers channel to drs channel
A5202_map_tmp: list[int] = [4, 2, 3, 1, 8, 6, 7, 5,
                            11, 9, 12, 10, 15, 13, 16, 14,
                            20, 18, 19, 17, 24, 22, 23, 21,
                            27, 25, 28, 26, 31, 29, 32, 30,
                            36, 34, 35, 33, 40, 38, 39, 37,
                            43, 41, 44, 42, 47, 45, 48, 46,
                            52, 50, 51, 49, 56, 54, 55, 53,
                            59, 57, 60, 58, 63, 61, 64, 62]
A5202_map = []
for i in range(0, 64):
    A5202_map.append(A5202_map_tmp[i] - 1)  # convert to zero-based index


def build_map_Cer_Sci():
    """
    Build a map for CER and SCI channels in FERS board 1
    """
    map_Cer_Sci = {}
    for iCer in [0, 1, 2, 3,
                 8, 9, 10, 11,
                 16, 17, 18, 19,
                 24, 25, 26, 27,
                 32, 33, 34, 35,
                 40, 41, 42, 43,
                 48, 49, 50, 51,
                 56, 57, 58, 59]:
        iSci = iCer + 4

        index_Cer = A5202_map.index(iCer)
        index_Sci = A5202_map.index(iSci)
        map_Cer_Sci[index_Cer] = index_Sci
    return map_Cer_Sci


def build_map_ixy_DRS():
    """
    Build a map for ixy and DRS channels.
    """
    map_ixy_DRS = {}
    for i in range(0, 64):
        iy = i // 8
        ix = i % 8
        if ix >= 4:
            ix -= 4
        # mirror, flip
        ix = 3 - ix
        iy = 7 - iy
        map_ixy_DRS[i] = (ix, iy)
    return map_ixy_DRS


def build_map_FERS1_ixy():
    """
    Build a map for ixy and FERS channels for board 1.
    """
    map_ixy_DRS = build_map_ixy_DRS()
    map_ixy_FERS = {}
    for ifers in range(0, 64):
        idrs = A5202_map[ifers]
        ix, iy = map_ixy_DRS[idrs]
        map_ixy_FERS[ifers] = (ix, iy)

    return map_ixy_FERS


def build_map_FERSs_ixy():
    """
    Build a map for ixy and FERS channels for both boards.
    """
    map_ixy_FERS1 = build_map_FERS1_ixy()
    map_ixy_FERS2 = {}
    map_ixy_FERS3 = {}
    map_ixy_FERS4 = {}
    map_ixy_FERS5 = {}
    for ifers, (ix, iy) in map_ixy_FERS1.items():
        map_ixy_FERS2[ifers] = (ix-4, iy)
        map_ixy_FERS3[ifers] = (ix-8, iy-4)
        map_ixy_FERS4[ifers] = (ix-12, iy)
        map_ixy_FERS5[ifers] = (ix-16, iy)

    map_ixy_FERSs = {}
    map_ixy_FERSs["Board1"] = map_ixy_FERS1
    map_ixy_FERSs["Board2"] = map_ixy_FERS2
    map_ixy_FERSs["Board3"] = map_ixy_FERS3
    map_ixy_FERSs["Board4"] = map_ixy_FERS4
    map_ixy_FERSs["Board5"] = map_ixy_FERS5
    return map_ixy_FERSs


def build_map_DRSVar():
    """
    Build a map for DRS channel variables.
    """
    map_DRSVar_Cer = {}
    map_DRSVar_Sci = {}
    for i in range(0, 64):
        if i < 32:
            iboard = 0
        else:
            iboard = 2
        igroup = i // 8
        if igroup >= 4:
            igroup -= 4
        ichan = i % 8
        isCer = (i % 8) < 4  # first 4 channels are CER
        varname = f"DRS_Board{iboard}_Group{igroup}_Channel{ichan}"
        if isCer:
            map_DRSVar_Cer[i] = varname
        else:
            map_DRSVar_Sci[i] = varname
    return map_DRSVar_Cer, map_DRSVar_Sci


def build_map_ixy_DRSVar():
    """
    Build a map for DRS channel variables and ixy
    """
    map_ixy_DRS = build_map_ixy_DRS()
    map_DRSVar_Cer, map_DRSVar_Sci = build_map_DRSVar()
    map_ixy_DRSVar_Cer = {}
    map_ixy_DRSVar_Sci = {}
    for i in range(0, 64):
        ix, iy = map_ixy_DRS[i]
        if i in map_DRSVar_Cer:
            varname = map_DRSVar_Cer[i]
            map_ixy_DRSVar_Cer[(ix, iy)] = varname
        if i in map_DRSVar_Sci:
            varname = map_DRSVar_Sci[i]
            map_ixy_DRSVar_Sci[(ix, iy)] = varname
    return map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci


def get_hodoscope_channels():
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
    # bottom z seems to have flipped left and right
    hodoscope_channels["BottomZ"] = [
        "DRS_Board1_Group1_Channel0",
        "DRS_Board1_Group0_Channel7",
    ]

    return hodoscope_channels


def hodoTS2iX(ts):
    """
    Convert hodoscope timestamp difference between right and left (right - left) 
    to iX coordinate.
    """
    # hodoscope slice
    ihodo = int((ts + 200.0) / 400 * 11.0) + 1
    ihodo = 12 - ihodo  # flip
    iXmin = (ihodo - 8) * 4 - 4
    iXmax = (ihodo - 8) * 4 - 1
    return iXmin, iXmax


if __name__ == "__main__":
    map_Cer_Sci = build_map_Cer_Sci()
    print("Map of CER to SCI channels in FERS1:")
    for cer, sci in map_Cer_Sci.items():
        print(f"CER {cer} -> SCI {sci}")

    map_FERS1_ixy = build_map_FERS1_ixy()
    print("\nMap of FERS1 channels to (ix, iy):")
    for ifers, (ix, iy) in map_FERS1_ixy.items():
        print(f"FERS1 channel {ifers} -> (ix={ix}, iy={iy})")

    map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()
    print("\nMap of DRS variable names to (ix, iy) (CER):")
    for (ix, iy), varname in map_ixy_DRSVar_Cer.items():
        print(f"(ix={ix}, iy={iy}) -> {varname}")
    print("\nMap of DRS variable names to (ix, iy) (SCI):")
    for (ix, iy), varname in map_ixy_DRSVar_Sci.items():
        print(f"(ix={ix}, iy={iy}) -> {varname}")

    print("\nHodoscope channels:")
    hodoscope_channels = get_hodoscope_channels()
    for group, channels in hodoscope_channels.items():
        print(f"{group}: {channels}")
