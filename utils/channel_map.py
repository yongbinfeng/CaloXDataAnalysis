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


def build_map_ixy_FERS1():
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
