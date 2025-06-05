# map of fers channel to drs channel
A5202_map: list[int] = [4, 2, 3, 1, 8, 6, 7, 5,
                        11, 9, 12, 10, 15, 13, 16, 14,
                        20, 18, 19, 17, 24, 22, 23, 21,
                        27, 25, 28, 26, 31, 29, 32, 30,
                        36, 34, 35, 33, 40, 38, 39, 37,
                        43, 41, 44, 42, 47, 45, 48, 46,
                        52, 50, 51, 49, 56, 54, 55, 53,
                        59, 57, 60, 58, 63, 61, 64, 62]


def build_map_Cer_Sci():
    """
    Build a map for CER and SCI channels.
    """
    map_Cer_Sci = {}
    for iCer in [1, 2, 3, 4,
                 9, 10, 11, 12,
                 17, 18, 19, 20,
                 25, 26, 27, 28,
                 33, 34, 35, 36,
                 41, 42, 43, 44,
                 49, 50, 51, 52,
                 57, 58, 59, 60]:
        iSci = iCer + 4

        index_Cer = A5202_map.index(iCer)
        index_Sci = A5202_map.index(iSci)
        map_Cer_Sci[index_Cer] = index_Sci
    return map_Cer_Sci
