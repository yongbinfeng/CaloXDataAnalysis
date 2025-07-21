def getRangesForFERSEnergySums(subtractPedestal=False, calibrate=False):
    suffix = ""
    xmin_board = 0
    xmin_total = 0
    xmax_board = 15000
    xmax_total = 200000
    if subtractPedestal:
        suffix = "_subtracted"
        xmin_board = -2000
        xmin_total = -5000
        xmax_board = 10000
        xmax_total = 5e4
    if calibrate:
        suffix += "_calibrated"
        xmin_board = -50
        xmin_total = -100
        xmax_board = 200
        xmax_total = 1e3
    return suffix, xmin_board, xmax_board, xmin_total, xmax_total
