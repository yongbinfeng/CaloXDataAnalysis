def getRangesForFERSEnergySums(subtractPedestal=False, calibrate=False, clip=False):
    suffix = ""
    xmin_board = 0
    xmin_total = 0
    xmax_board = 15000
    xmax_total = 200000
    title = "Energy (Raw)"
    if subtractPedestal:
        suffix = "_subtracted"
        xmin_board = -2000
        xmin_total = -5000
        xmax_board = 10000
        xmax_total = 5e4
        title = "Energy (PD subtracted)"
    if calibrate:
        suffix += "_calibrated"
        xmin_board = -50
        xmin_total = -100
        xmax_board = 200
        xmax_total = 1e3
        title = "# p.e."
    if clip:
        suffix += "_clipped"
        xmin_board = -50
        xmin_total = -100
        xmax_board = 200
        xmax_total = 1e3
        title = "# p.e. (clipped)"
    return suffix, xmin_board, xmax_board, xmin_total, xmax_total, title


def getBoardEnergyFitParameters(runNumber, is3mm=False, isCer=False):
    args = {}
    args["scanruns"] = {}
    args["cosmicruns"] = {}

    args["scanruns"]["3mm_Cer"] = {
        "addMIP": False,
        "addHE": False,
        "xmin": -10, "xmax": 30,
        "xfitmin": -10, "xfitmax": 30,
        "xgausmean": 5, "xgausmin": 0, "xgausmax": 10,
        "wgausmean": 5, "wgausmin": 1, "wgausmax": 10,
    }
    # for the scan runs, no difference between Cer and Sci
    args["scanruns"]["3mm_Sci"] = args["scanruns"]["3mm_Cer"].copy()
    args["scanruns"]["6mm_Cer"] = {
        "addMIP": False,
        "addHE": False,
        "xmin": -30, "xmax": 60,
        "xfitmin": -30, "xfitmax": 60,
        "xgausmean": 10, "xgausmin": 0, "xgausmax": 20,
        "wgausmean": 10, "wgausmin": 5, "wgausmax": 50,
    }
    args["scanruns"]["6mm_Sci"] = args["scanruns"]["6mm_Cer"].copy()

    args["cosmicruns"]["3mm_Cer"] = {
        "addMIP": False,
        "addHE": False,
        "xmin": -10, "xmax": 50,
        "xfitmin": -10, "xfitmax": 6,
        "xgausmean": 5, "xgausmin": 2, "xgausmax": 10,
        "wgausmean": 5, "wgausmin": 1, "wgausmax": 10,
    }
    args["cosmicruns"]["3mm_Sci"] = {
        "addMIP": True,
        "addHE": False,
        "xmin": -10, "xmax": 100,
        "xfitmin": -10, "xfitmax": 60,
        "xgausmean": 5, "xgausmin": 1, "xgausmax": 12,
        "wgausmean": 4, "wgausmin": 2, "wgausmax": 10,
        "xmipmean": 30, "xmipmin": 12, "xmipmax": 35,
        "wmipmean": 6, "wmipmin": 1, "wmipmax": 30,
    }
    args["cosmicruns"]["6mm_Cer"] = {
        "addMIP": False,
        "addHE": False,
        "xmin": -30, "xmax": 100,
        "xfitmin": -30, "xfitmax": 20,
        "xgausmean": 10, "xgausmin": 1, "xgausmax": 20,
        "wgausmean": 10, "wgausmin": 5, "wgausmax": 50,
    }
    args["cosmicruns"]["6mm_Sci"] = {
        "addMIP": True,
        "addHE": False,
        "xmin": -30, "xmax": 200,
        "xfitmin": -30, "xfitmax": 160,
        "xgausmean": 10, "xgausmin": 1, "xgausmax": 20,
        "wgausmean": 20, "wgausmin": 10, "wgausmax": 50,
        "xmipmean": 100, "xmipmin": 20, "xmipmax": 140,
        "wmipmean": 60, "wmipmin": 10, "wmipmax": 80,
    }

    from utils.utils import IsScanRun
    isscanrun = IsScanRun(runNumber)
    runtype = "scanruns" if isscanrun else "cosmicruns"

    boardtype = "3mm" if is3mm else "6mm"
    channeltype = "Cer" if isCer else "Sci"

    keystring = boardtype + "_" + channeltype
    if keystring not in args[runtype]:
        raise ValueError(
            f"Run {runNumber} with board type {boardtype} and channel type {channeltype} not found in energy fit parameters.")

    results = args[runtype][keystring].copy()
    # Update the results with the run number
    results["runNumber"] = runNumber

    return results


def getEventEnergyFitParameters(runNumber, isCer=False, clip=False):
    args = {}
    args["scanruns"] = {}
    args["cosmicruns"] = {}
    args["scanruns_clipped"] = {}
    args["cosmicruns_clipped"] = {}

    args["scanruns"]["Cer"] = {
        "addMIP": False,
        "addHE": False,
        "xmin": -50, "xmax": 500,
        "xfitmin": -50, "xfitmax": 450,
        "xgausmean": 110, "xgausmin": 50, "xgausmax": 250,
        "wgausmean": 40, "wgausmin": 20, "wgausmax": 70,
    }
    # for the scan runs, no difference between Cer and Sci
    args["scanruns"]["Sci"] = args["scanruns"]["Cer"].copy()

    args["scanruns_clipped"]["Cer"] = args["scanruns"]["Cer"].copy()
    args["scanruns_clipped"]["Sci"] = args["scanruns"]["Sci"].copy()

    args["cosmicruns"]["Cer"] = {
        "addMIP": True,
        "addHE": False,
        "xmin": -50, "xmax": 1000,
        "xfitmin": -50, "xfitmax": 250,
        "xgausmean": 105, "xgausmin": 50, "xgausmax": 120,
        "wgausmean": 40, "wgausmin": 20, "wgausmax": 50,
        "xmipmean": 170, "xmipmin": 120, "xmipmax": 200,
        "wmipmean": 60, "wmipmin": 30, "wmipmax": 80,
    }

    args["cosmicruns"]["Sci"] = {
        "addMIP": True,
        "addHE": False,
        "xmin": -50, "xmax": 1000,
        "xfitmin": -50, "xfitmax": 450,
        "xgausmean": 110, "xgausmin": 50, "xgausmax": 150,
        "wgausmean": 40, "wgausmin": 20, "wgausmax": 60,
        "xmipmean": 300, "xmipmin": 200, "xmipmax": 350,
        "wmipmean": 60, "wmipmin": 30, "wmipmax": 80,
    }

    args["cosmicruns_clipped"]["Cer"] = args["cosmicruns"]["Cer"].copy()
    args["cosmicruns_clipped"]["Cer"].update({
        "xgausmean": 150, "xgausmin": 100, "xgausmax": 180,
        "wgausmean": 20, "wgausmin": 10, "wgausmax": 30,
        "xmipmean": 200, "xmipmin": 150, "xmipmax": 250,
        "wmipmean": 40, "wmipmin": 20, "wmipmax": 80,
    })
    args["cosmicruns_clipped"]["Sci"] = args["cosmicruns"]["Sci"].copy()
    args["cosmicruns_clipped"]["Sci"].update({
        "xgausmean": 140, "xgausmin": 100, "xgausmax": 180,
        "wgausmean": 20, "wgausmin": 10, "wgausmax": 30,
        "xmipmean": 400, "xmipmin": 300, "xmipmax": 450,
        "wmipmean": 40, "wmipmin": 20, "wmipmax": 80,
    })

    from utils.utils import IsScanRun
    isscanrun = IsScanRun(runNumber)
    runtype = "scanruns" if isscanrun else "cosmicruns"
    if clip:
        runtype += "_clipped"

    channeltype = "Cer" if isCer else "Sci"

    if channeltype not in args[runtype]:
        raise ValueError(
            f"Run {runNumber} with channel type {channeltype} not found in energy fit parameters.")

    results = args[runtype][channeltype].copy()
    results["runNumber"] = runNumber

    # Return a copy to avoid modifying the original dictionary
    return results
