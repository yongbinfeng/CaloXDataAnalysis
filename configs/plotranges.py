def getRangesForFERSEnergySums(pdsub=False, calib=False, clip=False, HE=False, runNumber=1230) -> dict:
    scale = 1.0
    if runNumber >= 1350:
        scaleHG = 0.6
        scaleLG = 0.2
    configs = {}
    configs["raw"] = {
        "suffix": "raw",
        "title_LG": "ADC (Raw)",
        "title_HG": "ADC (Raw)",
        "title_Mix": "GeV",
        "xmin_board": {
            "HG_cer": 0,
            "LG_cer": 0,
            "HG_sci": 0,
            "LG_sci": 0,
            "Mix_cer": 0,
            "Mix_sci": 0,
        },
        "xmax_board": {
            "HG_cer": 1e5 * scaleHG,
            "LG_cer": 2e4 * scaleLG,
            "HG_sci": 3e5 * scaleHG,
            "LG_sci": 1e5 * scaleLG,
            "Mix_cer": 80,
            "Mix_sci": 80,
        },
        "xmin_total": {
            "HG_cer": 0,
            "LG_cer": 0,
            "HG_sci": 0,
            "LG_sci": 0,
            "Mix_cer": 0,
            "Mix_sci": 0,
        },
        "xmax_total": {
            "HG_cer": 2.5e5 * scaleHG,
            "LG_cer": 9e4 * scaleLG,
            "HG_sci": 1e6 * scaleHG,
            "LG_sci": 5e5 * scaleLG,
            "Mix_cer": 90,
            "Mix_sci": 90
        },
    }
    configs["pbsub"] = {
        "suffix": "pbsub",
        "title_LG": "ADC (PD sub)",
        "title_HG": "ADC (PD sub)",
        "title_Mix": "GeV",
        "xmin_board": {
            "HG_cer": -1e3,
            "LG_cer": -1e3,
            "HG_sci": -1e3,
            "LG_sci": -1e3,
            "Mix_cer": -10,
            "Mix_sci": -10,
        },
        "xmax_board": {
            "HG_cer": 9e4 * scaleHG,
            "LG_cer": 3e4 * scaleLG,
            "HG_sci": 4e5 * scaleHG,
            "LG_sci": 2e5 * scaleLG,
            "Mix_cer": 40,
            "Mix_sci": 40,
        },
        "xmin_total": {
            "HG_cer": -5e3,
            "LG_cer": -5e3,
            "HG_sci": -5e3,
            "LG_sci": -5e3,
            "Mix_cer": -10,
            "Mix_sci": -10,
        },
        "xmax_total": {
            "HG_cer": 2e5 * scaleHG,
            "LG_cer": 8e4 * scaleLG,
            "HG_sci": 8e5 * scaleHG,
            "LG_sci": 4e5 * scaleLG,
            "Mix_cer": 90,
            "Mix_sci": 90
        }
    }
    config = configs["raw"]
    if pdsub:
        config = configs["pbsub"]
    if HE:
        config["xmax_total"]["HG_cer"] = 2e5 * scaleHG
        config["xmax_total"]["LG_cer"] = 5e4 * scaleLG
        config["xmax_total"]["HG_sci"] = 1e6 * scaleHG
        config["xmax_total"]["LG_sci"] = 1e5 * scaleLG
        config["xmax_total"]["Mix_cer"] = 150
        config["xmax_total"]["Mix_sci"] = 150
    return config


def getFERSSaturationValue():
    # This is the saturation value for the FERS channels
    # return 8191
    # return 7800
    return 7500


def getDRSPlotRanges(subtractMedian=False, isAmplified=False, is6mm=False):
    xmin = -50
    xmax = 50
    if subtractMedian:
        xmin = -20
        xmax = 40
    if isAmplified:
        xmin = -1000
        xmax = 2500
    if is6mm and isAmplified:
        xmin = -200
        xmax = 500
    return xmin, xmax


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

    args["positron"] = {}
    args["positron"]["Cer"] = {
        "addMIP": False, "addHE": False,
        "xmin": -50, "xmax": 2000,
        "xfitmin": -50, "xfitmax": 1500,
        "xgausmean": 1000, "xgausmin": 100, "xgausmax": 1500,
        "wgausmean": 50, "wgausmin": 10, "wgausmax": 700,
    }
    args["positron"]["Sci"] = {
        "addMIP": False, "addHE": False,
        "xmin": -50, "xmax": 2e4,
        "xfitmin": 2000, "xfitmax": None,
        "xgausmean": 5e3, "xgausmin": 1e3, "xgausmax": 1e4,
        "wgausmean": 200, "wgausmin": 10, "wgausmax": 550
    }

    from utils.utils import IsScanRun
    isscanrun = IsScanRun(runNumber)
    runtype = "scanruns" if isscanrun else "cosmicruns"
    if runNumber >= 1150:
        runtype = "positron"
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


def getServiceDRSPlotRanges(channel, subtractMedian=True):

    service_drs_ranges = {
        "DRS_Board7_Group0_Channel0": (2000, -3000),
        "DRS_Board7_Group0_Channel1": (2000, -3000),
        "DRS_Board7_Group0_Channel2": (2000, -3000),
        "DRS_Board7_Group0_Channel3": (2000, -3000),
        "DRS_Board7_Group0_Channel4": (1000, -3000),
        "DRS_Board7_Group0_Channel5": (1000, -3000),
        "DRS_Board7_Group0_Channel6": (1000, -3000),
        "DRS_Board7_Group0_Channel7": (1000, -3000),
        "DRS_Board7_Group1_Channel0": (200, -500),
        "DRS_Board7_Group1_Channel1": (200, -3000),
        "DRS_Board7_Group1_Channel2": (200, -500),
        "DRS_Board7_Group1_Channel3": (100, -200),
        "DRS_Board7_Group1_Channel4": (100, -200),
        "DRS_Board7_Group1_Channel5": (100, -200),
        "DRS_Board7_Group1_Channel6": (500, -3000),
        "DRS_Board7_Group1_Channel7": (100, -200),
        "DRS_Board7_Group2_Channel0": (500, -3000),
        "DRS_Board7_Group2_Channel1": (500, -3000),
        "DRS_Board7_Group2_Channel2": (500, -3000),
        "DRS_Board7_Group2_Channel3": (500, -3000),
        "DRS_Board7_Group2_Channel4": (500, -3000),
        "DRS_Board7_Group2_Channel5": (500, -3000),
        "DRS_Board7_Group2_Channel6": (500, -3000),
        "DRS_Board7_Group2_Channel7": (500, -3000),
    }
    xmin = -50
    xmax = 50
    if subtractMedian:
        if channel in service_drs_ranges:
            xmax, xmin = service_drs_ranges[channel]
    return xmin, xmax


def getServiceDRSProcessedInfoRanges(channel, cat):
    service_drs_ranges = {
        "muon": {
            "peak_value": (-500, 50),
            "sum": (-2e3, 500)
        },
        "preshower": {
            # "peak_value": (-3000, 50),
            # "sum": (-7e4, 1e4)
            "peak_value": (-1000, 50),
            "sum": (-1e4, 500)
        },
        "Cerenkov1": {
            "peak_value": (-1500, 50),
            "sum": (-1e4, 2e3)
        },
        "Cerenkov2": {
            "peak_value": (-1500, 50),
            "sum": (-1e4, 2e3)
        },
        "Cerenkov3": {
            "peak_value": (-1500, 50),
            "sum": (-1e4, 2e3)
        },
    }
    if channel in service_drs_ranges:
        if cat in service_drs_ranges[channel]:
            return service_drs_ranges[channel][cat]
        else:
            print(f"Category {cat} not found for channel {channel}.")
            return -100, 100
    return -100, 100
