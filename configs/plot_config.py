from channels.channel_map import _DRS_6MM_RUN


def getRangesForFERSEnergySums(pdsub=False, calib=False, clip=False, HE=False, run_number=1230, beam_energy=80) -> dict:
    if run_number >= 1350:
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
            "LG_cer": 2e4 * scaleLG,
            "HG_sci": 1e6 * scaleHG,
            "LG_sci": 2e5 * scaleLG,
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
            "LG_cer": 3e4 * scaleLG,
            "HG_sci": 8e5 * scaleHG,
            "LG_sci": 2e5 * scaleLG,
            "Mix_cer": 90,
            "Mix_sci": 90
        }
    }
    config = configs["raw"]
    if pdsub:
        config = configs["pbsub"]
    if HE:
        config["xmax_total"]["HG_cer"] = 1e5 * scaleHG
        config["xmax_total"]["LG_cer"] = 1e4 * scaleLG
        config["xmax_total"]["HG_sci"] = 5e5 * scaleHG
        config["xmax_total"]["LG_sci"] = 3e4 * scaleLG
        config["xmax_total"]["Mix_cer"] = 80
        config["xmax_total"]["Mix_sci"] = 80
    # temporary fix; these can/should probably be improved later on...
    config["xmax_total"]["HG_cer"] *= beam_energy / 50
    config["xmax_total"]["LG_cer"] *= beam_energy / 50
    config["xmax_total"]["HG_sci"] *= beam_energy / 50
    config["xmax_total"]["LG_sci"] *= beam_energy / 50
    config["xmax_total"]["Mix_cer"] *= beam_energy / 50
    config["xmax_total"]["Mix_sci"] *= beam_energy / 50
    return config


def get_ttu_hodo_ranges():
    xmin = -0.032
    xmax = 4.064
    nbins = 64
    return xmin, xmax, nbins


def get_fers_saturation_value():
    # This is the saturation value for the FERS channels
    # return 8191
    # return 7800
    return 7500


def get_drs_plot_ranges(subtractMedian=False, is_amplified=False, is6mm=False, is_reference=False, is_mcp=False, run_number=None):
    xmin = -50
    xmax = 50
    if subtractMedian:
        xmin = -20
        xmax = 40
    # tb2026 calo layout (run_number > 1828): wider amplified range, and a
    # dedicated range for 6mm channels.
    tb2026_calo = run_number is not None and run_number > 1828
    if is_amplified:
        if tb2026_calo:
            xmin, xmax = -500, 3500
        else:
            xmin, xmax = -1500, 2500
    if is6mm and is_amplified:
        if tb2026_calo:
            return -100, 1000
        if run_number is not None and run_number >= _DRS_6MM_RUN:
            return -20, 40
        xmin = -100
        xmax = 300
    if is_reference:
        xmin = -500
        xmax = 2500
    if is_mcp:
        xmin = -500
        xmax = 2000
    return xmin, xmax


def get_drs_energy_range(is6mm=False):
    """DRS CFD energy integral histogram range: (nbins, xmin, xmax)."""
    if is6mm:
        return (180, -100, 800)
    return (200, -500, 6000)


def get_drs_peak_value_range():
    """DRS CFD peak value histogram range: (nbins, xmin, xmax)."""
    return (200, 0, 800)


def get_drs_noise_range():
    """DRS per-event baseline noise (RMS) histogram range: (nbins, xmin, xmax)."""
    return (200, 0, 30)


def get_drs_sum_vs_fers_ranges():
    """DRS energy sum vs FERS 2D ranges by FERS gain.
    Values are {var: (DRS_ymax, FERS_xmax)}.
    """
    return {
        "HG":  {"Cer": (20000, 8500),  "Sci": (30000, 8500)},
        "LG":  {"Cer": (20000, 2000),  "Sci": (30000, 4000)},
        "Mix": {"Cer": (40000, 40000), "Sci": (60000, 40000)},
    }


def get_fers_vs_drs_xmax():
    """FERS energy axis xmax for DRS-vs-FERS 2D correlation plots."""
    return {"FERS": 9000, "FERSLG": 3000, "FERSMix": 40000}


def get_fers_1d_range():
    """FERS per-channel 1D histogram range: {gain: (nbins, xmin, xmax)}."""
    return {"HG": (800, 0, 800), "LG": (300, 0, 300)}


def get_fers_max_range():
    """FERS max readout value histogram range: (nbins, xmin, xmax)."""
    return (100, 7500, 8500)


def get_fers_2d_hg_lg_range():
    """FERS 2D HG-vs-LG correlation range: (nbins_x, 0, hg_max, nbins_y, 0, lg_max)."""
    return (300, 0, 9000, 300, 0, 3000)


def get_fers_2d_mix_lg_range():
    """FERS 2D Mix-vs-LG correlation range: (nbins_x, 0, mix_max, nbins_y, 0, lg_max)."""
    return (300, 0, 40000, 300, 0, 2000)


def get_drs_cfd_finebins_range(is_cer):
    return (410, 450) if is_cer else (420, 490)


def get_drs_time_ns_range():
    """Broad range [ns] for scalar time_mcp / time_cfd_mcp timing histograms."""
    return (0, 20)


def get_drs_time_ns_finebins_range(is_cer):
    """Fine-binned range [ns] for time_cfd_mcp_finebins histograms."""
    return (6.5, 12.5) if is_cer else (7.5, 16.5)


def get_drs_time_arr_ns_range():
    """Waveform window [ns] for time_mcp 2D hist / profile."""
    return (0, 100)


def get_drs_prof_plot_ranges(subtractMedian=False, is_amplified=False, is6mm=False, is_reference=False, is_cer=False, run_number=None):
    xmin, xmax = get_drs_plot_ranges(
        subtractMedian, is_amplified, is6mm, is_reference)
    if not is_reference:
        xmax = xmax / 15.0
        xmin = xmin / 25.0
    if is_cer:
        xmax = xmax / 2.0
        xmin = xmin / 2.0
    # tb2026 (run_number > 1828): larger pulses -> double the upper range.
    if run_number is not None and run_number > 1828:
        xmax = xmax * 2.0
    return xmin, xmax


def getBoardEnergyFitParameters(run_number, is3mm=False, isCer=False):
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

    from utils.utils import is_scan_run
    isscanrun = is_scan_run(run_number)
    runtype = "scanruns" if isscanrun else "cosmicruns"

    boardtype = "3mm" if is3mm else "6mm"
    channeltype = "Cer" if isCer else "Sci"

    keystring = boardtype + "_" + channeltype
    if keystring not in args[runtype]:
        raise ValueError(
            f"Run {run_number} with board type {boardtype} and channel type {channeltype} not found in energy fit parameters.")

    results = args[runtype][keystring].copy()
    # Update the results with the run number
    results["run_number"] = run_number

    return results


def getEventEnergyFitParameters(run_number, isCer=False, clip=False):
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

    from utils.utils import is_scan_run
    isscanrun = is_scan_run(run_number)
    runtype = "scanruns" if isscanrun else "cosmicruns"
    if run_number >= 1150:
        runtype = "positron"
    if clip:
        runtype += "_clipped"

    channeltype = "Cer" if isCer else "Sci"

    if channeltype not in args[runtype]:
        raise ValueError(
            f"Run {run_number} with channel type {channeltype} not found in energy fit parameters.")

    results = args[runtype][channeltype].copy()
    results["run_number"] = run_number

    # Return a copy to avoid modifying the original dictionary
    return results


def get_service_drs_plot_ranges(channel, subtractMedian=True):

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


def get_service_drs_processed_info_ranges(channel, cat):
    service_drs_ranges = {
        "HoleVeto": {
            "peak_value": (-50, 1500),
            "sum": (-2e3, 1e4)
        },
        "TTUMuonVeto": {
            "peak_value": (-50, 1000),
            "sum": (-2e3, 1e4)
        },
        "PSD": {
            # "peak_value": (-50, 3000),
            # "sum": (-1e4, 7e4)
            "peak_value": (-50, 1000),
            "sum": (-500, 2e4)
        },
        "Cer474": {
            "peak_value": (-50, 1500),
            "sum": (-2e3, 3.5e4)
        },
        "Cer519": {
            "peak_value": (-50, 1500),
            "sum": (-2e3, 3e4)
        },
        "Cer537": {
            "peak_value": (-50, 1500),
            "sum": (-2e3, 2.5e4)
        },
        "KT1": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "KT2": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "T3": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "T4": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "ST1": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "ST3": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4)
        },
        "MCP": {
            "peak_value": (0, 400),
            "sum": (-500, 1000),
            "waveform": (-500, 1000)
        },
        "TailCatcher": {
            "peak_value": (-50, 1200),
            "sum": (-2e3, 1.2e4),
            "waveform": (-200, 3500)
        },
    }
    waveform_default = (-150, 3000)
    if channel.startswith("MCP"):
        return service_drs_ranges["MCP"].get(cat, (0, 400))
    if channel in service_drs_ranges:
        if cat in service_drs_ranges[channel]:
            return service_drs_ranges[channel][cat]
        if cat == "waveform":
            return waveform_default
        print(f"Category {cat} not found for channel {channel}.")
    if cat == "waveform":
        return waveform_default
    return -100, 100
