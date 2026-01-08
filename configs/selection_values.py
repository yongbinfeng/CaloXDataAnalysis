def get_service_drs_cut(service_drs: str) -> tuple:
    # Returns (ts_min, ts_max, value_cut, method)
    cut_default = (0, 1000, -5e4, "Sum")
    cuts = {
        "HoleVeto": (100, 350, -1000.0, "Min"),
        "PSD": (100, 400, -9000.0, "Sum"),
        "TTUMuonVeto": (100, 400, -2e3, "Sum"),
        "Cer474": (800, 950, -3000.0, "Sum"),
        "Cer519": (400, 600, -3000.0, "Sum"),
        "Cer537": (400, 600, -3000.0, "Sum"),
    }
    return cuts.get(service_drs, cut_default)
