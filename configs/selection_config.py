def get_service_drs_cut(service_drs: str, run_number: int = None) -> tuple:
    # format:
    # (ts_begin, ts_end, window_pre, window_post, cut_value, method)
    # where:
    # (ts_begin, ts_end): window to look for pulse/peak
    # (window_pre, window_post): integration window with respect to the pulse/peak time
    # method: "Sum" or "PeakValue"
    cut_default = (0, 1000, 5, 45, 5e4, "Sum")
    # Cer474 pulse window moved earlier for tb2026 (run_number >= 1700).
    cer474 = (800, 1000, 20, 150, 2000.0, "Sum") if (run_number is None or run_number < 1700) \
        else (400, 800, 20, 150, 2000.0, "Sum")
    cuts = {
        "HoleVeto":    (100, 350,  20, 150, 2e3,    "Sum"),
        "PSD":         (100, 400,  20, 150, 3500.0,  "Sum"),
        "TTUMuonVeto": (200, 400,  20, 150, 5e3,    "Sum"),
        "Cer474":      cer474,
        "Cer519":      (400, 800,  5, 100, 1000.0,  "Sum"),
        "Cer537":      (400, 800,  5, 80, 500.0,   "Sum"),
        "KT1":         (0,   1000, 5, 45, 1e3,    "Sum"),
        "KT2":         (0,   1000, 5, 45, 1e3,    "Sum"),
        "T3":          (0,   1000, 5, 45, 1e3,    "Sum"),
        "T4":          (0,   1000, 5, 45, 1e3,    "Sum"),
        "MCP_US_0":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_US_1":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_US_2":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_US_3":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_DS_0":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_DS_1":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_DS_2":    (500, 650,  5, 45, 20,     "PeakValue"),
        "MCP_DS_3":    (500, 650,  5, 45, 20,     "PeakValue"),
    }
    return cuts.get(service_drs, cut_default)


def get_particle_selection(particle_type: str) -> dict:
    """
    Returns a dictionary of detector requirements for a given particle.
    Requirement: True means 'Fired' (Pass), False means 'Vetoed' (Not Fired).
    """
    selections = {
        # Muon: PSD veto + TTUMuon passing
        "muon": {
            "TTUMuonVeto": True,
            "PSD": False,
        },
        # Pions: PSD veto + TTUMuon veto + Cherenkov passing
        "pion": {
            "TTUMuonVeto": False,
            "PSD": False,
            "Cer474": True, "Cer519": True, "Cer537": True
        },
        # Electrons: PSD pass + TTUMuon veto + Cherenkov passing
        "electron": {
            "TTUMuonVeto": False,
            "PSD": True,
            "Cer474": True, "Cer519": True, "Cer537": True
        },
        # Protons: PSD veto + TTUMuon veto + Cherenkov veto
        "proton": {
            "TTUMuonVeto": False,
            "PSD": False,
            "Cer474": False, "Cer519": False, "Cer537": False
        },
    }
    return selections.get(particle_type.lower(), {})
