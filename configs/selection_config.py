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
        "MCP_1":       (150, 250,  5, 45, 50,     "PeakValue"),
        "MCP_2":       (150, 250,  5, 45, 50,     "PeakValue"),
        "TailCatcher":  (0,   1000, 5, 45, 1e3,    "Sum"),
    }
    return cuts.get(service_drs, cut_default)


# Detector requirements per particle type.
# Requirement: True means 'Fired' (Pass), False means 'Vetoed' (Not Fired).
#
# Selections are split by detector era so each setup reads top-to-bottom
# without any patching logic. Pick the right table with get_particle_selection().

# tb2025 setup (run_number <= 1700): PSD and three Cherenkovs available.
PARTICLE_SELECTIONS_TB2025 = {
    # Muon: PSD veto + TTUMuon passing
    "muon": {
        "TTUMuonVeto": True,
        "PSD": False,
    },
    # Pion: PSD veto + TTUMuon veto + Cherenkovs passing
    "pion": {
        "TTUMuonVeto": False,
        "PSD": False,
        "Cer474": True, "Cer519": True, "Cer537": True,
    },
    # Electron: PSD pass + TTUMuon veto + Cherenkovs passing
    "electron": {
        "TTUMuonVeto": False,
        "PSD": True,
        "Cer474": True, "Cer519": True, "Cer537": True,
    },
    # Proton: PSD veto + TTUMuon veto + Cherenkovs veto
    "proton": {
        "TTUMuonVeto": False,
        "PSD": False,
        "Cer474": False, "Cer519": False, "Cer537": False,
    },
    # Timing studies: MCPs fire
    "mcp_clean": {
        "MCP_1": True, "MCP_2": True,
    },
}

# tb2026 setup (run_number > 1700): no PSD, no Cer537. Only Cer474/Cer519 left.
PARTICLE_SELECTIONS_TB2026 = {
    # Muon: TTUMuon passing
    "muon": {
        "TTUMuonVeto": True,
    },
    # Pion: TTUMuon veto + Cer519 fired but Cer474 NOT fired
    "pion": {
        "TTUMuonVeto": False,
        "Cer474": False, "Cer519": True,
    },
    # Electron: TTUMuon veto + both Cherenkovs fired
    "electron": {
        "TTUMuonVeto": False,
        "Cer474": True, "Cer519": True,
    },
    # Proton: TTUMuon veto + both Cherenkovs veto
    "proton": {
        "TTUMuonVeto": False,
        "Cer474": False, "Cer519": False,
    },
    # Timing studies: MCPs fire
    "mcp_clean": {
        "MCP_1": True, "MCP_2": True,
    },
}


def get_particle_selection(particle_type: str, run_number: int = None) -> dict:
    """
    Returns a dictionary of detector requirements for a given particle.
    Requirement: True means 'Fired' (Pass), False means 'Vetoed' (Not Fired).

    tb2026 runs (run_number > 1700) use PARTICLE_SELECTIONS_TB2026; all other
    runs use PARTICLE_SELECTIONS_TB2025.
    """
    selections = (PARTICLE_SELECTIONS_TB2026
                  if (run_number is not None and run_number > 1700)
                  else PARTICLE_SELECTIONS_TB2025)
    return dict(selections.get(particle_type.lower(), {}))
