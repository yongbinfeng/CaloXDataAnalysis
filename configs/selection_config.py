# Per-detector service-DRS cuts, format:
#   (ts_begin, ts_end, window_pre, window_post, cut_value, method)
#   (ts_begin, ts_end): window to look for the pulse/peak
#   (window_pre, window_post): integration window wrt the pulse/peak time
#   method: "Sum" or "PeakValue"
_SERVICE_DRS_CUT_DEFAULT = (0, 1000, 5, 45, 5e4, "Sum")

# Base table — applies to all runs unless overridden by an era table below.
_SERVICE_DRS_CUTS = {
    "HoleVeto":    (100, 350,  20, 150, 2e3,    "Sum"),
    "PSD":         (100, 400,  20, 150, 3500.0, "Sum"),
    "TTUMuonVeto": (200, 400,  20, 150, 5e3,    "Sum"),
    "Cer474":      (800, 1000, 20, 150, 2000.0, "Sum"),
    "Cer519":      (400, 800,  5,  100, 1000.0, "Sum"),
    "Cer537":      (400, 800,  5,  80,  500.0,  "Sum"),
    "KT1":         (0,   1000, 5,  45,  1e3,    "Sum"),
    "KT2":         (0,   1000, 5,  45,  1e3,    "Sum"),
    "T3":          (0,   1000, 5,  45,  1e3,    "Sum"),
    "T4":          (0,   1000, 5,  45,  1e3,    "Sum"),
    "MCP_US_0":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_US_1":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_US_2":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_US_3":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_DS_0":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_DS_1":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_DS_2":    (500, 650,  5,  45,  20,     "PeakValue"),
    "MCP_DS_3":    (500, 650,  5,  45,  20,     "PeakValue"),
}

# tb2026 (run_number >= 1700): Cer474 window moved earlier; service detectors
# (single MCPs, TailCatcher, scintillator telescopes) only present this era.
_SERVICE_DRS_CUTS_TB2026 = {
    "Cer474":      (400, 800,  20, 150, 2000.0, "Sum"),
    "MCP_1":       (150, 250,  5,  45,  50,     "PeakValue"),
    "MCP_2":       (150, 250,  5,  45,  50,     "PeakValue"),
    "TailCatcher": (0,   1000, 5,  45,  1e3,    "Sum"),
    "ST1":         (0,   1000, 5,  45,  5e3,    "Sum"),
    "ST3":         (0,   1000, 5,  45,  5e3,    "Sum"),
}

# tb2026 calo-MCP layout (run_number > 1828): earlier, lower-threshold MCP window.
_SERVICE_DRS_CUTS_CALO_MCP = {
    "MCP_US_0":    (400, 550,  5,  45,  50,     "PeakValue"),
    "MCP_US_1":    (400, 550,  5,  45,  200,     "PeakValue"),
    "MCP_DS_0":    (400, 550,  5,  45,  50,     "PeakValue"),
    "MCP_DS_1":    (400, 550,  5,  45,  200,     "PeakValue"),
}


def get_service_drs_cut(service_drs: str, run_number: int = None) -> tuple:
    """Return the (ts_begin, ts_end, window_pre, window_post, cut, method) tuple
    for a service-DRS detector, layering era overrides on the base table."""
    cuts = dict(_SERVICE_DRS_CUTS)
    if run_number is not None and run_number >= 1700:
        cuts.update(_SERVICE_DRS_CUTS_TB2026)
    if run_number is not None and run_number > 1828:
        cuts.update(_SERVICE_DRS_CUTS_CALO_MCP)
    return cuts.get(service_drs, _SERVICE_DRS_CUT_DEFAULT)


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

# tb2026 calo-MCP layout (run_number > 1828): per-particle overrides layered on
# top of the era table (e.g. mcp_clean uses the US/DS MCPs instead of MCP_1/2).
PARTICLE_SELECTIONS_CALO_MCP = {
    "mcp_clean": {
        "MCP_US_0": True, "MCP_DS_0": True,
    },
}


def get_particle_selection(particle_type: str, run_number: int = None) -> dict:
    """
    Returns a dictionary of detector requirements for a given particle.
    Requirement: True means 'Fired' (Pass), False means 'Vetoed' (Not Fired).

    tb2026 runs (run_number > 1700) use PARTICLE_SELECTIONS_TB2026; all other
    runs use PARTICLE_SELECTIONS_TB2025. For run_number > 1828, per-particle
    overrides from PARTICLE_SELECTIONS_CALO_MCP replace the era entry.
    """
    particle_type = particle_type.lower()
    selections = (PARTICLE_SELECTIONS_TB2026
                  if (run_number is not None and run_number > 1700)
                  else PARTICLE_SELECTIONS_TB2025)
    if run_number is not None and run_number > 1828 and particle_type in PARTICLE_SELECTIONS_CALO_MCP:
        selections = PARTICLE_SELECTIONS_CALO_MCP
    return dict(selections.get(particle_type, {}))
