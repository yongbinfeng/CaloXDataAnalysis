def get_service_drs_cut(service_drs: str) -> tuple:
    cut_default = (0, 1000, -5e4, "Sum")
    cuts = {
        "HoleVeto": (100, 350, -2e3, "Sum"),
        "PSD": (100, 400, -3500.0, "Sum"),
        "TTUMuonVeto": (200, 400, -2e3, "Sum"),
        "Cer474": (800, 900, -2000.0, "Sum"),
        "Cer519": (450, 550, -1000.0, "Sum"),
        "Cer537": (400, 500, -500.0, "Sum"),
        "KT1": (0, 1000, -1e3, "Sum"),
        "KT2": (0, 1000, -1e3, "Sum"),
        "T3": (0, 1000, -1e3, "Sum"),
        "T4": (0, 1000, -1e3, "Sum"),
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
