def get_service_drs_cut_values(service_drs: str) -> float:
    values = {
        "HoleVeto": -1000.0,
        "PSD": -9000.0,
        "TTUMuonVeto": -100.0,
        "Cer474": -1000.0,
        "Cer519": -1000.0,
        "Cer537": -1000.0,
    }
    return values.get(service_drs, -1000.0)
