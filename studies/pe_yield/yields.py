"""p.e. per GeV of beam energy for each fibre type: Mix ADC (pedestal-subtracted HG, LG scaled to HG where HG saturates; recovered from
the calibrated dump with FERS_response) divided by the gain of the channel's board (gain.json per_board; boards with too few lit
channels for a fit: the median of the fitted boards without resolved peaks), summed over all channels of a type, averaged over the
positron events. Also with one gain for all boards (the resolved-board value) for comparison.
Events: beam PID (positron), 1-3 hodoscope fibres per axis, calibrated (S+C)/E above SEL_MIN_FRACTION.
Usage (repository root): python3 studies/pe_yield/yields.py      Output: $PE_YIELD_DATA/yields.json"""
import json, sys
sys.path.insert(0, ".")
import numpy as np
from configs.selection_config import get_particle_selection
from config import DATA, RESPONSE, RUNS, BAD_RUNS, MAP_RUN, GAIN_BOARDS, SEL_MIN_FRACTION, FAR

TYPES = ("Sci", "CerQ", "CerP")
SCE_BINS = np.linspace(-0.2, 3.0, 161)        # (S+C)/E histogram


def load(run):
    d = np.load(f"{DATA}/fers/fers_{run}.npz"); R = json.load(open(RESPONSE))
    meta = d["meta"].astype(int); E = float(d["energy"])
    resp = np.array([R.get(n.replace("_calib", ""), {}).get("response", np.nan) for n in d["names"]])
    resp = np.where(resp > 0, resp, np.nan)       # missing / dead channels: left out (the calibrated value is 0)
    pid = np.ones(len(d["vals"]), bool)
    for det, fire in get_particle_selection("electron", run).items():
        pid &= d[f"is_{det}_fired" if fire else f"is_{det}_vetoed"].astype(bool)
    nx, ny = d["TTU_Hodo_nHitX"], d["TTU_Hodo_nHitY"]
    pid &= (nx >= 1) & (nx <= 3) & (ny >= 1) & (ny <= 3)
    sce = np.nansum(d["vals"], 1) / E
    return d, meta, E, d["vals"] / resp, pid, sce


def run_yield(run, gain, common):
    d, meta, E, adc, pid, sce = load(run)
    typ = np.where(meta[:, 2] == 0, "Sci", np.where(meta[:, 3] == 1, "CerQ", "CerP"))
    sel = pid & (sce > SEL_MIN_FRACTION)
    mu = np.nanmean(adc[sel], 0)
    w = np.nan_to_num(mu) * (typ == "Sci"); cx, cy = (w * meta[:, 4]).sum() / w.sum(), (w * meta[:, 5]).sum() / w.sum()
    r = np.hypot(meta[:, 4] - cx, meta[:, 5] - cy)
    out = dict(run=run, energy=E, n_pid=int(pid.sum()), n_sel=int(sel.sum()), centre=[float(cx), float(cy)],
               sce_median=float(np.median(sce[pid])), sce_hist=np.histogram(sce[pid], SCE_BINS)[0].tolist(), types={})
    for t in TYPES:
        k = typ == t; s = np.nansum(adc[sel][:, k] / gain[k], 1)
        out["types"][t] = dict(pe_per_gev=float(s.mean() / E), stat=float(s.std() / np.sqrt(len(s)) / E), rms_over_mean=float(s.std() / s.mean()),
                               pe_per_gev_common=float(np.nansum(mu[k]) / common / E),
                               share_ref=float(np.nansum(mu[k & np.isin(meta[:, 0], GAIN_BOARDS)]) / np.nansum(mu[k])),
                               far=float(np.nansum(mu[k & (r >= FAR)]) / np.nansum(mu[k])))
    out["types"]["Cer"] = {k: out["types"]["CerQ"][k] + out["types"]["CerP"][k] for k in ("pe_per_gev", "pe_per_gev_common")}
    if run == MAP_RUN:
        out["map"] = dict(meta=meta.tolist(), pe_per_event=(mu / gain).tolist())
    return out


def main():
    GJ = json.load(open(f"{DATA}/gain.json")); common = GJ["summary"]["gain"]; PB = {int(b): v for b, v in GJ["per_board"].items()}
    meta = np.load(f"{DATA}/fers/fers_{next(iter(RUNS))}.npz")["meta"].astype(int)
    other = float(np.median([v["gain"] for v in PB.values() if not v["resolved"]]))
    gain = np.array([PB[b]["gain"] if b in PB else other for b in meta[:, 0]])
    res = dict(gain=GJ["summary"], per_board={b: v["gain"] for b, v in PB.items()}, unfitted_boards=other,
               selection=dict(min_fraction=SEL_MIN_FRACTION), sce_bins=SCE_BINS.tolist(), runs=[], bad_runs=[])
    for run in RUNS:
        res["runs"].append(run_yield(run, gain, common)); y = res["runs"][-1]
        print(run, y["energy"], y["n_sel"], {t: round(v["pe_per_gev"], 2) for t, v in y["types"].items()})
    for run in BAD_RUNS:
        _, _, E, _, pid, sce = load(run)
        res["bad_runs"].append(dict(run=run, energy=E, n_pid=int(pid.sum()), sce_median=float(np.median(sce[pid])),
                                    sce_hist=np.histogram(sce[pid], SCE_BINS)[0].tolist()))
    json.dump(res, open(f"{DATA}/yields.json", "w"))


if __name__ == "__main__":
    main()
