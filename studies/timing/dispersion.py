"""Dispersion in the fibres from the timing resolution vs position in the z scans (studies/zscan datasets).
Per channel and run, the width of the in-time CFD time (relative to the MCP) in bins of the number of photoelectrons of the hit
(peak amplitude / single-p.e. amplitude, one calibration per channel from the out-of-time pulses of all runs, spe_of). If dispersion
spreads the arrival time of every photon by k d (d: light path to the readout), a hit of N p.e. has
    sigma^2 = sigma_c^2 + (sigma_1^2 + k^2 d^2) / N,
so at fixed N, sigma^2 is linear in X over the scan, with d(sigma^2)/dX = -2 k^2 d_mean / N, and N x d(sigma^2)/dX (the per-photon
slope) is the same in every p.e. bin. Fit per channel and bin: sigma^2 = a_stretch + s X; errors scaled by sqrt(chi2/ndf) where above 1
(run-to-run scatter). Per channel: the per-photon slope, weighted over the bins, with the chi2 of its constancy (the 1/N test).
Control: the width of MCP_US_1 - MCP_DS_1 (both fixed upstream of the table, on the same DRS chip) vs X: it must not change.
Hits: CFD time within W_T of the arrival peak (wider than the z-scan in-time window, which cuts the 1-p.e. distributions of the
dispersing fibres at ~2 sigma far from the readout); random hits of the same p.e. bin from the sidebands, subtracted in the width.
Width: gaussian core (iterative 2-sigma clipping with the flat random background removed, corrected for the clipping).
Usage (repository root): python3 studies/timing/dispersion.py DATASET     -> $TIMING_DATA/dispersion_<dataset>.json"""
import json, os, sys
import numpy as np
from config import DATA as OUT, ZSCAN
sys.path.insert(0, ZSCAN)
from common import Dataset, DATA, TS_NS, A_HIT, SIDEBANDS, load_drs, intime, spe_of, pid, mcp_ok
from datasets import drs_dump
from channels.maps.services import get_mcp_channels

PE_BINS = ((0.5, 1.5), (1.5, 2.5), (2.5, 4.5), (4.5, 8.5), (8.5, 1e9))
MIN_HITS = 30                 # hits in a run and p.e. bin for a width
MIN_RUNS = 4                  # runs with a width for a slope (2025 has five)
MODE = {"testfibres": "loose"}  # event selection per dataset (default "spot", the z-scan selection)
MCP_CTRL = ("MCP_US_1", "MCP_DS_1")
W_T = {"Sci": 40.0, "Cer": 15.0}  # [TS] hit window around the arrival peak: 8 ns scintillating (decay tail), 3 ns Cherenkov
SPAN_NS = sum(hi - lo for lo, hi in SIDEBANDS) * TS_NS


def core_sigma(t, bkg=0.0, nit=8):
    """gaussian core width and its error: iterative 2-sigma clipped rms, the flat background (bkg hits per ns, uniform within the
    clipping window: variance (4 s)^2 / 12) removed, x 1.137 for the clipping"""
    if len(t) < MIN_HITS: return np.nan, np.nan
    m, s = np.median(t), 1.4826 * np.median(np.abs(t - np.median(t)))
    for _ in range(nit):
        k = np.abs(t - m) < 2 * s; nb = bkg * 4 * s; ns = k.sum() - nb
        if ns < MIN_HITS or s <= 0: return np.nan, np.nan
        m = t[k].mean(); v = (np.sum((t[k] - m) ** 2) - nb * 4 * s ** 2 / 3) / ns
        if v <= 0: return np.nan, np.nan
        s = np.sqrt(v) * 1.137
    return float(s), float(s / np.sqrt(2 * ns))


def line_fit(x, gi, y, e):
    """y = a_group + s x (weighted); returns s, its error scaled by sqrt(chi2/ndf) where above 1, chi2/ndf, intercepts"""
    G = sorted(set(gi)); A = np.column_stack([(gi == g).astype(float) for g in G] + [x]); w = 1 / e ** 2
    C = np.linalg.inv(A.T @ (A * w[:, None])); c = C @ (A.T @ (w * y)); ndf = len(y) - A.shape[1]
    chi2 = float(np.sum(w * (y - A @ c) ** 2) / ndf) if ndf > 0 else np.nan
    return float(c[-1]), float(np.sqrt(C[-1, -1] * max(chi2, 1.0))), chi2, {str(g): float(v) for g, v in zip(G, c[:-1])}


def wmean(v, e):
    v, e = np.asarray(v), np.asarray(e); w = 1 / e ** 2; m = float(np.sum(w * v) / w.sum())
    chi2 = float(np.sum(w * (v - m) ** 2) / (len(v) - 1)) if len(v) > 1 else np.nan
    return m, float(1 / np.sqrt(w.sum())), chi2


def analyse(name):
    ds = Dataset(name); mode = MODE.get(name, "spot")
    runs = ds.runs; x = np.array([ds.x[r] for r in runs]); G = sorted(set(ds.group.values())); gi = np.array([G.index(ds.group[r]) for r in runs])
    if name == "testfibres":     # partial-overlap runs (testfibres.py) left out, as in the light fits
        J = json.load(open(f"{DATA}/zscan_testfibres.json"))["selections"]["loose"]; full = dict(zip(J["runs"], J["full_overlap"]))
    else:
        full = {r: True for r in runs}
    hits = {}; aout = {}; eout = {}; ctrl = []; amp_us = []
    for i, r in enumerate(runs):
        D = load_drs(ds, r, mode)
        if D is None: ctrl.append((np.nan, np.nan)); amp_us.append(np.nan); continue
        for j, ch in enumerate(D["chans"]):
            mt = ds.meta[ch]; T, A = D["T"][:, j], D["A"][:, j]; it = intime(D["E"][:, j], T, A, W_T[mt["type"]] if mt["type"] in W_T else W_T["Cer"], A_HIT[mt["is6mm"]])
            if it is None: continue
            side = (np.abs(T) < 3000) & (A > A_HIT[mt["is6mm"]]) & np.any([(T >= it["t0"] + lo) & (T < it["t0"] + hi) for lo, hi in SIDEBANDS], 0)
            hits.setdefault(ch, {})[r] = (T[it["inn"]] * TS_NS, A[it["inn"]], A[side])
            out = (np.abs(D["T"][:, j]) < 3000) & (np.abs(D["T"][:, j] - it["t0"]) > 60)
            aout.setdefault(ch, []).append(D["A"][out, j]); eout.setdefault(ch, []).append(D["E"][out, j])
        d = np.load(drs_dump(name, r)); daq = ds.daq(r); a, b = (d[f"{m}_TS_cfd"] for m in MCP_CTRL)
        s = pid(d, daq) & mcp_ok(d, daq) & (a > 0) & (b > 0) & (d[f"{MCP_CTRL[0]}_peak_value"] > 50)
        ctrl.append(core_sigma((a[s] - b[s]) * TS_NS)); amp_us.append(float(np.median(d[f"{MCP_CTRL[0]}_peak_value"][pid(d, daq) & mcp_ok(d, daq)])))
    chans = {}
    for ch, hr in hits.items():
        mt = ds.meta[ch]; a1, _ = spe_of(np.concatenate(aout[ch]), np.concatenate(eout[ch]), mt["is6mm"])
        if not np.isfinite(a1): continue
        rec = dict(type=mt["type"], label=mt.get("label", ch), is6mm=mt["is6mm"], spe_amp=a1, bins=[])
        for lo, hi in PE_BINS:
            S, N = [], []
            for r in runs:
                if r not in hr or not full.get(r): S.append((np.nan, np.nan)); N.append(np.nan); continue
                t, amp, aside = hr[r]; m = (amp / a1 >= lo) & (amp / a1 < hi); bkg = np.sum((aside / a1 >= lo) & (aside / a1 < hi)) / SPAN_NS
                S.append(core_sigma(t[m], bkg)); N.append(float(np.mean(amp[m] / a1)) if m.any() else np.nan)
            S = np.array(S); ok = np.isfinite(S[:, 0]); b = dict(pe=[lo, hi], sigma=S[:, 0].tolist(), sigma_err=S[:, 1].tolist(), n_pe=N)
            if ok.sum() >= MIN_RUNS and len(set(gi[ok])) < ok.sum() - 1:
                s, es, chi2, icpt = line_fit(x[ok], gi[ok], S[ok, 0] ** 2, 2 * S[ok, 0] * S[ok, 1])
                nm = float(np.nanmean(np.array(N)[ok]))
                b.update(slope=s, slope_err=es, chi2=chi2, icpt=icpt, n_mean=nm, per_photon=s * nm, per_photon_err=es * nm,
                         sigma_mean=float(np.nanmean(S[ok, 0])))
            rec["bins"].append(b)
        fb = [b for b in rec["bins"] if "slope" in b]
        if fb:
            rec["per_photon"], rec["per_photon_err"], rec["per_photon_chi2"] = wmean([b["per_photon"] for b in fb], [b["per_photon_err"] for b in fb])
        chans[ch] = rec
    c = np.array(ctrl); ok = np.isfinite(c[:, 0])
    control = dict(sigma=c[:, 0].tolist(), sigma_err=c[:, 1].tolist(), amp_us=amp_us)
    if ok.sum() >= MIN_RUNS:
        s, es, chi2, _ = line_fit(x[ok], gi[ok], c[ok, 0] ** 2, 2 * c[ok, 0] * c[ok, 1]); control.update(slope=s, slope_err=es, chi2=chi2)
    types = {}
    for t in sorted({v["type"] for v in chans.values()}):
        v = [r for r in chans.values() if r["type"] == t and "per_photon" in r]
        if not v: continue
        m, e, chi2 = wmean([r["per_photon"] for r in v], [r["per_photon_err"] for r in v])
        per_bin = []
        for k in range(len(PE_BINS)):
            q = [(r["bins"][k]["per_photon"], r["bins"][k]["per_photon_err"]) for r in v if "per_photon" in r["bins"][k]]
            per_bin.append(dict(zip(("value", "err", "chi2"), wmean(*zip(*q)))) | dict(n=len(q)) if q else None)
        types[t] = dict(per_photon=m, err=float(e * np.sqrt(max(chi2, 1.0)) if np.isfinite(chi2) else e), chi2=chi2, n=len(v), per_bin=per_bin)
    res = dict(dataset=name, title=ds.cfg["title"], mode=mode, runs=runs, x_m=x.tolist(), groups=[ds.group[r] for r in runs],
               mcp_ref=get_mcp_channels(ds.daq(runs[0]))["MCP_DS_1"],
               full=[bool(full.get(r)) for r in runs], pe_bins=[list(b) for b in PE_BINS], channels=chans, types=types, control=control)
    os.makedirs(OUT, exist_ok=True); json.dump(res, open(f"{OUT}/dispersion_{name}.json", "w"))
    print(f"{name} ({mode}): {len(chans)} channels; control {MCP_CTRL[0]}-{MCP_CTRL[1]}: d(sigma^2)/dX "
          f"{1e3 * control.get('slope', np.nan):+.0f} +- {1e3 * control.get('slope_err', np.nan):.0f} ps^2/mm")
    for t, v in types.items():
        print(f"  {t:12s} {v['n']:3d} ch: per-photon N d(sigma^2)/dX = {1e3 * v['per_photon']:+6.0f} +- {1e3 * v['err']:3.0f} ps^2/mm (chi2 between channels {v['chi2']:.1f});"
              " per p.e. bin: " + "  ".join(f"{1e3 * b['value']:+.0f}±{1e3 * b['err']:.0f}" if b else "-" for b in v["per_bin"]))
    return res


if __name__ == "__main__":
    analyse(sys.argv[1])
