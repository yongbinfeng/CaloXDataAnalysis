"""Which DRS integral gives the best attenuation measurement? Same light-per-hit analysis for each variant of
dump_drs.py: default long / short integral, linear baseline long / short, linear baseline with the peak searched
only in the in-time region, long / short.
Light per hit as in analysis/intime_analysis.py: in time = |T - peak| < 6 TS (Cherenkov) or 20 TS (scintillating), peak above
40 ADC; randoms from sidebands 100-300 TS on both sides (whole-record search) or inside the search window (in-time search). Photoelectrons: single-p.e. amplitude peak of the out-of-time pulses
(searched in 80-200 ADC) and the short integral of pulses within that peak.
Same event selection as the analysis (common.py). Usage: python3 studies/zscan/checks/variants.py DATASET
Output: $ZSCAN_DATA/checks/variants_<dataset>.json and a printed comparison."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import Dataset, DATA, TYPES, MIN_SIGNAL, pid_hodo, mcp_ok, beam_window, fit_channels, wmean


VARIANTS = {"default long": ("E", "T", "A"), "default short": ("Es", "T", "A"), "linear bl. long": ("El", "T", "A"),
            "linear bl. short": ("Esl", "T", "A"), "lin. + in-time search, long": ("Eli", "Tli", "Ali"),
            "lin. + in-time search, short": ("Esli", "Tli", "Ali")}
W_CER, W_SCI, SB = 6.0, 20.0, (100.0, 300.0)
# With the peak searched only in the in-time region, every event has a "pulse" there (in a noise-only window, its largest noise
# excursion, 13-42 ADC typically), and the 100-300 TS sidebands lie outside the search window. Hits therefore need a peak above
# A_HIT (about half a photoelectron) in all variants, and the in-time-search randoms come from sidebands inside the search window.
A_HIT = 40.0
SB_IN = ((-90.0, -40.0), (60.0, 100.0))
PE_AMP = (80.0, 200.0)


def peak_in(x, lo, hi, step):
    h, ed = np.histogram(x, bins=np.arange(lo, hi + step, step)); k = int(np.argmax(h)); return ed[k] + step / 2 if h[k] >= 20 else np.nan


def per_run(run, chans, variants=VARIANTS):
    d = np.load(f"{DATA}/drs/drs_{run}.npz"); names = list(d["channels"]); daq = ds.daq(run)
    s1 = pid_hodo(d, daq) & mcp_ok(d, daq) & (np.abs(d["T"]) < 3000).any(1)
    if s1.sum() < ds.cfg["min_events"]: return None
    sel = s1 & beam_window(d, s1, ds, run, "spot")
    # read each array once (indexing an NpzFile re-reads and decompresses the array every time)
    keys = {k for v in variants.values() for k in v}
    arr = {k: np.ascontiguousarray(d[k][sel].T) for k in keys}          # channel-major
    out = {}
    for ch in chans:
        if ch not in names: continue
        i = names.index(ch); w = W_SCI if ds.meta[ch]["type"] == "Sci" else W_CER; res = {}
        for vname, (ek, tk, ak) in variants.items():
            E, T, A = arr[ek][i], arr[tk][i], arr[ak][i]
            ok = np.abs(T) < 3000; E, T, A = E[ok], T[ok], A[ok]; n = len(T)
            win = T[(T > 380) & (T < 720)]
            if len(win) < 100: res[vname] = None; continue
            h, ed = np.histogram(win, bins=np.arange(380, 724, 2.0)); k = int(np.argmax(h)); t0 = float(np.median(win[np.abs(win - ed[k] - 1) < 8]))
            dt = T - t0; big = A > A_HIT; inn = (np.abs(dt) < w) & big
            if tk == "T":                                     # whole-record search: sidebands far from the peak
                side = (np.abs(dt) >= SB[0]) & (np.abs(dt) < SB[1]) & big; span = 2 * (SB[1] - SB[0])
            else:                                             # in-time search: sidebands inside the search window
                side = np.zeros(len(dt), bool); span = 0.0
                for lo, hi in SB_IN: side |= (dt >= lo) & (dt < hi); span += hi - lo
                side &= big
            n_rand = side.sum() * (2 * w) / span; n_hit = inn.sum() - n_rand
            if n_hit < 50: res[vname] = None; continue
            e_side = E[side].mean() if side.sum() else 0.0
            lph = (E[inn].sum() - n_rand * e_side) / n_hit; err = E[inn].std() / np.sqrt(n_hit)
            # single photoelectron from the out-of-time pulses (|dt| beyond 100 TS): amplitude peak, then this integral at that amplitude
            out_t = np.abs(dt) >= SB[0]; a1 = peak_in(A[out_t], *PE_AMP, 4.0)
            spe = float(np.median(E[out_t & (np.abs(A - a1) < 0.2 * a1)])) if np.isfinite(a1) and (out_t & (np.abs(A - a1) < 0.2 * a1)).sum() > 30 else np.nan
            res[vname] = (lph, err, n_hit / n, float((A[inn] > 1500).mean()), spe, a1, n_hit)
        out[ch] = res
    return out

ds = Dataset(sys.argv[1]); runs = list(ds.runs); chans = list(ds.meta)
summ = {}
for r in list(runs):
    s_ = per_run(r, chans)
    if s_ is None: runs.remove(r); continue
    summ[r] = s_
xs = np.array([ds.x[r] for r in runs]); groups = [ds.group[r] for r in runs]
results = {"runs": runs, "x_m": xs.tolist(), "selection": "mcp", "variants": {}}
print(f"{ds.name}: {len(runs)} runs")
for vname in VARIANTS:
    lph, err, frac, sat, spe = {}, {}, {}, {}, {}
    for ch in chans:
        vals = [summ[r].get(ch, {}).get(vname) for r in runs]
        if any(v is None for v in vals): continue
        lph[ch] = np.array([v[0] for v in vals]); err[ch] = np.array([v[1] for v in vals]); frac[ch] = np.array([v[2] for v in vals])
        sat[ch] = max(v[3] for v in vals); spe[ch] = np.nanmedian([v[4] for v in vals])
    good = [ch for ch in lph if np.all(lph[ch] > MIN_SIGNAL) and np.all(err[ch] > 0) and sat[ch] < 0.01]
    att, _ = fit_channels(xs, groups, lph, err, good, per_chip=ds.cfg["offsets"] == "chip")
    # restricting the peak search leaves no out-of-time pulses: take the single p.e. from the matching full-search variant
    if "in-time search" in vname:
        ref = results["variants"]["linear bl. " + ("short" if "short" in vname else "long")]
        spe = {ch: (ref[ch]["spe"] if ch in ref else np.nan) for ch in spe}
    per = {ch: dict(type=ds.meta[ch]["type"], inv_lambda=float(att[ch][0]), err=float(att[ch][1]), resid=float(att[ch][2]),
                    stat=float(np.median(err[ch] / lph[ch])), hit_fraction=float(np.median(frac[ch])), spe=float(spe[ch]),
                    pe_per_hit=float(np.median(lph[ch]) / spe[ch]) if spe[ch] > 0 else float("nan"), light_per_hit=lph[ch].tolist())
           for ch in good}
    results["variants"][vname] = per
    line = []
    for t in TYPES:
        ps = [p for p in per.values() if p["type"] == t and p["resid"] < 0.06]
        if not ps: line.append(f"{t}: -"); continue
        m, e, chi2 = wmean([p["inv_lambda"] for p in ps], [p["err"] for p in ps])
        line.append(f"{t}: n {len(ps)} stat {100 * np.median([p['stat'] for p in ps]):.1f}% -> {m:+.3f}+-{e:.3f}/m chi2 {chi2:.1f} "
                    f"hits {np.median([p['hit_fraction'] for p in ps]):.2f} pe/hit {np.nanmedian([p['pe_per_hit'] for p in ps]):.1f}")
    print(f"   {vname:30s} | " + " | ".join(line))
os.makedirs(f"{DATA}/checks", exist_ok=True)
json.dump(results, open(f"{DATA}/checks/variants_{ds.name}.json", "w"), default=float)
