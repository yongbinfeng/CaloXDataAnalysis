"""Baseline drift along the DRS record, from the integral-variant dumps: per event, the same pulse integrated with the constant
baseline (median of samples 0-200) and with the two-anchor linear baseline (line through samples 0-120 and 800-1000); their
difference is the integrated drift under the integration window. Per channel and run: its median over in-time hits, relative to
the light per hit; spread over channels and boards; and whether it changes along the scan (only that would bias an attenuation slope).
Usage: python3 studies/zscan/checks/baseline_drift.py DATASET     Output: $ZSCAN_DATA/checks/baseline_drift_<dataset>.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from common import Dataset, DATA, chip, pid_hodo, mcp_ok

ds = Dataset(sys.argv[1]); runs = ds.runs
rows = {}                                         # ch -> list over runs of (dL, dS, lph_long, lph_short, n)
for r in runs:
    d = np.load(f"{DATA}/drs/drs_{r}.npz"); names = list(d["channels"]); daq = ds.daq(r)
    sel = pid_hodo(d, daq) & mcp_ok(d, daq) & (np.abs(d["T"]) < 3000).any(1)
    if sel.sum() < ds.cfg["min_events"]: continue
    E, El, Es, Esl, T = (np.ascontiguousarray(d[k][sel].T) for k in ("E", "El", "Es", "Esl", "T"))
    for i, ch in enumerate(names):
        if ch not in ds.meta: continue
        t = T[i]; ok = np.abs(t) < 3000; w = t[ok & (t > 380) & (t < 720)]
        if len(w) < 100: continue
        h, ed = np.histogram(w, bins=np.arange(380, 724, 2.0)); t0 = ed[int(np.argmax(h))] + 1
        win = 20.0 if ds.meta[ch]["type"] == "Sci" else 6.0
        inn = ok & (np.abs(t - t0) < win)
        if inn.sum() < 50: continue
        rows.setdefault(ch, []).append((r, float(np.median(E[i][inn] - El[i][inn])), float(np.median(Es[i][inn] - Esl[i][inn])),
                                        float(np.mean(E[i][inn])), float(np.mean(Es[i][inn])), int(inn.sum())))
xs = ds.x
per = {}
for ch, v in rows.items():
    if len(v) < 0.8 * len(runs): continue
    r_, dL, dS, lL, lS, n = map(np.array, zip(*v)); x = np.array([xs[q] for q in r_])
    relL, relS = dL / lL, dS / lS
    sL = np.polyfit(x, relL, 1)[0] if np.ptp(x) > 0 else np.nan; sS = np.polyfit(x, relS, 1)[0] if np.ptp(x) > 0 else np.nan
    per[ch] = dict(type=ds.meta[ch]["type"], chip=chip(ch), drift_long_adc=float(np.median(dL)), drift_short_adc=float(np.median(dS)),
                   rel_long=float(np.median(relL)), rel_short=float(np.median(relS)), rel_long_rms_runs=float(relL.std()),
                   rel_short_rms_runs=float(relS.std()), slope_rel_long_per_m=float(sL), slope_rel_short_per_m=float(sS))
print(f"{ds.name}: {len(runs)} runs, {len(per)} channels")
for t in ("Sci", "CerQuartz", "CerPlastic"):
    ps = [p for p in per.values() if p["type"] == t]
    if not ps: continue
    q = lambda k, f=1: np.percentile([f * p[k] for p in ps], [16, 50, 84])
    print(f"   {t:10s} n {len(ps):3d}: drift under the long integral {q('drift_long_adc')[1]:+.0f} ADC (16-84% {q('drift_long_adc')[0]:+.0f}..{q('drift_long_adc')[2]:+.0f}), "
          f"= {q('rel_long', 100)[1]:+.2f}% of the light per hit (16-84% {q('rel_long', 100)[0]:+.2f}..{q('rel_long', 100)[2]:+.2f}); short: {q('drift_short_adc')[1]:+.0f} ADC, "
          f"{q('rel_short', 100)[1]:+.2f}%")
    print(f"   {'':10s}        run-to-run rms of that fraction: long {100 * np.median([p['rel_long_rms_runs'] for p in ps]):.2f}%, short {100 * np.median([p['rel_short_rms_runs'] for p in ps]):.2f}%; "
          f"change along the scan: long {100 * np.median([p['slope_rel_long_per_m'] for p in ps]):+.2f}%/m (16-84% {np.percentile([100 * p['slope_rel_long_per_m'] for p in ps], 16):+.2f}..{np.percentile([100 * p['slope_rel_long_per_m'] for p in ps], 84):+.2f}), "
          f"short {100 * np.median([p['slope_rel_short_per_m'] for p in ps]):+.2f}%/m")
# by DRS board: is the drift a board property?
brd = {}
for ch, p in per.items(): brd.setdefault("_".join(ch.split("_")[1:3]), []).append(p["drift_long_adc"])
print("   drift under the long integral by DRS board [ADC]: " + ", ".join(f"{b}: {np.median(v):+.0f} ({len(v)})" for b, v in sorted(brd.items())))
os.makedirs(f"{DATA}/checks", exist_ok=True)
json.dump(per, open(f"{DATA}/checks/baseline_drift_{ds.name}.json", "w"), indent=1)
