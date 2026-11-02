"""Average waveform vs time slice away from the pulse (drs_profiles.root, blsub vs raw TS, all events): shape of the baseline
along the 1024-sample record, its offset under the in-time pulse (interpolated from samples 100-380 and 650-1000 with a quadratic),
what that offset adds to the long (50-sample) and short (15-sample) integrals, and whether it changes along the scan.
Usage: python3 studies/zscan/checks/baseline_shape.py DATASET     Output: $ZSCAN_DATA/checks/baseline_shape_<dataset>.json"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ROOT
ROOT.gROOT.SetBatch(True)
from common import Dataset, DATA
ds = Dataset(sys.argv[1]); J = json.load(open(f"{DATA}/zscan_{ds.name}.json")); D = J["selections"][J["main_selection"]]["drs"]
FIT_RANGES = ((100, 380), (650, 1000)); T_PULSE = 480   # drs_profiles.root: framework per-run DRS profiles (results/root/Run<N>)
runs = D["runs"]; xs = np.array(D["x_m"])
prof = {}
for r in runs:
    f = ROOT.TFile(f"results/root/Run{r}/drs_profiles.root")
    for ch in ds.meta:
        h = f.Get(f"prof_{ch}_blsub_VS_ts")
        if h: prof.setdefault(ch, []).append(np.array([h.GetBinContent(b) for b in range(1, 1025)]))
    f.Close()
ts = np.arange(1024) + 0.5; mfit = np.zeros(1024, bool)
for lo, hi in FIT_RANGES: mfit |= (ts >= lo) & (ts < hi)
res = {}
for ch, P in prof.items():
    P = np.array(P)
    if len(P) != len(runs): continue
    mean = P.mean(0)
    c = np.polyfit(ts[mfit], mean[mfit], 2); at_pulse = np.polyval(c, T_PULSE)
    per_run = np.array([np.polyval(np.polyfit(ts[mfit], p[mfit], 2), T_PULSE) for p in P])
    win = lambda lo, hi: float(mean[(ts >= lo) & (ts < hi)].mean())
    lph = np.array(D["channels"][ch]["lph"]) if ch in D["channels"] else None
    res[ch] = dict(type=ds.meta[ch]["type"], w0_200=win(0, 200), w300_380=win(300, 380), w650_750=win(650, 750), w800_1000=win(800, 1000),
                   at_pulse=float(at_pulse), slope_at_pulse_per_m=float(np.polyfit(xs, per_run, 1)[0]), rms_runs=float(per_run.std()),
                   long_bias=float(50 * at_pulse), short_bias=float(15 * at_pulse),
                   long_bias_rel=float(50 * at_pulse / np.median(lph)) if lph is not None else float("nan"),
                   mean_signal=float(np.median(D["channels"][ch]["all"])) if ch in D["channels"] else np.nan)
print(f"{ds.name}: {len(res)} channels; average baseline (blsub) in windows of time slices [ADC], median over channels (16-84%):")
for k in ("w0_200", "w300_380", "at_pulse", "w650_750", "w800_1000"):
    v = np.array([p[k] for p in res.values()]); print(f"   {k:10s}: {np.median(v):+.2f} ({np.percentile(v, 16):+.2f}..{np.percentile(v, 84):+.2f})")
for t in ("Sci", "CerQuartz", "CerPlastic"):
    ps = [p for p in res.values() if p["type"] == t]
    print(f"   {t:10s}: offset under the pulse {np.median([p['at_pulse'] for p in ps]):+.2f} ADC/sample -> long integral {np.median([p['long_bias'] for p in ps]):+.0f} ADC "
          f"({100 * np.nanmedian([p['long_bias_rel'] for p in ps]):+.1f}% of the light per hit), short {np.median([p['short_bias'] for p in ps]):+.0f} ADC; "
          f"run-to-run rms {np.median([p['rms_runs'] for p in ps]):.2f} ADC/sample; change along the scan {np.median([p['slope_at_pulse_per_m'] for p in ps]):+.2f} ADC/sample/m "
          f"(16-84% {np.percentile([p['slope_at_pulse_per_m'] for p in ps], 16):+.2f}..{np.percentile([p['slope_at_pulse_per_m'] for p in ps], 84):+.2f})")
dim = sorted(res, key=lambda c: res[c]["mean_signal"] if np.isfinite(res[c]["mean_signal"]) else 0)[:5]
print("   dimmest channels (closest to 'not hit'): " + "; ".join(f"{c.replace('DRS_', '')}: 0-200 {res[c]['w0_200']:+.1f}, 300-380 {res[c]['w300_380']:+.1f}, "
                                                               f"650-750 {res[c]['w650_750']:+.1f}, 800-1000 {res[c]['w800_1000']:+.1f}" for c in dim))
json.dump(res, open(f"{DATA}/checks/baseline_shape_{ds.name}.json", "w"), indent=1)
