"""Shape of the per-channel CFD-integral distribution at several scan points (positrons, same selection as the analysis):
overlays normalised to unit area, a Gaussian fitted to the core of each, and skewness / core-fit quality per point.
Usage: python3 studies/zscan/checks/integral_shapes.py DATASET     Plots: results/plots/ZScan/<dataset>/checks/
"""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ROOT
from common import Dataset, DATA, load_drs, SAT_FRAC
from datasets import PLOTS
from plotting.my_function import DrawHistos
ds = Dataset(sys.argv[1]); J = json.load(open(f"{DATA}/zscan_{ds.name}.json")); MODE = J["main_selection"]; D = J["selections"][MODE]["drs"]; C = D["channels"]
OUTD = f"{PLOTS}/{ds.name}/checks"; os.makedirs(OUTD, exist_ok=True)


def pick(t, bright):
    """brightest (or a dim but well-measured) unsaturated channel of a type"""
    c = sorted((np.median(p["all"]), ch) for ch, p in C.items() if p["type"] == t and p["high_amp"]["1500"] < SAT_FRAC
               and p["fit"]["lph"]["main"][J["main_fit"]][2] < 0.05)
    return c[-1][1] if bright else c[len(c) // 4][1]


runs = D["runs"]; idx = np.linspace(0, len(runs) - 1, 5).astype(int); use = [runs[i] for i in idx]
data = {}
for r in use:
    L = load_drs(ds, r, MODE); data[r] = (L["chans"], L["E"])
COLS = [ROOT.kBlue + 1, ROOT.kAzure + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kRed + 1]
names = []
for t, bright, tag in (("Sci", True, "bright"), ("Sci", False, "dim"), ("CerQuartz", True, "bright"), ("CerPlastic", True, "bright")):
    ch = pick(t, bright); p = C[ch]
    allv = np.concatenate([data[r][1][:, data[r][0].index(ch)] for r in use]); lo, hi = np.percentile(allv, [0.2, 99.8])
    lo = min(lo, -0.05 * hi)
    for logy in (False, True):
        hs, lg = [], []
        print(f"{ds.name} {t} {tag} {ch} (all-event mean {np.median(p['all']):.0f} ADC)")
        for r, col in zip(use, COLS):
            e = data[r][1][:, data[r][0].index(ch)]
            h = ROOT.TH1F(f"h_{ch}_{r}_{logy}", "", 80, lo, hi)
            for v in e: h.Fill(float(v))
            h.Scale(1.0 / h.Integral()); h.SetDirectory(0); hs.append(h)
            # Gaussian on the core (+-1.5 sigma around the peak, iterated)
            m, s = np.median(e), 1.4826 * np.median(np.abs(e - np.median(e)))
            for _ in range(5):
                k = np.abs(e - m) < 1.5 * s; m, s = e[k].mean(), e[k].std() / 0.7368    # 0.7368: rms of a Gaussian truncated at 1.5 sigma
            sk = ((e - e.mean()) ** 3).mean() / e.std() ** 3
            frac_tail = float((e > m + 3 * s).mean()); frac_low = float((e < m - 3 * s).mean())
            x = 1000 * ds.x[r]
            lg.append(f"X = {x:+.0f} mm, mean {e.mean():.0f}")
            if not logy:
                print(f"   X {x:+7.1f}: mean {e.mean():7.1f}  median {np.median(e):7.1f}  core Gaussian mu {m:7.1f} sigma {s:6.1f} (sigma/mu {s / m:.2f})  "
                      f"skew {sk:+.2f}  beyond +3 sigma {100 * frac_tail:.1f}%  below -3 sigma {100 * frac_low:.1f}%")
        nm = f"integral_{t}_{tag}" + ("_log" if logy else "")
        DrawHistos(hs, lg, lo, hi, "CFD integral [ADC]", 1e-4 if logy else 0, max(h.GetMaximum() for h in hs) * (30 if logy else 1.7),
                   "Fraction of events", nm, dology=logy, drawoptions="HIST", legendoptions="L", mycolors=COLS,
                   legendPos=[0.50, 0.66, 0.93, 0.89], outdir=OUTD, lumitext=f"{ds.name}, {t} ({tag})", usePDF=True)
        names.append(nm)
json.dump(names, open(f"{DATA}/checks/integral_shape_plots_{ds.name}.json", "w"))
