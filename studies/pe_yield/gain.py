"""FERS HG gain (ADC per p.e.) from the single-p.e. structure of the beam spectra.
Model of a channel's HG spectrum in the gain run: sum over n of a_n x T(x - off - n g), smeared by sqrt(n) s1, where T is the channel's
own pedestal-run spectrum (pedestal + dark counts, which the FERS hold samples at random times and so smears) and a_n >= 0 free
(the in-time light of a shower is not Poisson). For fixed (g, off, s1) the model is linear in a_n (non-negative least squares);
g is scanned and the chi2 profiled over off and s1. Channels with the single-p.e. peaks resolved give a sharp minimum; where the dark
counts smear the peaks the profile is flat and the fit says nothing.
Fitted: every channel with a little light (MIN_LIGHT..MAX_LIGHT above the pedestal run). On GAIN_BOARDS the peaks are resolved and each
channel gives its gain; on the other boards the single profiles are shallow, so the gain is taken per board: median of the channel
minima, error from their spread, and the sum of the channel profiles as a check.
Usage (repository root): python3 studies/pe_yield/gain.py      Output: $PE_YIELD_DATA/gain.json"""
import json
import numpy as np
from config import DATA, PED_RUN, GAIN_RUN, GAIN_BOARDS

EDGES = np.arange(-80, 451, 1.0)       # [ADC] above the pedestal-run median
XC = 0.5 * (EDGES[1:] + EDGES[:-1])
NMAX = 9                               # p.e. components
GAINS = np.arange(45, 100.01, 0.5)
OFFS = np.arange(-4, 4.01, 2.0)
S1S = (0.0, 4.0, 8.0)
MAX_LIGHT = 400                        # [ADC] mean light above the pedestal run: more and the comb is washed out
MIN_LIGHT = 30                         # [ADC] less: mostly zero p.e., no handle on the spacing
MIN_DARK = 30                          # [ADC] 99th percentile of the pedestal run above its median: less = no SiPM signal (dead)
GOOD = dict(err=2.0, chi2ndf=1.5)      # good fit: profile minimum narrower than this and chi2/ndf below


def nnls(A, b):
    """least squares with non-negative coefficients (active set, few columns)"""
    free = np.ones(A.shape[1], bool)
    while True:
        x = np.zeros(A.shape[1]); x[free] = np.linalg.lstsq(A[:, free], b, rcond=None)[0]
        if (x >= 0).all() or free.sum() == 1: return np.clip(x, 0, None)
        free[np.argmin(np.where(free, x, np.inf))] = False


def components(tx, t, off, g, s1):
    """zero-p.e. template shifted by off + n g and smeared by sqrt(n) s1, n = 0..NMAX-1, on the EDGES bins"""
    out = []
    for n in range(NMAX):
        y = np.interp(XC - off - n * g, tx, t, left=0, right=0); sig = np.sqrt(n) * s1
        if sig > 0.5:
            k = np.exp(-0.5 * np.arange(-int(4 * sig), int(4 * sig) + 1) ** 2 / sig ** 2); y = np.convolve(y, k / k.sum(), "same")
        out.append(y)
    return np.stack(out, 1)


def fit(beam, ped):
    p0 = np.median(ped); t, e = np.histogram(ped - p0, np.arange(-200, 1500, 1.0)); tx = e[:-1] + 0.5; t = t / t.sum()
    h = np.histogram(beam - p0, EDGES)[0].astype(float); w = 1 / np.sqrt(h + 1)
    chi = np.full((len(S1S), len(OFFS), len(GAINS)), np.inf); best = None
    for i, s1 in enumerate(S1S):
        for j, off in enumerate(OFFS):
            for k, g in enumerate(GAINS):
                A = components(tx, t, off, g, s1); a = nnls(A * w[:, None], h * w); chi[i, j, k] = np.sum(((A @ a - h) * w) ** 2)
                if best is None or chi[i, j, k] < best[0]: best = (chi[i, j, k], g, off, s1, A @ a, a)
    c2, g, off, s1, model, a = best
    prof = chi.min((0, 1)); ok = prof < c2 + 1
    return dict(gain=float(g), err=float(0.5 * (GAINS[ok].max() - GAINS[ok].min() + GAINS[1] - GAINS[0])), off=float(off), s1=float(s1),
                chi2ndf=float(c2 / (len(h) - NMAX - 3)), edge=bool(g in (GAINS[0], GAINS[-1])), frac=(a / a.sum()).tolist(),
                profile=(prof - c2).tolist(), hist=h.tolist(), model=model.tolist())


def main():
    P = np.load(f"{DATA}/raw/raw_{PED_RUN}.npz"); B = np.load(f"{DATA}/raw/raw_{GAIN_RUN}.npz")
    meta = P["meta"].astype(int); assert (meta == B["meta"]).all()
    ped, beam = P["hg"].astype(float), B["hg"].astype(float)
    dark = np.percentile(ped, 99, 0) - np.median(ped, 0); light = beam.mean(0) - ped.mean(0)
    typ = np.where(meta[:, 2] == 0, "Sci", np.where(meta[:, 3] == 1, "CerQ", "CerP"))
    todo = [i for i in range(len(meta)) if dark[i] > MIN_DARK and MIN_LIGHT < light[i] < MAX_LIGHT]
    res = []
    for n, i in enumerate(todo):
        r = fit(beam[:, i], ped[:, i])
        r.update(board=int(meta[i, 0]), channel=int(meta[i, 1]), type=str(typ[i]), light=float(light[i]), dark_p99=float(dark[i]),
                 good=bool(r["err"] < GOOD["err"] and r["chi2ndf"] < GOOD["chi2ndf"] and not r["edge"]))
        res.append(r)
        print(f"{n + 1}/{len(todo)} B{r['board']} ch{r['channel']} {r['type']}: g {r['gain']:.1f} +- {r['err']:.1f} chi2/ndf {r['chi2ndf']:.2f}"
              f" light {r['light']:.0f}{' good' if r['good'] else ''}", flush=True)
    ref = [r for r in res if r["board"] in GAIN_BOARDS and r["good"]]
    use = np.array([r["gain"] for r in ref])
    summary = dict(gain=float(np.median(use)), spread=float(0.5 * np.diff(np.percentile(use, [16, 84]))[0]), n=int(len(use)),
                   per_type={t: dict(median=float(np.median(v)), n=len(v)) for t in ("Sci", "CerP", "CerQ")
                             if (v := [r["gain"] for r in ref if r["type"] == t])})
    boards = {}
    for b in sorted(set(meta[:, 0])):
        rb = [r for r in res if r["board"] == b and r["chi2ndf"] < GOOD["chi2ndf"] and not r["edge"]]
        if b in GAIN_BOARDS: rb = [r for r in rb if r["good"]]
        if len(rb) < 3: continue
        g = np.array([r["gain"] for r in rb]); prof = np.sum([r["profile"] for r in rb], 0)
        boards[int(b)] = dict(gain=float(np.median(g)), err=float(1.2533 * 0.5 * np.diff(np.percentile(g, [16, 84]))[0] / np.sqrt(len(g))),
                              spread=float(0.5 * np.diff(np.percentile(g, [16, 84]))[0]), n=len(g), profile_sum=prof.tolist(),
                              profile_min=float(GAINS[np.argmin(prof)]), resolved=b in GAIN_BOARDS)
    json.dump(dict(ped_run=PED_RUN, gain_run=GAIN_RUN, boards=list(GAIN_BOARDS), edges=EDGES.tolist(), gains=GAINS.tolist(),
                   summary=summary, per_board=boards, channels=res), open(f"{DATA}/gain.json", "w"))
    print("boards 9/10 (resolved peaks)", summary)
    for b, v in boards.items():
        print(f"  board {b:2d}: {v['gain']:.1f} +- {v['err']:.1f} (spread {v['spread']:.1f}, {v['n']} ch), summed-profile minimum {v['profile_min']:.1f}")


if __name__ == "__main__":
    main()
