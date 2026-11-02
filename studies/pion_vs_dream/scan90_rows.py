"""tb2025 90/93.3-degree 40 GeV positrons: per Y row, ln(type-summed signal) vs position along the fibres, all runs in the row.
Position: table X from the config and, independently, from the DRS arrival time.
15120 / 15130 are the second files of runs 1512 / 1513 (own config entries), analysed with the tb2025 setup in scratchpad/mirror."""
import json, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
from utils.utils import get_hist_mpv
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
rl = json.load(open("data/Runlist.json"))
X = {r: float(rl[str(r)]["Table position X"].replace("mm", "")) for r in list(range(1499, 1514)) + [15120, 15130]}
ROWS = {"Y=1128": [1499, 1506, 1510], "Y~1064": [1501, 1507, 1511, 1513, 15130], "Y=1000": [1502, 1508, 1512], "Y=936": [1504, 1509, 15120]}
MIRROR = SCR + "/mirror"
rootdir = lambda run: f"{MIRROR}/results/root/Run{run}" if run >= 10000 else f"results/root/Run{run}"
EL = ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
TIME_SLOPE = -6.0          # ns per m of table X (2026 scan: -5.4 quartz .. -6.2 Sci)
order = P._build_channel_order(P.build_drs_boards(run_number=1499))
vtype = {o[0]: o[4] for o in order}
TYPES = ("Sci", "CerQuartz", "CerPlastic")

def drs(run):
    p = f"{rootdir(run)}/drs_energy_stats.json"
    if not os.path.exists(p): return None
    st = json.load(open(p)); sums = {t: sum(v[0] for ch, v in st.items() if vtype.get(ch) == t) for t in TYPES}
    f = ROOT.TFile(f"{rootdir(run)}/drs_stats.root"); mpv = {}
    for ch in st:
        h = f.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
        if h and h.GetEntries() >= 50 and st[ch][0] > 100: mpv[ch] = get_hist_mpv(h)[0]
    f.Close()
    return sums, mpv

def fers(run):
    p = f"{SCR}/scan90ch_{run}.npz"
    if not os.path.exists(p): return None
    d = np.load(p); sel = np.ones(len(d["vals"]), bool)
    for f in EL: sel &= d[f].astype(bool)
    v, m = d["vals"][sel], d["meta"]
    mean = v.mean(0); w = np.clip(mean * (m[:, 2] == 0), 0, None)
    sums = {"Sci": mean[m[:, 2] == 0].sum(), "CerQuartz": mean[(m[:, 2] == 1) & (m[:, 3] == 1)].sum(), "CerPlastic": mean[(m[:, 2] == 1) & (m[:, 3] == 0)].sum()}
    return sums, (np.average(m[:, 4], weights=w), np.average(m[:, 5], weights=w))

for row, runs in ROWS.items():
    D = {r: drs(r) for r in runs}; F = {r: fers(r) for r in runs}
    have = [r for r in runs if D[r] or F[r]]
    if len(have) < 2: print(f"== {row}: fewer than 2 runs ready"); continue
    # arrival-time position: median channel MPV shift vs the first run with DRS, converted with TIME_SLOPE
    rd = [r for r in have if D[r]]; xt = {}
    if rd:
        ref = rd[0]
        for r in rd:
            common = [ch for ch in D[r][1] if ch in D[ref][1]]
            dt = np.median([D[r][1][ch] - D[ref][1][ch] for ch in common]) * 0.2 if common else np.nan
            xt[r] = X[ref] + dt / TIME_SLOPE * 1000
    print(f"== {row}")
    print(f"   {'run':>5s} {'deg':>5s} {'X cfg':>7s} {'X time':>7s} | {'ix':>6s} {'iy':>6s} | " + " ".join(f"{'DRS ' + t:>14s}" for t in TYPES) + " | " + " ".join(f"{'FERS ' + t:>15s}" for t in TYPES))
    for r in have:
        print(f"   {r:5d} {rl[str(r)]['Table Degree']:>5s} {X[r]:7.1f} {xt.get(r, np.nan):7.1f} | "
              + (f"{F[r][1][0]:+6.2f} {F[r][1][1]:+6.2f}" if F[r] else f"{'':6s} {'':6s}") + " | "
              + " ".join(f"{D[r][0][t]:14.0f}" if D[r] else f"{'':>14s}" for t in TYPES) + " | "
              + " ".join(f"{F[r][0][t]:15.2f}" if F[r] else f"{'':>15s}" for t in TYPES))
    # runs whose shower lands in the same towers (FERS energy-weighted ix within 0.3 tower of the row median)
    ixs = {r: F[r][1][0] for r in have if F[r]}
    med = np.median(list(ixs.values())) if ixs else np.nan
    same = [r for r in have if r in ixs and abs(ixs[r] - med) < 0.3]
    for src, S, subset in (("DRS", D, have), ("FERS", F, have), ("DRS same-ix", D, same), ("FERS same-ix", F, same)):
        rr = [r for r in subset if S[r]]
        if len(set(X[r] for r in rr)) < 2: continue
        for xlab, xs in (("X cfg", {r: X[r] for r in rr}), ("X time", {r: xt.get(r, np.nan) for r in rr})):
            x = np.array([xs[r] for r in rr]) / 1000
            if np.isnan(x).any(): continue
            out = []
            for t in TYPES:
                y = np.log([S[r][0][t] for r in rr]); p = np.polyfit(x, y, 1); res = 100 * (y - np.polyval(p, x))
                out.append(f"{t} {p[0]:+.3f}/m (lambda {abs(1 / p[0]):5.2f} m) res% [" + " ".join(f"{e:+.1f}" for e in res) + "]")
            print(f"   {src:12s} vs {xlab:6s} runs {rr}: " + " ; ".join(out))
