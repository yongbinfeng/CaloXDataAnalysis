"""2026 90-degree DRS scan: attenuation per type, split by how often a channel's pulse exceeds 800 ADC (clipping risk)."""
import json, os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
sys.path.insert(0, "."); sys.path.insert(0, "scripts/compare")
from utils.root_setup import setup_root
setup_root(batch_mode=True, load_functions=False)
import plot_drs_mpv_scan as P
rl = json.load(open("data/Runlist.json"))
runs = [r for r in P.SCAN_RUNS if os.path.exists(f"results/root/Run{r}/drs_energy_stats.json")]
X = {r: float(rl[str(r)]["Table position X"].replace("mm", "")) for r in runs}
st = {r: json.load(open(f"results/root/Run{r}/drs_energy_stats.json")) for r in runs}
ovf = {}
for r in runs:
    f = ROOT.TFile(f"results/root/Run{r}/drs_stats.root")
    for ch in st[r]:
        h = f.Get(f"hist_{ch}_peak_value")
        if h and h.GetEntries() > 0:
            ovf[ch] = max(ovf.get(ch, 0.0), h.GetBinContent(h.GetNbinsX() + 1) / h.GetEntries())
    f.Close()
order = P._build_channel_order(P.build_drs_boards(run_number=P.REFERENCE_RUN))
rows = []
for ch, b, g, c, var, tx, ty, is6 in order:
    pts = [(X[r], st[r][ch][0]) for r in runs if ch in st[r]]
    if len(pts) < 15: continue
    x, m = np.array(pts).T
    if np.median(m) < 100: continue
    rows.append((var, ovf.get(ch, 0.0), np.polyfit(x, np.log(m), 1)[0] * 1000, np.median(m)))
for lab, lo, hi in (("never > 800 ADC (max ovf < 0.1%)", -1, 0.001), ("0.1-5% of events > 800 ADC", 0.001, 0.05), ("> 5% of events > 800 ADC", 0.05, 2)):
    print(lab)
    for var in P._TYPE_ORDER:
        s = np.array([r[2] for r in rows if r[0] == var and lo < r[1] <= hi]); w = np.array([r[3] for r in rows if r[0] == var and lo < r[1] <= hi])
        if len(s) == 0: print(f"   {var:10s}: none"); continue
        sw = np.average(s, weights=w)
        print(f"   {var:10s}: n {len(s):3d}  slope median {np.median(s):+.3f}/m (lambda {abs(1/np.median(s)):5.1f} m), weighted {sw:+.3f}/m (lambda {abs(1/sw):5.1f} m)")
