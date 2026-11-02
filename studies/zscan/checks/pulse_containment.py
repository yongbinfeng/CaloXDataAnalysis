"""Is the short DRS integral (15 samples = 3 ns, -4..+11 from the CFD crossing) long enough for scintillation pulses?
Average in-time pulse per fibre type (main-detector channels, positron runs), aligned at the CFD anchor sample, and the fraction
of the full pulse area (to +200 samples = 40 ns) inside the short (-4..+11) and long (-5..+45) windows, for small (1-2 p.e.) and
large pulses. Usage (repository root): python3 studies/zscan/checks/pulse_containment.py RUN DATASET"""
import json, os, sys
run = int(sys.argv[1]) if len(sys.argv) > 1 else 1897; dataset = sys.argv[2] if len(sys.argv) > 2 else "2026A"
sys.argv = ["x", "--run", str(run), "--last-event", "30000"]
sys.path.insert(0, "."); sys.path.insert(0, "studies/zscan")
from utils.root_setup import setup_root
setup_root(n_threads=10, batch_mode=True, load_functions=True)
import ROOT
import numpy as np
from core.analysis_manager import CaloXAnalysisManager
from utils.parser import get_args
from common import Dataset, DATA
META = Dataset(dataset).meta

PRE, POST = 20, 200
ROOT.gInterpreter.Declare(f"""
ROOT::RVec<float> pc_seg(const ROOT::RVec<float>& w, float ts) {{
    ROOT::RVec<float> o({PRE + POST}, -99999.f);
    if (ts < 0) return o;
    const int a = int(ts) + 1;                          // anchor = first sample above 20%, as compute_cfd_integral
    for (int i = 0; i < {PRE + POST}; ++i) {{ int k = a - {PRE} + i; if (k >= 0 && k < (int)w.size()) o[i] = w[k]; }}
    return o;
}}""")
m = CaloXAnalysisManager(get_args()).prepare()
m.apply_beam_pid_selection(flag_only=True, particle="electron")
rdf = m.get_rdf()
cols = set(map(str, rdf.GetColumnNames()))
chans = [c for c in META if f"{c}_blsub" in cols]
for c in chans:
    rdf = rdf.Define(f"seg_{c}", f"pc_seg({c}_blsub, {c}_TS_cfd)")
d = rdf.AsNumpy([f"seg_{c}" for c in chans] + [f"{c}_TS_cfd_mcp" for c in chans] + [f"{c}_peak_value" for c in chans])
x = np.arange(-PRE, POST)
SHORT, LONG = (x >= -4) & (x < 11), (x >= -5) & (x < 45)
res = {}
for c in chans:
    T = np.asarray(d[f"{c}_TS_cfd_mcp"]); A = np.asarray(d[f"{c}_peak_value"])
    win = T[(T > 380) & (T < 720)]
    if len(win) < 300: continue
    h, ed = np.histogram(win, bins=np.arange(380, 724, 2.0)); t0 = ed[np.argmax(h)] + 1
    w = 6.0 if META[c]["type"] != "Sci" else 20.0
    S = np.stack([np.asarray(s, dtype=np.float32) for s in d[f"seg_{c}"]])
    ok = (np.abs(T - t0) < w) & (S.min(1) > -9e4) & (A < 1500)
    for cls, sel in (("small", ok & (A > 40) & (A < 150)), ("large", ok & (A > 400))):
        if sel.sum() < 50: continue
        p = S[sel].mean(0); full = p[x >= -5].sum()
        res.setdefault(c, {"type": META[c]["type"]})[cls] = dict(n=int(sel.sum()), short=float(p[SHORT].sum() / full),
                                                                  long=float(p[LONG].sum() / full), shape=p.tolist())
json.dump(res, open(f"{DATA}/checks/pulse_containment_{run}.json", "w"))
for t in sorted({v["type"] for v in res.values()}):
    for cls in ("small", "large"):
        v = [q[cls] for q in res.values() if q["type"] == t and cls in q]
        if v: print(f"{t:10s} {cls:5s} ({len(v):2d} ch, {np.median([q['n'] for q in v]):6.0f} hits/ch): short {np.median([q['short'] for q in v]):.2f}"
                    f"  long {np.median([q['long'] for q in v]):.2f}  of the area to +40 ns")
