"""Resolutions for the pion_response_final correction (positron depth, measured attenuation, quartz/plastic separately,
own electron scale), both dY options, next to the in-sample (section 4) correction on the same events, and DREAM."""
import json, os, sys
import numpy as np
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
src = open(code("pion_final.py")).read()
# same correction as pion_final.py, but also hand back the per-event signals (on the electron scale) and the in-sample version
src = src.replace("""    Cx = Qx + Px; sE, cE = e_at(E)""", """    Cx = Qx + Px; sE, cE = e_at(E)
    def insample(sig):
        y = np.log(np.maximum(sig, 1e-3)); m = ok.copy()
        Am = A.copy(); Am[:, 1] = dX - np.median(dX); Am[:, 2] = dY - np.median(dY)
        for _ in range(5):
            c, *_ = np.linalg.lstsq(Am[m], y[m], rcond=None)
            r = y - Am @ c; s = 1.4826*np.median(np.abs(r[m]-np.median(r[m]))); m = m & (np.abs(r) < 2.5*s)
        return sig*np.exp(-(Am[:, 1:] @ c[1:]))
    ARR["final"] = (Sx/sE, Cx/cE, Sox/sE); ARR["insample"] = (insample(S)/sE, insample(C)/cE, insample(So)/sE); ARR["raw"] = (S/sE, C/cE, So/sE)""")
exec("ARR = {}\n" + src[:src.index("RUNS = [")])
L = {t: (ATT[t]["lam"], ATT[t]["elam"]) for t in ("Sci", "CerQuartz", "CerPlastic")}
d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
dres = {int(e): (s, c, q) for e, s, c, q in zip(d32["energy_GeV"], d32["S_res_pct"], d32["C_res_pct"], d32["QS_res_pct"])}
def core(x):
    m, s = np.median(x), 1.4826*np.median(np.abs(x - np.median(x)))
    for _ in range(5):
        k = np.abs(x - m) < 2*s; m, s = x[k].mean(), x[k].std()
    return m, s
def resol(S, C, So):
    sig = {"S": S, "C": C, "DR": (S - 0.29*C)/0.71, "QS": qsfun(S, C, So)}
    return {k: (100*v.std()/v.mean(), 100*core(v)[1]/core(v)[0]) for k, v in sig.items()}
out = {}
print("resolution [%]: rms/mean | Gaussian core sigma/mean.  versions: raw, in-sample (section 4), final dY not depth, final dY = depth")
for run in (1441, 1439, 1438, 1437, 1429, 1432, 1433, 1434):
    vers = {}
    for yd in (0, 1):
        ARR.clear(); E, _ = pion(run, L["Sci"][0], L["CerQuartz"][0], L["CerPlastic"][0], yd)
        vers[f"final_dY{yd}"] = resol(*ARR["final"])
        if yd == 0: vers["raw"] = resol(*ARR["raw"]); vers["insample"] = resol(*ARR["insample"])
    out[str(int(E))] = vers
    dd = dres.get(int(E))
    for k in ("S", "C", "DR", "QS"):
        print(f"{E:4.0f} {k:3s} | " + " | ".join(f"{v}: {vers[v][k][0]:5.1f} / {vers[v][k][1]:5.1f}" for v in ("raw", "insample", "final_dY0", "final_dY1"))
              + (f" | DREAM {dd[{'S': 0, 'C': 1, 'QS': 2}[k]]:.1f}" if dd and k != "DR" else ""))
json.dump(out, open(f"{SCR}/pion_final_res.json", "w"), indent=1)
