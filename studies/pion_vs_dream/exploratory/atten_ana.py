import numpy as np
PION = ["is_TTUMuonVeto_vetoed", "is_PSD_vetoed", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired", "is_HoleVeto_vetoed"]
NS_PER_CM_SCI = 1 / 16.9 - 1 / 29.98          # light in scintillating fibre (~16.9 cm/ns) vs particle at c
dz_all, S_all, Cq_all, Cp_all, E_all = [], [], [], [], []
for run in (1429, 1432, 1433):
    d = np.load(f"time_{run}.npz"); E = float(d["energy"])
    sel = d["MCP_DS_1_peak_value"] > 50
    for f in PION: sel &= d[f].astype(bool)
    T, A = d["Sci_t"][sel], d["Sci_a"][sel]
    good = (A > 100) & (A < 1400)
    ref = np.array([np.median(T[good[:, j], j]) if good[:, j].sum() > 20 else np.nan for j in range(T.shape[1])])
    dt = T - ref; w = A * good * np.isfinite(dt) * (np.abs(dt) < 20)
    tbar = np.nansum(np.where(w > 0, w * dt, 0), 1) / np.maximum(w.sum(1), 1e-6)
    ok = w.sum(1) > 0
    dz = -tbar * 0.2 / NS_PER_CM_SCI           # earlier light = deeper shower; cm, relative to the run median
    dz_all.append(dz[ok] - np.median(dz[ok]))
    S_all.append(d["FERS_energyMix_Sci_calib_sum"][sel][ok] / E); Cq_all.append(d["FERS_Cq"][sel][ok] / E); Cp_all.append(d["FERS_Cp"][sel][ok] / E)
    E_all.append(np.full(ok.sum(), E))
dz, S, Cq, Cp = map(np.concatenate, (dz_all, S_all, Cq_all, Cp_all))
print(f"events with a valid MCP, 80+100+120 GeV pions: {len(dz)};  depth proxy spread (16-84%): {np.percentile(dz, 16):+.0f} .. {np.percentile(dz, 84):+.0f} cm")
edges = np.percentile(dz, [0, 10, 25, 50, 75, 90, 97, 100])
print(f"{'depth bin [cm]':>18s} {'N':>5s} {'S/E':>7s} {'Cq/E':>7s} {'Cp/E':>7s} {'Cq/Cp':>7s}")
cent, mS, mQ, mP = [], [], [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (dz >= a) & (dz < b)
    cent.append(np.median(dz[m])); mS.append(np.median(S[m])); mQ.append(np.median(Cq[m])); mP.append(np.median(Cp[m]))
    print(f"{a:+8.0f} .. {b:+5.0f}  {m.sum():5d} {mS[-1]:7.3f} {mQ[-1]:7.3f} {mP[-1]:7.3f} {mQ[-1] / mP[-1]:7.3f}")
core = (dz > np.percentile(dz, 5)) & (dz < np.percentile(dz, 95))
for name, x in (("S", S), ("C quartz", Cq), ("C plastic", Cp)):
    good = core & (x > 0.05)
    slope = np.polyfit(dz[good], np.log(x[good]), 1)[0]
    print(f"{name:10s}: d ln(signal)/dz = {slope * 100:+.2f} % per cm  ->  effective attenuation length {1 / slope / 100 if slope > 0 else float('inf'):.2f} m")
