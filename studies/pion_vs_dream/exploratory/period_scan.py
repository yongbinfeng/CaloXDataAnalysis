import numpy as np
SEL = ["is_TTUMuonVeto_vetoed", "is_PSD_fired", "is_Cer474_fired", "is_Cer519_fired", "is_Cer537_fired"]
HX, HY, S_, C_ = [], [], [], []
for run in (1526, 1527):
    d = np.load(f"hodo_{run}.npz"); E = float(d["energy"])
    sel = np.ones(len(d["TTU_Hodo_X"]), bool)
    for f in SEL: sel &= d[f].astype(bool)
    sel &= (d["TTU_Hodo_nHitX"] >= 1) & (d["TTU_Hodo_nHitX"] <= 2) & (d["TTU_Hodo_nHitY"] >= 1) & (d["TTU_Hodo_nHitY"] <= 2)
    S, C = d["FERS_energyMix_Sci_calib_sum"][sel] / E, d["FERS_energyMix_Cer_calib_sum"][sel] / E
    k = (S > 0.8) & (S < 1.3)
    HX.append(d["hodoX_c"][sel][k] * 10); HY.append(d["hodoY_c"][sel][k] * 10)        # mm, hodoscope frame
    S_.append(S[k] / np.median(S[k])); C_.append(C[k] / np.median(C[k]))
hx, hy, S, C = map(np.concatenate, (HX, HY, S_, C_))
# remove a smooth trend in (x, y): quadratic, so cm-scale effects do not leak into the periodogram
A = np.column_stack([np.ones(len(hx)), hx, hy, hx * hx, hy * hy, hx * hy])
def detrend(v):
    c, *_ = np.linalg.lstsq(A, v, rcond=None); return v - A @ c
rS, rC = detrend(S), detrend(C)
periods = np.arange(1.5, 12.01, 0.05)
rng = np.random.default_rng(1)
def scan(coord, res):
    amps = []
    for P in periods:
        ph = 2 * np.pi * coord / P
        M = np.column_stack([np.cos(ph), np.sin(ph)])
        a, *_ = np.linalg.lstsq(M, res, rcond=None); amps.append(np.hypot(*a))
    amps = np.array(amps)
    # noise level: same scan with the residuals shuffled between events
    noise = []
    for _ in range(20):
        sh = rng.permutation(res); nm = []
        for P in periods[::10]:
            ph = 2 * np.pi * coord / P; M = np.column_stack([np.cos(ph), np.sin(ph)])
            a, *_ = np.linalg.lstsq(M, sh, rcond=None); nm.append(np.hypot(*a))
        noise.append(np.max(nm))
    return amps, np.percentile(noise, 95)
print(f"{len(S)} positron events (60 + 100 GeV, peak region, 1-2 hodoscope fibres per axis); response normalised to each run's median")
for name, coord in (("hodoscope X", hx), ("hodoscope Y", hy)):
    for lab, res in (("S", rS), ("C", rC)):
        amps, n95 = scan(coord, res)
        top = np.argsort(amps)[::-1]
        best = []
        for i in top:
            if all(abs(periods[i] - periods[j]) > 0.4 for j in best): best.append(i)
            if len(best) == 3: break
        print(f"   {lab} vs {name}: largest sinusoid amplitudes " + ", ".join(f"{amps[i] * 100:.2f}% at {periods[i]:.2f} mm" for i in best)
              + f"   (95% noise level for the largest of a scan: {n95 * 100:.2f}%)")
np.savez("period_scan_data.npz", hx=hx, hy=hy, rS=rS, rC=rC)
