"""Positron energy resolution of the 2025 FERS (runs of config.RUNS), and its parts: photostatistics, electronic noise, position.
Per event: S = sum of the calibrated scintillating channels, C = sum of the Cherenkov channels (quartz + plastic), each calibrated to the
beam energy, and (S+C)/2. Events: positron PID, 1-3 hodoscope fibres per axis; two selections:
  spot  impact point within BEAM_SPOT of the run median (hodoscope), both axes
  all   the whole beam spot
Sums: over the whole module, and over the towers within REGION towers of the shower centre of the run (fewer channels, less noise).
Position corrections: the response depends on where the shower sits in the module (fibre pattern, quartz / plastic towers); each signal
is corrected with a quadratic surface in (x, y), fitted with outliers clipped on the even events and applied to the odd ones, with
(x, y) = the centre of gravity of the shower (energy-weighted S of the towers within 2.5 towers of the hottest one) or the impact
point (hodoscope).
Width: gaussian core (iterative 2-sigma clipping, x 1.137); resolution sigma/mean.
Photostatistics: sigma/E = 1/sqrt(N_pe/GeV x E) with the p.e. yields of yields.json (20-100 GeV mean).
Noise: width of the calibrated sums in the pedestal run PED_RUN (pedestal and dark counts of the summed channels), in GeV; the channels
are correlated (the width of the sum is larger than the sum of the channel variances).
S - C: for positrons both measure E; their difference removes the fluctuations common to both (leakage, containment).
Fits over FIT_EMIN-100 GeV: (sigma/E)^2 = a^2/E + b^2 + (c/E)^2 with c the measured noise.
Usage (repository root): python3 studies/pe_yield/resolution.py      Output: $PE_YIELD_DATA/resolution.json"""
import json, sys
sys.path.insert(0, ".")
import numpy as np
from configs.selection_config import get_particle_selection
from config import DATA, RUNS, RESPONSE, PED_RUN, MAP_RUN

BEAM_SPOT = 1.0               # [cm] "spot" selection: impact point within this of the run median, both axes
FIT_EMIN = 20                 # [GeV] lowest energy in the fits (10 GeV: offsets of the channels far from the shower)
SIGNALS = ("S", "C", "SC")    # SC = (S + C) / 2
HBINS = np.linspace(0.6, 1.4, 81)
PROFILE_BINS = np.arange(-1.5, 1.51, 0.25)   # [towers] shower centre for the response profiles
REGION = 4.5                  # [towers] radius of the region sum around the shower centre of the run
COG_NEAR = 2.5                # [towers] towers around the hottest one in the centre of gravity
CORRECTIONS = ("cog", "impact")


def core(x, nit=8):
    """gaussian core mean, width and error of the width"""
    m, s = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
    for _ in range(nit):
        k = np.abs(x - m) < 2 * s; m, s = x[k].mean(), x[k].std() * 1.137
    return float(m), float(s), float(s / np.sqrt(2 * k.sum()))


def quad_surface(x, y, z, nit=5):
    """z = c0 + c1 x + c2 y + c3 x^2 + c4 y^2 + c5 x y, least squares with the points beyond 2.5 sigma of the residuals clipped
    (the shower tails); returns a function of (x, y)"""
    A = lambda x, y: np.column_stack([np.ones_like(x), x, y, x * x, y * y, x * y])
    k = np.ones(len(z), bool)
    for _ in range(nit):
        c = np.linalg.lstsq(A(x[k], y[k]), z[k], rcond=None)[0]; r = z - A(x, y) @ c
        s = 1.4826 * np.median(np.abs(r[k])); k = np.abs(r) < 2.5 * s
    return lambda x, y: A(x, y) @ c


def events(run, E):
    d = np.load(f"{DATA}/fers/fers_{run}.npz"); meta = d["meta"].astype(int); v = d["vals"]
    sel = np.ones(len(v), bool)
    for det, fire in get_particle_selection("electron", run).items():
        sel &= d[f"is_{det}_fired" if fire else f"is_{det}_vetoed"].astype(bool)
    nx, ny, hx, hy = d["TTU_Hodo_nHitX"], d["TTU_Hodo_nHitY"], d["TTU_Hodo_X"], d["TTU_Hodo_Y"]
    sel &= (nx >= 1) & (nx <= 3) & (ny >= 1) & (ny <= 3)
    sci = meta[:, 2] == 0; cer = meta[:, 2] == 1; v = np.nan_to_num(v[sel]); ix, iy = meta[:, 4].astype(float), meta[:, 5].astype(float)
    w = np.clip(v[:, sci], 0, None); imax = np.argmax(w, 1)
    near = np.hypot(ix[sci][None, :] - ix[sci][imax][:, None], iy[sci][None, :] - iy[sci][imax][:, None]) <= COG_NEAR
    wn = w * near; cx = (wn * ix[sci]).sum(1) / wn.sum(1); cy = (wn * iy[sci]).sum(1) / wn.sum(1)
    x0, y0 = np.median(cx), np.median(cy); region = np.hypot(ix - x0, iy - y0) <= REGION
    out = dict(S=v[:, sci].sum(1) / E, C=v[:, cer].sum(1) / E, S_reg=v[:, sci & region].sum(1) / E, C_reg=v[:, cer & region].sum(1) / E,
               cog_x=cx - x0, cog_y=cy - y0, impact_x=hx[sel] - np.median(hx[sel]), impact_y=hy[sel] - np.median(hy[sel]))
    out["spot"] = (np.abs(out["impact_x"]) < BEAM_SPOT) & (np.abs(out["impact_y"]) < BEAM_SPOT)
    ok = (out["S"] > 0.5) & (out["C"] > 0.3)       # removes the few events with a partial shower
    return {k: x[ok] for k, x in out.items()}, (meta, region)


def resolutions(ev, reg=""):
    """raw and position-corrected (centre of gravity, impact point) resolutions of S, C, (S+C)/2 of the module (reg="") or region
    (reg="_reg") sum; S - C and the S-C correlation after the impact-point correction"""
    S, C = ev["S" + reg], ev["C" + reg]; sig = {"S": S, "C": C, "SC": (S + C) / 2}
    even = np.arange(len(S)) % 2 == 0; odd = ~even; out = dict(n=int(len(S)), signals={}); corr = {}
    for k, x in sig.items():
        m0, s0, e0 = core(x[odd]); rec = dict(mean=m0, res=s0 / m0, err=e0 / m0, hist=np.histogram(x[odd] / m0, HBINS)[0].tolist())
        for c in CORRECTIONS:
            a, b = ev[c + "_x"], ev[c + "_y"]
            f = quad_surface(a[even], b[even], x[even] / np.median(x[even])); xc = x[odd] / f(a[odd], b[odd]); m1, s1, e1 = core(xc)
            rec.update({f"res_{c}": s1 / m1, f"err_{c}": e1 / m1, f"hist_{c}": np.histogram(xc / m1, HBINS)[0].tolist()})
            if c == "impact": corr[k] = xc / m1
        out["signals"][k] = rec
    out["s_minus_c"] = float(core(corr["S"] - corr["C"])[1]); out["corr_sc"] = float(np.corrcoef(corr["S"], corr["C"])[0, 1])
    return out


def profiles(ev):
    """response of S and C (normalised to the median) vs the shower centre in x and y"""
    out = {}
    for k in ("S", "C"):
        x = ev[k] / np.median(ev[k])
        for ax in ("impact_x", "impact_y", "cog_x", "cog_y"):
            b = np.digitize(ev[ax], PROFILE_BINS) - 1
            out[f"{k}_{ax}"] = [float(np.median(x[b == i])) if np.sum(b == i) > 30 else None for i in range(len(PROFILE_BINS) - 1)]
    return out


def noise(region_of):
    """width of the calibrated S and C sums in the pedestal run [GeV], whole module and the region of each run (region_of: (board,
    channel) -> in region)"""
    P = np.load(f"{DATA}/raw/raw_{PED_RUN}.npz"); meta = P["meta"].astype(int); hg = P["hg"].astype(float)
    R = json.load(open(RESPONSE))
    resp = np.array([R.get(f"FERS_Board{b}_energyMix_{c}", {}).get("response", np.nan) for b, c in meta[:, :2]])
    x = (hg - hg.mean(0)) * np.where(np.isfinite(resp) & (resp > 0), resp, 0)
    reg = np.array([region_of.get((b, c), False) for b, c in meta[:, :2]])
    out = {}
    for k, m in (("S", meta[:, 2] == 0), ("C", meta[:, 2] == 1)):
        out[k] = float(np.std(x[:, m].sum(1))); out[k + "_reg"] = float(np.std(x[:, m & reg].sum(1)))
        out[k + "_uncorrelated"] = float(np.sqrt(np.var(x[:, m], 0).sum()))
    for t in ("", "_reg"): out["SC" + t] = 0.5 * np.sqrt(out["S" + t] ** 2 + out["C" + t] ** 2)
    return out


def fit_res(E, r, e, c=0.0):
    """(sigma/E)^2 = a^2/E + b^2 + (c/E)^2, c fixed; weighted least squares in (a^2, b^2); errors scaled by sqrt(chi2/ndf) where above 1"""
    E, r, e = map(np.asarray, (E, r, e)); y = r ** 2 - (c / E) ** 2; w = 1 / (2 * r * e) ** 2; A = np.column_stack([1 / E, np.ones_like(E)])
    C = np.linalg.inv(A.T @ (A * w[:, None])); p = C @ (A.T @ (w * y)); chi2 = float(np.sum(w * (y - A @ p) ** 2) / max(len(E) - 2, 1))
    C = C * max(chi2, 1.0)
    a = np.sqrt(max(p[0], 0)); ea = np.sqrt(C[0, 0]) / (2 * max(a, 1e-6))
    return dict(a=float(a), a_err=float(ea), b2=float(p[1]), b2_err=float(np.sqrt(C[1, 1])), b=float(np.sqrt(max(p[1], 0))),
                b_upper=float(np.sqrt(max(p[1] + 1.645 * np.sqrt(C[1, 1]), 0))), c=float(c), chi2ndf=chi2)


def main():
    Y = json.load(open(f"{DATA}/yields.json"))
    pe = {t: float(np.mean([y["types"][t]["pe_per_gev"] for y in Y["runs"] if y["energy"] >= FIT_EMIN])) for t in ("Sci", "Cer")}
    res = dict(beam_spot_cm=BEAM_SPOT, region_towers=REGION, fit_emin=FIT_EMIN, hbins=HBINS.tolist(), profile_bins=PROFILE_BINS.tolist(), pe_per_gev=pe,
               photostat={"S": 1 / np.sqrt(pe["Sci"]), "C": 1 / np.sqrt(pe["Cer"]), "SC": 0.5 * np.sqrt(1 / pe["Sci"] + 1 / pe["Cer"])}, runs=[], fits={})
    region_of = None
    for run, E in RUNS.items():
        ev, (meta, region) = events(run, E); r = dict(run=run, energy=E)
        if run == MAP_RUN: region_of = {(b, c): bool(g) for (b, c), g in zip(meta[:, :2], region)}; r["profiles"] = profiles(ev)
        for selname, m in (("all", np.ones(len(ev["S"]), bool)), ("spot", ev["spot"])):
            sub = {k: x[m] for k, x in ev.items()}
            r[selname] = resolutions(sub); r[selname + "_reg"] = resolutions(sub, "_reg")
        res["runs"].append(r)
    nz = noise(region_of); res["noise_gev"] = nz
    use = [r for r in res["runs"] if r["energy"] >= FIT_EMIN]; En = [r["energy"] for r in use]
    for selname in ("spot", "all", "spot_reg", "all_reg"):
        reg = "_reg" if selname.endswith("_reg") else ""
        for k in SIGNALS:
            for tag in ("",) + tuple("_" + c for c in CORRECTIONS):
                res["fits"][f"{selname}_{k}{tag}"] = fit_res(En, [r[selname]["signals"][k]["res" + tag] for r in use],
                                                             [r[selname]["signals"][k]["err" + tag] for r in use], nz[k + reg])
        res["fits"][f"{selname}_s_minus_c"] = fit_res(En, [r[selname]["s_minus_c"] for r in use], [r[selname]["s_minus_c"] / np.sqrt(r[selname]["n"]) for r in use])
    json.dump(res, open(f"{DATA}/resolution.json", "w"))
    print(f"p.e./GeV S {pe['Sci']:.1f}, C {pe['Cer']:.2f} -> photostatistics S {100 * res['photostat']['S']:.1f}%/sqrt(E), C {100 * res['photostat']['C']:.1f}%/sqrt(E)")
    print("noise [GeV]: " + ", ".join(f"{k} {v:.2f}" for k, v in nz.items()))
    for selname in ("spot", "all", "spot_reg", "all_reg"):
        print(f"== {selname}: run E events | resolution % raw / CoG / impact (S; C; (S+C)/2) | S-C | corr")
        for r in res["runs"]:
            x = r[selname]
            print(f"   {r['run']} {r['energy']:4.0f} {x['n']:6d} | " + "; ".join("/".join(f"{100 * x['signals'][k][t]:5.2f}" for t in ("res", "res_cog", "res_impact")) for k in SIGNALS)
                  + f" | {100 * x['s_minus_c']:5.2f} | {x['corr_sc']:+.2f}")
    for k, f in res["fits"].items():
        print(f"  fit {k:22s}: a = {100 * f['a']:5.1f} +- {100 * f['a_err']:.1f} %/sqrt(E), b = {100 * f['b']:4.2f} % (< {100 * f['b_upper']:.2f}), c = {f['c']:.2f} GeV (chi2/ndf {f['chi2ndf']:.1f})")


if __name__ == "__main__":
    main()
