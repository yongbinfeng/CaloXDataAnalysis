"""Z-scan analysis of one dataset (datasets.py): arrival time and light vs table X per DRS channel, sums over channels, and the FERS
whole-module sums where FERS was recorded.

DRS, per channel and run (events with an MCP pulse and a valid time reference, common.load_drs): in-time hits (common.intime): light per
hit, hit fraction, light per event, arrival time (median of the in-time hits); the all-event mean of all positrons (no MCP needed; used
for channels hit in more than half of the events); fractions of events above 1000 / 1500 / 2000 ADC peak (saturation).
Fits, per run set of the dataset:  ln(signal) = a_group + X / lambda (+ run offsets: per DRS chip or common to all channels, see
common.fit_channels);  time = t_group + slope X  (light speed = 1 / |slope|).
Sums: per event over the unsaturated channels of a category ([3, 6 mm SiPM] x fibre type), in-time integrals (events without a hit count
as zero), randoms subtracted; also in photoelectrons (3 mm channels, single p.e. from the out-of-time pulses); common run offsets of all
sums. Variants: randoms from the sideband before the pulse only; 6 mm hit threshold 15 / 20 / 30 ADC instead of 25.
Lateral drift of the shower: 1/lambda_ch = 1/lambda_type + ky gy_ch (common.drift_fit).
FERS (calorimeter FERS recorded): calibrated value per channel, mean over events = light per event, all positrons and only those with
an MCP pulse; per-channel fits and whole-module type sums; position matching (common.matching).
Test fibres (dataset testfibres): each fibre on its own, testfibres.py.
Usage (repository root): python3 studies/zscan/analysis.py DATASET     Output: $ZSCAN_DATA/zscan_<DATASET>.json"""
import json, sys
import numpy as np
from datasets import DATASETS
from common import (Dataset, DATA, TYPES, TS_NS, SIDEBANDS, A_HIT, MIN_SIGNAL, AMP_LEVELS, SAT_FRAC, MAX_RESID, HF_ALL, MIN_FERS,
                    MCP_SELECT, load_drs, intime, spe_of, gradients, fit_channels, fit_sum, fit_time, detrend, wmean, drift_fit, matching, load_fers)

THRESHOLDS = {"base": A_HIT, "thr6_15": {False: 40.0, True: 15.0}, "thr6_20": {False: 40.0, True: 20.0}, "thr6_30": {False: 40.0, True: 30.0}}
SUM_VARIANTS = {"base": ("base", SIDEBANDS), "sbpre": ("base", SIDEBANDS[:1]), "thr6_15": ("thr6_15", SIDEBANDS),
                "thr6_20": ("thr6_20", SIDEBANDS), "thr6_30": ("thr6_30", SIDEBANDS)}
EST = ("lph", "lpe", "all")
CNAME = {"Sci": "scintillating", "CerQuartz": "quartz", "CerPlastic": "plastic"}


def run_drs(ds, run, mode):
    """per-channel quantities of one run, and what the sums need: per threshold variant the mean, covariance and hodoscope regression of
    the per-event in-time integrals of all channels; per sum variant the randoms per event"""
    L = load_drs(ds, run, mode)
    if L is None: return None
    chans, E, A, T, N = L["chans"], L["E"], L["A"], L["T"], L["n"]; nch = len(chans)
    X = {v: np.zeros((N, nch)) for v in THRESHOLDS}; rand = {sv: np.zeros(nch) for sv in SUM_VARIANTS}
    out, spe = {}, {}
    for i, ch in enumerate(chans):
        is6 = ds.meta[ch]["is6mm"]; win = ds.w_in(ch); e, t, a = E[:, i], T[:, i], A[:, i]
        o = {f"amp{int(l)}": float((a > l).mean()) for l in AMP_LEVELS}
        o.update(all=float(L["Eall"][:, i].mean()), all_err=float(L["Eall"][:, i].std() / np.sqrt(len(L["Eall"]))), t0=np.nan)
        out[ch] = o
        b = intime(e, t, a, win, A_HIT[is6])
        if b is None: continue
        inn = b["inn"]; n_hit = inn.sum() - b["n_rand"]
        if n_hit < 50: continue                                            # no usable in-time signal: not fitted, not summed
        s = e[inn].sum() - b["n_rand"] * b["e_side"]
        o.update(t0=b["t0"], lph=s / n_hit, lph_err=float(e[inn].std() / np.sqrt(n_hit)), hf=n_hit / N, lpe=s / N,
                 lpe_err=float(np.where(inn, e, 0).std() / np.sqrt(N)),
                 time=float(np.median(t[inn])) * TS_NS, time_err=1.2533 * float(np.std(t[inn])) * TS_NS / np.sqrt(inn.sum()))
        oot = (np.abs(t) < 3000) & (np.abs(t - b["t0"]) > 100)               # out-of-time pulses: single-p.e. calibration
        spe[ch] = (a[oot], e[oot])
        for v, thr in THRESHOLDS.items():
            for sv, (tv, sb) in SUM_VARIANTS.items():
                if tv != v: continue
                r = intime(e, t, a, win, thr[is6], sb)
                rand[sv][i] = r["n_rand"] * r["e_side"] / N
                if sb == SIDEBANDS: X[v][:, i] = np.where(r["inn"], e, 0.0)
    C = gradients(X["base"], L["hx"], L["hy"])
    with np.errstate(divide="ignore", invalid="ignore"):
        for i, ch in enumerate(chans): out[ch]["gx"], out[ch]["gy"] = float(C[1, i] / C[0, i]), float(C[2, i] / C[0, i])
    stats = {v: dict(mean=X[v].mean(0), cov=np.cov(X[v], rowvar=False)) for v in THRESHOLDS}
    return dict(chans=chans, n=N, n_pos=L["n_pos"], ch=out, spe=spe, stats=stats, rand=rand, grad=C)


def series(R, runs, ch, key):
    return np.array([R[r]["ch"][ch].get(key, np.nan) for r in runs], float)


def channel_fits(ds, runs, xs, groups, sig, err, good):
    """per run set: plain fit and fit with run offsets (per chip or global, ds.cfg["offsets"]); {rs: {kind: {ch: fit}}}, offsets"""
    fits, offs = {}, {}
    for rs, rr in ds.runsets.items():
        idx = [runs.index(r) for r in rr if r in runs]
        s = {ch: sig[ch][idx] for ch in good}; e = {ch: err[ch][idx] for ch in good}; x_, g_ = xs[idx], [groups[i] for i in idx]
        plain, _ = fit_channels(x_, g_, s, e, good, n_iter=1, per_chip=False)
        withoff, dl = fit_channels(x_, g_, s, e, good, n_iter=4, per_chip=ds.cfg["offsets"] == "chip")
        fits[rs] = dict(plain=plain, offsets=withoff); offs[rs] = {k: v.tolist() for k, v in dl.items()}
    return fits, offs


def summarise(per, ests, used, runsets):
    """weighted means of the per-channel 1/lambda per fibre type: {est: {rs: {kind: {type: ...}}}}"""
    S = {}
    for est in ests:
        for rs in runsets:
            for kind in ("plain", "offsets"):
                for t in TYPES:
                    ps = [p for p in per.values() if p["type"] == t and est in p["fit"] and used(p, est, rs, kind)]
                    v = [p["fit"][est][rs][kind][0] for p in ps]; m, e, chi2 = wmean(v, [p["fit"][est][rs][kind][1] for p in ps])
                    S.setdefault(est, {}).setdefault(rs, {}).setdefault(kind, {})[t] = dict(n=len(ps), inv_lambda=m, err=e, chi2ndf=chi2,
                                                                                           median=float(np.median(v)) if v else np.nan)
    return S


def analyse_drs(ds, mode):
    R = {}
    for r in ds.runs:
        o = run_drs(ds, r, mode)
        if o is None: print(f"   run {r}: too few events with an MCP pulse, left out"); continue
        R[r] = o
    runs = [r for r in ds.runs if r in R]; xs = np.array([ds.x[r] for r in runs]); groups = [ds.group[r] for r in runs]
    chans = R[runs[0]]["chans"]; assert all(R[r]["chans"] == chans for r in runs)
    S = {k: {ch: series(R, runs, ch, k) for ch in chans} for k in ("lph", "lph_err", "hf", "lpe", "lpe_err", "all", "all_err", "t0",
                                                                       "time", "time_err", "gx", "gy", "amp1000", "amp1500", "amp2000")}
    fin = lambda a: np.all(np.isfinite(a))
    good = {"lph": [ch for ch in chans if fin(S["lph"][ch]) and np.all(S["lph"][ch] > MIN_SIGNAL) and np.all(S["lph_err"][ch] > 0)
                    and np.ptp(S["t0"][ch]) < 40]}                     # the in-time peak found consistently (no jumps between runs)
    good["lpe"] = [ch for ch in good["lph"] if np.all(S["lpe"][ch] > 0) and np.all(S["lpe_err"][ch] > 0)]
    good["all"] = [ch for ch in chans if fin(S["all"][ch]) and np.all(S["all"][ch] > 0) and fin(S["hf"][ch]) and np.all(S["hf"][ch] > HF_ALL)]
    F, OFF = {}, {}
    for est in EST:
        F[est], OFF[est] = channel_fits(ds, runs, xs, groups, S[est], S[est + "_err"], good[est])
    main_rs, main_fit = "main", ds.cfg["main_fit"]
    per = {}
    for ch in good["lph"]:
        ts, ets, tr = fit_time(xs, groups, S["time"][ch], S["time_err"][ch])
        key = (lambda c: c.rsplit("_", 1)[0]) if ds.cfg["offsets"] == "chip" else (lambda c: "all")
        per[ch] = dict(ds.meta[ch], high_amp={str(int(l)): float(np.nanmax(S[f"amp{int(l)}"][ch])) for l in AMP_LEVELS},
                       gx=float(np.nanmean(S["gx"][ch])), gy=float(np.nanmean(S["gy"][ch])),
                       **{k: S[k][ch].tolist() for k in ("lph", "lph_err", "hf", "lpe", "lpe_err", "all", "all_err", "t0", "time", "time_err")},
                       time_slope=ts, time_slope_err=ets, time_resid=tr,
                       fit={est: {rs: {kind: F[est][rs][kind][ch] for kind in ("plain", "offsets")} for rs in ds.runsets}
                            for est in EST if ch in good[est]},
                       # run offsets of the main fit (main run set), for display
                       offsets={est: OFF[est]["main"].get(key(ch), [0.0] * len(ds.runsets["main"])) for est in EST if ch in good[est]})
    unsat = lambda p: p["high_amp"]["1500"] < SAT_FRAC
    used = lambda p, est, rs, kind: unsat(p) and p["fit"][est][rs][kind][2] < MAX_RESID
    summary = summarise(per, EST, used, ds.runsets)
    drift = {}
    for est in EST:
        for rs in ds.runsets:
            chs = [ch for ch, p in per.items() if est in p["fit"] and used(p, est, rs, main_fit) and np.isfinite(p["gy"])]
            drift.setdefault(est, {})[rs] = drift_fit(chs, [per[c]["type"] for c in chs], [per[c]["fit"][est][rs][main_fit][0] for c in chs],
                                                      [per[c]["fit"][est][rs][main_fit][1] for c in chs], [per[c]["gy"] for c in chs])
    sums = sums_drs(ds, R, runs, xs, groups, chans, S, drift)
    # all-event means of the bright channels (no MCP needed), summed per fibre type
    sums_all = {}
    for rs, rr in ds.runsets.items():
        idx = [runs.index(r) for r in rr if r in runs]
        for t in TYPES:
            cs = [ch for ch in good["all"] if ds.meta[ch]["type"] == t and ch in per and unsat(per[ch])]
            if not cs: continue
            v = np.sum([S["all"][c][idx] for c in cs], 0); e = np.sqrt(np.sum([S["all_err"][c][idx] ** 2 for c in cs], 0))
            sums_all.setdefault(rs, {})[t] = dict(n=len(cs), values=v.tolist(), errs=e.tolist(), fit=fit_sum(xs[idx], [groups[i] for i in idx], v, e))
    # position check: arrival-time residual per run (median over well-timed channels) -> X offset with the median slope
    tch = [p for p in per.values() if np.isfinite(p["time_slope"]) and p["time_resid"] < 0.1]
    pos = {}
    if tch:
        rr_ = [np.array(p["time"]) - np.polyval(np.polyfit(xs, p["time"], 1), xs) for p in tch if fin(p["time"])]
        slope = float(np.median([p["time_slope"] for p in tch]))
        pos = dict(n=len(rr_), slope_ns_per_m=slope, resid_ns=np.median(rr_, axis=0).tolist(), dx_mm=(1000 * np.median(rr_, axis=0) / slope).tolist())
    return dict(runs=runs, x_m=xs.tolist(), groups=groups, n_sel=[R[r]["n"] for r in runs], n_positrons=[R[r]["n_pos"] for r in runs],
                channels=per, summary=summary, drift=drift, sums=sums, sums_all=sums_all, position_check=pos,
                n_saturated={t: int(sum(1 for p in per.values() if p["type"] == t and not unsat(p))) for t in TYPES},
                offsets={est: OFF[est] for est in EST})


def sums_drs(ds, R, runs, xs, groups, chans, S, drift):
    """per-event sums per category from the per-run means, covariances and randoms; fits per run set and variant, common run offsets"""
    J0 = [j for j, ch in enumerate(chans) if np.nanmax(S["amp1500"][ch]) < SAT_FRAC and np.all(np.isfinite(S["t0"][ch]))]
    spe = {}
    for j in J0:
        ch = chans[j]; is6 = ds.meta[ch]["is6mm"]
        a = np.concatenate([R[r]["spe"][ch][0] for r in runs if ch in R[r]["spe"]]); e = np.concatenate([R[r]["spe"][ch][1] for r in runs if ch in R[r]["spe"]])
        spe[ch] = spe_of(a, e, is6) if len(a) else (np.nan, np.nan)
    cats = {}
    for t in TYPES:
        for s6 in (False, True):
            J = [j for j in J0 if ds.meta[chans[j]]["type"] == t and ds.meta[chans[j]]["is6mm"] == s6]
            if not J: continue
            Jp = [j for j in J if not s6 and np.isfinite(spe[chans[j]][1]) and spe[chans[j]][1] > 0]
            q = np.array([spe[chans[j]][1] for j in Jp]); a1 = [spe[chans[j]][0] for j in Jp]
            cats[f"{CNAME[t]} {'6' if s6 else '3'} mm"] = dict(type=t, is6mm=s6, J=J, Jp=Jp, q=q, channels=[chans[j] for j in J],
                                                             channels_pe=[chans[j] for j in Jp], spe_amp=float(np.median(a1)) if a1 else np.nan,
                                                             spe_int=float(np.median(q)) if len(q) else np.nan)
    out = {}
    for rs, rr in ds.runsets.items():
        idx = [runs.index(r) for r in rr if r in runs]; x_, g_ = xs[idx], [groups[i] for i in idx]
        G = sorted(set(g_)); gi = np.array([G.index(g) for g in g_])
        res = {}
        for name, c in cats.items():
            J, Jp, q = c["J"], c["Jp"], c["q"]; vals = {}
            for sv, (tv, _) in SUM_VARIANTS.items():
                v = {k: [] for k in ("adc", "adc_err", "pe", "pe_err")}
                for r in [runs[i] for i in idx]:
                    st, rd, N = R[r]["stats"][tv], R[r]["rand"][sv], R[r]["n"]
                    v["adc"].append(float(st["mean"][J].sum() - rd[J].sum())); v["adc_err"].append(float(np.sqrt(st["cov"][np.ix_(J, J)].sum() / N)))
                    if Jp:
                        v["pe"].append(float((st["mean"][Jp] / q).sum() - (rd[Jp] / q).sum()))
                        v["pe_err"].append(float(np.sqrt((st["cov"][np.ix_(Jp, Jp)] / np.outer(q, q)).sum() / N)))
                vals[sv] = v
            with np.errstate(divide="ignore", invalid="ignore"):
                gy = float(np.median([R[runs[i]]["grad"][2, J].sum() / R[runs[i]]["grad"][0, J].sum() for i in idx]))
            fits = {sv: {u: fit_sum(x_, g_, np.array(vals[sv][u]), np.array(vals[sv][u + "_err"])) for u in ("adc", "pe")
                         if vals[sv][u] and np.all(np.array(vals[sv][u]) > 0)} for sv in SUM_VARIANTS}
            res[name] = dict({k: c[k] for k in ("type", "is6mm", "channels", "channels_pe", "spe_amp", "spe_int")}, values=vals, fits=fits, gy=gy)
        # common run offsets: median residual of all ADC sums (base variant) per run, detrended within each group; every sum refitted
        rsd = np.array([c["fits"]["base"]["adc"]["resid"] for c in res.values() if "adc" in c["fits"]["base"]])
        delta = detrend(np.median(rsd, axis=0), x_, gi) if len(rsd) else np.zeros(len(idx))
        dr = drift["lph"].get(rs); ky, eky = (dr["ky"], dr["eky"]) if dr else (np.nan, np.nan)
        for c in res.values():
            for sv, fs in c["fits"].items():
                for u, f in fs.items():
                    y, ey = np.array(c["values"][sv][u]) * np.exp(-delta), np.array(c["values"][sv][u + "_err"]) * np.exp(-delta)
                    fc = fit_sum(x_, g_, y, ey)
                    f.update(common=dict(inv_lambda=fc["inv_lambda"], err=fc["err"], chi2ndf=fc["chi2ndf"], resid_rms=fc["resid_rms"]),
                             drift_corr=dict(inv_lambda=f["inv_lambda"] - ky * c["gy"], err=float(np.hypot(f["err"], eky * c["gy"])),
                                             common=fc["inv_lambda"] - ky * c["gy"]))
        out[rs] = dict(categories=res, common_offsets=delta.tolist(), ky=ky)
    return out


def analyse_fers(ds, mode, with_mcp):
    L = {r: load_fers(ds, r, mode, with_mcp) for r in ds.runs}
    runs = ds.runs; xs = np.array([ds.x[r] for r in runs]); groups = [ds.group[r] for r in runs]
    names, meta = L[runs[0]]["names"], L[runs[0]]["meta"]
    sig = {n: np.array([L[r]["mean"][i] for r in runs]) for i, n in enumerate(names)}
    err = {n: np.array([max(L[r]["err"][i], 1e-6) for r in runs]) for i, n in enumerate(names)}
    chans = [n for n in names if np.all(sig[n] > MIN_FERS)]
    fits, offs = {}, {}
    for rs, rr in ds.runsets.items():
        idx = [runs.index(r) for r in rr]
        s = {n: sig[n][idx] for n in chans}; e = {n: err[n][idx] for n in chans}; x_, g_ = xs[idx], [groups[i] for i in idx]
        plain, _ = fit_channels(x_, g_, s, e, chans, n_iter=1, per_chip=False)
        withoff, dl = fit_channels(x_, g_, s, e, chans, n_iter=4, per_chip=False)
        fits[rs] = dict(plain=plain, offsets=withoff); offs[rs] = dl.get("all", np.zeros(len(idx))).tolist()
    per = {}
    for i, n in enumerate(names):
        if n not in chans: continue
        per[n] = dict(meta[n], signal=sig[n].tolist(), err=err[n].tolist(), gx=float(np.nanmean([L[r]["gx"][i] for r in runs])),
                      gy=float(np.nanmean([L[r]["gy"][i] for r in runs])), lg_max=float(max(L[r]["lg_max"][i] for r in runs)),
                      fit={"lpe": {rs: {k: fits[rs][k][n] for k in ("plain", "offsets")} for rs in ds.runsets}})
    used = lambda p, est, rs, kind: p["fit"][est][rs][kind][2] < MAX_RESID
    sums, drift = {}, {}
    for rs, rr in ds.runsets.items():
        idx = [runs.index(r) for r in rr]
        for t in TYPES:
            ks = [i for i, n in enumerate(names) if meta[n]["type"] == t]
            v = np.array([L[runs[j]]["mean"][ks].sum() for j in idx]); e = np.array([np.sqrt((L[runs[j]]["err"][ks] ** 2).sum()) for j in idx])
            sums.setdefault(rs, {})[t] = dict(values=v.tolist(), errs=e.tolist(), fit=fit_sum(xs[idx], [groups[j] for j in idx], v, e))
        chs = [n for n, p in per.items() if used(p, "lpe", rs, "plain") and np.isfinite(p["gy"])]
        drift[rs] = drift_fit(chs, [per[c]["type"] for c in chs], [per[c]["fit"]["lpe"][rs]["plain"][0] for c in chs],
                              [per[c]["fit"]["lpe"][rs]["plain"][1] for c in chs], [per[c]["gy"] for c in chs])
    return dict(n_sel=[L[r]["n"] for r in runs], centroid=[L[r]["cy"] for r in runs], channels=per, offsets=offs,
                summary=summarise(per, ("lpe",), used, ds.runsets), sums=sums, drift=drift, lg_max=float(max(p["lg_max"] for p in per.values())))


def report(ds, res):
    main = res["selections"][res["main_selection"]]; D = main["drs"]
    print(f"{ds.name}: {len(D['runs'])} runs, {len(D['channels'])} DRS channels, saturated {D['n_saturated']}")
    for t in TYPES:
        v = [p["time_slope"] for p in D["channels"].values() if p["type"] == t and np.isfinite(p["time_slope"]) and p["time_resid"] < 0.1]
        if v: print(f"   speed {t:10s}: dt/dX median {np.median(v):+.2f} ns/m ({len(v)} ch.) -> {1 / abs(np.median(v)) / 0.2998:.3f} c")
    for name, c in D["sums"]["main"]["categories"].items():
        f = c["fits"]["base"].get("adc")
        if f: print(f"   DRS sum {name:20s} {len(c['channels']):3d} ch: 1/lambda {f['inv_lambda']:+.3f} +- {f['err']:.3f}, "
                    f"common offsets {f['common']['inv_lambda']:+.3f} +- {f['common']['err']:.3f} (chi2/ndf {f['common']['chi2ndf']:.1f})")
    for est in ("lph", "lpe"):
        s = D["summary"][est]["main"][ds.cfg["main_fit"]]
        print(f"   per channel {est}: " + ", ".join(f"{t} {s[t]['inv_lambda']:+.3f} +- {s[t]['err']:.3f} ({s[t]['n']})" for t in TYPES))
    if "fers" in main:
        print("   FERS whole-module sums: " + ", ".join(f"{t} {main['fers']['sums']['main'][t]['fit']['inv_lambda']:+.3f} +- "
                                                    f"{main['fers']['sums']['main'][t]['fit']['err']:.3f}" for t in TYPES))


if __name__ == "__main__":
    if DATASETS[sys.argv[1]].get("channel_map"):
        import testfibres; testfibres.main_analysis(); sys.exit()
    ds = Dataset(sys.argv[1])
    res = dict(dataset=ds.name, title=ds.cfg["title"], desc=ds.cfg["desc"], runs=ds.runs, x_m=[ds.x[r] for r in ds.runs],
               groups=[ds.group[r] for r in ds.runs], runsets=ds.runsets, main_fit=ds.cfg["main_fit"], offsets=ds.cfg["offsets"],
               mcp_select=MCP_SELECT, fers=ds.cfg["fers"], selections={})
    modes = ["spot"]
    if ds.cfg["fers"]:
        res["matching"] = matching(ds); modes = ["matched", "spot"]
        M = res["matching"]
        print(f"shower drift: {M['centroid_slope_rows_per_m']:+.3f} tower rows/m = {M['hodo_cm_per_m']:+.2f} hodoscope cm/m")
    res["main_selection"] = modes[0]
    for mode in modes:
        sel = dict(drs=analyse_drs(ds, mode))
        if ds.cfg["fers"]:
            sel["fers"] = analyse_fers(ds, mode, False); sel["fers_mcp"] = analyse_fers(ds, mode, True)
        res["selections"][mode] = sel
    report(ds, res)
    json.dump(res, open(f"{DATA}/zscan_{ds.name}.json", "w"), default=float)
