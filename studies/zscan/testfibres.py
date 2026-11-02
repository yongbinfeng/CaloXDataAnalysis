"""Test fibres (dataset testfibres, datasets.py): light speed and attenuation of each test fibre with the z-scan method (common.py).
Events: positrons with a pulse in MCP_DS_0 above its cut (= the test-fibre mcp_clean, the same MCP requirement as the z scans) and a valid
time reference, no hodoscope requirement and no beam window ("loose": each bundle sees the shower averaged over the beam profile, and the
hodoscope 1-3 fibre cut alone removes ~3/4 of the positrons). Same in-time hits as the z scans (+-20 TS scintillating / +-6 TS others
around the arrival peak, peak above 40 ADC, randoms from the sidebands before and after the pulse subtracted), arrival time = median CFD
time of the in-time hits (relative to the MCP). Long integral for the scintillating bundles, short for the others.
Second selection, "spot": the z-scan events (in addition 1-3 hodoscope fibres per axis and beam within 1 cm of the run median; about 7x
fewer events), shown for comparison. The framework mcp_clean keeps all MCP events, without particle ID.
Test-fibre specifics (as in the framework scan, studies/testing_fiber/plot_drs_mpv_scan_testfibers.py):
  - single small bundles: each fibre on its own, no sums
  - light per crossing = light per hit x (1 - exp(-mu)), mu from the p.e. per hit (zero-truncated Poisson: p.e. per hit = mu / (1 - exp(-mu))):
    the bundles are small, so most events without a hit are events in which no shower particle crossed them; the light per event follows
    that geometry (upper limit on the attenuation), the light per hit misses the crossings with zero p.e. Single p.e. from the out-of-time
    pulses of all runs together (one calibration per fibre)
  - one normalisation per continuous stretch (40 GeV, 40 GeV after the lead block moved, 60 GeV); time: one intercept per stretch
  - light fits leave out the runs where the scintillating bundles are hit less than 90% of their plateau (partial overlap: edges of the
    holder, metal blocks of the module surface)
Usage (repository root): python3 studies/zscan/analysis.py testfibres; python3 studies/zscan/page.py testfibres  (both call this module)
Output: $ZSCAN_DATA/zscan_testfibres.json, results/html/ZScan/ZScan_testfibres.html"""
import json, os, types
import numpy as np
from common import Dataset, DATA, TS_NS, A_HIT, load_drs, intime, spe_of, fit_sum, fit_time, wmean
from datasets import PLOTS, HTML, DATASETS

NAME = "testfibres"
OVERLAP_FRACTION = 0.9            # light fits: runs where the scintillating bundles are hit less than this fraction of their plateau are left out
PLATEAU = (-0.220, 0.050)         # [m] X range of the hit-fraction plateau
MIN_SPAN = 0.100                  # [m] travel a stretch needs for its own lambda
SAT = 1500.0                      # [ADC] peak above which the DRS clips
FRAMEWORK = "results/plots/PositionScan"   # framework scan, shown for comparison
EST = {"lpc": "light per crossing", "lph": "light per hit", "lpe": "light per event"}
SELECTIONS = {"spot": "z-scan selection (hodoscope 1-3 fibres per axis, beam within 1 cm)", "loose": "loose (positrons with an MCP pulse, no hodoscope requirement)"}
MAIN_SELECTION = "loose"
KEYS = ("t0", "n_hit", "hf", "hf_err", "lph", "lph_err", "lpe", "lpe_err", "sat", "time", "time_err")


def run_fibres(ds, run, mode):
    """per fibre in one run: in-time hits, hit fraction, light per hit / per event, clipped fraction of the hits, arrival time; the
    out-of-time pulses for the single-p.e. calibration"""
    L = load_drs(ds, run, mode)
    if L is None: return None
    N = L["n"]; out, oot = {}, {}
    for i, ch in enumerate(L["chans"]):
        e, t, a = L["E"][:, i], L["T"][:, i], L["A"][:, i]
        b = intime(e, t, a, ds.w_in(ch), A_HIT[False])
        if b is None: continue
        inn = b["inn"]; n_hit = inn.sum() - b["n_rand"]
        if n_hit < 50: continue
        s = e[inn].sum() - b["n_rand"] * b["e_side"]
        out[ch] = dict(t0=b["t0"], n_hit=float(n_hit), hf=n_hit / N, hf_err=float(np.sqrt(inn.sum() + b["n_rand"]) / N),
                       lph=s / n_hit, lph_err=float(e[inn].std() / np.sqrt(n_hit)), lpe=s / N, lpe_err=float(np.where(inn, e, 0).std() / np.sqrt(N)),
                       sat=float((a[inn] > SAT).mean()),
                       time=float(np.median(t[inn])) * TS_NS, time_err=1.2533 * float(np.std(t[inn])) * TS_NS / np.sqrt(inn.sum()))
        o = (np.abs(t) < 3000) & (np.abs(t - b["t0"]) > 100)
        oot[ch] = (a[o], e[o])
    return dict(n=N, n_pos=L["n_pos"], ch=out, oot=oot)


def mu_of(m):
    """mu with mu / (1 - exp(-mu)) = m (mean of a zero-truncated Poisson); nan for m <= 1"""
    if not m > 1.001: return np.nan
    lo, hi = 1e-6, 2.0 * m + 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if mid / -np.expm1(-mid) < m else (lo, mid)
    return 0.5 * (lo + hi)


def fit_light(xs, groups, y, ey, ok):
    """ln(light) = a_stretch + X / lambda over the runs ok (stretches with at least 3 of them), and per stretch with >= MIN_SPAN of travel"""
    ok = ok & np.isfinite(y) & (y > 0) & np.isfinite(ey) & (ey > 0)
    g = np.array(groups)
    for k in set(groups):
        if (ok & (g == k)).sum() < 3: ok &= g != k
    if ok.sum() < 4: return None
    f = fit_sum(xs[ok], list(g[ok]), y[ok], ey[ok])
    f.update(groups=sorted(set(g[ok])), used=ok.tolist(), per_group={})
    for k in f["groups"]:
        m = ok & (g == k)
        if np.ptp(xs[m]) >= MIN_SPAN:
            fk = fit_sum(xs[m], [k] * int(m.sum()), y[m], ey[m])
            f["per_group"][k] = dict(inv_lambda=fk["inv_lambda"], err=fk["err"], chi2ndf=fk["chi2ndf"])
    return f


def analyse(ds, mode):
    R = {}
    for r in ds.runs:
        o = run_fibres(ds, r, mode)
        if o is None: print(f"   run {r}: too few events with an MCP pulse, left out"); continue
        R[r] = o
    runs = [r for r in ds.runs if r in R]; xs = np.array([ds.x[r] for r in runs]); groups = [ds.group[r] for r in runs]
    ser = lambda ch, k: np.array([R[r]["ch"].get(ch, {}).get(k, np.nan) for r in runs], float)
    # overlap of the shower with the bundles: hit fraction of the scintillating bundles relative to their plateau
    mon = [ch for ch, m in ds.meta.items() if m["type"] == "Sci"]
    inpl = (xs >= PLATEAU[0]) & (xs <= PLATEAU[1])
    ratio = np.nanmedian([ser(ch, "hf") / np.nanmedian(ser(ch, "hf")[inpl]) for ch in mon], axis=0)
    full = ratio >= OVERLAP_FRACTION
    fibres = {}
    for ch, m in ds.meta.items():
        S = {k: ser(ch, k) for k in KEYS}
        hit = [r for r in runs if ch in R[r]["oot"]]
        fb = dict(label=m["label"], type=m["type"], channel=ch, n_runs=len(hit), **{k: v.tolist() for k, v in S.items()})
        fibres[m["label"]] = fb
        if len(hit) < 5: continue
        a1, q = spe_of(np.concatenate([R[r]["oot"][ch][0] for r in hit]), np.concatenate([R[r]["oot"][ch][1] for r in hit]), False)
        S["pe"] = S["lph"] / q; S["mu"] = np.array([mu_of(v) for v in S["pe"]]); p = -np.expm1(-S["mu"])
        S["lpc"], S["lpc_err"] = S["lph"] * p, S["lph_err"] * p
        ts, ets, tr = fit_time(xs, groups, S["time"], S["time_err"])
        t1, et1, tr1 = fit_time(xs, ["all"] * len(xs), S["time"], S["time_err"])
        ok = np.isfinite(S["time"]); icpt = {}
        for g in sorted(set(groups)):
            mg = ok & (np.array(groups) == g)
            if mg.any(): icpt[g] = wmean(S["time"][mg] - ts * xs[mg], S["time_err"][mg])[0]
        fb.update(spe_amp=a1, spe_int=q, pe=S["pe"].tolist(), mu=S["mu"].tolist(), lpc=S["lpc"].tolist(), lpc_err=S["lpc_err"].tolist(),
                  time_slope=ts, time_slope_err=ets, time_resid=tr, time_icpt=icpt, time_slope_1line=t1, time_slope_1line_err=et1,
                  time_resid_1line=tr1, fits={est: f for est in EST if (f := fit_light(xs, groups, S[est], S[est + "_err"], full)) is not None})
    return dict(runs=runs, x_m=xs.tolist(), x_table_m=[ds.x[r] - ds.cfg["x_offset"](r) for r in runs], groups=groups, n_sel=[R[r]["n"] for r in runs],
                n_positrons=[R[r]["n_pos"] for r in runs], overlap_ratio=ratio.tolist(), full_overlap=full.tolist(),
                overlap_fraction=OVERLAP_FRACTION, fibres=fibres)


def report(res, mode):
    print(f"{NAME}, {mode}: {len(res['runs'])} runs, median {np.median(res['n_sel']):.0f} events; partial overlap (left out of the light fits) at X = "
          f"{[round(1000 * x, 1) for x, f in zip(res['x_m'], res['full_overlap']) if not f]} mm")
    for l, f in sorted(res["fibres"].items()):
        if "time_slope" not in f: print(f"   {l:16s} no usable signal ({f['n_runs']} runs with hits)"); continue
        c = f["fits"].get("lpc")
        print(f"   {l:16s} dt/dX {f['time_slope']:+.2f} +- {f['time_slope_err']:.2f} ns/m ({1 / abs(f['time_slope']) / 0.2998:.3f} c)   "
              + (f"1/lambda (per crossing) {c['inv_lambda']:+.3f} +- {c['err']:.3f} /m" if c else "no light fit"))


def main_analysis():
    ds = Dataset(NAME)
    res = dict(dataset=NAME, title=ds.cfg["title"], desc=ds.cfg["desc"], main_selection=MAIN_SELECTION, selections={})
    for mode in SELECTIONS:
        res["selections"][mode] = analyse(ds, mode); report(res["selections"][mode], mode)
    json.dump(res, open(f"{DATA}/zscan_{NAME}.json", "w"), default=float)


# ---------------- page ----------------
def page():
    import ROOT
    from plotting.my_function import DrawHistos
    from plotting.calox_plot_helper import plot_fiber_summary, fiber_base_type
    from utils.html_generator import generate_html
    JJ = json.load(open(f"{DATA}/zscan_{NAME}.json")); MS = JJ["main_selection"]; J = dict(JJ["selections"][MS], title=JJ["title"], desc=JJ["desc"])
    OUTD = f"{PLOTS}/{NAME}"; os.makedirs(OUTD, exist_ok=True); OUT = f"{HTML}/ZScan_{NAME}.html"
    pm = types.SimpleNamespace(use_jsroot=False, get_output_dir=lambda: OUTD)
    XM = np.array(J["x_m"]) * 1000; GRP = J["groups"]; GN = sorted(set(GRP)); gi = np.array([GN.index(g) for g in GRP])
    full = np.array(J["full_overlap"]); F = J["fibres"]
    labels = sorted(F, key=lambda l: (fiber_base_type(l), l)); good = [l for l in labels if "time_slope" in F[l]]
    GCOL = [ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kRed + 1]; GMK = [20, 21, 22]
    fw_pe = f"{FRAMEWORK}/DRS_PE_Scan/drs_pe_scan_channels.json"; fw_t = f"{FRAMEWORK}/DRS_CFD_MPV_Scan/drs_cfd_mpv_scan_channels.json"
    FW_PE = json.load(open(fw_pe)) if os.path.exists(fw_pe) else {}; FW_T = json.load(open(fw_t)) if os.path.exists(fw_t) else {}

    def fw_att(l):
        a = FW_PE.get(l, {}).get("attenuation", {}).get("light_per_crossing")
        return (1000 * a["slope_per_mm"], 1000 * a["slope_err_per_mm"]) if a else None

    def fw_speed(l):     # framework: cm/ns -> mm/ns
        s = FW_T.get(l)
        return (10 * s["speed_cmns"], 10 * s["e_speed_cmns"]) if s and s.get("speed_cmns") else None

    def graph(xs, ys, es=None):
        g = ROOT.TGraphErrors(len(xs))
        for i, (x, y) in enumerate(zip(xs, ys)):
            g.SetPoint(i, float(x), float(y)); g.SetPointError(i, 0, float(es[i]) if es is not None else 0)
        return g

    TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
    def table(cols, rows):
        head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
        body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
        return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
    fmt = lambda m, e: f"{m:+.3f} &plusmn; {e:.3f}" if m is not None and np.isfinite(m) else "&ndash;"
    def lam(m, e):
        if not np.isfinite(m): return "&ndash;"
        return f"{1 / m:.1f} ({1 / (m + e):.1f}&ndash;{1 / (m - e):.1f})" if m - e > 0 else (f"&gt; {1 / (m + e):.1f}" if m + e > 0 else "none seen")
    details = lambda title, body: f'<details style="margin:10px 0"><summary><b>{title}</b></summary>{body}</details>'

    # summary plots: speed and 1/lambda per fibre (x: the framework scan)
    speed = {l: (1000 / abs(F[l]["time_slope"]), 1000 * F[l]["time_slope_err"] / F[l]["time_slope"] ** 2) for l in good}
    main = [plot_fiber_summary(speed, "Light speed [mm/ns]", 120, 230, "speed_fibres", pm, "test fibres, z-scan method",
                               secondary={l: v for l in good if (v := fw_speed(l))}, note="#bullet z-scan method   #times framework scan")]
    att = {l: (F[l]["fits"]["lpc"]["inv_lambda"], F[l]["fits"]["lpc"]["err"]) for l in good if "lpc" in F[l]["fits"]}
    main.append(plot_fiber_summary(att, "1/#lambda [m^{-1}] (light per crossing)", -1.0, 2.0, "inv_lambda_fibres", pm, "test fibres, z-scan method",
                                   secondary={l: v for l in att if (v := fw_att(l))}, note="#bullet z-scan method   #times framework scan"))

    def lam_pt(a):       # 1/lambda +- err [/m] -> lambda [m] with its 1-sigma interval (open upwards when 1/lambda - err <= 0); None if 1/lambda <= 0
        m, e = a
        if not m > 0: return None
        return (1 / m, 1 / m - 1 / (m + e), 1 / (m - e) - 1 / m if m - e > 0 else 1e4)
    main.append(plot_fiber_summary({l: v for l in att if (v := lam_pt(att[l]))}, "#lambda [m] (light per crossing)", 0.3, 3000, "lambda_fibres", pm,
                                   "test fibres, z-scan method", secondary={l: v for l in att if (a := fw_att(l)) and (v := lam_pt(a))}, logy=True,
                                   note="#bullet z-scan method   #times framework scan   (1/#lambda #leq 0: not shown)"))
    mon = [l for l in good if F[l]["type"] == "Sci"]
    G = [graph(XM, J["overlap_ratio"])] + [graph(XM, np.array(F[l]["hf"]) / np.nanmedian(np.array(F[l]["hf"])[(XM >= 1000 * PLATEAU[0]) & (XM <= 1000 * PLATEAU[1])])) for l in mon]
    cut = ROOT.TLine(XM.min() - 20, OVERLAP_FRACTION, XM.max() + 20, OVERLAP_FRACTION); cut.SetLineStyle(2)
    DrawHistos(G, ["median (used)"] + mon, XM.min() - 20, XM.max() + 20, "X along the test fibres [mm] (readout #rightarrow)", 0, 1.3,
               "Hit fraction / plateau", "overlap_monitor", dology=False, drawoptions="P", legendoptions="P", markerstyles=[20, 24, 25],
               mycolors=[ROOT.kBlack, ROOT.kGreen + 2, ROOT.kBlue + 1], legendPos=[0.20, 0.18, 0.60, 0.36], extraToDraw=[cut], outdir=OUTD,
               lumitext="scintillating bundles", usePDF=False)
    main.append("overlap_monitor")

    # per fibre: arrival time, light per crossing, hit fraction vs X
    per_fibre = {}
    for l in good:
        f = F[l]; nm = l.replace(" ", "_"); per_fibre[l] = []
        t, et = np.array(f["time"], float), np.array(f["time_err"], float)
        G, L_, C, M, ex = [], [], [], [], []
        for k, g in enumerate(GN):
            m = (gi == k) & np.isfinite(t)
            if not m.any(): continue
            G.append(graph(XM[m], t[m], et[m])); L_.append(g); C.append(GCOL[k]); M.append(GMK[k])
            if g in f["time_icpt"]:
                fl = ROOT.TF1(f"ft_{nm}_{k}", f"{f['time_icpt'][g]} + {f['time_slope']}*x/1000.", XM[m].min() - 10, XM[m].max() + 10)
                fl.SetLineColor(GCOL[k]); fl.SetLineWidth(2); ex.append((fl, "SAME"))
        lo, hi = np.nanmin(t), np.nanmax(t); pad = 0.2 * (hi - lo) + 0.1
        txt = ROOT.TLatex(0.20, 0.25, f"dt/dX = {f['time_slope']:+.2f} #pm {f['time_slope_err']:.2f} ns/m, v = {1 / abs(f['time_slope']) / 0.2998:.3f} c")
        txt.SetNDC(); txt.SetTextSize(0.035); ex.append(txt)
        DrawHistos(G, L_, XM.min() - 20, XM.max() + 20, "X along the test fibres [mm]", lo - pad, hi + 2 * pad, "Arrival time - MCP [ns]", f"time_{nm}",
                   dology=False, drawoptions="P", legendoptions="P", markerstyles=M, mycolors=C, legendPos=[0.55, 0.74, 0.92, 0.89], extraToDraw=ex,
                   outdir=OUTD, lumitext=l, usePDF=False)
        per_fibre[l].append(f"time_{nm}")
        for est, ytit in (("lpc", "Light per crossing [ADC]"), ("hf", "Hit fraction")):
            y = np.array(f[est], float); ey = np.array(f[est + "_err"], float)
            G, L_, C, M, ex = [], [], [], [], []
            fit = f["fits"].get(est)
            for k, g in enumerate(GN):
                m = (gi == k) & np.isfinite(y) & full
                if not m.any(): continue
                G.append(graph(XM[m], y[m], ey[m])); C.append(GCOL[k]); M.append(GMK[k]); L_.append(g)
                if fit and g in fit["groups"]:
                    fl = ROOT.TF1(f"fl_{est}_{nm}_{k}", f"exp({fit['coef'][fit['groups'].index(g)]} + {fit['inv_lambda']}*x/1000.)", XM[m].min() - 10, XM[m].max() + 10)
                    fl.SetLineColor(GCOL[k]); fl.SetLineWidth(2); ex.append((fl, "SAME"))
            m = np.isfinite(y) & ~full
            if m.any():
                g1 = graph(XM[m], y[m], ey[m]); g1.SetMarkerStyle(24); g1.SetMarkerColor(ROOT.kGray + 1); g1.SetLineColor(ROOT.kGray + 1); ex.append((g1, "P SAME"))
            if fit:
                txt = ROOT.TLatex(0.20, 0.25, f"1/#lambda = {fit['inv_lambda']:+.2f} #pm {fit['err']:.2f} /m (#chi^{{2}}/ndf {fit['chi2ndf']:.1f})"); txt.SetNDC(); txt.SetTextSize(0.035); ex.append(txt)
            ok = np.isfinite(y)
            if not ok.any(): continue
            lo, hi = np.nanmin((y - ey)[ok]), np.nanmax((y + ey)[ok]); pad = 0.25 * (hi - lo) + 1e-3 * abs(hi)
            DrawHistos(G, L_, XM.min() - 20, XM.max() + 20, "X along the test fibres [mm]", max(lo - pad, 0), hi + 2 * pad, ytit, f"{est}_{nm}",
                       dology=False, drawoptions="P", legendoptions="P", markerstyles=M, mycolors=C, legendPos=[0.55, 0.74, 0.92, 0.89], extraToDraw=ex,
                       outdir=OUTD, lumitext=l, usePDF=False)
            per_fibre[l].append(f"{est}_{nm}")

    # tables
    rows = []
    for l in good:
        f = F[l]; s = 1 / abs(f["time_slope"]); fs = fw_speed(l)
        rows.append([l, int(np.isfinite(f["time"]).sum()), fmt(f["time_slope"], f["time_slope_err"]), f"{s / 0.2998:.3f} c", f"{1000 * s:.0f}",
                     f"{0.2998 / s:.2f}", f"{f['time_resid']:.3f}", fmt(f["time_slope_1line"], f["time_slope_1line_err"]), f"{fs[0]:.0f}" if fs else "&ndash;"])
    shift = {g: np.median([1000 * (F[l]["time_icpt"][g] - F[l]["time_icpt"][GN[0]]) / F[l]["time_slope"] for l in good if g in F[l]["time_icpt"] and GN[0] in F[l]["time_icpt"]])
             for g in GN[1:]}
    speed_tab = table(["Fibre", "runs", "dt/dX [ns/m]", "speed", "[mm/ns]", "n<sub>eff</sub> = c/v", "residual rms [ns]",
                       "one line through all stretches [ns/m]", "framework scan [mm/ns]"], rows)
    rows = []
    for l in good:
        f = F[l]; c = f["fits"].get("lpc"); a = fw_att(l)
        pe = np.array(f["pe"], float)[full]; hf = np.array(f["hf"], float)[full]
        per_g = [fmt(c["per_group"][g]["inv_lambda"], c["per_group"][g]["err"]) if c and g in c["per_group"] else "&ndash;" for g in GN]
        rows.append([l, f"{np.nanmedian(pe):.1f}", f"{np.nanmedian(hf):.2f}", f"{100 * np.nanmax(f['sat']):.0f}%",
                     f"<b>{fmt(c['inv_lambda'], c['err'])}</b> ({c['chi2ndf']:.1f})" if c else "&ndash;", f"<b>{lam(c['inv_lambda'], c['err'])}</b>" if c else "&ndash;",
                     *per_g, *[fmt(f["fits"][e]["inv_lambda"], f["fits"][e]["err"]) if e in f["fits"] else "&ndash;" for e in ("lph", "lpe")],
                     fmt(*a) if a else "&ndash;"])
    att_tab = table(["Fibre", "p.e. per hit", "hit fraction", "clipped hits (max)", "1/&lambda; [/m], light per crossing (&chi;<sup>2</sup>/ndf)", "&lambda; [m]",
                     *[f"{g} only" for g in GN], "light per hit", "light per event", "framework scan (per crossing)"], rows)
    rows = []
    for l in labels:
        row = [l]
        for m_ in SELECTIONS:
            f = JJ["selections"][m_]["fibres"][l]; c = f.get("fits", {}).get("lpc")
            row += [fmt(f["time_slope"], f["time_slope_err"]) if "time_slope" in f else "&ndash;", fmt(c["inv_lambda"], c["err"]) if c else "&ndash;"]
        rows.append(row)
    nsel = {m_: np.median(JJ["selections"][m_]["n_sel"]) for m_ in SELECTIONS}
    cmp_tab = table(["Fibre"] + [f"{h}, {m_}" for m_ in SELECTIONS for h in ("dt/dX [ns/m]", "1/&lambda; [/m] per crossing")], rows)
    cmp_text = ("<p>Selections: " + "; ".join(f"<b>{m_}</b>{' (shown above)' if m_ == MS else ''}: {SELECTIONS[m_]}, median {nsel[m_]:.0f} events per run"
                                             for m_ in SELECTIONS) + ".</p>")
    bad = [l for l in labels if l not in good]
    run_tab = table(["Run", "table X [mm]", "X along the fibres [mm]", "stretch", "events (MCP)", "MCP fraction", "overlap (hit fraction / plateau)", "light fits"],
                    [[r, f"{1000 * xt:+.1f}", f"{XM[i]:+.1f}", GRP[i], J["n_sel"][i], f"{100 * J['n_sel'][i] / max(J['n_positrons'][i], 1):.0f}%",
                      f"{J['overlap_ratio'][i]:.2f}", "" if full[i] else "left out"] for i, (r, xt) in enumerate(zip(J["runs"], J["x_table_m"]))])
    spe_tab = table(["Fibre", "channel", "single p.e. height / integral [ADC]", "integral"],
                    [[l, F[l]["channel"], f"{F[l]['spe_amp']:.0f} / {F[l]['spe_int']:.0f}", "long (10 ns)" if F[l]["type"] == "Sci" else "short (3 ns)"] for l in good])
    nav = " &middot; ".join([f'<a href="ZScan_overview.html">overview</a>'] + [f'<a href="ZScan_{d}.html">{d}</a>' for d in DATASETS if d != NAME])
    intro = f"""
<p>{J['desc']}</p>
<p>Other pages: {nav}. Framework scan: <a href="../PositionScan/">PositionScan</a>.</p>
<h3 style="margin:14px 0 4px 0">Light speed (arrival time vs X along the fibre)</h3>
{speed_tab}
<p>One intercept per stretch. Relative to the first stretch, the arrival times of the later ones sit as if X were
{", ".join(f"{v:+.0f} mm ({g})" for g, v in shift.items())} further along (median over the fibres): a step of about 0.1 ns, or a block move
somewhat smaller than 482 mm. A single line through all stretches (as in the framework scan) takes this step into the slope and gives speeds
2&ndash;5% higher (sapphire 7%); within the stretches the slopes agree with the calorimeter fibres of the z scans. n<sub>eff</sub> = c/v is an effective index.</p>
<h3 style="margin:14px 0 4px 0">Attenuation</h3>
{att_tab}
<p>Light per crossing = light per hit &times; (1 &minus; e<sup>&minus;&mu;</sup>), &mu; from the p.e. per hit. Clipped hits: fraction of the
in-time hits above {SAT:.0f} ADC peak (largest run); the light of a clipped bundle is a lower limit. Fits: one normalisation per stretch, runs with
partial overlap left out; a stretch gets its own &lambda; with at least {1000 * MIN_SPAN:.0f} mm of travel.</p>
{"<p>No usable in-time signal: " + ", ".join(bad) + ".</p>" if bad else ""}
{details("Event selection: loose vs z-scan (spot)", cmp_text + cmp_tab)}
{details("Runs", run_tab)}
{details("Method", f'''<p><b>Events</b> ({SELECTIONS[MS]}): positrons (beam Cherenkov counters fired, muon veto not), a pulse in MCP_DS_0 above its cut (the
test-fibre <code>mcp_clean</code>, as in the z scans) and a valid time reference (MCP_DS_1). The z scans in addition require 1&ndash;3 hodoscope
fibres per axis and the beam within 1 cm of the run median; for single channels of the module this keeps the shower in place, while a test
bundle sees the shower averaged over the beam profile, and the hodoscope cut alone removes about 3/4 of the positrons. Both selections are
compared under &ldquo;Event selection&rdquo;.</p>
<p><b>Hits</b>: CFD time within &plusmn;4 ns (scintillating) / &plusmn;1.2 ns (others) of the arrival peak and peak above 40 ADC; randoms
from the sidebands before and after the pulse subtracted. Light per hit / per event (events without a hit count as zero) / per crossing.
Arrival time: median of the in-time hits. Single p.e.: amplitude peak of the out-of-time pulses of all runs, median integral of the
pulses near it.</p>
<p><b>Overlap</b>: hit fraction of the scintillating bundles (about 40 p.e. per hit, so a hit means a crossing) relative to its median in
X = {1000 * PLATEAU[0]:.0f}&hellip;{1000 * PLATEAU[1]:+.0f} mm; runs below {OVERLAP_FRACTION:.0%} are left out of the light fits.</p>''' + spe_tab)}
"""
    # one grid per section, as the other z-scan pages (page.py)
    sections = [("Summary: light speed, attenuation, overlap", main)] + [(f"{l}: arrival time, light per crossing, hit fraction vs X", ns)
                                                                         for l, ns in per_fibre.items()]
    pl = sum(([f"{n}.png" for n in ns] + ["NEWLINE"] for _, ns in sections), [])[:-1]
    generate_html(pl, OUTD, plots_per_row=3, output_html=OUT, title=f"Z scan: {J['title']}", intro_text=intro)
    html = open(OUT).read(); parts = html.split('    <div class="grid">\n')
    html = parts[0] + "".join(f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{t}</h3>\n    <div class="grid">\n' + part
                              for (t, _), part in zip(sections, parts[1:]))
    open(OUT, "w").write(html)
    print(OUT)
