"""Does cleaning the DRS noise help? Pulse features before and after.

Follows drs_noise_fft.py, which found two removable noise components: a DRS4
fixed pattern in cell space and the amplifier resonance at 200-600 MHz. This
script applies the corrections to raw waveforms and compares pulse features
on the same events:

    raw       baseline-subtracted waveform, as the analysis sees it
    cell      minus the per-cell offset pattern (from a disjoint set of events)
    cell+LP   cell, then a zero-phase low-pass at --lowpass MHz

Features per channel: pre-pulse noise sigma, peak amplitude, integral, CFD
time against the MCP, 10-90% rise time, and the MCP-aligned mean pulse.
Event selection is made once on the raw waveform and reused for every
variant, so differences are the cleaning's doing.

Outputs follow the analysis conventions:
    results/root/Run<N>/drs_cleaning.root, drs_cleaning_summary.json
    results/plots/Run<N>/DRS_Cleaning/*.png
    results/html/Run<N>/DRS/DRS_Cleaning.html, DRS_Cleaning_Channels.html

    python3 studies/drs_noise_fft/drs_noise_cleaning.py --run 1994 \\
        --channels data/channel_maps/testingfibers.json [--jsroot]
"""
import argparse
import json
import os
import sys

import numpy as np
import ROOT

from channels.channel_map import get_mcp_channels
from configs.plot_style import PlotStyle
from core.plot_manager import PlotManager
from utils.plot_helper import get_run_paths, save_hists_to_file

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DT_NS = 0.2
FS_GHZ = 5.0
NOISE_WIN = (20, 380)
SIGNAL_WIN = (380, 990)
CFD_FRACTION = 0.20
INTEGRAL_WIN = (-10, 50)          # samples around the peak, as in the analysis energy
MIN_AMP_RAW = 160.0               # timing/shape features need a real pulse
CORE_NS = 2.0
VARIANTS = ["raw", "cell", "cell+LP"]
COLOURS = [1, ROOT.kAzure + 1, ROOT.kRed + 1]

STYLE_LIN = PlotStyle(dology=False, drawoptions="HIST", mycolors=COLOURS,
                      addOverflow=False, addUnderflow=False,
                      legendPos=[0.55, 0.72, 0.90, 0.90], legendoptions="L")
STYLE_LOG = PlotStyle(dology=True, drawoptions="HIST", mycolors=COLOURS,
                      addOverflow=False, addUnderflow=False,
                      legendPos=[0.55, 0.72, 0.90, 0.90], legendoptions="L")
STYLE_BAR = PlotStyle(dology=False, drawoptions="HIST", mycolors=COLOURS,
                      addOverflow=False, addUnderflow=False,
                      legendPos=[0.60, 0.75, 0.90, 0.90], legendoptions="L",
                      W_ref=1200, H_ref=500)
STYLE_WAVE = PlotStyle(dology=False, drawoptions="HIST", mycolors=COLOURS,
                       addOverflow=False, addUnderflow=False,
                       legendPos=[0.60, 0.72, 0.90, 0.90], legendoptions="L",
                       W_ref=900, H_ref=500)


def parse_args():
    from configs.run_config import run_number, jsonFile
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", type=int, default=run_number)
    p.add_argument("--channels", required=True, metavar="FILE")
    p.add_argument("--json-file", default=jsonFile)
    p.add_argument("--nevents", type=int, default=6000,
                   help="first half derives the cell pattern, second half is evaluated")
    p.add_argument("--mcp", default=None)
    p.add_argument("--lowpass", type=float, default=500.0, metavar="MHz",
                   help="cut-off of the zero-phase low-pass in the cell+LP variant")
    p.add_argument("--jsroot", action="store_true")
    return p.parse_args()


# ----------------------------------------------------------------- waveforms
def load_waveforms(args, cols):
    with open(args.json_file) as f:
        files = json.load(f)[str(args.run)]
    files = files if isinstance(files, list) else [files]
    chain = ROOT.TChain("EventTree")
    for fn in files:
        if not os.path.exists(fn):
            sys.exit(f"ERROR: {fn} not found")
        chain.Add(fn)
    raw = ROOT.RDataFrame(chain).Range(args.nevents).AsNumpy(columns=cols)
    return {c: (np.stack([np.asarray(v, dtype=np.float32) for v in raw[c]])
                if raw[c].dtype == object else np.asarray(raw[c]))
            for c in cols}


def rebaseline(W):
    return W - np.median(W[:, NOISE_WIN[0]:NOISE_WIN[1]], axis=1, keepdims=True)


def cell_pattern(W, start):
    """Per-cell mean offset from the pre-pulse samples of the given events."""
    idx = np.arange(*NOISE_WIN)
    seg = W[:, NOISE_WIN[0]:NOISE_WIN[1]]
    seg = seg - seg.mean(axis=1, keepdims=True)
    ok = np.abs(seg).max(axis=1) < 6 * seg.std(axis=1).mean()
    cells = (start[ok][:, None] + idx[None, :]) % 1024
    cnt = np.bincount(cells.ravel(), minlength=1024)
    acc = np.bincount(cells.ravel(), weights=seg[ok].ravel().astype(np.float64), minlength=1024)
    return np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0).astype(np.float32)


def subtract_pattern(W, start, pattern):
    cells = (start[:, None] + np.arange(1024)[None, :]) % 1024
    return rebaseline(W - pattern[cells])


def lowpass_kernel(fc_mhz, ntaps=41):
    """Windowed-sinc FIR, symmetric so it adds no delay."""
    fc = fc_mhz * 1e6 / (FS_GHZ * 1e9)
    n = np.arange(ntaps) - (ntaps - 1) / 2
    k = 2 * fc * np.sinc(2 * fc * n) * np.hamming(ntaps)
    return (k / k.sum()).astype(np.float32)


def lowpass(W, kernel):
    return np.stack([np.convolve(w, kernel, mode="same") for w in W])


# ------------------------------------------------------------------ features
def cfd_time(w, lo, hi):
    seg = w[lo:hi]
    p = int(seg.argmax())
    a = float(seg[p])
    if a <= 0:
        return np.nan, a, p + lo
    thr = CFD_FRACTION * a
    i = p
    while i > 0 and seg[i] > thr:
        i -= 1
    if i >= p:
        return np.nan, a, p + lo
    y0, y1 = seg[i], seg[i + 1]
    return lo + i + ((thr - y0) / (y1 - y0) if y1 != y0 else 0.0), a, p + lo


def rise_time(w, peak_idx, amp):
    """10% to 90% of the peak on the leading edge, in ns."""
    seg = w[:peak_idx + 1]
    below10 = np.where(seg < 0.1 * amp)[0]
    below90 = np.where(seg < 0.9 * amp)[0]
    if len(below10) == 0 or len(below90) == 0:
        return np.nan
    return (below90[-1] - below10[-1]) * DT_NS


def features(W, t_mcp, sel):
    """Per-event features for the selected events. Noise uses every event."""
    out = {"noise": W[:, NOISE_WIN[0]:NOISE_WIN[1]].std(axis=1)}
    peak, integ, dt, rise = [], [], [], []
    for e in np.where(sel)[0]:
        t, a, p = cfd_time(W[e], *SIGNAL_WIN)
        lo, hi = max(p + INTEGRAL_WIN[0], 0), min(p + INTEGRAL_WIN[1], 1024)
        peak.append(a)
        integ.append(float(W[e, lo:hi].sum()))
        dt.append((t - t_mcp[e]) * DT_NS)
        rise.append(rise_time(W[e], p, a))
    out.update(peak=np.array(peak), integral=np.array(integ),
               dt=np.array(dt), rise=np.array(rise))
    return out


def robust_sigma(x):
    x = x[np.isfinite(x)]
    if len(x) < 5:
        return np.nan
    return 0.7413 * (np.percentile(x, 75) - np.percentile(x, 25))


def core_stats(dt):
    dt = dt[np.isfinite(dt)]
    if len(dt) < 15:
        return np.nan, np.nan
    core = np.abs(dt - np.median(dt)) < CORE_NS
    return robust_sigma(dt[core]), float(core.mean())


# ------------------------------------------------- can the ringing be subtracted?
COH_FIT, COH_PRED = slice(20, 380), slice(620, 990)
COH_GRID = np.arange(240.0, 500.0, 4.0)      # MHz
_T = np.arange(1024) * DT_NS * 1e-9


def _design(f_mhz, sl):
    w = 2 * np.pi * f_mhz * 1e6
    return np.column_stack([np.cos(w * _T[sl]), np.sin(w * _T[sl])])


def ringing_coherence(W, max_events=200):
    """Fit a sinusoid to the pre-pulse samples and extrapolate it past the pulse.

    If the ~300 MHz ringing kept its phase across the record it could be fitted
    where there is no pulse and subtracted underneath it. Returns the median
    fraction of post-pulse variance that the extrapolation removes: positive
    means coherent and subtractable, zero or negative means random phase.
    """
    quiet = np.abs(W[:, SIGNAL_WIN[0]:SIGNAL_WIN[1]]).max(axis=1) < 40   # no pulse at all
    gains, freqs = [], []
    for w in W[quiet][:max_events].astype(np.float64):
        w = w - w[COH_FIT].mean()
        best = None
        for f in COH_GRID:
            X = _design(f, COH_FIT)
            coef, *_ = np.linalg.lstsq(X, w[COH_FIT], rcond=None)
            res = ((w[COH_FIT] - X @ coef) ** 2).mean()
            if best is None or res < best[0]:
                best = (res, f, coef)
        _, f, coef = best
        post = w[COH_PRED] - w[COH_PRED].mean()
        pred = _design(f, COH_PRED) @ coef
        gains.append(1 - ((post - pred) ** 2).mean() / (post ** 2).mean())
        freqs.append(f)
    if not gains:
        return np.nan, np.nan, np.nan
    return (float(np.median(gains)), float(np.median(freqs)),
            float(np.percentile(freqs, 75) - np.percentile(freqs, 25)))


# ---------------------------------------------------------------- histograms
FEATURES = [  # key, axis title, log y, range rule
    ("noise", "pre-pulse noise #sigma [ADC]", False, (0, 15)),
    ("peak", "peak amplitude [ADC]", True, "p99"),
    ("integral", "integral, -10..+50 samples [ADC]", True, "p99"),
    ("dt", "t_{CFD} - t_{MCP} [ns]", False, "core"),
    ("rise", "rise time 10-90% [ns]", False, (0, 8)),
]


def make_hist(name, title, values, lo, hi, nbins=60):
    h = ROOT.TH1D(name, title, nbins, lo, hi)
    h.SetDirectory(0)
    for v in values[np.isfinite(values)]:
        h.Fill(float(v))
    return h


def feature_range(rule, raw_values):
    if isinstance(rule, tuple):
        return rule
    v = raw_values[np.isfinite(raw_values)]
    if len(v) == 0:
        return 0.0, 1.0
    if rule == "p99":
        return 0.0, float(np.percentile(v, 99) * 1.1 or 1.0)
    if rule == "core":
        m = float(np.median(v))
        return m - 6.0, m + 6.0
    return float(v.min()), float(v.max())


def mean_pulse(W, t_mcp, sel, halfwidth=100):
    """MCP-aligned mean waveform around the median pulse position."""
    ev = np.where(sel)[0]
    if len(ev) == 0:
        return None
    shift = np.round(t_mcp[ev] - np.median(t_mcp[ev])).astype(int)
    centre = int(np.median([cfd_time(W[e], *SIGNAL_WIN)[2] for e in ev[:200]]))
    acc = np.zeros(2 * halfwidth)
    n = 0
    for e, s in zip(ev, shift):
        lo = centre + s - halfwidth
        if lo < 0 or lo + 2 * halfwidth > 1024:
            continue
        acc += W[e, lo:lo + 2 * halfwidth]
        n += 1
    return (acc / max(n, 1), centre - halfwidth)


# ------------------------------------------------------------------------
def main():
    args = parse_args()
    paths = get_run_paths(args.run)
    for p in paths.values():
        os.makedirs(p, exist_ok=True)
    with open(args.channels) as f:
        chmap = json.load(f)
    labels = {name: lab for lab, name in chmap.items()}
    mcp = args.mcp or get_mcp_channels(args.run)["MCP_DS_0"]
    start_branch = {c: c.rsplit("_", 1)[0] + "_StartIndexCell" for c in labels}

    data = load_waveforms(args, list(labels) + [mcp] + sorted(set(start_branch.values())))
    nev = data[mcp].shape[0]
    half = nev // 2
    print(f"run {args.run}: {nev} events; cell pattern from the first {half}, "
          f"features from the last {nev - half}")

    # MCP reference from the raw, unamplified MCP: same for every variant
    Wm = rebaseline(data[mcp])
    t_mcp, a_mcp = np.array([cfd_time(-Wm[e], *SIGNAL_WIN)[:2] for e in range(nev)]).T
    has_mcp = (a_mcp > 50) & np.isfinite(t_mcp)
    kernel = lowpass_kernel(args.lowpass)

    hists, summary, examples, pulses = [], {}, {}, {}
    ev_mask = np.zeros(nev, bool); ev_mask[half:] = True
    for c, lab in labels.items():
        start = data[start_branch[c]].astype(int)
        W_raw = rebaseline(data[c])
        pattern = cell_pattern(W_raw[:half], start[:half])
        variants = {"raw": W_raw}
        variants["cell"] = subtract_pattern(data[c], start, pattern)
        variants["cell+LP"] = rebaseline(lowpass(variants["cell"], kernel))

        # one selection, on the raw waveform, applied to every variant
        amp_raw = np.array([cfd_time(W_raw[e], *SIGNAL_WIN)[1] for e in range(nev)])
        sel = ev_mask & has_mcp & (amp_raw > MIN_AMP_RAW) & (amp_raw < 2500)
        feats = {v: features(W, t_mcp, sel) for v, W in variants.items()}
        for v in feats:                      # noise only on the evaluation half
            feats[v]["noise"] = feats[v]["noise"][half:]

        gain, f_med, f_iqr = ringing_coherence(W_raw[half:])
        summary[lab] = {"n_selected": int(sel.sum()),
                        "ringing_extrapolation_gain": gain,
                        "ringing_fit_MHz": f_med, "ringing_fit_iqr_MHz": f_iqr}
        for v in VARIANTS:
            f = feats[v]
            cs, cf = core_stats(f["dt"])
            summary[lab][v] = {
                "noise_sigma": float(np.median(f["noise"])),
                "peak_median": float(np.nanmedian(f["peak"])) if len(f["peak"]) else np.nan,
                "integral_median": float(np.nanmedian(f["integral"])) if len(f["integral"]) else np.nan,
                "rise_median_ns": float(np.nanmedian(f["rise"])) if len(f["rise"]) else np.nan,
                "dt_core_sigma_ns": float(cs), "dt_core_fraction": float(cf),
            }
            r = summary[lab][v]
            r["peak_snr"] = (r["peak_median"] / r["noise_sigma"]
                             if r["noise_sigma"] and np.isfinite(r["peak_median"]) else np.nan)
        for key, title, _, rule in FEATURES:
            lo, hi = feature_range(rule, feats["raw"][key])
            for v in VARIANTS:
                hists.append(make_hist(f"h_{lab}_{v.replace('+', '_')}_{key}",
                                       f"{lab} {v};{title};events", feats[v][key], lo, hi))
        pulses[lab] = {v: mean_pulse(W, t_mcp, sel) for v, W in variants.items()}
        ex = np.where(sel)[0]
        if len(ex):
            e = ex[len(ex) // 2]
            examples[lab] = {v: W[e].copy() for v, W in variants.items()}
        print(f"  {lab:16s} n={sel.sum():4d}  noise " +
              " -> ".join(f"{summary[lab][v]['noise_sigma']:.2f}" for v in VARIANTS) +
              "   dt core sigma " +
              " -> ".join(f"{summary[lab][v]['dt_core_sigma_ns']:.2f}" for v in VARIANTS))

    save_hists_to_file(hists, os.path.join(paths["root"], "drs_cleaning.root"))
    with open(os.path.join(paths["root"], "drs_cleaning_summary.json"), "w") as f:
        json.dump({"run": args.run, "nevents": nev, "lowpass_MHz": args.lowpass,
                   "min_amp_raw": MIN_AMP_RAW, "channels": summary}, f, indent=1)

    # ------------------------------------------------------------ plots
    pm = PlotManager(paths["root"], paths["plots"], paths["html"], args.run,
                     use_jsroot=args.jsroot,
                     selection_text=f"**Cleaning study.** Cell pattern from events 0-{half - 1}, "
                                    f"features from events {half}-{nev - 1}. Timing and shape "
                                    f"features use events with an MCP pulse and raw amplitude "
                                    f"> {MIN_AMP_RAW:.0f} ADC (the same events for every variant). "
                                    f"Low-pass cut-off {args.lowpass:.0f} MHz.")
    pm.set_output_dir("DRS_Cleaning")
    labs = list(labels.values())

    def bar_plot(name, ytitle, yrange, series, legends):
        """Per-channel bars drawn directly: DrawHistos would rebuild the axis
        and drop the channel names, so this registers its own PNG/JSON."""
        c = ROOT.TCanvas(f"c_{name}", "", 1200, 500)
        c.SetBottomMargin(0.30); c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.08)
        hs = []
        for i, vals in enumerate(series):
            h = ROOT.TH1D(f"{name}_{i}", "", len(labs), 0, len(labs)); h.SetDirectory(0)
            for j, (lab, val) in enumerate(zip(labs, vals)):
                h.SetBinContent(j + 1, val if np.isfinite(val) else 0.0)
                h.GetXaxis().SetBinLabel(j + 1, lab)
            h.SetLineColor(COLOURS[i % len(COLOURS)]); h.SetLineWidth(2)
            h.GetXaxis().LabelsOption("v"); h.GetXaxis().SetLabelSize(0.045)
            h.GetYaxis().SetTitle(ytitle); h.GetYaxis().SetTitleOffset(0.8)
            h.SetMinimum(yrange[0]); h.SetMaximum(yrange[1])
            hs.append(h)
        hs[0].Draw("hist")
        for h in hs[1:]:
            h.Draw("hist same")
        lg = ROOT.TLegend(0.70, 0.74, 0.96, 0.91); lg.SetBorderSize(0); lg.SetFillStyle(0)
        for h, l in zip(hs, legends):
            lg.AddEntry(h, l, "l")
        lg.Draw()
        tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.05); tx.SetTextFont(62)
        tx.DrawLatex(0.09, 0.94, "CaloX")
        tx.SetTextFont(42); tx.SetTextAlign(31)
        tx.DrawLatex(0.97, 0.94, f"Run {args.run}")
        os.makedirs(pm.get_output_dir(), exist_ok=True)
        if pm.use_jsroot:
            pm._canvas_jsons[name] = ROOT.TBufferJSON.ToJSON(c).Data()
        else:
            c.SaveAs(os.path.join(pm.get_output_dir(), f"{name}.png"))
        pm.add_plot(name)
        keep.extend(hs + [c, lg])

    keep = []
    col = lambda v, key: [summary[lab][v][key] for lab in labs]
    bar_plot("summary_noise_sigma", "pre-pulse noise #sigma [ADC]", (0, 10),
             [col(v, "noise_sigma") for v in VARIANTS], VARIANTS)
    bar_plot("summary_peak_snr", "median peak / noise #sigma", (0, 220),
             [col(v, "peak_snr") for v in VARIANTS], VARIANTS)
    bar_plot("summary_peak_ratio", "median peak relative to raw", (0.5, 1.25),
             [[summary[lab][v]["peak_median"] / summary[lab]["raw"]["peak_median"]
               if summary[lab]["raw"]["peak_median"] else np.nan for lab in labs] for v in VARIANTS],
             VARIANTS)
    bar_plot("summary_integral_ratio", "median integral relative to raw", (0.5, 1.25),
             [[summary[lab][v]["integral_median"] / summary[lab]["raw"]["integral_median"]
               if summary[lab]["raw"]["integral_median"] else np.nan for lab in labs] for v in VARIANTS],
             VARIANTS)
    bar_plot("summary_dt_core_sigma", "core #sigma(t_{CFD} - t_{MCP}) [ns]", (0, 2.0),
             [col(v, "dt_core_sigma_ns") for v in VARIANTS], VARIANTS)
    bar_plot("summary_dt_core_fraction", "fraction of events within #pm2 ns", (0, 1.05),
             [col(v, "dt_core_fraction") for v in VARIANTS], VARIANTS)
    bar_plot("summary_rise_time", "median rise time 10-90% [ns]", (0, 4),
             [col(v, "rise_median_ns") for v in VARIANTS], VARIANTS)
    bar_plot("summary_ringing_coherence",
             "post-pulse variance removed by extrapolating the pre-pulse sinusoid", (-0.5, 0.5),
             [[summary[lab]["ringing_extrapolation_gain"] for lab in labs]],
             ["single tone fitted on samples 20-380, tested on 620-990"])
    pm.add_newline()
    # example waveforms: one Cherenkov-type and one scintillating channel
    for lab in [l for l in labs if l in examples][:4]:
        hs = []
        for i, v in enumerate(VARIANTS):
            h = ROOT.TH1D(f"wave_{lab}_{i}", "", 1024, 0, 1024); h.SetDirectory(0)
            for b, val in enumerate(examples[lab][v]):
                h.SetBinContent(b + 1, float(val))
            hs.append(h)
        ymax = float(max(examples[lab]["raw"].max() * 1.2, 50))
        pm.plot_1d(hs, f"example_waveform_{lab}", "sample", (SIGNAL_WIN[0] - 40, SIGNAL_WIN[1] - 250),
                   "ADC (baseline subtracted)", (-0.25 * ymax, ymax),
                   legends=VARIANTS, style=STYLE_WAVE, extra_text=lab)
    summary_html = pm.generate_html("DRS/DRS_Cleaning.html", plots_per_row=2,
                                    title=f"DRS noise cleaning, run {args.run}")

    # per-channel pages: one row per channel
    pm.set_output_dir("DRS_Cleaning")
    hist_by_name = {h.GetName(): h for h in hists}
    for lab in labs:
        for key, title, logy, rule in FEATURES:
            hs = [hist_by_name[f"h_{lab}_{v.replace('+', '_')}_{key}"] for v in VARIANTS]
            lo, hi = hs[0].GetXaxis().GetXmin(), hs[0].GetXaxis().GetXmax()
            ymax = max(h.GetMaximum() for h in hs) * (8 if logy else 1.3) or 1
            pm.plot_1d(hs, f"{lab}_{key}", title, (lo, hi), "events",
                       (0.5 if logy else 0, ymax), legends=VARIANTS,
                       style=STYLE_LOG if logy else STYLE_LIN, extra_text=lab)
        # mean pulse
        hs = []
        for i, v in enumerate(VARIANTS):
            mp = pulses[lab][v]
            h = ROOT.TH1D(f"pulse_{lab}_{i}", "", 200, 0, 200); h.SetDirectory(0)
            if mp is not None:
                for b, val in enumerate(mp[0]):
                    h.SetBinContent(b + 1, float(val))
            hs.append(h)
        ymax = max(h.GetMaximum() for h in hs) * 1.2 or 1
        pm.plot_1d(hs, f"{lab}_meanpulse", "sample (MCP-aligned window)", (0, 200),
                   "mean ADC", (-0.2 * ymax, ymax), legends=VARIANTS, style=STYLE_LIN, extra_text=lab)
        pm.add_newline()
    channels_html = pm.generate_html("DRS/DRS_Cleaning_Channels.html", plots_per_row=6,
                                     title=f"DRS noise cleaning per channel, run {args.run}")
    print(f"\n{summary_html}\n{channels_html}")


if __name__ == "__main__":
    main()
