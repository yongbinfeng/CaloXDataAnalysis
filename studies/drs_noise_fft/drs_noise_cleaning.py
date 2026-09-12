"""Does cleaning the DRS noise help? Pulse features before and after.

Follows drs_noise_fft.py, which found two removable noise components: a DRS4
fixed pattern in cell space and the amplifier resonance at 200-600 MHz. This
script applies the corrections to raw waveforms and compares pulse features
on the same events:

    raw       baseline-subtracted waveform, as the analysis sees it
    cell      minus the per-cell offset pattern (from a disjoint set of events)
    cell+LP   cell, then a zero-phase low-pass at --lowpass MHz

Features per channel: pre-pulse noise sigma, peak amplitude, integral, CFD
time against the MCP, 10-90% rise time, and the mean pulse. The mean pulse
is aligned exactly as the analysis aligns its _VS_ts_mcp profiles (group
reference channel, plus the MCP that get_mcp_reference names for the run,
if any) and averaged over every event, so the raw curve reproduces the
profile in drs_profiles.root, which is overlaid as a check. Timing and rise
time use events with an MCP pulse and an in-time pulse above 5 sigma in the
raw waveform; peak and integral use every event with an MCP pulse. Each
selection is made once, on the raw waveform, and reused for every variant.

The whole run is processed, one channel at a time, so the mean pulses are
over exactly the events the analysis's profiles are over. The cell pattern
is derived from the first half of the run and applied to all of it; the
noise-sigma reduction is reported for both halves so the held-out one can
be checked against the in-sample one.

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

from channels.channel_map import get_mcp_channels, get_mcp_reference
from channels.channel_map import get_service_drs_channels
from configs.selection_config import get_particle_selection, get_service_drs_cut
from configs.plot_style import PlotStyle
from core.plot_manager import PlotManager
from utils.plot_helper import get_run_paths, save_hists_to_file

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DT_NS = 0.2
FS_GHZ = 5.0
BASELINE_WIN = (0, 200)           # analysis: median of these samples
NOISE_WIN = (20, 380)             # pre-pulse samples used for the noise sigma
SIGNAL_WIN = (380, 990)
CFD_FRACTION = 0.20
INTEGRAL_WIN = (-10, 50)          # samples around the peak, as in the analysis energy
MIN_SNR_RAW = 5.0                 # timing/shape features need a real pulse: > 5 sigma
IN_TIME_NS = 4.0                  # ...and it has to sit where this channel's pulses sit
CORE_NS = 2.0
REF_MIN_AMP, REF_FRACTION, REF_OFFSET = 500.0, 0.5, 790   # process_dynamic_led + update_ts
MCP_OFFSET = 500
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
    p.add_argument("--nevents", type=int, default=None,
                   help="process only the first N events (default: the whole run)")
    p.add_argument("--mcp", default=None)
    p.add_argument("--lowpass", type=float, default=500.0, metavar="MHz",
                   help="cut-off of the zero-phase low-pass in the cell+LP variant")
    p.add_argument("--jsroot", action="store_true")
    return p.parse_args()


# ----------------------------------------------------------------- waveforms
class WaveformLoader:
    """Read one branch of the run at a time: a full run of one DRS channel is
    ~100 MB as float32, all of them at once is not."""

    def __init__(self, args):
        with open(args.json_file) as f:
            files = json.load(f)[str(args.run)]
        files = files if isinstance(files, list) else [files]
        self.chain = ROOT.TChain("EventTree")
        for fn in files:
            if not os.path.exists(fn):
                sys.exit(f"ERROR: {fn} not found")
            self.chain.Add(fn)
        self.nevents = args.nevents
        self._cache = {}

    def __call__(self, col, keep=False):
        if col in self._cache:
            return self._cache[col]
        rdf = ROOT.RDataFrame(self.chain)
        if self.nevents:
            rdf = rdf.Range(self.nevents)
        a = rdf.AsNumpy(columns=[col])[col]
        out = (np.stack([np.asarray(v, dtype=np.float32) for v in a])
               if a.dtype == object else np.asarray(a))
        if keep:
            self._cache[col] = out
        return out


def pipeline_flip_list(run_number):
    """Channels the analysis inverts (get_drs_branches_to_flip), so the study
    sees the same polarity the analysis does."""
    from channels.channel_map import build_drs_boards
    from variables.drs import get_drs_branches_to_flip
    boards = build_drs_boards(run_number)
    refs = [c.get_channel_name() for b in boards.values() for c in b if c.is_reference]
    return set(get_drs_branches_to_flip(run_number, drs_channels_ref=refs, drsboards=boards))


def mcp_clean_mask(load, run_number, flipped):
    """The analysis's 'mcp_clean' selection, evaluated the way SelectionManager
    does: for each detector it names, the peak of the (flipped) baseline-
    subtracted waveform in the detector's window compared with its cut.
    Returns the mask and a description."""
    reqs = get_particle_selection("mcp_clean", run_number)
    channels = get_service_drs_channels(run_number)
    mask, parts = None, []
    for det, required in reqs.items():
        ch = channels[det]
        ts_min, ts_max, _, _, cut, method = get_service_drs_cut(det, run_number)
        W = rebaseline(load(ch))
        if ch in flipped:
            W = -W
        seg = W[:, ts_min:ts_max]
        val = seg.max(axis=1) if method == "PeakValue" else seg.sum(axis=1)
        fired = val > cut
        ok = fired if required else ~fired
        mask = ok if mask is None else (mask & ok)
        parts.append(f"{det} {'fired' if required else 'vetoed'} ({method} in [{ts_min},{ts_max}) > {cut})")
    return mask, " and ".join(parts)


def rebaseline(W):
    return W - np.median(W[:, BASELINE_WIN[0]:BASELINE_WIN[1]], axis=1, keepdims=True)


def led_time(w_flipped):
    """process_dynamic_led: 50% crossing below the peak, integer slice, or nan."""
    p = int(w_flipped.argmax())
    a = w_flipped[p]
    if a < REF_MIN_AMP:
        return np.nan
    i = p
    while i > 0 and w_flipped[i] > REF_FRACTION * a:
        i -= 1
    return float(i + 1)


def pipeline_shift(ref_ts, mcp_ts_cfd_ref):
    """Per-event sample shift that turns a raw sample index into the
    analysis's ts_mcp: ts - ref_TS + 790, minus the MCP's ref-corrected CFD
    plus 500 when the analysis has an MCP reference for this run."""
    shift = REF_OFFSET - ref_ts
    if mcp_ts_cfd_ref is not None:
        shift = shift - mcp_ts_cfd_ref + MCP_OFFSET
    return shift


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


def cfd_pipeline(w, i_start, i_end, min_amp=10.0, frac=CFD_FRACTION):
    """compute_cfd_integral as the analysis runs it on the MCP: peak searched
    in [i_start, i_end), min_amplitude 10, walk back stops at i_start, and the
    crossing is always interpolated between scan_idx and scan_idx+1 -- so a
    noise peak yields a (garbage) time rather than an invalid one, exactly as
    the analysis's ts_mcp does for events without an MCP pulse."""
    seg = w[i_start:i_end]
    p = int(seg.argmax()) + i_start
    a = w[p]
    if a < min_amp:
        return np.nan, a
    thr = frac * a
    i = p
    while i > i_start and w[i] > thr:
        i -= 1
    y0, y1 = w[i], w[i + 1]
    return (i + (thr - y0) / (y1 - y0)) if y1 != y0 else float(i + 1), a


def rise_time(w, peak_idx, amp):
    """10% to 90% of the peak on the leading edge, in ns."""
    seg = w[:peak_idx + 1]
    below10 = np.where(seg < 0.1 * amp)[0]
    below90 = np.where(seg < 0.9 * amp)[0]
    if len(below10) == 0 or len(below90) == 0:
        return np.nan
    return (below90[-1] - below10[-1]) * DT_NS


def features(W, t_mcp, sel, ref=None):
    """Per-event features for the selected events. Noise uses every event.

    t_mcp is the MCP time already corrected by its group reference; ref is
    this channel's group-reference time per event, so dt mirrors the
    analysis's _TS_cfd_mcp and is free of the jitter between DRS boards."""
    out = {"noise": W[:, NOISE_WIN[0]:NOISE_WIN[1]].std(axis=1)}
    if ref is None:
        ref = np.zeros(W.shape[0])
    peak, integ, dt, rise = [], [], [], []
    for e in np.where(sel)[0]:
        t, a, p = cfd_time(W[e], *SIGNAL_WIN)
        lo, hi = max(p + INTEGRAL_WIN[0], 0), min(p + INTEGRAL_WIN[1], 1024)
        peak.append(a)
        integ.append(float(W[e, lo:hi].sum()))
        dt.append((t - ref[e] - t_mcp[e]) * DT_NS)
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
    ("dt", "t_{CFD} - t_{MCP} [ns] (reference-corrected)", False, "core"),
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


def mean_pulse(W, shift):
    """Mean waveform in the analysis's ts_mcp coordinate, as its TProfile does:
    every event with a valid reference enters, shifted by its own offset."""
    acc = np.zeros(1024)
    cnt = np.zeros(1024)
    for w, sh in zip(W, shift):
        if not np.isfinite(sh):
            continue
        sh = int(np.floor(sh))                         # Profile1D bins floor(ts_mcp)
        lo, hi = max(0, sh), min(1024, 1024 + sh)      # destination range
        if hi <= lo:
            continue
        acc[lo:hi] += w[lo - sh:hi - sh]
        cnt[lo:hi] += 1
    return np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0), cnt


# ------------------------------------------------------------------------
def main():
    args = parse_args()
    paths = get_run_paths(args.run)
    for p in paths.values():
        os.makedirs(p, exist_ok=True)
    with open(args.channels) as f:
        chmap = json.load(f)
    labels = {name: lab for lab, name in chmap.items()}
    mcp_channels = get_mcp_channels(args.run)
    mcp_name = get_mcp_reference(args.run)       # what the analysis aligns to, if anything
    mcp_det = mcp_channels.get(mcp_name)
    if args.mcp:
        mcp = args.mcp
    elif mcp_det:
        mcp = mcp_det
    else:
        sys.exit(f"run {args.run} has no MCP reference; pass --mcp <channel>")
    start_branch = {c: c.rsplit("_", 1)[0] + "_StartIndexCell" for c in labels}
    ref_branch = {c: c.rsplit("_", 1)[0] + "_Channel8" for c in list(labels) + [mcp]}
    print(f"analysis alignment for run {args.run}: group reference channel"
          + (f" + {mcp_name} ({mcp_det})" if mcp_det else " only (no MCP reference for this run)"))
    print(f"timing reference for the features: {mcp}")

    if mcp_det and mcp_det not in ref_branch:
        ref_branch[mcp_det] = mcp_det.rsplit("_", 1)[0] + "_Channel8"
    load = WaveformLoader(args)
    flipped = pipeline_flip_list(args.run)

    # group reference times (process_dynamic_led on the flipped Channel8)
    ref_ts = {}
    for rb in sorted(set(ref_branch.values())):
        Wr = -rebaseline(load(rb))                      # reference channels are flipped
        ref_ts[rb] = np.array([led_time(Wr[e]) for e in range(Wr.shape[0])])
        del Wr
    nev = len(next(iter(ref_ts.values())))
    half = nev // 2
    print(f"run {args.run}: {nev} events; cell pattern from the first {half}, "
          f"everything evaluated on all {nev} (noise sigma also on the held-out half)")

    # MCP reference for the timing features, from the raw unamplified MCP and
    # corrected by its own group reference, as the analysis does: the same
    # numbers for every variant
    Wm = -rebaseline(load(mcp))                          # MCPs are flipped for run >= 1839
    mcp_win = get_service_drs_cut(mcp_name, args.run)[:2] if mcp_name else SIGNAL_WIN
    t_mcp_raw, a_mcp = np.array([cfd_pipeline(Wm[e], *mcp_win) for e in range(nev)]).T
    del Wm
    t_mcp = t_mcp_raw - ref_ts[ref_branch[mcp]]
    mcp_clean, mcp_clean_desc = mcp_clean_mask(load, args.run, flipped)
    has_mcp = mcp_clean & np.isfinite(t_mcp)              # the analysis's --mcp-clean events
    print(f"MCP CFD window {list(mcp_win)}: {int(np.isfinite(t_mcp).sum())} events get a time; "
          f"mcp_clean [{mcp_clean_desc}] keeps {int(has_mcp.sum())}")
    kernel = lowpass_kernel(args.lowpass)
    # the analysis's MCP term for ts_mcp: same CFD, same window, garbage times
    # from noise included, because that is what enters its profiles
    mcp_cfd_ref = None
    if mcp_det:
        if mcp_det == mcp:
            mcp_cfd_ref = REF_OFFSET + t_mcp_raw - ref_ts[ref_branch[mcp_det]]
        else:
            Wd = -rebaseline(load(mcp_det))
            t = np.array([cfd_pipeline(Wd[e], *mcp_win)[0] for e in range(nev)])
            del Wd
            mcp_cfd_ref = REF_OFFSET + t - ref_ts[ref_branch[mcp_det]]
    profiles = ROOT.TFile(os.path.join(paths["root"], "drs_profiles.root"), "READ") \
        if os.path.exists(os.path.join(paths["root"], "drs_profiles.root")) else None

    hists, summary, examples, pulses = [], {}, {}, {}
    ev_mask = np.ones(nev, bool)
    for c, lab in labels.items():
        start = load(start_branch[c], keep=True).astype(int)
        W_raw = load(c)
        if c in flipped:                                  # same polarity as the analysis
            W_raw = -W_raw
            print(f"  {lab}: inverted, as the analysis does (it is in the flip list)")
        W_raw = rebaseline(W_raw)
        pattern = cell_pattern(W_raw[:half], start[:half])
        variants = {"raw": W_raw}
        variants["cell"] = subtract_pattern(W_raw, start, pattern)
        variants["cell+LP"] = rebaseline(lowpass(variants["cell"], kernel))

        # selections, on the raw waveform only, applied to every variant
        raw_cfd = np.array([cfd_time(W_raw[e], *SIGNAL_WIN) for e in range(nev)])
        amp_raw, tpk_raw = raw_cfd[:, 1], raw_cfd[:, 2]
        sigma_raw = float(np.median(W_raw[:, NOISE_WIN[0]:NOISE_WIN[1]].std(axis=1)))
        ref_c = ref_ts[ref_branch[c]]
        sel_mcp = ev_mask & has_mcp & np.isfinite(ref_c)         # peak, integral
        bright = sel_mcp & (amp_raw > MIN_SNR_RAW * sigma_raw) & (amp_raw < 2500)
        dt_raw = (tpk_raw - ref_c - t_mcp) * DT_NS
        # where this channel's pulses sit relative to the MCP: the mode, so
        # out-of-time pulses on dim channels do not drag it around
        if bright.sum() >= 20:
            hcount, edges = np.histogram(dt_raw[bright], bins=np.arange(-60, 60.1, 2.0))
            centre = 0.5 * (edges[hcount.argmax()] + edges[hcount.argmax() + 1])
        else:
            centre = float(np.nanmedian(dt_raw[bright])) if bright.any() else 0.0
        sel = bright & (np.abs(dt_raw - centre) < IN_TIME_NS)   # timing, rise, examples
        feats = {v: features(W, t_mcp, sel, ref_c) for v, W in variants.items()}
        for v, W in variants.items():
            noise_all = feats[v]["noise"]
            feats[v]["noise_in_sample"] = noise_all[:half]       # pattern derived here
            feats[v]["noise"] = noise_all[half:]                 # held-out half
            f_all = features(W, t_mcp, sel_mcp, ref_c)           # every mcp_clean event
            feats[v]["peak"], feats[v]["integral"] = f_all["peak"], f_all["integral"]

        gain, f_med, f_iqr = ringing_coherence(W_raw[half:])
        summary[lab] = {"n_selected": int(sel.sum()), "n_mcp": int(sel_mcp.sum()),
                        "sigma_raw": sigma_raw, "in_time_centre_ns": float(centre),
                        "in_time_fraction_of_bright": float(sel.sum() / max(bright.sum(), 1)),
                        "ringing_extrapolation_gain": gain,
                        "ringing_fit_MHz": f_med, "ringing_fit_iqr_MHz": f_iqr}
        for v in VARIANTS:
            f = feats[v]
            cs, cf = core_stats(f["dt"])
            summary[lab][v] = {
                "noise_sigma": float(np.median(f["noise"])),
                "noise_sigma_in_sample": float(np.median(f["noise_in_sample"])),
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
        shift = pipeline_shift(ref_ts[ref_branch[c]], mcp_cfd_ref)
        # The three treatments are compared on the analysis's mcp_clean events
        # (what check_drs_mcp --mcp-clean books), so the panel compares like
        # with like. "raw, all events" is kept for reference: it is what a
        # profile booked without --mcp-clean looks like, diluted by the events
        # whose MCP time is noise.
        hit = np.where(has_mcp, shift, np.nan)
        pulses[lab] = {v: mean_pulse(W, hit) for v, W in variants.items()}
        raw_all = mean_pulse(W_raw, shift)                # not plotted: only to tell which
        summary[lab]["mean_pulse_peak_all"] = float(raw_all[0].max())   # profile the file holds
        summary[lab]["mcp_clean_fraction"] = float(has_mcp.mean())
        summary[lab]["mean_pulse_peak_mcp_clean"] = float(pulses[lab]["raw"][0].max())
        if profiles:
            hp = profiles.Get(f"prof_{c}_blsub_VS_ts_mcp")
            if hp:
                prof_y = np.array([hp.GetBinContent(i) for i in range(1, 1025)])
                pulses[lab]["analysis profile"] = (prof_y, None)
                summary[lab]["analysis_profile_peak"] = float(prof_y.max())
                pk = int(prof_y.argmax())
                lo_b, hi_b = max(pk - 60, 0), min(pk + 150, 1024)
                # the profile was booked either over every event or with
                # --mcp-clean; say which of the study's curves it agrees with
                d_all = float(np.abs(raw_all[0][lo_b:hi_b] - prof_y[lo_b:hi_b]).max())
                d_clean = float(np.abs(pulses[lab]["raw"][0][lo_b:hi_b] - prof_y[lo_b:hi_b]).max())
                summary[lab]["profile_matches"] = "raw, mcp_clean" if d_clean <= d_all else "raw"
                summary[lab]["raw_vs_profile_max_abs_diff"] = d_all
                summary[lab]["mcp_clean_vs_profile_max_abs_diff"] = d_clean
        ex = np.where(sel)[0]
        if len(ex):
            e = ex[len(ex) // 2]
            examples[lab] = {v: W[e].copy() for v, W in variants.items()}
        print(f"  {lab:16s} mcp_clean={sel_mcp.sum():4d} in-time={sel.sum():4d} "
              f"({summary[lab]['in_time_fraction_of_bright']:.0%} of bright)  noise " +
              " -> ".join(f"{summary[lab][v]['noise_sigma']:.2f}" for v in VARIANTS) +
              f" (in-sample {summary[lab]['cell']['noise_sigma_in_sample']:.2f})" +
              "   dt core sigma " +
              " -> ".join(f"{summary[lab][v]['dt_core_sigma_ns']:.2f}" for v in VARIANTS))

    save_hists_to_file(hists, os.path.join(paths["root"], "drs_cleaning.root"))
    with open(os.path.join(paths["root"], "drs_cleaning_summary.json"), "w") as f:
        json.dump({"run": args.run, "nevents": nev, "lowpass_MHz": args.lowpass,
                   "min_snr_raw": MIN_SNR_RAW, "in_time_ns": IN_TIME_NS,
                   "alignment": "group reference" + (f" + {mcp_name}" if mcp_det else " only"),
                   "timing_reference": mcp,
                   "channels": summary}, f, indent=1)

    # ------------------------------------------------------------ plots
    pm = PlotManager(paths["root"], paths["plots"], paths["html"], args.run,
                     use_jsroot=args.jsroot,
                     selection_text=f"**Cleaning study.** Cell pattern from events 0-{half - 1}; "
                                    f"noise sigma quoted on the held-out events {half}-{nev - 1}; "
                                    f"everything else on all {nev} events. Peak and integral: "
                                    f"every mcp_clean event ({mcp_clean_desc}). Timing and rise time: "
                                    f"mcp_clean events whose raw pulse is above {MIN_SNR_RAW:.0f} sigma "
                                    f"and within {IN_TIME_NS:.0f} ns of where this channel's pulses "
                                    f"sit. Mean pulse: every event, aligned as the analysis aligns "
                                    f"its _VS_ts_mcp profiles, with that profile overlaid. Selections "
                                    f"are made on the raw waveform and reused for every variant. "
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
    treatments_text = (
        "\n**What the curves are.** Every plot overlays the same events under three treatments of "
        "the waveform:\n"
        f"* **raw** — the waveform as the analysis sees it: baseline subtracted (median of samples "
        f"{BASELINE_WIN[0]}–{BASELINE_WIN[1] - 1}), polarity as in the analysis. Nothing else.\n"
        f"* **cell** — raw minus the DRS4 fixed pattern. Each of the 1024 sampling cells has a "
        f"small constant offset (2–4 ADC RMS) left over from the voltage calibration; the trigger "
        f"cell (StartIndexCell) moves it around in sample number every event, so it looks like "
        f"noise unless events are realigned by cell. The per-cell offsets are measured on the "
        f"pre-pulse samples of the first half of the run, subtracted from every sample according "
        f"to its cell, and the baseline is re-taken. It removes 10–30% of the noise variance and "
        f"cannot change the pulse.\n"
        f"* **cell+LP** — cell, then a zero-phase low-pass filter: a symmetric windowed-sinc FIR "
        f"({lowpass_kernel.__defaults__[0]} taps, cut-off {args.lowpass:.0f} MHz) aimed at the "
        f"amplifier noise, which peaks at 200–600 MHz. Symmetric means it adds no time shift. It "
        f"lowers the noise sigma a further ~12% but also clips the 2 ns Cherenkov pulses, costing "
        f"6–10% of their peak; the integral and the rise time are unchanged.\n"
        "* **analysis profile** (mean-pulse panel only) — prof_..._VS_ts_mcp read straight from "
        "drs_profiles.root. The raw curve reproduces it; anything between raw and the other "
        "curves is the cleaning.")
    summary_html = pm.generate_html("DRS/DRS_Cleaning.html", plots_per_row=2,
                                    intro_text=treatments_text,
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
        # mean pulse in the analysis's ts_mcp coordinate, zoomed on the peak
        curves = list(VARIANTS) + (["analysis profile"] if "analysis profile" in pulses[lab] else [])
        ref_curve = pulses[lab]["analysis profile"][0] if "analysis profile" in pulses[lab] else pulses[lab]["raw"][0]
        pk = int(ref_curve.argmax())
        lo_s, hi_s = max(pk - 60, 0), min(pk + 140, 1024)
        hs = []
        for i, v in enumerate(curves):
            h = ROOT.TH1D(f"pulse_{lab}_{i}", "", 1024, 0, 1024); h.SetDirectory(0)
            for b, val in enumerate(pulses[lab][v][0]):
                h.SetBinContent(b + 1, float(val))
            hs.append(h)
        ymax = max(h.GetMaximum() for h in hs) * 1.25 or 1
        ncol = len(curves)
        colours = (COLOURS + [ROOT.kGray + 2])[:ncol]
        styles = ([1, 1, 1, 2])[:ncol]
        pm.plot_1d(hs, f"{lab}_meanpulse", "ts_{mcp} (mcp_clean events)", (lo_s, hi_s),
                   "mean ADC", (-0.2 * ymax, ymax), legends=curves,
                   style=PlotStyle(dology=False, drawoptions="HIST",
                                   mycolors=colours, linestyles=styles,
                                   addOverflow=False, addUnderflow=False,
                                   legendPos=[0.48, 0.62, 0.90, 0.90], legendoptions="L"),
                   extra_text=lab)
        pm.add_newline()
    channels_html = pm.generate_html("DRS/DRS_Cleaning_Channels.html", plots_per_row=6,
                                     intro_text=treatments_text,
                                     title=f"DRS noise cleaning per channel, run {args.run}")
    print(f"\n{summary_html}\n{channels_html}")


if __name__ == "__main__":
    main()
