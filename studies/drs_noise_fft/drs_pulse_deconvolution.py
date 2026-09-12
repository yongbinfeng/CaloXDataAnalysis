"""Recover the photon arrival-time profile n(t) from the mean DRS pulses.

The measured pulse is n(t) convolved with the response h(t) of SiPM,
amplifier and DRS. The Cherenkov fibres see prompt light, so their mean
pulse *is* h(t). It is built per DRS board, from that board's Cherenkov
channels, and used to unfold the other channels of the same board: a
response taken from another board carries the residual timing jitter
between boards (~0.3 ns after the reference correction), and that alone
turns a chi2/ndf of 0.4 into 50 for a 6 ns scintillator pulse. Two ways:

  fit   a parametric n(t) -- exponential rise and one or two exponential
        decays -- convolved with h(t) and fitted to the profile (chi2 on the
        TProfile errors). Nothing is divided, so noise is never amplified.
  NNLS  non-negative least squares by projected gradient: n(t) as a free
        histogram of arrival times, constrained only to be >= 0, stopped
        when chi2/ndf reaches 1. Non-parametric: it shows structure the
        exponential model would miss.

Both work in a -6..+12 ns window around the peak, with a linear baseline as
a nuisance. Beyond that the bright scintillating pulses droop (to -13% of
the peak at +40 ns) in a way the 60 ADC Cherenkov pulses do not, so the
response measured on the latter does not apply there; the decay time is
fixed by the first two or three decay times anyway.

On each board one Cherenkov channel is held out of h(t) and unfolded as a
control: it should come back as a spike, whose width is the method's
resolution floor.

Time and scale are kept absolute. The profiles are already on one time
axis (group reference + MCP alignment, ts_mcp), so n(t) stays on it and
differences between channels are light-propagation and path differences.
h(t) has unit area, so n(t) is in ADC per sample and its integral is the
pulse area: with a constant ADC-per-photon, the areas compare the number
of photons captured.

Reads prof_<ch>_blsub_VS_ts_mcp from results/root/Run<N>/drs_profiles.root
(book them with scripts/check_drs_mcp.py --mcp-clean). Writes
results/root/Run<N>/drs_deconvolution.json and .root (the unfolded n(t) as
histograms on the ts_mcp axis), plots under
results/plots/Run<N>/DRS_Deconvolution/, results/html/Run<N>/DRS/
DRS_Deconvolution.html, and the overlay pages DRS_Unfolded_overlay.html
(NNLS) and DRS_Unfolded_fit_overlay.html with absolute and peak = 1 views.

    python3 studies/drs_noise_fft/drs_pulse_deconvolution.py --run 1994 \\
        --channels data/channel_maps/testingfibers.json
"""
import argparse
import json
import os
import sys

import numpy as np
import ROOT

from configs.plot_style import PlotStyle
from core.plot_manager import PlotManager
from utils.overlay import build_channel_overlay
from utils.plot_helper import get_run_paths, save_hists_to_file

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DT = 0.2                                  # ns per sample
CHERENKOV = ("Multiclad", "Singleclad", "Quartz700", "Quartz")   # prompt light: response
EXCLUDE = ("Quartz_1",)                   # that channel is MCP_US_1
MIN_RESPONSE_PEAK = 40.0                  # ADC; a Cherenkov profile below this is too dim for h(t)
MAX_TAU_SLOW = 40.0                       # ns; a "slow component" beyond this is a baseline, not light
WIN_PRE, WIN_POST = 30, 60                # fit window: -6 .. +12 ns around the peak
H_PRE, H_POST = 40, 200                   # window used to build h(t)
H_ONSET_FRACTION = 0.02                   # h(t) starts at the last sample below this before its peak
H_PEDESTAL = (0, 30)                      # samples of the h window (-8..-2 ns) whose median is the local pedestal
H_LENGTH = 150                            # samples of h(t) kept after the onset (30 ns)
SHAPE_SYS = 0.02                          # channel-to-channel spread of h(t), as a fraction of the pulse
MAX_NNLS_ITER = 20000
NBINS = 1024                              # the ts_mcp axis of the profiles
OVERLAY_COLOURS = [633, 601, 418, 617, 807, 434, 881, 829, 861, 909, 843, 403]   # predefined ROOT


def parse_args():
    from configs.run_config import run_number
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", type=int, default=run_number)
    p.add_argument("--channels", required=True, metavar="FILE")
    p.add_argument("--jsroot", action="store_true")
    return p.parse_args()


# ------------------------------------------------------------- profiles
def read_profile(f, channel):
    h = f.Get(f"prof_{channel}_blsub_VS_ts_mcp")
    if not h:
        return None
    y = np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])
    e = np.array([h.GetBinError(i) for i in range(1, h.GetNbinsX() + 1)])
    base = np.median(y[200:400])
    noise = np.std(y[200:400] - base)
    return y - base, e, noise


def subsample_peak(y, pk):
    """Parabolic refinement of the peak position, in samples."""
    if 0 < pk < len(y) - 1:
        y0, y1, y2 = y[pk - 1], y[pk], y[pk + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            return pk + 0.5 * (y0 - y2) / den
    return float(pk)


def shift_fractional(y, delta):
    """Shift a sampled curve by delta samples (linear interpolation)."""
    x = np.arange(len(y))
    return np.interp(x - delta, x, y, left=0.0, right=0.0)


def build_response(profiles, labels):
    """Unit-area response from the peak-aligned mean of the Cherenkov channels."""
    stack, weights = [], []
    for lab in labels:
        y, e, noise = profiles[lab]
        pk = int(y.argmax())
        seg = y[pk - H_PRE:pk + H_POST].copy()
        # the profiles sit on a shallow pedestal of 1-5% of the peak just
        # before the pulse (the far baseline was taken 40+ ns earlier); take
        # it out here or it defines the onset instead of the pulse
        seg -= np.median(seg[H_PEDESTAL[0]:H_PEDESTAL[1]])
        seg /= seg[H_PRE]
        seg = shift_fractional(seg, subsample_peak(seg, H_PRE) - H_PRE)
        stack.append(seg)
        weights.append(y[pk])                       # bright channels carry less relative noise
    stack, weights = np.array(stack), np.array(weights) / np.sum(weights)
    h = (stack * weights[:, None]).sum(axis=0)
    spread = np.std(stack, axis=0)
    # onset: walk back from the peak to the last sample below a threshold set
    # above h's own pre-pulse noise. (Searching forward from the start, or
    # using a fixed 1%, finds noise instead; one dim channel in the average
    # then moves the onset by nanoseconds and nothing fits.)
    pk = int(h.argmax())
    pre_noise = np.std(h[:max(pk - 25, 1)])
    thr = max(H_ONSET_FRACTION * h[pk], 3 * pre_noise)
    below = np.where(h[:pk] < thr)[0]
    onset = int(below[-1]) if len(below) else 0
    h = h[onset:onset + H_LENGTH]                   # causal from its onset; the
    h /= h[h > 0].sum()                             # tail keeps its sign
    return h, spread, onset


# ------------------------------------------------------------- forward fit
def model_n(t, t0, tau_r, tau_d, frac2=0.0, tau_d2=1.0):
    """Exponential rise, one or two exponential decays; unit-normalised area."""
    u = np.clip(t - t0, 0, None)
    n = (1 - np.exp(-u / tau_r)) * ((1 - frac2) * np.exp(-u / tau_d) / tau_d
                                    + frac2 * np.exp(-u / tau_d2) / tau_d2)
    n[t < t0] = 0.0
    s = n.sum()
    return n / s if s > 0 else n


def convolve(n, h):
    return np.convolve(n, h, mode="full")[:len(n)]


def nelder_mead(fun, x0, steps, maxiter=4000, tol=1e-8):
    """Small Nelder-Mead: no scipy on the analysis machines."""
    n = len(x0)
    simplex = [np.array(x0, float)]
    for i in range(n):
        p = np.array(x0, float); p[i] += steps[i]; simplex.append(p)
    vals = [fun(p) for p in simplex]
    for _ in range(maxiter):
        order = np.argsort(vals); simplex = [simplex[i] for i in order]; vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol * (1 + abs(vals[0])):
            break
        c = np.mean(simplex[:-1], axis=0)
        xr = c + (c - simplex[-1]); fr = fun(xr)
        if fr < vals[0]:
            xe = c + 2 * (c - simplex[-1]); fe = fun(xe)
            simplex[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            simplex[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (simplex[-1] - c); fc = fun(xc)
            if fc < vals[-1]:
                simplex[-1], vals[-1] = xc, fc
            else:
                for i in range(1, n + 1):
                    simplex[i] = simplex[0] + 0.5 * (simplex[i] - simplex[0]); vals[i] = fun(simplex[i])
    i = int(np.argmin(vals))
    return simplex[i], vals[i]


def error_model(e, y, noise):
    """Profile statistical error, plus the response's shape systematic and
    the profile's baseline noise, in quadrature. Without the systematic a
    2% response uncertainty on a 600 ADC pulse is 12 ADC against a 1 ADC
    statistical error, and no model can reach chi2/ndf = 1."""
    return np.sqrt(e ** 2 + (SHAPE_SYS * np.abs(y)) ** 2 + noise ** 2)


def fit_profile(t, y, sigma, h, two_components):
    """Fit A * (n(t) conv h) to y. Returns parameters, errors, chi2, ndf, model."""
    # amplitude is linear: solve it analytically for every shape trial
    basis = np.column_stack([np.ones_like(t), t - t.mean()])   # linear baseline nuisance

    def linear_part(shape):
        """Weighted LS for amplitude and baseline given the pulse shape."""
        X = np.column_stack([shape, basis]) / sigma[:, None]
        coef, *_ = np.linalg.lstsq(X, y / sigma, rcond=None)
        return coef

    def chi2_of(theta):
        t0, ltr, ltd = theta[:3]
        extra = ((1 / (1 + np.exp(-theta[3])), min(np.exp(theta[4]), MAX_TAU_SLOW))
                 if two_components else ())
        shape = convolve(model_n(t, t0, np.exp(ltr), np.exp(ltd), *extra), h)
        coef = linear_part(shape)
        model = coef[0] * shape + basis @ coef[1:]
        return (((y - model) / sigma) ** 2).sum(), coef

    t_delta = t[int(y.argmax())] - h.argmax() * DT     # where a delta n(t) would have to sit
    best = None
    for ltr in np.log([0.2, 0.6, 1.5]):
        for ltd in np.log([0.1, 0.5, 1.5, 4.0, 10.0]):
            for dt0 in (-1.5, 0.0, 1.5):
                x0 = [t_delta + dt0, ltr, ltd] + ([0.0, np.log(20.0)] if two_components else [])
                steps = [0.5, 0.3, 0.3] + ([0.5, 0.3] if two_components else [])
                xb, fb = nelder_mead(lambda th: chi2_of(th)[0], x0, steps)
                if best is None or fb < best[1]:
                    best = (xb, fb)
    theta, chi2 = best
    chi2v, coef = chi2_of(theta)
    A = coef[0]
    # errors from the numerical Hessian of chi2 (Delta chi2 = 1)
    k = len(theta); H = np.zeros((k, k)); eps = 1e-3
    for i in range(k):
        for j in range(k):
            ei = np.zeros(k); ej = np.zeros(k); ei[i] = eps; ej[j] = eps
            H[i, j] = (chi2_of(theta + ei + ej)[0] - chi2_of(theta + ei - ej)[0]
                       - chi2_of(theta - ei + ej)[0] + chi2_of(theta - ei - ej)[0]) / (4 * eps * eps)
    try:
        cov = 2 * np.linalg.inv(H)
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(k, np.nan)
    t0, ltr, ltd = theta[:3]
    out = {"t0_ns": float(t0), "t0_err": float(err[0]),
           "tau_rise_ns": float(np.exp(ltr)), "tau_rise_err": float(np.exp(ltr) * err[1]),
           "tau_decay_ns": float(np.exp(ltd)), "tau_decay_err": float(np.exp(ltd) * err[2]),
           "amplitude": float(A), "baseline_adc": float(coef[1]), "baseline_slope_adc_per_ns": float(coef[2]),
           "chi2": float(chi2v), "ndf": int(len(y) - k - 3)}
    extra = ()
    if two_components:
        f2 = 1 / (1 + np.exp(-theta[3])); td2 = min(np.exp(theta[4]), MAX_TAU_SLOW)
        out.update({"slow_fraction": float(f2), "slow_fraction_err": float(f2 * (1 - f2) * err[3]),
                    "tau_slow_ns": float(td2), "tau_slow_err": float(td2 * err[4])})
        extra = (f2, td2)
    n = model_n(t, t0, np.exp(ltr), np.exp(ltd), *extra)
    return out, A * convolve(n, h) + basis @ coef[1:], A * n


# ------------------------------------------------------------- NNLS
def nnls_unfold(t, y, sigma, h, max_iter=MAX_NNLS_ITER):
    """n >= 0 minimising chi2 of y = n conv h, by projected gradient
    (Landweber) with the step from the largest eigenvalue, which makes the
    chi2 decrease monotonically. The baseline is taken once from the
    pre-pulse samples; re-solving it every step fights the gradient and
    stalls. Stops at chi2/ndf = 1 or after max_iter. Returns n, iterations
    and chi2/ndf."""
    w = 1 / sigma ** 2
    H = lambda n: convolve(n, h)
    HT = lambda r: np.convolve(r, h[::-1], mode="full")[len(h) - 1:len(h) - 1 + len(r)]
    v = np.random.default_rng(0).normal(size=len(y))
    for _ in range(50):
        v = HT(w * H(v)); v /= np.linalg.norm(v)
    alpha = 1.0 / float(v @ HT(w * H(v)))
    n = np.zeros_like(y)
    ndf = len(y) - 1
    chi2 = np.inf
    for k in range(1, max_iter + 1):
        r = y - H(n)
        chi2 = float((w * r * r).sum())
        if chi2 / ndf <= 1.0:
            return n, k, chi2 / ndf
        n = np.clip(n + alpha * HT(w * r), 0, None)
    return n, max_iter, chi2 / ndf


def decay_time_from_curve(t, n):
    """Exponential fit to the tail of n(t): from 60% of the peak on the way down
    to 5%. Returns tau, or nan if there is no tail to fit."""
    pk = int(n.argmax()); a = n[pk]
    if a <= 0:
        return np.nan
    after = n[pk:]
    lo = np.where(after < 0.6 * a)[0]; hi = np.where(after < 0.05 * a)[0]
    if len(lo) == 0 or len(hi) == 0 or hi[0] - lo[0] < 4:
        return np.nan
    seg = after[lo[0]:hi[0]]; tt = t[pk + lo[0]:pk + hi[0]]
    ok = seg > 0
    if ok.sum() < 4:
        return np.nan
    slope = np.polyfit(tt[ok], np.log(seg[ok]), 1)[0]
    return float(-1 / slope) if slope < 0 else np.nan


def width_rms(t, n):
    n = np.clip(n, 0, None); s = n.sum()
    if s <= 0:
        return np.nan
    m = (t * n).sum() / s
    return float(np.sqrt(((t - m) ** 2 * n).sum() / s))


# ------------------------------------------------------------- main
def main():
    args = parse_args()
    paths = get_run_paths(args.run)
    for p in paths.values():
        os.makedirs(p, exist_ok=True)
    with open(args.channels) as f:
        chmap = json.load(f)
    fpath = os.path.join(paths["root"], "drs_profiles.root")
    if not os.path.exists(fpath):
        sys.exit(f"ERROR: {fpath} not found; run scripts/check_drs_mcp.py --run {args.run} "
                 f"--channels {args.channels} --mcp-clean first")
    f = ROOT.TFile(fpath, "READ")
    profiles = {lab: read_profile(f, ch) for lab, ch in chmap.items()}
    profiles = {k: v for k, v in profiles.items() if v is not None and k not in EXCLUDE}

    board_of = {lab: chmap[lab].split("_Group")[0] for lab in profiles}
    boards = sorted(set(board_of.values()))
    responses, controls = {}, {}
    for b in boards:
        cher = [l for l in profiles if board_of[l] == b and l.split("_")[0] in CHERENKOV
                and profiles[l][0].max() > MIN_RESPONSE_PEAK]
        if len(cher) < 2:
            print(f"{b}: only {len(cher)} usable Cherenkov channel(s); cannot build a response, "
                  f"its channels are skipped")
            continue
        controls[b] = cher[-1]                        # held out, unfolded as a check
        resp_labels = cher[:-1]
        h, h_spread, onset = build_response(profiles, resp_labels)
        responses[b] = (h, resp_labels)
        print(f"{b}: h(t) from {', '.join(resp_labels)}; control {controls[b]}; "
              f"FWHM {(h > h.max() / 2).sum() * DT:.2f} ns")

    results = {"run": args.run,
               "responses": {b: {"from": r[1], "control": controls[b],
                                 "fwhm_ns": float((r[0] > r[0].max() / 2).sum() * DT)}
                             for b, r in responses.items()},
               "channels": {}}
    curves = {}
    print(f"\n{'channel':16s} {'kind':8s} {'tau_rise':>9s} {'tau_decay':>13s} {'slow comp.':>22s} "
          f"{'chi2/ndf':>9s}  {'NNLS: it':>8s} {'chi2':>5s} {'tau':>6s} {'fwhm':>5s}")
    for lab in profiles:
        b = board_of[lab]
        if b not in responses or lab in responses[b][1] or profiles[lab][0].max() < 10:
            continue
        h = responses[b][0]
        y_full, e_full, noise = profiles[lab]
        pk = int(y_full.argmax())
        lo, hi = pk - WIN_PRE, pk + WIN_POST
        y, e = y_full[lo:hi].copy(), error_model(e_full[lo:hi], y_full[lo:hi], noise)
        y -= np.median(y[:12])                        # baseline: the -6..-3.6 ns samples
        t = (np.arange(len(y)) - WIN_PRE) * DT                # 0 = profile peak (fit coordinate)
        t_abs = t + pk * DT                                   # the same samples on the ts_mcp axis
        kind = "control" if lab == controls[b] else "target"
        fit1, m1, n1 = fit_profile(t, y, e, h, two_components=False)
        fit2, m2, n2 = fit_profile(t, y, e, h, two_components=True)
        use2 = ((fit1["chi2"] - fit2["chi2"]) > 25 and 0.02 < fit2.get("slow_fraction", 0) < 0.9
                and fit2.get("tau_slow_ns", 0) < 0.95 * MAX_TAU_SLOW)
        fit, model, n_fit = (fit2, m2, n2) if use2 else (fit1, m1, n1)
        n_rl, iters, chi2ndf_rl = nnls_unfold(t, y, e, h)
        tau_rl = decay_time_from_curve(t, n_rl)
        w_rl = width_rms(t, n_rl)
        fit["light_onset_ns"] = float(fit["t0_ns"] + pk * DT)          # absolute, ts_mcp axis
        fit["area_adc_ns"] = float(fit["amplitude"] * DT)             # sum of n(t) = pulse area
        res = {"kind": kind, "board": b, "peak_adc": float(y.max()),
               "peak_time_ns": float(pk * DT), "fit": fit, "two_components": bool(use2),
               "nnls": {"iterations": int(iters), "chi2_ndf": float(chi2ndf_rl),
                        "tau_decay_ns": tau_rl, "rms_width_ns": w_rl,
                        "fwhm_ns": float((n_rl > n_rl.max() / 2).sum() * DT),
                        "area_adc_ns": float(n_rl.sum() * DT),
                        "peak_time_ns": float(t_abs[int(n_rl.argmax())]),
                        "mean_time_ns": float((t_abs * n_rl).sum() / n_rl.sum()) if n_rl.sum() > 0 else None}}
        results["channels"][lab] = res
        curves[lab] = (t_abs, y, e, model, n_fit, n_rl)
        slow = (f"{fit['slow_fraction']:.0%} at {fit['tau_slow_ns']:.1f}±{fit['tau_slow_err']:.1f}"
                if use2 else "-")
        print(f"{lab:16s} {kind:8s} {fit['tau_rise_ns']:5.2f}±{min(fit['tau_rise_err'], 99):4.2f} "
              f"{fit['tau_decay_ns']:6.2f}±{min(fit['tau_decay_err'], 99):5.2f} {slow:>22s} "
              f"{fit['chi2'] / fit['ndf']:9.2f}  {iters:8d} {chi2ndf_rl:5.2f} {tau_rl:6.2f} "
              f"{res['nnls']['fwhm_ns']:5.2f}")
    print("  (tau in ns; NNLS tau from an exponential fit to the unfolded tail; controls should be spikes)")
    print(f"\n{'channel':16s} {'light onset':>11s} {'NNLS peak':>9s} {'NNLS mean':>9s}   {'area fit':>9s} {'area NNLS':>9s}   [ns on the ts_mcp axis; area in ADC·ns]")
    for lab, r in results["channels"].items():
        print(f"{lab:16s} {r['fit']['light_onset_ns']:11.2f} {r['nnls']['peak_time_ns']:9.2f} "
              f"{r['nnls']['mean_time_ns'] or float('nan'):9.2f}   {r['fit']['area_adc_ns']:9.0f} {r['nnls']['area_adc_ns']:9.0f}")

    # the unfolded n(t) on the ts_mcp axis, as histograms: the overlay pages and
    # anything else that compares channels read these
    unfolded = {"fit": [], "nnls": []}
    for lab, (t_abs, y, e, model, n_fit, n_rl) in curves.items():
        for key, n in (("fit", n_fit), ("nnls", n_rl)):
            hh = ROOT.TH1D(f"unfolded_{key}_{chmap[lab]}", f"{lab} unfolded ({key});t_{{mcp}} [ns];n(t) [ADC / 0.2 ns]",
                           NBINS, 0.0, NBINS * DT)
            hh.SetDirectory(0)
            for i, v in enumerate(n):
                hh.SetBinContent(int(round(t_abs[i] / DT)) + 1, float(v))
            unfolded[key].append((hh, lab))
    save_hists_to_file([hh for key in unfolded for hh, _ in unfolded[key]],
                       os.path.join(paths["root"], "drs_deconvolution.root"))

    with open(os.path.join(paths["root"], "drs_deconvolution.json"), "w") as fo:
        json.dump(results, fo, indent=1)

    # ---------------------------------------------------------- plots
    resp_desc = "; ".join(f"{b}: {', '.join(r[1])} (control {controls[b]})" for b, r in responses.items())
    pm = PlotManager(paths["root"], paths["plots"], paths["html"], args.run, use_jsroot=args.jsroot,
                     selection_text=(
                         f"**Pulse unfolding.** h(t) is the peak-aligned mean of the mcp_clean profiles "
                         f"of the Cherenkov channels on the same DRS board (prompt light, so the pulse "
                         f"is the response) — {resp_desc}. A response from another board carries the "
                         f"residual inter-board jitter and does not fit. Each other channel's "
                         f"profile is unfolded two ways. **fit**: exponential rise and decay (a "
                         f"second decay only where it lowers chi2 by more than 25) convolved with "
                         f"h(t) and fitted with chi2 on the profile errors. **NNLS**: n(t) as a free "
                         f"non-negative histogram, projected-gradient least squares stopped at "
                         f"chi2/ndf = 1. Both in a -6..+12 ns window (baseline from the pre-pulse "
                         f"samples; the fit also carries a linear-baseline nuisance), because the "
                         f"bright scintillating pulses droop beyond that in a way the Cherenkov "
                         f"response does not describe. Time 0 is the profile peak. One Cherenkov "
                         f"channel per board is kept out of h(t) as a control: it should unfold to "
                         f"a spike."))
    pm.set_output_dir("DRS_Deconvolution")
    STYLE = PlotStyle(W_ref=900, H_ref=600, dology=False, drawoptions="HIST", mycolors=[1, ROOT.kRed + 1, ROOT.kAzure + 1],
                      addOverflow=False, addUnderflow=False, legendPos=[0.50, 0.70, 0.90, 0.90],
                      legendoptions="L")
    keep = []
    # the responses, one per board, overlaid
    hs = []
    for i, (b, (h, labs_b)) in enumerate(responses.items()):
        hh = ROOT.TH1D(f"resp_{i}", "", len(h), -0.5 * DT, (len(h) - 0.5) * DT); hh.SetDirectory(0)
        for j, v in enumerate(h):
            hh.SetBinContent(j + 1, v)
        hs.append(hh)
    keep += hs
    pm.plot_1d(hs, "response_h", "t since onset [ns]", (0, 20), "h(t), unit area",
               (min(h.min() for h, _ in responses.values()) * 1.5, max(h.max() for h, _ in responses.values()) * 1.2),
               legends=[f"{b}: FWHM {results['responses'][b]['fwhm_ns']:.2f} ns" for b in responses], style=STYLE)
    pm.add_newline()
    for lab, (t, y, e, model, n_fit, n_rl) in curves.items():
        nb = len(t); x0, x1 = t[0] - DT / 2, t[-1] + DT / 2      # t is absolute here
        hy = ROOT.TH1D(f"y_{lab}", "", nb, x0, x1); hm = ROOT.TH1D(f"m_{lab}", "", nb, x0, x1)
        for i in range(nb):
            hy.SetBinContent(i + 1, y[i]); hy.SetBinError(i + 1, e[i]); hm.SetBinContent(i + 1, model[i])
        for hh_ in (hy, hm):
            hh_.SetDirectory(0)
        keep += [hy, hm]
        ymax = max(y.max(), model.max()) * 1.25
        r = results["channels"][lab]
        pm.plot_1d([hy, hm], f"{lab}_pulse_fit", "t_{mcp} [ns]", (x0, x1), "mean ADC",
                   (-0.1 * ymax, ymax), legends=["profile", "fit: n(t) #otimes h(t)"],
                   style=STYLE, extra_text=f"{lab} ({r['kind']}, #chi^{{2}}/ndf {r['fit']['chi2'] / r['fit']['ndf']:.1f})")
        hf = ROOT.TH1D(f"nf_{lab}", "", nb, x0, x1); hr = ROOT.TH1D(f"nr_{lab}", "", nb, x0, x1)
        for i in range(nb):
            hf.SetBinContent(i + 1, n_fit[i]); hr.SetBinContent(i + 1, n_rl[i])
        nmax = max(n_fit.max(), n_rl.max()) or 1
        for hh_ in (hf, hr):
            hh_.SetDirectory(0)
        keep += [hf, hr]
        r = results["channels"][lab]
        txt = (f"#tau_{{d}} = {r['fit']['tau_decay_ns']:.2f} ns (fit), {r['nnls']['tau_decay_ns']:.2f} (NNLS)"
               if np.isfinite(r["nnls"]["tau_decay_ns"]) else f"#tau_{{d}} = {r['fit']['tau_decay_ns']:.2f} ns (fit)")
        pm.plot_1d([hf, hr], f"{lab}_unfolded", "t_{mcp} [ns]", (x0, x1), "n(t) [ADC / 0.2 ns]",
                   (-0.1 * nmax, 1.25 * nmax), legends=["fit model", "NNLS"],
                   style=PlotStyle(W_ref=900, H_ref=600, dology=False, drawoptions="HIST", mycolors=[ROOT.kRed + 1, ROOT.kAzure + 1],
                                   addOverflow=False, addUnderflow=False,
                                   legendPos=[0.50, 0.74, 0.90, 0.90], legendoptions="L"),
                   extra_text=lab)
        pm.add_newline()
    # summary: every channel's n(t) on the common time axis, absolute and peak = 1
    xr = {}
    for key in ("nnls", "fit"):
        hs = [hh for hh, _ in unfolded[key]]
        labs = [lab for _, lab in unfolded[key]]
        lo = min(hh.GetBinLowEdge(hh.FindFirstBinAbove(0)) for hh in hs) - 2
        hi = max(hh.GetBinLowEdge(hh.FindLastBinAbove(0) + 1) for hh in hs) + 2
        xr[key] = (lo, hi)
        pal = OVERLAY_COLOURS[:len(hs)]
        style = PlotStyle(dology=False, drawoptions="HIST", mycolors=pal, addOverflow=False,
                          addUnderflow=False, legendPos=[0.55, 0.50, 0.90, 0.90], legendoptions="L",
                          W_ref=900, H_ref=600)
        ymax = max(hh.GetMaximum() for hh in hs) * 1.2
        pm.plot_1d(hs, f"summary_unfolded_{key}_absolute", "t_{mcp} [ns]", (lo, hi),
                   "n(t) [ADC / 0.2 ns]", (-0.05 * ymax, ymax), legends=labs, style=style,
                   extra_text=f"{key.upper()}, absolute scale", prepend=True)
        hn = []
        for hh in hs:
            c = hh.Clone(hh.GetName() + "_norm"); c.SetDirectory(0)
            if c.GetMaximum() > 0:
                c.Scale(1.0 / c.GetMaximum())
            hn.append(c)
        keep += hn
        pm.plot_1d(hn, f"summary_unfolded_{key}_peak1", "t_{mcp} [ns]", (lo, hi),
                   "n(t) / peak", (-0.05, 1.25), legends=labs, style=style,
                   extra_text=f"{key.upper()}, peak = 1", prepend=True)
    html = pm.generate_html("DRS/DRS_Deconvolution.html", plots_per_row=2,
                            title=f"DRS pulse unfolding, run {args.run}")
    print(f"\n{html}")
    for key, page in (("nnls", "DRS/DRS_Unfolded_overlay.html"), ("fit", "DRS/DRS_Unfolded_fit_overlay.html")):
        out = build_channel_overlay(
            unfolded[key], output_html=os.path.join(paths["html"], page),
            xlabel="t_mcp [ns]", ylabel="n(t) [ADC / 0.2 ns]",
            title=f"Unfolded light n(t), {key.upper()} — Run {args.run}",
            intro_text=("Photon arrival profiles unfolded from the mcp_clean profiles with the "
                        "Cherenkov response of the same DRS board. Time is the analysis's ts_mcp "
                        "axis (group reference + MCP aligned), so offsets between channels are "
                        "light-propagation and path differences; the scale is ADC per 0.2 ns with "
                        "h(t) of unit area, so areas compare photons captured. Use Peak = 1 to "
                        "compare shapes."),
            filename=f"DRS_unfolded_{key}_Run{args.run}", xrange=xr[key])
        if out:
            print(out)


if __name__ == "__main__":
    main()
