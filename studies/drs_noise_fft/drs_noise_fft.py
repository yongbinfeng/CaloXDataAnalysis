"""Frequency-domain look at DRS waveforms: noise spectrum and matched-filter timing.

Reads raw waveforms for the channels in a --channels file plus one MCP, and
produces three things from the pre-pulse (noise-only) samples:

  1. the noise power spectrum per channel, its split by frequency band, and
     the coherence between channels (common-mode pickup would be coherent);
  2. the DRS4 fixed pattern in *cell* space, found by realigning every event
     by its StartIndexCell, and how much noise it accounts for;
  3. a matched-filter time estimate versus the 20% CFD, using the MCP-aligned
     profiles from drs_profiles.root as templates (skipped if that file is
     missing -- run scripts/check_drs_mcp.py first).

Writes PNG plots and results.json next to this script by default, and prints
the tables. Nothing here changes the analysis; it is a study.

    python3 studies/drs_noise_fft/drs_noise_fft.py --run 1994 \
        --channels data/channel_maps/testingfibers.json
"""
import argparse
import json
import os
import sys

import numpy as np
import ROOT

from channels.channel_map import get_mcp_channels, get_mcp_reference
from utils.plot_helper import get_run_paths

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DT_NS = 0.2                         # DRS at 5 GSPS
NOISE_WIN = (20, 380)               # pre-pulse samples: baseline only
SIGNAL_WIN = (380, 990)             # where the pulses are; 990 avoids the last-cell dip
CFD_FRACTION = 0.20                 # as in compute_cfd_integral
BANDS = [(0, 100), (100, 200), (200, 400), (400, 600), (600, 1000), (1000, 2500)]
AMP_BINS = [(20, 40), (40, 80), (80, 160), (160, 400), (400, 2500)]
CORE_NS = 2.0                       # |dt - median| below this counts as "in core"


def parse_args():
    from configs.run_config import run_number, jsonFile
    here = os.path.dirname(os.path.abspath(__file__))
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run", type=int, default=run_number)
    p.add_argument("--channels", required=True, metavar="FILE",
                   help="JSON mapping labels to DRS channel names")
    p.add_argument("--json-file", default=jsonFile, help="run -> ROOT file map")
    p.add_argument("--nevents", type=int, default=4000)
    p.add_argument("--mcp", default=None,
                   help="MCP channel for the timing reference (default: the one "
                        "the analysis uses, channels.maps.services.get_mcp_reference)")
    p.add_argument("--outdir", default=here)
    return p.parse_args()


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_waveforms(args, cols):
    with open(args.json_file) as f:
        files = json.load(f)[str(args.run)]
    files = files if isinstance(files, list) else [files]
    chain = ROOT.TChain("EventTree")
    for fn in files:
        if not os.path.exists(fn):
            sys.exit(f"ERROR: {fn} not found")
        chain.Add(fn)
    rdf = ROOT.RDataFrame(chain).Range(args.nevents)
    raw = rdf.AsNumpy(columns=cols)
    out = {}
    for c in cols:
        a = raw[c]
        if a.dtype == object:
            out[c] = np.stack([np.asarray(v, dtype=np.float64) for v in a])
        else:
            out[c] = np.asarray(a)
    return out


def pipeline_flip_list(run_number):
    """Channels the analysis inverts (get_drs_branches_to_flip), so the study
    sees the same polarity the analysis does."""
    from channels.channel_map import build_drs_boards
    from variables.drs import get_drs_branches_to_flip
    boards = build_drs_boards(run_number)
    refs = [c.get_channel_name() for b in boards.values() for c in b if c.is_reference]
    return set(get_drs_branches_to_flip(run_number, drs_channels_ref=refs, drsboards=boards))


def baseline_subtract(W):
    return W - np.median(W[:, NOISE_WIN[0]:NOISE_WIN[1]], axis=1, keepdims=True)


def clean_noise_segments(W):
    """Pre-pulse samples, mean removed, events with a stray pulse dropped."""
    seg = W[:, NOISE_WIN[0]:NOISE_WIN[1]]
    seg = seg - seg.mean(axis=1, keepdims=True)
    ok = np.abs(seg).max(axis=1) < 6 * seg.std(axis=1).mean()
    return seg, ok


# --------------------------------------------------------------------------
# 1. noise spectrum
# --------------------------------------------------------------------------
def noise_spectra(W, labels):
    n = NOISE_WIN[1] - NOISE_WIN[0]
    win = np.hanning(n)
    freq = np.fft.rfftfreq(n, DT_NS * 1e-9) / 1e6          # MHz
    res, spectra = {}, {}
    for c, lab in labels.items():
        seg, ok = clean_noise_segments(W[c])
        X = np.fft.rfft(seg[ok] * win, axis=1)
        P = (np.abs(X) ** 2).mean(axis=0)
        P[0] = 0.0
        spectra[c] = X
        tot = P.sum()
        res[lab] = {
            "sigma": float(seg[ok].std()),
            # ignore the Nyquist bin: alternate-cell offsets put a spike there
            "peak_MHz": float(freq[:-1][P[:-1].argmax()]),
            "band_fraction": {f"{a}-{b}": float(P[(freq >= a) & (freq < b)].sum() / tot)
                              for a, b in BANDS},
            "psd": (P / P.max()).tolist(),
        }
    return freq, res, spectra


def coherence_summary(spectra, labels, freq):
    """Magnitude-squared coherence per pair, grouped by hardware relation."""
    def board_group(c):
        parts = c.split("_")
        return "_".join(parts[1:3]), parts[3]
    band = (freq > 10) & (freq < 500)
    groups = {"same group": [], "same board": [], "other board": []}
    cols = list(labels)
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            Xa, Xb = spectra[a], spectra[b]
            m = min(len(Xa), len(Xb))
            Xa, Xb = Xa[:m], Xb[:m]
            coh = (np.abs((Xa * np.conj(Xb)).mean(0)) ** 2
                   / ((np.abs(Xa) ** 2).mean(0) * (np.abs(Xb) ** 2).mean(0)))
            (ba, ga), (bb, gb) = board_group(a), board_group(b)
            key = ("same group" if (ba, ga) == (bb, gb)
                   else "same board" if ba == bb else "other board")
            groups[key].append(float(coh[band].mean()))
    return {k: {"pairs": len(v), "mean_coherence": float(np.mean(v)) if v else None}
            for k, v in groups.items()}


# --------------------------------------------------------------------------
# 2. DRS cell pattern
# --------------------------------------------------------------------------
def cell_patterns(W, start_cells, labels):
    idx = np.arange(NOISE_WIN[0], NOISE_WIN[1])
    res, patterns = {}, {}
    for c, lab in labels.items():
        seg, ok = clean_noise_segments(W[c])
        seg, start = seg[ok], start_cells[c][ok]
        cells = (start[:, None] + idx[None, :]) % 1024
        cnt = np.bincount(cells.ravel(), minlength=1024)
        acc = np.bincount(cells.ravel(), weights=seg.ravel(), minlength=1024)
        pat = np.where(cnt > 0, acc / np.maximum(cnt, 1), 0.0)
        resid = seg - pat[cells]
        s0, s1 = seg.std(), resid.std()
        patterns[c] = pat
        res[lab] = {"sigma": float(s0), "pattern_rms": float(pat.std()),
                    "sigma_after": float(s1),
                    "variance_removed": float(1 - s1 ** 2 / s0 ** 2),
                    "min_entries_per_cell": int(cnt.min())}
    return res, patterns


# --------------------------------------------------------------------------
# 3. timing: CFD versus matched filter
# --------------------------------------------------------------------------
def cfd_time(w, lo, hi, frac=CFD_FRACTION):
    """Leading-edge crossing at frac*peak, linear interpolation. (time, peak)"""
    seg = w[lo:hi]
    p = int(seg.argmax())
    a = seg[p]
    if a <= 0:
        return np.nan, a
    thr = frac * a
    i = p
    while i > 0 and seg[i] > thr:
        i -= 1
    if i >= p:
        return np.nan, a
    y0, y1 = seg[i], seg[i + 1]
    return lo + i + ((thr - y0) / (y1 - y0) if y1 != y0 else 0.0), a


def make_template(profile, pre=40, post=120):
    """Unit-norm template around the profile peak, or None if there is no pulse.

    A template without a real pulse makes the matched filter lock onto a fixed
    lag regardless of the data; the peak-to-RMS check catches that.
    """
    y = np.array([profile.GetBinContent(i) for i in range(1, profile.GetNbinsX() + 1)])
    y -= np.median(y)
    pp = int(y.argmax())
    if y[pp] < 5 * y[:max(pp - 60, 1)].std() or pp - pre < 0 or pp + post > len(y):
        return None
    t = y[pp - pre:pp + post]
    t -= t.mean()
    return t / np.linalg.norm(t), pre


def mf_time(w, tmpl, pre, lo, hi):
    """Lag of maximum correlation with the template, parabolic sub-sample."""
    L = len(tmpl)
    seg = w[lo - pre:hi + (L - pre)]
    corr = np.correlate(seg, tmpl, mode="valid")
    k = int(corr.argmax())
    d = 0.0
    if 0 < k < len(corr) - 1:
        y0, y1, y2 = corr[k - 1], corr[k], corr[k + 1]
        den = y0 - 2 * y1 + y2
        d = 0.5 * (y0 - y2) / den if den else 0.0
    edge = k == 0 or k == len(corr) - 1
    return lo + k + d, edge


def robust_sigma(x):
    return 0.7413 * (np.percentile(x, 75) - np.percentile(x, 25))


def led_time(w_flipped, min_amp=500.0, frac=0.5):
    """process_dynamic_led on a flipped reference channel: 50% crossing, integer."""
    p = int(w_flipped.argmax())
    a = w_flipped[p]
    if a < min_amp:
        return np.nan
    i = p
    while i > 0 and w_flipped[i] > frac * a:
        i -= 1
    return float(i + 1)


def timing_study(W, labels, mcp, profiles_path, ref_ts, ref_of):
    if not os.path.exists(profiles_path):
        print(f"\n{profiles_path} not found: skipping the timing study "
              f"(run scripts/check_drs_mcp.py --run <run> first).")
        return None
    prof = ROOT.TFile(profiles_path, "READ")
    nev = W[mcp].shape[0]
    t_mcp, a_mcp = np.array([cfd_time(-W[mcp][e], *SIGNAL_WIN) for e in range(nev)]).T
    # correct every time by its group reference, as the analysis does: this is
    # what removes the ~1.3 ns jitter between DRS boards
    t_mcp = t_mcp - ref_ts[ref_of[mcp]]
    has = (a_mcp > 50) & np.isfinite(t_mcp)
    ev = np.where(has)[0]
    print(f"\nevents with an MCP pulse (>50 ADC): {has.sum()} of {nev}; "
          f"MCP arrival spread {robust_sigma(t_mcp[has]) * DT_NS:.2f} ns")
    res = {"_mcp_events": int(has.sum()),
           "_mcp_spread_ns": float(robust_sigma(t_mcp[has]) * DT_NS)}
    for c, lab in labels.items():
        h = prof.Get(f"prof_{c}_blsub_VS_ts_mcp")
        tmpl = make_template(h) if h else None
        if tmpl is None:
            res[lab] = {"skipped": "no usable template (no pulse in the profile)"}
            continue
        t, pre = tmpl
        tc, amp = np.array([cfd_time(W[c][e], *SIGNAL_WIN) for e in ev]).T
        tm, edge = zip(*(mf_time(W[c][e], t, pre, *SIGNAL_WIN) for e in ev))
        tm, edge = np.array(tm), np.array(edge)
        ref_c = ref_ts[ref_of[c]][ev]
        dc = (tc - ref_c - t_mcp[ev]) * DT_NS
        dm = (tm - ref_c - t_mcp[ev]) * DT_NS
        bins = {}
        for lo, hi in AMP_BINS:
            s = (amp >= lo) & (amp < hi) & np.isfinite(dc) & ~edge
            if s.sum() < 15:
                continue
            row = {"n": int(s.sum())}
            for name, d in (("cfd", dc[s]), ("mf", dm[s])):
                core = np.abs(d - np.median(d)) < CORE_NS
                row[name] = {"core_sigma_ns": float(robust_sigma(d[core])) if core.sum() > 5 else None,
                             "core_fraction": float(core.mean())}
            bins[f"{lo}-{hi}"] = row
        res[lab] = {"edge_lock_fraction": float(edge.mean()), "bins": bins}
    return res


# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------
def _colour(i, n):
    return ROOT.TColor.GetColorPalette(int(i * 250 / max(n, 1)))


def plot_psd(freq, spec, mcp_label, outdir, keep):
    cv = ROOT.TCanvas("psd", "", 1400, 600)
    cv.Divide(2, 1)
    labs = list(spec)
    # left: normalised PSD, every channel
    cv.cd(1).SetLogx(); cv.cd(1).SetLogy()
    mg = ROOT.TMultiGraph()
    for i, lab in enumerate(labs):
        p = np.array(spec[lab]["psd"])[1:]
        g = ROOT.TGraph(len(p), freq[1:], p)
        special = lab == mcp_label
        g.SetLineColor(ROOT.kBlack if special else _colour(i, len(labs)))
        g.SetLineWidth(3 if special else 1)
        mg.Add(g, "l"); keep.append(g)
    mg.SetTitle("Noise PSD, pre-pulse samples, each channel scaled to its own maximum;"
                "frequency [MHz];PSD (a.u.)")
    mg.Draw("a"); mg.GetXaxis().SetLimits(10, 2500); mg.SetMinimum(2e-3); mg.SetMaximum(2)
    tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.035)
    tx.DrawLatex(0.15, 0.20, "thick black: MCP (unamplified)"); keep.append(tx)
    # right: variance by band, stacked
    cv.cd(2)
    nb = len(BANDS)
    hs = ROOT.THStack("bands", "Share of noise variance per frequency band;;fraction")
    cols = [ROOT.kBlue - 9, ROOT.kAzure + 1, ROOT.kRed + 1, ROOT.kOrange + 7,
            ROOT.kGreen + 2, ROOT.kGray + 1]
    for j, (a, b) in enumerate(BANDS):
        h = ROOT.TH1D(f"band{j}", f"{a}-{b} MHz", len(labs), 0, len(labs))
        for i, lab in enumerate(labs):
            h.SetBinContent(i + 1, spec[lab]["band_fraction"][f"{a}-{b}"])
            h.GetXaxis().SetBinLabel(i + 1, lab)
        h.SetFillColor(cols[j]); h.SetLineColor(ROOT.kWhite)
        hs.Add(h); keep.append(h)
    hs.Draw("hist"); hs.GetXaxis().LabelsOption("v"); hs.SetMaximum(1.0)
    cv.cd(2).SetBottomMargin(0.28)
    lg = cv.cd(2).BuildLegend(0.12, 0.72, 0.42, 0.90); lg.SetNColumns(2); keep.append(lg)
    keep.append(hs)
    cv.SaveAs(os.path.join(outdir, "noise_psd.png"))


def plot_cell_pattern(patterns, cellres, labels, outdir, keep):
    cv = ROOT.TCanvas("cell", "", 1400, 600)
    cv.Divide(2, 1)
    cv.cd(1)
    lab0 = next(iter(labels.values()))
    c0 = next(iter(labels))
    g = ROOT.TGraph(1024, np.arange(1024, dtype=float), patterns[c0])
    g.SetTitle(f"DRS fixed pattern in cell space, {lab0};cell;mean offset [ADC]")
    g.SetLineColor(ROOT.kBlue + 1); g.Draw("al"); keep.append(g)
    cv.cd(2)
    labs = list(cellres)
    h0 = ROOT.TH1D("s0", "Noise sigma before / after subtracting the cell pattern;;sigma [ADC]",
                   len(labs), 0, len(labs))
    h1 = ROOT.TH1D("s1", "", len(labs), 0, len(labs))
    for i, lab in enumerate(labs):
        h0.SetBinContent(i + 1, cellres[lab]["sigma"]); h0.GetXaxis().SetBinLabel(i + 1, lab)
        h1.SetBinContent(i + 1, cellres[lab]["sigma_after"])
    h0.SetFillColor(ROOT.kGray + 1); h0.SetLineColor(ROOT.kGray + 1)
    h1.SetFillColor(ROOT.kAzure + 1); h1.SetLineColor(ROOT.kAzure + 1)
    h0.SetMinimum(0); h0.SetMaximum(10)
    h0.Draw("hist"); h1.Draw("hist same"); h0.GetXaxis().LabelsOption("v")
    cv.cd(2).SetBottomMargin(0.28)
    lg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
    lg.AddEntry(h0, "raw", "f"); lg.AddEntry(h1, "cell pattern subtracted", "f"); lg.Draw()
    keep += [h0, h1, lg]
    cv.SaveAs(os.path.join(outdir, "cell_pattern.png"))


def plot_timing(timing, outdir, keep):
    cv = ROOT.TCanvas("tim", "", 1400, 600)
    cv.Divide(2, 1)
    # one channel per fibre type, the one with the cleanest high-amplitude
    # timing, so the plot spans the types rather than the first six in the file
    top = f"{AMP_BINS[-1][0]}-{AMP_BINS[-1][1]}"
    best = {}
    for lab, r in timing.items():
        if lab.startswith("_"):
            continue
        b = r.get("bins", {}).get(top)
        if not b:
            continue
        kind = lab.split("_")[0]
        if kind not in best or b["cfd"]["core_fraction"] > best[kind][1]:
            best[kind] = (lab, b["cfd"]["core_fraction"])
    show = [v[0] for v in sorted(best.values(), key=lambda x: -x[1])][:6]
    MIN_CORE = 0.20                     # sigma is meaningless below this
    for pad, (title, key) in enumerate([("core sigma", "core_sigma_ns"),
                                        ("fraction of events in core (|dt| < 2 ns)", "core_fraction")]):
        cv.cd(pad + 1).SetLogx()
        if key == "core_sigma_ns":
            cv.cd(pad + 1).SetLogy()
        mg = ROOT.TMultiGraph()
        for j, lab in enumerate(show):
            for style, meth, ls in ((20, "cfd", 1), (24, "mf", 2)):
                xs, ys = [], []
                for (lo, hi) in AMP_BINS:
                    b = timing[lab]["bins"].get(f"{lo}-{hi}")
                    if not b or b[meth][key] is None:
                        continue
                    if key == "core_sigma_ns" and b[meth]["core_fraction"] < MIN_CORE:
                        continue
                    xs.append(np.sqrt(lo * hi)); ys.append(b[meth][key])
                if not xs:
                    continue
                g = ROOT.TGraph(len(xs), np.array(xs), np.array(ys, dtype=float))
                g.SetMarkerStyle(style); g.SetMarkerColor(j + 1); g.SetLineColor(j + 1)
                g.SetLineStyle(ls); g.SetTitle(f"{lab} ({'CFD' if meth == 'cfd' else 'matched'})")
                mg.Add(g, "lp"); keep.append(g)
        if key == "core_sigma_ns":
            title += f" (bins with >{MIN_CORE:.0%} of events in core)"
        mg.SetTitle(f"{title};pulse amplitude [ADC];" +
                    ("sigma [ns]" if key == "core_sigma_ns" else "fraction"))
        mg.Draw("a"); mg.GetXaxis().SetLimits(15, 3000)
        if key == "core_sigma_ns":
            mg.SetMinimum(0.1); mg.SetMaximum(5)
        else:
            mg.SetMinimum(0); mg.SetMaximum(1.05)
        lg = cv.cd(pad + 1).BuildLegend(0.12, 0.60, 0.50, 0.88); lg.SetTextSize(0.022)
        keep += [mg, lg]
    cv.SaveAs(os.path.join(outdir, "timing_cfd_vs_matched.png"))


# --------------------------------------------------------------------------
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    with open(args.channels) as f:
        chmap = json.load(f)
    labels = {name: lab for lab, name in chmap.items()}
    mcp_label = get_mcp_reference(args.run)
    if args.mcp:
        mcp, mcp_label = args.mcp, args.mcp
    elif mcp_label is None:
        sys.exit(f"run {args.run} has no MCP reference; pass --mcp <channel>")
    else:
        mcp = get_mcp_channels(args.run)[mcp_label]
    labels_all = dict(labels); labels_all[mcp] = mcp_label
    start_branches = {c: c.rsplit("_", 1)[0] + "_StartIndexCell" for c in labels}
    ref_of = {c: c.rsplit("_", 1)[0] + "_Channel8" for c in labels_all}

    cols = list(labels_all) + sorted(set(start_branches.values()) | set(ref_of.values()))
    data = load_waveforms(args, cols)
    flipped = pipeline_flip_list(args.run)
    for c in labels:                                  # same polarity as the analysis
        if c in flipped:
            data[c] = -data[c]
            print(f"  {labels[c]}: inverted, as the analysis does (it is in the flip list)")
    W = {c: baseline_subtract(data[c]) for c in labels_all}
    starts = {c: data[start_branches[c]].astype(int) for c in labels}
    ref_ts = {}
    for rb in set(ref_of.values()):
        Wr = -baseline_subtract(data[rb])                 # reference channels are flipped
        ref_ts[rb] = np.array([led_time(Wr[e]) for e in range(Wr.shape[0])])
    print(f"run {args.run}: {W[mcp].shape[0]} events, {len(labels)} channels + MCP")

    # 1. spectra
    freq, spec, spectra = noise_spectra(W, labels_all)
    coh = coherence_summary({c: spectra[c] for c in labels_all}, labels_all, freq)
    print(f"\n{'channel':18s} {'sigma':>5s}  " + " ".join(f"{a}-{b}".rjust(9) for a, b in BANDS) + "   peak")
    for lab, r in spec.items():
        print(f"{lab:18s} {r['sigma']:5.2f}  " +
              " ".join(f"{r['band_fraction'][f'{a}-{b}']:8.0%} " for a, b in BANDS) +
              f"  {r['peak_MHz']:4.0f} MHz")
    print("\ncoherence 10-500 MHz (1/N_events is the floor for unrelated noise):")
    for k, v in coh.items():
        print(f"  {k:12s} pairs={v['pairs']:3d}  mean={v['mean_coherence']:.4f}")

    # 2. cell pattern
    cellres, patterns = cell_patterns(W, starts, labels)
    print(f"\n{'channel':18s} {'sigma':>5s} {'pattern':>8s} {'after':>6s} {'removed':>8s}")
    for lab, r in cellres.items():
        print(f"{lab:18s} {r['sigma']:5.2f} {r['pattern_rms']:8.2f} {r['sigma_after']:6.2f} "
              f"{r['variance_removed']:7.0%}")

    # 3. timing
    profiles = os.path.join(get_run_paths(args.run)["root"], "drs_profiles.root")
    timing = timing_study(W, labels, mcp, profiles, ref_ts, ref_of)
    if timing:
        print(f"\n{'channel':16s} " + " ".join(f"{'%d-%d'%b:>21s}" for b in AMP_BINS))
        for lab, r in timing.items():
            if lab.startswith("_"):
                continue
            if "bins" not in r:
                print(f"{lab:16s} skipped: {r['skipped']}"); continue
            row = f"{lab:16s} "
            for lo, hi in AMP_BINS:
                b = r["bins"].get(f"{lo}-{hi}")
                if not b:
                    row += f"{'-':>21s} "; continue
                cell = []
                for m in ("cfd", "mf"):
                    s, fr = b[m]["core_sigma_ns"], b[m]["core_fraction"]
                    cell.append(f"{s if s is not None else float('nan'):4.2f}({fr:3.0%})")
                row += f"{cell[0]}/{cell[1]:>10s} "
            print(row)
        print("  (core sigma [ns] and core fraction: CFD / matched filter)")

    # plots + json
    keep = []
    plot_psd(freq, spec, mcp_label, args.outdir, keep)
    plot_cell_pattern(patterns, cellres, labels, args.outdir, keep)
    if timing:
        plot_timing(timing, args.outdir, keep)
    out = {"run": args.run, "nevents": int(W[mcp].shape[0]), "channels": args.channels,
           "mcp": mcp, "mcp_label": mcp_label, "noise": {k: {kk: vv for kk, vv in v.items() if kk != "psd"}
                                  for k, v in spec.items()},
           "coherence": coh, "cell_pattern": cellres, "timing": timing}
    with open(os.path.join(args.outdir, "results.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nwrote plots and results.json to {args.outdir}")


if __name__ == "__main__":
    main()
