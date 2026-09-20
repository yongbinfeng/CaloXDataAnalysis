"""
In-time vs out-of-time analysis of the test fibres, one run at a time.

For every channel in the --channels map, events are split by the MCP-aligned
CFD time: "in time" means within IN_TIME_WINDOW_TS of the channel's timing MPV
(IN_TIME_WINDOW_SCINT_TS for the scintillating fibres, whose slow pulses have a
wider timing distribution)
(read from results/root/Run{run}/drs_stats.root, so the DRS stats of the run
must exist), everything else is "out of time". The out-of-time sample holds
the noise pedestal of the integral and the SiPM dark counts, the in-time
sample the real hits on top of that; the dark counts also give the
single-photoelectron calibration of every fibre.

Entry points:
  run_intime_analysis(manager, args)   from scripts/check_drs_mcp.py --intime,
                                       after the DRS stats have been written
  scripts/compare/plot_drs_intime_energy_testfibers.py   standalone, same args

Outputs, per run (results/{root,plots,html}/Run{run}/):
  root/drs_intime.root             in/out/all histograms and profiles per channel
  plots/DRS_InTime/                PNGs, intime_summary.json, intime_summary.md
  html/DRS/DRS_InTime.html         static page: summary plots, tables, per-fibre plots
  html/DRS/DRS_InTime_Prof_overlay.html      interactive overlays (tick curves on/off):
  html/DRS/DRS_InTime_Prof_aligned_overlay.html  same, every fibre shifted so its
                                             in-time peak sits at TS = 0 (shape comparison)
  html/DRS/DRS_InTime_Amp_overlay.html       mean waveform / amplitude / energy, with
  html/DRS/DRS_InTime_Energy_overlay.html    in-time, out-of-time and all-event curves
  html/DRS/DRS_InTime_*_Sci_overlay.html     for every fibre; scintillating fibres on
                                             their own pages (different binning)
"""


from html import escape as html_escape
import json
import os
import numpy as np
import ROOT

from core.plot_manager import PlotManager, PlotStyle
from utils.data_loader import getRunInfo
from utils.overlay import build_channel_overlay
from plotting.calox_plot_helper import fiber_base_type, plot_fiber_summary
from utils.plot_helper import get_run_paths
from utils.utils import get_hist_mpv
from variables.drs import get_arr_name

PLOT_SUBDIR = "DRS_InTime"          # under results/plots/Run{run}/
HTML_PAGE = "DRS/DRS_InTime.html"     # under results/html/Run{run}/
ROOT_FILE = "drs_intime.root"         # under results/root/Run{run}/

IN_TIME_WINDOW_TS = 6.0          # |TS_cfd_mcp - MPV| below this counts as in time (Cherenkov fibres)
# The scintillating fibres' CFD times have shoulders out to about +-20 TS (time
# walk on the slow scintillation rise, mostly for the smaller pulses); with 6 TS
# about a quarter of their real hits would land in the out-of-time sample.
IN_TIME_WINDOW_SCINT_TS = 20.0


# Sideband for the random coincidences (dark counts and noise maxima that land
# inside the in-time window by chance): |TS_cfd_mcp - MPV| in [lo, hi) on both
# sides of the peak, far enough from it to contain no real hits, no late
# scintillation pulses and no afterpulses. The random density per TS measured
# there, times the window width, is the expected number of randoms in time.
# Showers of a second particle a few to 20 ns off the trigger time populate
# |dt| up to about 100 TS in the scintillating fibres, so the sideband starts
# beyond that.
SIDEBAND_TS = (100.0, 300.0)


def in_time_window(label):
    """Half-width [TS] of the in-time window for a fibre label."""
    return IN_TIME_WINDOW_SCINT_TS if "scint" in label.lower() else IN_TIME_WINDOW_TS
PE_AMP_SEARCH = (80.0, 200.0)    # 1 p.e. amplitude peak is looked for here [ADC]
PE_INT_SEARCH = (500.0, 1200.0)  # 1 p.e. integral peak is looked for here [ADC]
NOISE_AMP_SEARCH = (0.0, 60.0)   # 0 p.e. (noise) amplitude peak [ADC]

# Labels that are not real test fibres and are left out of the fibre-count
# and correlation study: GradedIndex_* carry no signal, Quartz_1 is the
# upstream MCP channel (MCP_US_1) in the service map.
EXCLUDE_FROM_CORR = {"GradedIndex_0", "GradedIndex_1", "Quartz_1"}
CORR_REFERENCE = "Scintillating_0"   # fibre used as the shower tag in the summary

_STYLE = PlotStyle(dology=True, drawoptions="HIST",
                   mycolors=[ROOT.kBlue + 1, ROOT.kRed + 1],
                   addOverflow=False, addUnderflow=False,
                   legendPos=[0.45, 0.72, 0.88, 0.88])
_STYLE_PROF = PlotStyle(dology=False, drawoptions="HIST",
                        mycolors=[ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kBlack],
                        addOverflow=False, addUnderflow=False,
                        legendPos=[0.45, 0.68, 0.88, 0.88])
PROFILE_TS_RANGE = (440.0, 620.0)   # MCP-aligned TS window shown for the mean waveform


def _ranges(label):
    """(energy range, amplitude range) as (nbins, xmin, xmax), by fibre label."""
    if "scint" in label.lower():
        return (600, -5000.0, 115000.0), (240, 0.0, 2400.0)
    return (360, -1000.0, 8000.0), (200, 0.0, 400.0)


def _timing_mpv(fstats, ch):
    """Timing MPV [TS] from the run's drs_stats.root (fine-binned, else coarse)."""
    for name in (f"hist_{ch}_TS_cfd_mcp_finebins", f"hist_{ch}_TS_cfd_mcp"):
        h = fstats.Get(name)
        if h and h.GetEntries() >= 5:
            return get_hist_mpv(h)[0]
    return None


def _peak_in(h, lo, hi):
    """Bin centre of the highest bin of h between lo and hi."""
    ax = h.GetXaxis()
    b_lo, b_hi = ax.FindBin(lo), ax.FindBin(hi)
    return ax.GetBinCenter(max(range(b_lo, b_hi + 1), key=h.GetBinContent))


def _shift_to_peak(h, b_peak, name):
    """TH1F copy of h (1024 one-TS bins) with bin b_peak moved to x = 0."""
    out = ROOT.TH1F(name, "", 1024, -512.5, 511.5)   # bin i centre = i - 513
    for b in range(1, h.GetNbinsX() + 1):
        i = b - b_peak + 513
        if 1 <= i <= 1024:
            out.SetBinContent(i, h.GetBinContent(b))
            out.SetBinError(i, h.GetBinError(b))
    out.SetDirectory(0)
    return out


def _poisson_binomial(ps):
    """P(N = k) for independent Bernoulli trials with success probabilities ps."""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1.0 - p, p])
    return dist


def _draw_conditional_matrix(labels, cond, outdir, name, lumi):
    """cond[i][j] = P(fibre j in time | fibre i in time); x = given, y = observed."""
    n = len(labels)
    h = ROOT.TH2F(f"h_{name}", "", n, 0.5, n + 0.5, n, 0.5, n + 0.5)
    for i, li in enumerate(labels):
        h.GetXaxis().SetBinLabel(i + 1, li)
        h.GetYaxis().SetBinLabel(i + 1, li)
        for j in range(n):
            h.SetBinContent(i + 1, j + 1, cond[i][j])
    c = ROOT.TCanvas(f"c_{name}", "", 900, 800)
    c.SetLeftMargin(0.20); c.SetBottomMargin(0.22); c.SetRightMargin(0.13); c.SetTopMargin(0.07)
    ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetPaintTextFormat(".2f")
    h.GetXaxis().LabelsOption("v"); h.GetXaxis().SetLabelSize(0.03); h.GetYaxis().SetLabelSize(0.03)
    h.GetXaxis().SetTitle("given: fibre in time"); h.GetXaxis().SetTitleOffset(2.6)
    h.GetYaxis().SetTitle("P(fibre in time | given)"); h.GetYaxis().SetTitleOffset(2.6)
    h.SetMinimum(0.0); h.SetMaximum(1.0); h.SetMarkerSize(0.75)
    h.Draw("COLZ TEXT")
    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.03); lat.SetTextAlign(31)
    lat.DrawLatex(0.87, 0.95, lumi)
    c.SaveAs(os.path.join(outdir, f"{name}.png"))
    c._keep = [h, lat]


def run_intime_analysis(manager, args):
    """Run the analysis on manager.rdf (DRS columns defined, selections applied)
    and write the per-run outputs. Returns the path of the static HTML page."""
    ch_map = json.load(open(args.channels))
    rdf = manager.rdf

    stats_path = f"results/root/Run{args.run}/drs_stats.root"
    if not os.path.exists(stats_path):
        raise SystemExit(f"{stats_path} not found; run check_drs_mcp.py for this run first")
    fstats = ROOT.TFile(stats_path, "READ")

    # ---- book everything lazily, one RDF pass -------------------------------
    booked = {}
    n_all = rdf.Count()
    for label, ch in ch_map.items():
        mpv = _timing_mpv(fstats, ch)
        if mpv is None:
            print(f"Warning: no timing histogram for {label} ({ch}); skipped")
            continue
        win = f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) < {in_time_window(label)}"
        rdf = rdf.Define(f"{ch}_intime", win)
        r_in, r_out = rdf.Filter(f"{ch}_intime"), rdf.Filter(f"!{ch}_intime")
        e_rng, a_rng = _ranges(label)
        booked[label] = dict(
            ch=ch, mpv=mpv,
            n_in=r_in.Count(),
            n_sb=rdf.Filter(f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) >= {SIDEBAND_TS[0]} && "
                            f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) < {SIDEBAND_TS[1]}").Count(),
            e_sb=rdf.Filter(f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) >= {SIDEBAND_TS[0]} && "
                            f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) < {SIDEBAND_TS[1]}")
                   .Histo1D((f"h_e_sb_{ch}", "", *e_rng), f"{ch}_energy"),
            # exact means and spreads over the events: the histograms below have
            # finite ranges and their means would be biased by the overflow
            e_in_mean=r_in.Mean(f"{ch}_energy"), e_in_std=r_in.StdDev(f"{ch}_energy"),
            e_out_mean=r_out.Mean(f"{ch}_energy"), e_out_std=r_out.StdDev(f"{ch}_energy"),
            e_sb_mean=rdf.Filter(f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) >= {SIDEBAND_TS[0]} && "
                                 f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) < {SIDEBAND_TS[1]}").Mean(f"{ch}_energy"),
            e_sb_std=rdf.Filter(f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) >= {SIDEBAND_TS[0]} && "
                                f"abs({ch}_TS_cfd_mcp - {mpv:.3f}) < {SIDEBAND_TS[1]}").StdDev(f"{ch}_energy"),
            a_in_mean=r_in.Mean(f"{ch}_peak_value"), a_out_mean=r_out.Mean(f"{ch}_peak_value"),
            e_in=r_in.Histo1D((f"h_e_in_{ch}", "", *e_rng), f"{ch}_energy"),
            e_out=r_out.Histo1D((f"h_e_out_{ch}", "", *e_rng), f"{ch}_energy"),
            a_in=r_in.Histo1D((f"h_a_in_{ch}", "", *a_rng), f"{ch}_peak_value"),
            a_out=r_out.Histo1D((f"h_a_out_{ch}", "", *a_rng), f"{ch}_peak_value"),
            # fine-binned out-of-time copies for the single-p.e. calibration; the
            # scintillating fibres' own ranges are too coarse to resolve the peaks
            a_out_fine=r_out.Histo1D((f"h_a_out_fine_{ch}", "", 200, 0.0, 400.0), f"{ch}_peak_value"),
            e_out_fine=r_out.Histo1D((f"h_e_out_fine_{ch}", "", 360, -1000.0, 8000.0), f"{ch}_energy"),
            # mean baseline-subtracted waveform vs MCP-aligned TS
            p_in=r_in.Profile1D((f"p_in_{ch}", "", 1024, 0, 1024), get_arr_name(ch, use_mcp=True), f"{ch}_blsub"),
            p_out=r_out.Profile1D((f"p_out_{ch}", "", 1024, 0, 1024), get_arr_name(ch, use_mcp=True), f"{ch}_blsub"),
            p_all=rdf.Profile1D((f"p_all_{ch}", "", 1024, 0, 1024), get_arr_name(ch, use_mcp=True), f"{ch}_blsub"),
        )
    # event-by-event in-time flags (this triggers the single event loop)
    flag_cols = [f"{b['ch']}_intime" for b in booked.values()]
    flags = rdf.AsNumpy(columns=flag_cols)
    F = np.column_stack([flags[c].astype(bool) for c in flag_cols])   # events x fibres

    # ---- per-fibre plots + summary -------------------------------------------
    btype, benergy = getRunInfo(args.run)
    lumi = f"{btype}, {benergy} GeV"
    paths = get_run_paths(args.run)
    pm = PlotManager(rootdir=paths["root"], plotdir=paths["plots"], htmldir=paths["html"],
                     use_jsroot=args.jsroot, run_number=args.run)
    pm.set_output_dir(PLOT_SUBDIR)
    outdir = pm.get_output_dir()
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(paths["html"], "DRS"), exist_ok=True)
    fout = ROOT.TFile(os.path.join(paths["root"], ROOT_FILE), "RECREATE")
    overlay = {"prof": [], "prof_aligned": [], "amp": [], "energy": [], "amp_sci": [], "energy_sci": []}

    n_total = n_all.GetValue()
    summary = {}
    fibre_plot_names = []   # per-fibre plots, re-added after the summaries
    for label, b in booked.items():
        n_in = b["n_in"].GetValue()
        e_in, e_out = b["e_in"].GetValue(), b["e_out"].GetValue()
        a_in, a_out = b["a_in"].GetValue(), b["a_out"].GetValue()
        p_in, p_out, p_all = b["p_in"].GetValue(), b["p_out"].GetValue(), b["p_all"].GetValue()
        frac_in = n_in / n_total if n_total else 0.0
        # exact means over events (the histogram means are clipped by their ranges)
        n_out = n_total - n_in
        E_in, E_out = b["e_in_mean"].GetValue(), b["e_out_mean"].GetValue()
        sem_in = b["e_in_std"].GetValue() / n_in ** 0.5 if n_in > 1 else 0.0
        sem_out = b["e_out_std"].GetValue() / n_out ** 0.5 if n_out > 1 else 0.0
        A_in, A_out = b["a_in_mean"].GetValue(), b["a_out_mean"].GetValue()
        # sideband subtraction: randoms expected in the window, and their mean energy
        n_sb, e_sb = b["n_sb"].GetValue(), b["e_sb"].GetValue()
        E_sb = b["e_sb_mean"].GetValue()
        sem_sb = b["e_sb_std"].GetValue() / n_sb ** 0.5 if n_sb > 1 else 0.0
        sb_width = 2.0 * (SIDEBAND_TS[1] - SIDEBAND_TS[0])
        n_rand = n_sb * 2.0 * in_time_window(label) / sb_width
        n_hit = max(n_in - n_rand, 0.0)
        if n_hit > 0:
            e_hit = (n_in * E_in - n_rand * E_sb) / n_hit
            e_hit_err = ((n_in * sem_in) ** 2 + (n_rand * sem_sb) ** 2) ** 0.5 / n_hit
        else:
            e_hit, e_hit_err = 0.0, 0.0
        frac_hit = n_hit / n_total if n_total else 0.0
        frac_hit_err = ((n_in + n_rand * 2.0 * in_time_window(label) / sb_width) ** 0.5 / n_total) if n_total else 0.0
        ax = p_all.GetXaxis()
        b_lo, b_hi = ax.FindBin(PROFILE_TS_RANGE[0]), ax.FindBin(PROFILE_TS_RANGE[1])
        prof_peak = {k: max(h.GetBinContent(i) for i in range(b_lo, b_hi + 1))
                     for k, h in (("in", p_in), ("out", p_out), ("all", p_all))}
        b_peak_in = max(range(b_lo, b_hi + 1), key=p_in.GetBinContent)
        peak_ts_in = ax.GetBinCenter(b_peak_in)
        aligned = {k: _shift_to_peak(h, b_peak_in, f"h_prof_{k}_aligned_{b['ch']}")
                   for k, h in (("in", p_in), ("out", p_out), ("all", p_all))}
        # in-time minus out-of-time mean energy (older estimate, kept for reference)
        e_diff = E_in - E_out
        e_diff_err = (sem_in ** 2 + sem_out ** 2) ** 0.5
        is_scint = "scint" in label.lower()
        # single-p.e. calibration from the out-of-time (dark count) sample, for
        # every fibre; note the scintillating fibres saturate in 11-19 % of
        # their in-time events, so their p.e. counts are lower limits
        a_out_fine, e_out_fine = b["a_out_fine"].GetValue(), b["e_out_fine"].GetValue()
        gain_a = _peak_in(a_out_fine, *PE_AMP_SEARCH)
        a0 = _peak_in(a_out_fine, *NOISE_AMP_SEARCH)
        gain_e = _peak_in(e_out_fine, *PE_INT_SEARCH)
        # amplitude: subtract the 0 p.e. noise level; integral: baseline is
        # already subtracted, and in-time windows sit on a real pulse
        pe_amp = (A_in - a0) / gain_a
        pe_int = E_in / gain_e
        summary[label] = dict(
            channel=b["ch"], timing_mpv_ts=round(b["mpv"], 2),
            n_events=n_total, n_in_time=n_in, in_time_fraction=round(frac_in, 4),
            mean_energy_in=round(E_in, 1), mean_energy_out=round(E_out, 1),
            sem_energy_in=round(sem_in, 1), sem_energy_out=round(sem_out, 1),
            n_out_time=n_total - n_in, out_time_fraction=round(1.0 - frac_in, 4),
            in_time_fraction_err=round((frac_in * (1 - frac_in) / n_total) ** 0.5, 4) if n_total else 0.0,
            mean_amp_in=round(A_in, 2), mean_amp_out=round(A_out, 2),
            noise_amp_adc=a0, one_pe_amp_adc=gain_a, one_pe_int_adc=gain_e,
            profile_peak_in=round(prof_peak["in"], 2), profile_peak_out=round(prof_peak["out"], 2),
            profile_peak_all=round(prof_peak["all"], 2), profile_peak_ts_in=peak_ts_in,
            energy_diff=round(e_diff, 1), sem_energy_diff=round(e_diff_err, 1),
            energy_diff_pe=round(e_diff / gain_e, 2) if gain_e else None,
            pe_per_hit_amp=round(pe_amp, 2) if pe_amp is not None else None,
            pe_per_hit_int=round(pe_int, 2) if pe_int is not None else None,
            pe_per_beam_particle=round(pe_int * frac_in, 3) if pe_int is not None else None,
            # sideband-subtracted: real hits only
            n_sideband=n_sb, sideband_ts=list(SIDEBAND_TS), n_random_in_window=round(n_rand, 1),
            mean_energy_sideband=round(E_sb, 1),
            n_hit=round(n_hit, 1), hit_fraction=round(frac_hit, 4), hit_fraction_err=round(frac_hit_err, 4),
            mean_energy_hit=round(e_hit, 1), sem_energy_hit=round(e_hit_err, 1),
            pe_per_hit=round(e_hit / gain_e, 2) if gain_e else None,
            pe_per_hit_err=round(e_hit_err / gain_e, 3) if gain_e else None,
        )

        # keep everything in drs_intime.root and for the overlay pages
        fout.cd()
        for kind, h in (("amp_in", a_in), ("amp_out", a_out), ("energy_in", e_in), ("energy_out", e_out),
                        ("energy_sideband", e_sb),
                        ("amp_out_fine", a_out_fine), ("energy_out_fine", e_out_fine),
                        ("prof_in", p_in), ("prof_out", p_out), ("prof_all", p_all)):
            h.SetName(f"h_{kind}_{b['ch']}")
            h.Write()
        for h in aligned.values():
            h.Write()
        for sample, hp, ha, he in (("in", p_in, a_in, e_in), ("out", p_out, a_out, e_out)):
            overlay["prof"].append((hp, label, sample))
            overlay["amp_sci" if is_scint else "amp"].append((ha, label, sample))
            overlay["energy_sci" if is_scint else "energy"].append((he, label, sample))
        overlay["prof"].append((p_all, label, "all"))
        for k in ("in", "out", "all"):
            overlay["prof_aligned"].append((aligned[k], label, k))

        e_rng, a_rng = _ranges(label)
        legends = [f"out of time ({int(a_out.GetEntries())})",
                   f"in time |t-MPV| < {in_time_window(label):.0f} TS ({n_in})"]
        pm.plot_1d([a_out, a_in], f"intime_amp_{label}_{b['ch']}",
                   "Peak amplitude [ADC]", (a_rng[1], a_rng[2]), "Events", (0.5, None),
                   legends=legends, style=_STYLE, extra_text=label, lumitext=lumi)
        pm.plot_1d([e_out, e_in], f"intime_energy_{label}_{b['ch']}",
                   "CFD energy [ADC]", (e_rng[1], e_rng[2]), "Events", (0.5, None),
                   legends=legends, style=_STYLE, extra_text=label, lumitext=lumi)
        pm.plot_1d([p_out, p_in, p_all], f"intime_profile_{label}_{b['ch']}",
                   "MCP-aligned TS", PROFILE_TS_RANGE, "Mean ADC (baseline subtracted)", (None, None),
                   legends=legends + ["all events"], style=_STYLE_PROF, extra_text=label, lumitext=lumi)
        fibre_plot_names += [f"intime_amp_{label}_{b['ch']}", f"intime_energy_{label}_{b['ch']}",
                             f"intime_profile_{label}_{b['ch']}"]
    # page order: summary plots first (their own rows), then three plots per fibre
    pm.reset_plots()

    # ---- per-fibre summaries: in-time fraction and mean energy ---------------
    frac_vals = {l: (v["in_time_fraction"], v["in_time_fraction_err"]) for l, v in summary.items()}
    plot_fiber_summary(frac_vals, "In-time fraction", 0.0, 1.0, "intime_fraction_summary", pm, lumi)
    e_in_vals = {l: (v["mean_energy_in"], v["sem_energy_in"]) for l, v in summary.items()}
    e_out_vals = {l: (v["mean_energy_out"], v["sem_energy_out"]) for l, v in summary.items()}
    plot_fiber_summary(e_in_vals, "Mean CFD energy [ADC]",
                       100.0, 200000.0, "intime_energy_summary", pm, lumi,
                       secondary=e_out_vals, logy=True,
                       note="markers: in time,  x: out of time")
    e_hit_vals = {l: (v["mean_energy_hit"], v["sem_energy_hit"]) for l, v in summary.items()
                  if v["mean_energy_hit"] > 0}
    plot_fiber_summary(e_hit_vals, "Light per hit, randoms subtracted [ADC]",
                       10.0, 200000.0, "intime_energy_diff_summary", pm, lumi, logy=True)

    # ---- correlation between fibres ------------------------------------------
    labels = list(booked)
    corr_labels = [l for l in labels if l not in EXCLUDE_FROM_CORR]
    idx = [labels.index(l) for l in corr_labels]
    Fc = F[:, idx]
    p_in = Fc.mean(axis=0)

    # number of in-time fibres per event, observed vs. independent fibres
    n_fib = len(corr_labels)
    n_obs = Fc.sum(axis=1)
    h_n = ROOT.TH1F("h_nfib_obs", "", n_fib + 1, -0.5, n_fib + 0.5)
    for k in n_obs:
        h_n.Fill(int(k))
    h_ind = ROOT.TH1F("h_nfib_ind", "", n_fib + 1, -0.5, n_fib + 0.5)
    for k, p in enumerate(_poisson_binomial(p_in)):
        h_ind.SetBinContent(k + 1, p * len(n_obs))
    pm.plot_1d([h_ind, h_n], "intime_nfibres",
               f"Number of in-time fibres per event ({n_fib} fibres)", (-0.5, n_fib + 0.5),
               "Events", (0.5, None),
               legends=["independent fibres (expected)", "observed"],
               style=_STYLE, extra_text="all test fibres", lumitext=lumi, prepend=True)
    # (the fibre-count plot is now first; the matrix goes in front of it below)

    # P(j in time | i in time) and the phi coefficient against the reference fibre
    cond = [[(Fc[Fc[:, i], j].mean() if Fc[:, i].any() else 0.0) for j in range(n_fib)]
            for i in range(n_fib)]
    _draw_conditional_matrix(corr_labels, cond, outdir, "intime_conditional_matrix", lumi)
    pm.add_plot("intime_conditional_matrix", prepend=True)
    for name in ("intime_fraction_summary", "intime_energy_summary", "intime_energy_diff_summary"):
        pm.add_plot(name)
    pm.add_newline()   # per-fibre rows start fresh: amplitude | energy | profile
    for name in fibre_plot_names:
        pm.add_plot(name)

    corr_summary = {}
    if CORR_REFERENCE in corr_labels:
        r = corr_labels.index(CORR_REFERENCE)
        ref_in = Fc[:, r]
        for j, l in enumerate(corr_labels):
            x = Fc[:, j]
            phi = float(np.corrcoef(ref_in.astype(float), x.astype(float))[0, 1]) if x.any() and (~x).any() else 0.0
            corr_summary[l] = dict(
                p_in=round(float(p_in[j]), 4),
                p_in_given_ref_in=round(float(x[ref_in].mean()), 4),
                p_in_given_ref_out=round(float(x[~ref_in].mean()), 4),
                phi_with_ref=round(phi, 3))

    # ---- summary tables (JSON, markdown, and HTML intro) ----------------------
    with open(os.path.join(outdir, "intime_summary.json"), "w") as f:
        json.dump({"per_fibre": summary, "correlation_reference": CORR_REFERENCE,
                   "correlation": corr_summary, "n_in_time_fibres_observed": h_n and
                   [int(h_n.GetBinContent(k + 1)) for k in range(n_fib + 1)],
                   "n_in_time_fibres_independent": [round(float(v), 2) for v in _poisson_binomial(p_in) * len(n_obs)]},
                  f, indent=2)

    def cell(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.3f}" if abs(v) < 10 else f"{v:.0f}"
        return str(v)

    cols = [("Fibre", "label"), ("In-time frac.", "in_time_fraction"),
            ("<E> in", "mean_energy_in"), ("<E> out", "mean_energy_out"),
            ("<A> in", "mean_amp_in"), ("<A> out", "mean_amp_out"),
            ("1 p.e. amp", "one_pe_amp_adc"), ("1 p.e. int", "one_pe_int_adc"),
            ("p.e./hit (amp)", "pe_per_hit_amp"), ("p.e./hit (int)", "pe_per_hit_int"),
            ("p.e./beam e", "pe_per_beam_particle"),
            ("prof. peak in", "profile_peak_in"), ("prof. peak out", "profile_peak_out"),
            ("prof. peak all", "profile_peak_all")]
    rows = [[cell(l if k == "label" else summary[l][k]) for _, k in cols] for l in sorted(summary)]
    fcols = [("Fibre", "label"), ("Type", "type"), ("N in time", "n_in_time"), ("N out of time", "n_out_time"),
             ("In-time frac.", "in_time_fraction"), ("Out-of-time frac.", "out_time_fraction"),
             ("randoms in window", "n_random_in_window"), ("hit frac. (randoms subtracted)", "hit_fraction")]
    ecols = [("Fibre", "label"), ("Type", "type"), ("<E> in [ADC]", "mean_energy_in"), ("+- ", "sem_energy_in"),
             ("<E> out [ADC]", "mean_energy_out"), ("+-", "sem_energy_out"),
             ("<E> sideband [ADC]", "mean_energy_sideband")]
    dcols = [("Fibre", "label"), ("Type", "type"), ("light per hit [ADC]", "mean_energy_hit"), ("+-", "sem_energy_hit"),
             ("p.e. per hit", "pe_per_hit"), ("+- ", "pe_per_hit_err"),
             ("<E> in - <E> out [ADC]", "energy_diff"), ("in p.e.", "energy_diff_pe")]
    def row_of(cols_, l):
        out = []
        for _, k in cols_:
            if k == "label":
                out.append(l)
            elif k == "type":
                out.append(fiber_base_type(l))
            elif k.startswith("sem_") or k == "n_random_in_window":   # one decimal
                out.append(f"{summary[l][k]:.1f}")
            else:
                out.append(cell(summary[l][k]))
        return out
    frows = [row_of(fcols, l) for l in sorted(summary, key=lambda l: (fiber_base_type(l), l))]
    erows = [row_of(ecols, l) for l in sorted(summary, key=lambda l: (fiber_base_type(l), l))]
    drows = [row_of(dcols, l) for l in sorted(summary, key=lambda l: (fiber_base_type(l), l))]
    ccols = [("Fibre", "label"), ("P(in)", "p_in"), (f"P(in | {CORR_REFERENCE} in)", "p_in_given_ref_in"),
             (f"P(in | {CORR_REFERENCE} out)", "p_in_given_ref_out"), (f"phi with {CORR_REFERENCE}", "phi_with_ref")]
    crows = [[cell(l if k == "label" else corr_summary[l][k]) for _, k in ccols] for l in sorted(corr_summary)]

    def md_table(cols_, rows_):
        return "\n".join(["| " + " | ".join(h for h, _ in cols_) + " |", "|" + "---|" * len(cols_)]
                         + ["| " + " | ".join(r) + " |" for r in rows_])

    def html_table(cols_, rows_):
        esc = html_escape   # headers such as "<E> in" must not become tags
        th = "".join(f"<th style='padding:2px 8px;border-bottom:1px solid #999'>{esc(h)}</th>" for h, _ in cols_)
        trs = "".join("<tr>" + "".join(f"<td style='padding:2px 8px;text-align:right'>{esc(c)}</td>" for c in r) + "</tr>"
                      for r in rows_)
        return (f"<table style='border-collapse:collapse;font-size:13px;margin:10px 0'>"
                f"<thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table>")

    nfib_md = md_table([("N in-time fibres", ""), ("observed", ""), ("independent", "")],
                       [[str(k), str(int(h_n.GetBinContent(k + 1))), f"{h_ind.GetBinContent(k + 1):.1f}"]
                        for k in range(n_fib + 1)])
    md = (f"# In-time vs out-of-time summary, run {args.run} ({lumi})\n\n"
          f"{n_total} events{' after the mcp_clean selection' if args.mcp_clean else ''}; "
          f"in time = |TS_cfd_mcp - MPV| < {IN_TIME_WINDOW_TS:.0f} TS ({IN_TIME_WINDOW_SCINT_TS:.0f} TS for the "
          f"scintillating fibres). ADC units throughout.\n\n"
          f"## In-time vs out-of-time fraction by fibre\n\n{md_table(fcols, frows)}\n\n"
          f"## Mean energy by fibre\n\n{md_table(ecols, erows)}\n\n"
          f"## Light per hit, randoms subtracted, by fibre\n\n{md_table(dcols, drows)}\n\n"
          f"## Per fibre, full\n\n{md_table(cols, rows)}\n\n"
          f"## Correlation with {CORR_REFERENCE}\n\n{md_table(ccols, crows)}\n\n"
          f"## Number of in-time fibres per event ({n_fib} fibres)\n\n{nfib_md}\n")
    with open(os.path.join(outdir, "intime_summary.md"), "w") as f:
        f.write(md)
    print(md)

    intro = (f"<p>Run {args.run}, {lumi}, {n_total} events"
             f"{' after the mcp_clean selection' if args.mcp_clean else ''}. "
             f"In time = |TS_cfd_mcp - MPV| &lt; {IN_TIME_WINDOW_TS:.0f} TS, or {IN_TIME_WINDOW_SCINT_TS:.0f} TS for "
             f"the scintillating fibres whose slow pulses have a wider timing distribution (red); out of time = the rest (blue). "
             f"Single-p.e. calibration from the out-of-time sample: 1 p.e. amp = peak of the amplitude "
             f"distribution in [{PE_AMP_SEARCH[0]:.0f}, {PE_AMP_SEARCH[1]:.0f}] ADC, 1 p.e. int = peak of the "
             f"energy distribution in [{PE_INT_SEARCH[0]:.0f}, {PE_INT_SEARCH[1]:.0f}] ADC. "
             f"p.e./hit (amp) = (&lt;A&gt; in - noise peak) / 1 p.e. amp; p.e./hit (int) = &lt;E&gt; in / 1 p.e. int; "
             f"p.e./beam e = p.e./hit (int) x in-time fraction. Sideband subtraction: the random coincidences "
             f"(dark counts, noise) inside the in-time window are estimated from |TS_cfd_mcp - MPV| in "
             f"[{SIDEBAND_TS[0]:.0f}, {SIDEBAND_TS[1]:.0f}] TS on both sides; light per hit = "
             f"(N in &lt;E&gt; in - N random &lt;E&gt; sideband) / (N in - N random), hit fraction = "
             f"(N in - N random) / N, p.e. per hit = light per hit / 1 p.e. int. "
             f"prof. peak = maximum of the mean baseline-subtracted waveform vs MCP-aligned TS in "
             f"[{PROFILE_TS_RANGE[0]:.0f}, {PROFILE_TS_RANGE[1]:.0f}] for the in-time, out-of-time and all-event samples.</p>"
             + "<p><b>In-time vs out-of-time fraction by fibre</b></p>" + html_table(fcols, frows)
             + "<p><b>Mean energy by fibre</b> (standard error of the mean)</p>" + html_table(ecols, erows)
             + "<p><b>Light per hit, randoms subtracted, by fibre</b> (the light left by the beam particle in a "
               "real hit; the last two columns are the older in-minus-out estimate, which over-subtracts the "
               "pedestal for the high-yield fibres; the scintillating fibres saturate in part of their hits, so "
               "their p.e. values are lower limits)</p>" + html_table(dcols, drows)
             + "<p><b>Per fibre, full</b></p>" + html_table(cols, rows)
             + f"<p>Correlation between fibres: {', '.join(sorted(EXCLUDE_FROM_CORR))} are left out. "
             f"phi is the correlation coefficient of the two in-time flags.</p>"
             + html_table(ccols, crows))
    html = pm.generate_html(HTML_PAGE, plots_per_row=3,
                            title=f"Test fibres: in-time vs out-of-time, run {args.run}",
                            intro_text=intro)
    print(f"HTML: {html}")

    # interactive overlays: every fibre's in-time / out-of-time (/ all) curve on one
    # set of axes, one colour per fibre and one line style per sample; tick them
    # on and off per curve or per sample; peak-normalisation is on the page
    overlay_pages = [
        ("prof", "DRS/DRS_InTime_Prof_overlay.html", "MCP-aligned TS", "Mean ADC (baseline subtracted)",
         "mean waveform", PROFILE_TS_RANGE),
        ("prof_aligned", "DRS/DRS_InTime_Prof_aligned_overlay.html", "TS relative to the in-time peak",
         "Mean ADC (baseline subtracted)", "mean waveform, peaks aligned at 0", (-40.0, 120.0)),
        ("amp", "DRS/DRS_InTime_Amp_overlay.html", "Peak amplitude [ADC]", "Events",
         "pulse amplitude, Cherenkov fibres", None),
        ("amp_sci", "DRS/DRS_InTime_Amp_Sci_overlay.html", "Peak amplitude [ADC]", "Events",
         "pulse amplitude, scintillating fibres", None),
        ("energy", "DRS/DRS_InTime_Energy_overlay.html", "CFD energy [ADC]", "Events",
         "CFD energy, Cherenkov fibres", None),
        ("energy_sci", "DRS/DRS_InTime_Energy_Sci_overlay.html", "CFD energy [ADC]", "Events",
         "CFD energy, scintillating fibres", None),
    ]
    for key, page, xl, yl, what, xr in overlay_pages:
        if not overlay[key]:
            continue
        out = build_channel_overlay(
            overlay[key], output_html=os.path.join(paths["html"], page),
            xlabel=xl, ylabel=yl,
            title=f"Test fibres: {what}, in/out of time — Run {args.run}",
            intro_text=(f"Tick the curves to superimpose, or click the legend. [in] = "
                        f"|TS_cfd_mcp - MPV| < {IN_TIME_WINDOW_TS:.0f} TS ({IN_TIME_WINDOW_SCINT_TS:.0f} TS for "
                        f"the scintillating fibres), [out] = the rest"
                        f"{', [all] = every event' if key == 'prof' else ''}. "
                        f"Histograms are in {ROOT_FILE}."),
            filename=f"DRS_InTime_{key}_Run{args.run}", xrange=xr,
            variants=["in", "out", "all"] if key.startswith("prof") else ["in", "out"])
        if out:
            print(f"Overlay: {out}")
    fout.Close()
    PlotManager.print_html_summary()

    return html
