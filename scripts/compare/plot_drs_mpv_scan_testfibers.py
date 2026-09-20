"""
DRS MPV / Energy position scan.

Combines two complementary position scans:
  - Time scan  : MPV of TS_cfd_mcp histogram vs table position X per channel.
  - Energy scan: DRS_SumMean (mean energy) vs table position X per channel.

Usage
-----
  python scripts/compare/plot_drs_mpv_scan.py [--mode {time,energy,both}] [--jsroot]
"""

import argparse
import math
import os
import json
import ROOT

from channels.channel_map import build_drs_boards
from core.plot_manager import PlotManager
from plotting.my_function import DrawHistos
from plotting.calox_plot_helper import create_pave_text, fiber_base_type, plot_fiber_summary
from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var, get_hist_mpv
auto_timer("Total Execution Time")

# Runs with no usable beam signal in the fibres are left out:
#   2028, 2033  pi+ runs (beam file 94) mixed into the positron scan
#   2040, 2042  beam mostly missing the table (beam profile issue); 2042 cut short
SKIP_RUNS = {2028, 2033, 2040, 2041, 2042}
SCAN_RUNS = [r for r in range(1998, 2043) if r not in SKIP_RUNS]
REFERENCE_RUN = 1998

# y ranges of the per-fiber summary plots
SPEED_SUMMARY_RANGE     = (15.0, 22.0)     # cm/ns
TIME_P0_SUMMARY_RANGE   = (484.0, 495.0)   # TS
ENERGY_P0_SUMMARY_RANGE = (500.0, 1700.0)  # ADC

# An MPV closer than this to the edge of the fine-binned histogram means the
# peak was clipped by the histogram window (get_hist_mpv searches +-6 TS).
MPV_EDGE_TS = 6.0

PLOTDIR = "results/plots/PositionScan"
HTMLDIR = "results/html/PositionScan"

# beam type spellings in data/Runlist.json -> ROOT latex
_BTYPES = {"pion": "#pi^{+}", "pions": "#pi^{+}", "pi+": "#pi^{+}",
           "positron": "e^{+}", "positrons": "e^{+}", "e+": "e^{+}",
           "mu+": "#mu^{+}"}
# marker style per beam energy [GeV]: filled circle, open square, open triangle
ENERGY_MARKERS = {40: 20, 60: 25, 80: 26}
LEGEND_POS = [0.56, 0.79, 0.88, 0.89]
X_BEAM = {}   # table position [mm] -> (beam type, beam energy [GeV]); filled in main()


def _graphs_by_beam(pts, name, with_errors=True):
    """Split scan points (x, y[, err]) by beam energy into one graph each, so
    the runs at different energies get their own marker. Returns
    (graphs, legend labels, marker styles), lowest energy first."""
    groups = {}
    for pt in pts:
        groups.setdefault(X_BEAM.get(pt[0], ("", 0)), []).append(pt)
    graphs, labels, markers = [], [], []
    for (btype, energy), gpts in sorted(groups.items(), key=lambda kv: kv[0][1]):
        gr = ROOT.TGraphErrors(len(gpts)) if with_errors else ROOT.TGraph(len(gpts))
        gr.SetName(f"{name}_{energy}GeV")
        gr.SetTitle("")
        for i, pt in enumerate(sorted(gpts)):
            gr.SetPoint(i, pt[0], pt[1])
            if with_errors:
                gr.SetPointError(i, 0.0, pt[2])
        gr.SetLineWidth(2)
        graphs.append(gr)
        labels.append(f"{_BTYPES.get(btype.lower(), btype)} {energy} GeV" if energy else "beam unknown")
        markers.append(ENERGY_MARKERS.get(energy, 21))
    return graphs, labels, markers



# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _build_channel_order(drsboards):
    """Return list of (ch, board_no, group_no, chan_no, var, tower_x, tower_y, is6mm)."""
    order = []
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            order.append((
                chan.get_channel_name(blsub=False),
                board.board_no, chan.group_no, chan.channel_no,
                get_channel_var(chan),
                chan.i_tower_x, chan.i_tower_y,
                chan.is6mm,
            ))
    return order


def _make_pm(outdir, use_jsroot=False):
    pm = PlotManager(
        rootdir=f"results/root/Run{REFERENCE_RUN}",
        plotdir=PLOTDIR,
        htmldir=HTMLDIR,
        use_jsroot=use_jsroot,
    )
    pm.set_output_dir(outdir)
    return pm


# ---------------------------------------------------------------------------
# Time scan — data collection
# ---------------------------------------------------------------------------

def _collect_mpv_data(drsboards, run_x):
    """Read MPV of the TS_cfd_mcp histogram for each run/channel.

    Prefers the fine-binned (0.1 TS) histogram, whose window (see
    configs/plot_config.py) also keeps out the pile-up spike near 1000 TS that
    can fool the MPV finder on the coarse 0-1024 TS histogram. If a run was
    produced with an older, narrower window the peak may be clipped: that shows
    up as an MPV within one search window of the histogram edge, in which case
    fall back to the coarse histogram and say so."""
    mpv_data = {}
    for run, x_val in sorted(run_x.items(), key=lambda kv: kv[1]):
        root_path = f"results/root/Run{run}/drs_stats.root"
        if not os.path.exists(root_path):
            print(f"Warning: {root_path} not found, skipping run {run}")
            continue
        infile = ROOT.TFile(root_path, "READ")
        for _, board in drsboards.items():
            for chan in board:
                if chan.is_reference:
                    continue
                ch = chan.get_channel_name(blsub=False)
                hist_fb = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
                mpv = None
                if hist_fb and hist_fb.Integral() >= 5:
                    mpv, mpv_err = get_hist_mpv(hist_fb)
                    xax = hist_fb.GetXaxis()
                    if not (xax.GetXmin() + MPV_EDGE_TS <= mpv <= xax.GetXmax() - MPV_EDGE_TS):
                        print(f"Warning: run {run} {ch}: fine-binned MPV {mpv:.1f} TS "
                              f"sits at the histogram edge [{xax.GetXmin():.0f}, "
                              f"{xax.GetXmax():.0f}]; using the coarse histogram")
                        mpv = None
                if mpv is None:
                    hist = infile.Get(f"hist_{ch}_TS_cfd_mcp")
                    if not hist or hist.GetEntries() < 5:
                        continue
                    mpv, mpv_err = get_hist_mpv(hist)
                    if hist_fb and not (hist_fb.GetXaxis().GetXmin() <= mpv
                                        <= hist_fb.GetXaxis().GetXmax()):
                        # The coarse finder landed outside the fine window too
                        # (typically on the pile-up spike near 1000 TS): no
                        # usable peak in this run, so leave the point out.
                        print(f"Warning: run {run} {ch}: no usable timing peak "
                              f"(coarse MPV {mpv:.1f} TS); point dropped")
                        continue
                mpv_data.setdefault(ch, []).append((x_val, mpv, mpv_err))
        infile.Close()
    return mpv_data


# ---------------------------------------------------------------------------
# Time scan — fitting
# ---------------------------------------------------------------------------

def _fit_channel(pts):
    """Weighted least-squares independent linear fit. Returns (p0, p1, e0, e1)."""
    w    = [1.0 / ye**2 for _, _, ye in pts]
    sw   = sum(w)
    swx  = sum(wi * x     for wi, (x, _, _) in zip(w, pts))
    swy  = sum(wi * y     for wi, (_, y, _) in zip(w, pts))
    swxx = sum(wi * x**2  for wi, (x, _, _) in zip(w, pts))
    swxy = sum(wi * x * y for wi, (x, y, _) in zip(w, pts))
    denom = sw * swxx - swx ** 2
    p1 = (sw * swxy - swx * swy) / denom
    p0 = (swy - p1 * swx) / sw
    e1 = (sw   / denom) ** 0.5
    e0 = (swxx / denom) ** 0.5
    return p0, p1, e0, e1


# ---------------------------------------------------------------------------
# Time scan — per-channel plots
# ---------------------------------------------------------------------------

def _plot_time_channels(channel_order, mpv_data, pm, lumi, xmin, xmax,
                        ch_label_map=None):
    """Draw one MPV-vs-X TGraph per channel. Returns (p0_map, p1_map, e0_map, e1_map, plots_by_type)."""
    plots_by_type = {}
    p0_map, p1_map, e0_map, e1_map = {}, {}, {}, {}

    for ch, b, g, c, var, tx, ty, _is6mm in channel_order:
        pts = sorted(mpv_data.get(ch, []))
        if not pts:
            continue

        graphs, labels, markers = _graphs_by_beam(pts, f"gr_mpv_{ch}")
        gr = graphs[0]

        p0, p1, e0, e1 = _fit_channel(pts)

        fit = ROOT.TF1(f"fit_{ch}", "pol1", xmin, xmax)
        fit.SetParameters(p0, p1)
        fit.SetParError(0, e0)
        fit.SetParError(1, e1)
        fit.SetLineColor(ROOT.kRed + 1)
        fit.SetLineWidth(2)
        fit.SetLineStyle(2)
        gr.GetListOfFunctions().Add(fit)

        ys = [y for _, y, _ in pts]
        ymid  = (max(ys) + min(ys)) / 2
        yhalf = max(30.0, (max(ys) - min(ys)) / 2 * 1.5 + 10)

        speed    = 0.5 / abs(p1) if p1 != 0 else float("nan")
        e_speed  = speed * e1 / abs(p1) if p1 != 0 else float("nan")

        pave_bgc = create_pave_text(0.14, 0.82, 0.52, 0.89)
        pave_bgc.AddText(f"B: {b}, G: {g}, C: {c}")

        pave_fit = create_pave_text(0.22, 0.24, 0.42, 0.38)
        pave_fit.AddText(f"p0: {p0:.1f} #pm {e0:.1f} TS")
        pave_fit.AddText(f"p1: {p1:.3f} #pm {e1:.3f} TS/mm")
        pave_fit.AddText(f"v = 1/|p_{1}|: {speed:.2f} #pm {e_speed:.2f} cm/ns")

        label     = (ch_label_map or {}).get(ch, "")
        plot_name = (f"{label}_drs_cfd_mpv_scan_{ch}" if label
                     else f"drs_cfd_mpv_scan_{ch}_{var}")
        DrawHistos(
            graphs, labels,
            xmin, xmax, "Table position X [mm]",
            ymid - yhalf, ymid + yhalf, "MPV [TS]",
            plot_name,
            outdir=pm.get_output_dir(),
            drawoptions="EPL",
            markerstyles=markers,
            mycolors=[ROOT.kBlue + 1] * len(graphs),
            legendPos=LEGEND_POS,
            dology=False,
            usePDF=False,
            lumitext=lumi,
            extra_text=label or var,
            extraToDraw=[pave_bgc, pave_fit],
        )
        plots_by_type.setdefault(var, []).append(plot_name)
        p0_map[ch] = p0
        p1_map[ch] = p1
        e0_map[ch] = e0
        e1_map[ch] = e1

    return p0_map, p1_map, e0_map, e1_map, plots_by_type


# ---------------------------------------------------------------------------
# Energy scan — data collection and plotting
# ---------------------------------------------------------------------------

def _collect_summean(run_x):
    """Read DRS_SumMean (mean energy, vals[0]) per channel from each run's JSON."""
    data = {}
    for run, x_val in sorted(run_x.items(), key=lambda kv: kv[1]):
        path = f"results/root/Run{run}/drs_energy_stats.json"
        if not os.path.exists(path):
            print(f"Warning: {path} not found, skipping run {run}")
            continue
        with open(path) as f:
            stats = json.load(f)
        for ch, vals in stats.items():
            data.setdefault(ch, []).append((x_val, vals[0]))
    return data


def _fit_constant(pts, xmin, xmax):
    """Hand-written constant (pol0) fit to (x, y) points with x in [xmin, xmax].

    Points carry no errors, so p0 is the plain mean and e0 the standard error
    of the mean. Falls back to all points if none lie inside the window.
    Returns (p0, e0, n_used)."""
    ys = [y for x, y in pts if xmin <= x <= xmax] or [y for _, y in pts]
    n  = len(ys)
    p0 = sum(ys) / n
    e0 = (sum((y - p0) ** 2 for y in ys) / (n - 1)) ** 0.5 / n ** 0.5 if n > 1 else 0.0
    return p0, e0, n


def _plot_energy_channels(channel_order, data, pm, lumi, xmin, xmax, ch_label_map=None):
    """Draw one SumMean-vs-position TGraph per channel with constant fit.
    Returns (p0_map, e0_map, plots_by_type)."""
    plots_by_type = {}
    p0_map, e0_map = {}, {}

    FIT_XMIN, FIT_XMAX = -220.0, 50.0

    for ch, b, g, c, var, tx, ty, _is6mm in channel_order:
        pts = sorted(data.get(ch, []))
        if not pts:
            continue

        graphs, labels, markers = _graphs_by_beam(pts, f"gr_summean_{ch}", with_errors=False)
        gr = graphs[0]

        # Constant (pol0) fit restricted to [-220, 50] mm, done by hand so it
        # does not depend on a ROOT minimizer being available.
        p0, e0, _n = _fit_constant(pts, FIT_XMIN, FIT_XMAX)

        fit = ROOT.TF1(f"fit_sum_{ch}", "pol0", FIT_XMIN, FIT_XMAX)
        fit.SetParameter(0, p0)
        fit.SetParError(0, e0)
        fit.SetLineColor(ROOT.kRed + 1)
        fit.SetLineWidth(2)
        fit.SetLineStyle(2)
        gr.GetListOfFunctions().Add(fit)
        p0_map[ch] = p0
        e0_map[ch] = e0

        ys = [m for _, m in pts]
        ymid = (max(ys) + min(ys)) / 2
        yhalf = max(50.0, (max(ys) - min(ys)) / 2 * 1.4 + 20)

        pave = create_pave_text(0.14, 0.82, 0.52, 0.89)
        pave.AddText(f"B: {b}, G: {g}, C: {c}")

        pave_fit = create_pave_text(0.55, 0.30, 0.88, 0.36)
        pave_fit.AddText(f"p0: {p0:.1f} #pm {e0:.1f} ADC")

        label     = (ch_label_map or {}).get(ch, "")
        plot_name = (f"{label}_drs_summean_scan_{ch}" if label
                     else f"drs_summean_scan_{ch}_{var}")
        DrawHistos(
            graphs, labels,
            xmin, xmax, "Table position X [mm]",
            ymid - yhalf, ymid + yhalf, "Mean DRS Sum [ADC]",
            plot_name,
            outdir=pm.get_output_dir(),
            drawoptions="EPL",
            markerstyles=markers,
            mycolors=[ROOT.kBlue + 1] * len(graphs),
            legendPos=LEGEND_POS,
            dology=False,
            usePDF=False,
            lumitext=lumi,
            extra_text=label or var,
            extraToDraw=[pave, pave_fit],
        )
        plots_by_type.setdefault(var, []).append(plot_name)
    return p0_map, e0_map, plots_by_type


# ---------------------------------------------------------------------------
# Speed summary plot
# ---------------------------------------------------------------------------

def _plot_speed_summary(ch_label_map, p1_map, e1_map, pm, lumi):
    """Speed v = 0.5/|p1| [cm/ns] per fiber, from the time-scan slope."""
    speed_data = {}   # label -> (speed, e_speed)
    for ch, label in ch_label_map.items():
        p1 = p1_map.get(ch)
        e1 = e1_map.get(ch)
        if not p1:
            continue
        s  = 0.5 / abs(p1)
        speed_data[label] = (s, s * abs(e1) / abs(p1))
    return plot_fiber_summary(speed_data, "Speed  v = 1/|p_{1}| [cm/ns]",
                               SPEED_SUMMARY_RANGE[0], SPEED_SUMMARY_RANGE[1],
                               "drs_cfd_speed_summary", pm, lumi)


def _plot_p0_summary_time(ch_label_map, p0_map, e0_map, pm, lumi):
    """Time-scan intercept p0 [TS] per fiber (MPV extrapolated to X = 0)."""
    p0_data = {label: (p0_map[ch], e0_map.get(ch, 0.0))
               for ch, label in ch_label_map.items() if p0_map.get(ch) is not None}
    return plot_fiber_summary(p0_data, "MPV at X = 0  p_{0} [TS]",
                               TIME_P0_SUMMARY_RANGE[0], TIME_P0_SUMMARY_RANGE[1],
                               "drs_cfd_p0_summary", pm, lumi)


# ---------------------------------------------------------------------------
# Energy p0 summary plot
# ---------------------------------------------------------------------------

def _plot_p0_summary_energy(ch_label_map, p0_map, e0_map, pm, lumi):
    """Mean DRS Sum p0 [ADC] per fiber, from the energy-scan constant fit."""
    p0_data = {label: (p0_map[ch], e0_map.get(ch) or 0.0)
               for ch, label in ch_label_map.items() if p0_map.get(ch) is not None}
    return plot_fiber_summary(p0_data, "Mean DRS Sum  p_{0} [ADC]",
                               ENERGY_P0_SUMMARY_RANGE[0], ENERGY_P0_SUMMARY_RANGE[1],
                               "drs_summean_p0_summary", pm, lumi)


# ---------------------------------------------------------------------------
# Channel fit-parameter table
# ---------------------------------------------------------------------------

def _print_channel_table(ch_label_map, p0_map, p1_map, e0_map, e1_map, outdir):
    """Print and save a table; return an HTML <table> string for embedding in the page."""
    _TS_NS = 0.2   # 1 TS = 200 ps = 0.2 ns
    rows = []
    for ch, p0 in sorted(p0_map.items(), key=lambda kv: ch_label_map.get(kv[0], kv[0])):
        label   = ch_label_map.get(ch, "")
        compact = ch.replace("DRS_", "")
        p1      = p1_map.get(ch, float("nan"))
        e0      = e0_map.get(ch, float("nan"))
        e1      = e1_map.get(ch, float("nan"))
        speed   = 0.5 / abs(p1)       if p1 != 0 else float("nan")  # cm/ns
        e_speed = speed * e1 / abs(p1) if p1 != 0 else float("nan")
        rows.append((label, compact, p0, e0, p1, e1, speed, e_speed))

    # --- print ---
    w_label   = max(len(r[0]) for r in rows) if rows else 10
    w_channel = max(len(r[1]) for r in rows) if rows else 30
    header = (f"{'Fiber type':<{w_label}}  {'Channel':<{w_channel}}  "
              f"{'p0 [TS]':>18}  {'p1 [TS/mm]':>20}  {'v [cm/ns]':>16}")
    sep = "-" * len(header)
    print(f"\n{sep}\n{header}\n{sep}")
    for label, ch_short, p0, e0, p1, e1, speed, e_speed in rows:
        print(f"{label:<{w_label}}  {ch_short:<{w_channel}}  "
              f"{p0:>7.2f} ± {e0:<7.2f}  "
              f"{p1:>8.4f} ± {e1:<8.4f}  "
              f"{speed:>6.2f} ± {e_speed:<6.2f}")
    print(sep)

    # --- save JSON alongside the DRS_CFD_MPV_Scan_channels.html page ---
    out = {
        label: {
            "channel": f"DRS_{ch_short}",
            "p0": round(p0, 4), "e0": round(e0, 4),
            "p1": round(p1, 6), "e1": round(e1, 6),
            "speed_cmns": round(speed, 4), "e_speed_cmns": round(e_speed, 4),
        }
        for label, ch_short, p0, e0, p1, e1, speed, e_speed in rows
    }
    out_path = os.path.join(outdir, "drs_cfd_mpv_scan_channels.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved to {out_path}\n")

    # --- build HTML table for embedding as intro_text ---
    cell  = 'style="padding:4px 10px;border:1px solid #ccc"'
    rcell = 'style="padding:4px 10px;border:1px solid #ccc;text-align:right"'
    th = "".join(f'<th {cell}>{h}</th>'
                 for h in ["Fiber type", "Channel",
                            "p0 [TS]", "±e0", "p1 [TS/mm]", "±e1",
                            "v [cm/ns]", "±ev"])
    trs = "".join(
        f'<tr><td {cell}>{lbl}</td><td {cell}>{ch}</td>'
        f'<td {rcell}>{p0:.2f}</td><td {rcell}>{e0:.2f}</td>'
        f'<td {rcell}>{p1:.4f}</td><td {rcell}>{e1:.4f}</td>'
        f'<td {rcell}>{speed:.2f}</td><td {rcell}>{e_speed:.2f}</td></tr>'
        for lbl, ch, p0, e0, p1, e1, speed, e_speed in rows
    )
    return (
        f'<table style="border-collapse:collapse;font-size:13px;margin:10px 0">'
        f'<thead><tr style="background:#f0f0f0">{th}</tr></thead>'
        f'<tbody>{trs}</tbody></table>'
    )


def _print_energy_table(ch_label_map, p0_map, e0_map, outdir):
    """Print and save energy p0 table; return HTML <table> string for embedding."""
    rows = []
    for ch, p0 in sorted(p0_map.items(), key=lambda kv: ch_label_map.get(kv[0], kv[0])):
        label   = ch_label_map.get(ch, "")
        compact = ch.replace("DRS_", "")
        e0      = e0_map.get(ch, float("nan"))
        rows.append((label, compact, p0, e0))

    w_label   = max(len(r[0]) for r in rows) if rows else 10
    w_channel = max(len(r[1]) for r in rows) if rows else 30
    header = (f"{'Fiber type':<{w_label}}  {'Channel':<{w_channel}}  "
              f"{'p0 [ADC]':>20}")
    sep = "-" * len(header)
    print(f"\n{sep}\n{header}\n{sep}")
    for label, ch_short, p0, e0 in rows:
        print(f"{label:<{w_label}}  {ch_short:<{w_channel}}  "
              f"{p0:>9.1f} ± {e0:<9.1f}")
    print(sep)

    out = {
        label: {
            "channel": f"DRS_{ch_short}",
            "p0_adc": round(p0, 2),
            "e0_adc": round(e0, 2),
        }
        for label, ch_short, p0, e0 in rows
    }
    out_path = os.path.join(outdir, "drs_summean_scan_channels.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved to {out_path}\n")

    cell  = 'style="padding:4px 10px;border:1px solid #ccc"'
    rcell = 'style="padding:4px 10px;border:1px solid #ccc;text-align:right"'
    th = "".join(f'<th {cell}>{h}</th>'
                 for h in ["Fiber type", "Channel", "p0 [ADC]", "±e0"])
    trs = "".join(
        f'<tr><td {cell}>{lbl}</td><td {cell}>{ch}</td>'
        f'<td {rcell}>{p0:.1f}</td><td {rcell}>{e0:.1f}</td></tr>'
        for lbl, ch, p0, e0 in rows
    )
    return (
        f'<table style="border-collapse:collapse;font-size:13px;margin:10px 0">'
        f'<thead><tr style="background:#f0f0f0">{th}</tr></thead>'
        f'<tbody>{trs}</tbody></table>'
    )


# ---------------------------------------------------------------------------
# Top-level scan runners
# ---------------------------------------------------------------------------

def run_time_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                  use_jsroot=False, single_page=False, ch_label_map=None):
    mpv_data = _collect_mpv_data(drsboards, run_x)

    pm = _make_pm("DRS_CFD_MPV_Scan", use_jsroot)
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    p0_map, p1_map, e0_map, e1_map, plots_by_type = _plot_time_channels(
        channel_order, mpv_data, pm, lumi, xmin, xmax,
        ch_label_map=ch_label_map)

    # Build and save the channel table before generating the HTML page so it can
    # be embedded as intro_text inside DRS_CFD_MPV_Scan_channels.html.
    table_html = ""
    if ch_label_map and p0_map:
        table_html = _print_channel_table(ch_label_map, p0_map, p1_map, e0_map, e1_map,
                                          outdir=pm.get_output_dir())

    if single_page:
        all_plots = [name for names in plots_by_type.values() for name in names]

        # Speed and p0 summaries — drawn first so they appear at the top of the page
        summary_plots = []
        if ch_label_map and p1_map:
            summary_plots.append(_plot_speed_summary(ch_label_map, p1_map, e1_map, pm, lumi))
            summary_plots.append(_plot_p0_summary_time(ch_label_map, p0_map, e0_map, pm, lumi))
        summary_plots = [name for name in summary_plots if name]

        pm.reset_plots()
        for name in summary_plots:
            pm.add_plot(name)
        if summary_plots:
            pm.add_newline()
        for name in all_plots:
            pm.add_plot(name)
        pm.generate_html(
            "DRS_CFD_MPV_Scan_channels.html",
            plots_per_row=4,
            title="DRS CFD MPV Scan — selected channels",
            intro_text=table_html,
        )
    else:
        for type_name, plot_names in plots_by_type.items():
            pm.reset_plots()
            for name in plot_names:
                pm.add_plot(name)
            pm.generate_html(
                f"DRS_CFD_MPV_Scan_{type_name}.html",
                plots_per_row=4,
                title=f"DRS CFD MPV Scan — {type_name}",
            )

    return p0_map, p1_map, e0_map, e1_map


def run_energy_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                    use_jsroot=False, single_page=False, ch_label_map=None):
    data = _collect_summean(run_x)
    if not any(data.values()):
        print("No drs_energy_stats.json data found for any scan run.")
        return

    pm = _make_pm("DRS_SumMean_Scan", use_jsroot)
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    p0_map, e0_map, plots_by_type = _plot_energy_channels(
        channel_order, data, pm, lumi, xmin, xmax, ch_label_map=ch_label_map)

    table_html = ""
    if ch_label_map and p0_map:
        table_html = _print_energy_table(ch_label_map, p0_map, e0_map,
                                         outdir=pm.get_output_dir())

    if single_page:
        all_plots = [name for names in plots_by_type.values() for name in names]

        p0_plot = None
        if ch_label_map and p0_map:
            p0_plot = _plot_p0_summary_energy(ch_label_map, p0_map, e0_map, pm, lumi)

        pm.reset_plots()
        if p0_plot:
            pm.add_plot(p0_plot)
            pm.add_newline()
        for name in all_plots:
            pm.add_plot(name)
        pm.generate_html(
            "DRS_SumMean_Scan_channels.html",
            plots_per_row=4,
            title="DRS SumMean Scan — selected channels",
            intro_text=table_html,
        )
    else:
        for type_name, plot_names in plots_by_type.items():
            pm.reset_plots()
            for name in plot_names:
                pm.add_plot(name)
            pm.generate_html(
                f"DRS_SumMean_Scan_{type_name}.html",
                plots_per_row=4,
                title=f"DRS SumMean Scan — {type_name}",
            )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Photoelectron scan — from the per-run in-time analysis (analysis/intime_analysis.py)
# ---------------------------------------------------------------------------

PE_FIT_XMIN, PE_FIT_XMAX = -220.0, 50.0    # plateau window for the constant fits [mm]

# Exponential attenuation fits, light = A_E * exp(x / lambda): one normalisation
# A_E per beam energy (the 60 GeV runs give more light) and one common lambda,
# plus a separate fit per energy. lambda > 0 means more light towards +X, which
# the timing scan identifies as the readout end. Runs whose in-time fraction is
# below EDGE_FRACTION of the fibre's plateau value are edge runs where the beam
# only partly overlaps the fibre; they are drawn but left out of the fits.
EXP_FIT_KEYS = ("mean_energy_hit", "pe_per_hit")
EDGE_FRACTION = 0.6
MIN_FIT_SPAN_MM = 100.0   # a group needs at least this much travel for its own lambda
# The runs at different beam energies give different amounts of light, so each
# beam energy gets its own normalisation; the slope within each energy group
# constrains the common lambda, and each group is also fitted on its own.


def _fit_exponential(groups):
    """groups: {energy: [(x, y, ey)]}. Weighted least squares of ln y = a_g + b x
    with a common slope b; errors scaled by sqrt(chi2/ndf) when the scatter
    exceeds the point errors. Returns (b, eb, {g: a_g}, {g: (b_g, eb_g)}) or
    None if there are too few points."""
    def _prep(pts):
        return [(x, math.log(y), (y / ey) ** 2) for x, y, ey in pts if y > 0 and ey > 0]
    def _regress(sets):
        """common slope over several groups, each with its own intercept"""
        sxx = sxy = 0.0
        stats = []
        for pts in sets:
            sw = sum(w for _, _, w in pts)
            xb = sum(w * x for x, _, w in pts) / sw
            lb = sum(w * l for _, l, w in pts) / sw
            sxx += sum(w * (x - xb) ** 2 for x, _, w in pts)
            sxy += sum(w * (x - xb) * (l - lb) for x, l, w in pts)
            stats.append((pts, xb, lb))
        if sxx <= 0:
            return None
        b = sxy / sxx
        chi2 = sum(w * (l - lb - b * (x - xb)) ** 2 for pts, xb, lb in stats for x, l, w in pts)
        ndf = sum(len(pts) for pts in sets) - len(sets) - 1
        eb = (1.0 / sxx) ** 0.5 * max(1.0, (chi2 / ndf) ** 0.5 if ndf > 0 else 1.0)
        intercepts = [lb - b * xb for _, xb, lb in stats]
        return b, eb, intercepts
    prepared = {g: _prep(pts) for g, pts in groups.items()}
    prepared = {g: pts for g, pts in prepared.items() if len(pts) >= 3}
    if not prepared or sum(len(v) for v in prepared.values()) < 4:
        return None
    keys = list(prepared)
    common = _regress([prepared[g] for g in keys])
    if common is None:
        return None
    b, eb, icpt = common
    per_group = {}
    for g in keys:
        xs = [x for x, _, _ in prepared[g]]
        if max(xs) - min(xs) < MIN_FIT_SPAN_MM:
            continue
        r = _regress([prepared[g]])
        if r is not None:
            per_group[g] = (r[0], r[1])
    return b, eb, dict(zip(keys, icpt)), per_group


def _lambda_m(b, eb):
    """slope of ln(light) per mm -> attenuation length [m] and its error"""
    if b == 0:
        return float("inf"), float("inf")
    return 1.0 / b / 1000.0, eb / b / b / 1000.0


def _fmt_lambda(lam, err, unit=" m"):
    """lambda for a label: 'flat' when the slope is consistent with zero, which
    would otherwise print as a huge or negative length with a huge error."""
    if abs(lam) > 20.0 and err > abs(lam):
        return "flat"
    return f"{lam:.1f} #pm {err:.1f}{unit}" if abs(lam) >= 1.0 else f"{lam:.2f} #pm {err:.2f}{unit}"

# quantity key in intime_summary.json -> (error key or None, plot suffix, y title, y minimum or None)
# Per-fibre scan plots use a linear y axis fitted to the data so that the
# position dependence is visible; a fixed minimum (0) is used where it makes sense.
# All per-hit quantities are sideband-subtracted (random coincidences removed),
# see analysis/intime_analysis.py.
_PE_QUANTITIES = {
    "mean_energy_hit": ("sem_energy_hit", "drs_energy_hit_scan", "Light per hit [ADC]", None),
    "mean_energy_in": ("sem_energy_in", "drs_energy_in_scan", "<E> in time, raw [ADC]", None),
    "pe_per_hit": ("pe_per_hit_err", "drs_pe_scan", "p.e. per hit", None),
    "hit_fraction": ("hit_fraction_err", "drs_intime_frac_scan", "Hit fraction", 0.0),
}


def _collect_intime(run_x):
    """label -> quantity -> [(x, value, error)] from each run's intime_summary.json."""
    data = {}
    for run, x_val in sorted(run_x.items(), key=lambda kv: kv[1]):
        path = f"results/plots/Run{run}/DRS_InTime/intime_summary.json"
        if not os.path.exists(path):
            print(f"Warning: {path} not found, skipping run {run}")
            continue
        with open(path) as f:
            per_fibre = json.load(f)["per_fibre"]
        for label, v in per_fibre.items():
            if v.get("pe_per_hit") is None or not v.get("one_pe_int_adc"):
                if "mean_energy_hit" not in v:
                    print(f"Warning: run {run} has no sideband-subtracted values; rerun the in-time analysis")
                continue
            d = data.setdefault(label, {k: [] for k in _PE_QUANTITIES})
            for key, (err_key, _, _, _) in _PE_QUANTITIES.items():
                d[key].append((x_val, v[key], v[err_key]))
    return data


def _plot_pe_channels(data, pm, lumi, xmin, xmax):
    """Scan plots per fibre. Returns (plot names in page order, plateau values,
    exponential-fit results {label: {key: {...}}})."""
    names, plateau, expfits = [], {}, {}
    # Edge runs, where the beam only partly overlaps the fibres: the in-time
    # fraction of the high-yield fibres (plateau fraction >= 0.15, so that it
    # measures overlap rather than attenuation) is well below its plateau value.
    ratios = {}
    for label, d in data.items():
        frac_pts = sorted(d["hit_fraction"])
        f_plateau = _fit_constant([(x, y) for x, y, _ in frac_pts], PE_FIT_XMIN, PE_FIT_XMAX)[0]
        if f_plateau >= 0.15:
            for x, f, _ in frac_pts:
                ratios.setdefault(x, []).append(f / f_plateau)
    edge_x = {x for x, r in ratios.items() if sorted(r)[len(r) // 2] < EDGE_FRACTION}
    if edge_x:
        print(f"Exponential fits leave out the edge runs at X = {sorted(edge_x)} mm")
    for label in sorted(data, key=lambda l: (fiber_base_type(l), l)):
        plateau[label], expfits[label] = {}, {}
        for key, (_, suffix, ytitle, y_floor) in _PE_QUANTITIES.items():
            pts = sorted(data[label][key])
            if not pts:
                continue
            graphs, labels, markers = _graphs_by_beam(pts, f"gr_{suffix}_{label}")
            gr = graphs[0]

            p0, e0, n_used = _fit_constant([(x, y) for x, y, _ in pts], PE_FIT_XMIN, PE_FIT_XMAX)
            plateau[label][key] = (p0, e0, n_used)
            pave_fit = create_pave_text(0.13, 0.71, 0.60, 0.85, size=0.033)
            if key in EXP_FIT_KEYS:
                # exponential fit, one normalisation per beam energy, common lambda
                groups = {}   # beam energy -> points (edge runs left out)
                for x, y, ey in pts:
                    if x not in edge_x:
                        groups.setdefault(X_BEAM.get(x, ("", 0))[1], []).append((x, y, ey))
                res = _fit_exponential(groups)
                if res is not None:
                    b, eb, icpt, per_group = res
                    lam, elam = _lambda_m(b, eb)
                    expfits[label][key] = {
                        "lambda_m": lam, "lambda_err_m": elam, "slope_per_mm": b, "slope_err_per_mm": eb,
                        "per_energy": {g: dict(zip(("lambda_m", "lambda_err_m"), _lambda_m(*v)))
                                       for g, v in sorted(per_group.items())},
                        "n_points": sum(len(v) for v in groups.values()), "edge_runs_x": sorted(edge_x)}
                    for g, a in sorted(icpt.items()):
                        xs_g = [x for x, _, _ in groups[g]]
                        f = ROOT.TF1(f"expfit_{suffix}_{label}_{g}", "[0]*exp([1]*x)",
                                     min(xs_g) - 15, max(xs_g) + 15)
                        f.SetParameters(math.exp(a), b)
                        f.SetLineColor(ROOT.kRed + 1)
                        f.SetLineWidth(2)
                        f.SetNpx(200)
                        # attach to the graph of that energy so it is drawn with it
                        i_gr = next((i for i, lab in enumerate(labels) if lab.endswith(f"{g} GeV")), 0)
                        graphs[i_gr].GetListOfFunctions().Add(f)
                    pave_fit.AddText(f"#lambda = {_fmt_lambda(lam, elam)}  (common)")
                    pave_fit.AddText("   ".join(f"{g} GeV: {_fmt_lambda(*_lambda_m(*v), unit='')}"
                                               for g, v in sorted(per_group.items())) + " m")
                else:
                    pave_fit.AddText("exponential fit: too few points")
            else:
                fit = ROOT.TF1(f"fit_{suffix}_{label}", "pol0", PE_FIT_XMIN, PE_FIT_XMAX)
                fit.SetParameter(0, p0)
                fit.SetLineColor(ROOT.kRed + 1)
                fit.SetLineWidth(2)
                fit.SetLineStyle(2)
                gr.GetListOfFunctions().Add(fit)
                pave_fit.AddText(f"plateau: {p0:.3g} #pm {e0:.2g}")
            plot_name = f"{label}_{suffix}"
            lo = min(y - ye for _, y, ye in pts)
            hi = max(y + ye for _, y, ye in pts)
            pad = 0.35 * (hi - lo) if hi > lo else 0.2 * abs(hi) + 1.0
            ymin_p = lo - pad if y_floor is None else y_floor
            ymax_p = hi + 2.0 * pad     # room for the legend and the fit text above the data
            DrawHistos(
                graphs, labels,
                xmin, xmax, "Table position X [mm]",
                ymin_p, ymax_p, ytitle,
                plot_name,
                outdir=pm.get_output_dir(),
                drawoptions="EPL",
                markerstyles=markers,
                mycolors=[ROOT.kBlue + 1] * len(graphs),
                legendPos=LEGEND_POS,
                dology=False,
                usePDF=False,
                lumitext=lumi,
                extra_text=label,
                extraToDraw=[pave_fit],
            )
            names.append(plot_name)
    return names, plateau, expfits


def _print_pe_table(plateau, expfits, outdir):
    """Plateau values and attenuation lengths per fibre; saves JSON and markdown,
    returns the HTML table."""
    cols = [("Fibre", None), ("Type", None), ("light per hit [ADC]", "mean_energy_hit"),
            ("<E> in, raw [ADC]", "mean_energy_in"), ("p.e. per hit", "pe_per_hit"),
            ("hit frac.", "hit_fraction"), ("runs", None),
            ("lambda [m]", None), ("40 GeV", None), ("60 GeV", None)]
    def lam_str(d):
        return f"{d['lambda_m']:.2g} &pm; {d['lambda_err_m']:.2g}" if d else ""
    rows, out = [], {}
    for label in sorted(plateau, key=lambda l: (fiber_base_type(l), l)):
        v = plateau[label]
        if "pe_per_hit" not in v:
            continue
        out[label] = {k: {"plateau": round(v[k][0], 4), "err": round(v[k][1], 4), "n_runs": v[k][2]}
                      for k in v}
        out[label]["attenuation"] = expfits.get(label, {})
        e_hit = expfits.get(label, {}).get("mean_energy_hit", {})
        rows.append([label, fiber_base_type(label)]
                    + [f"{v[k][0]:.3g} &pm; {v[k][1]:.2g}" for k in
                       ("mean_energy_hit", "mean_energy_in", "pe_per_hit", "hit_fraction")]
                    + [str(v["pe_per_hit"][2]), lam_str(e_hit),
                       lam_str(e_hit.get("per_energy", {}).get(40)), lam_str(e_hit.get("per_energy", {}).get(60))])
    with open(os.path.join(outdir, "drs_pe_scan_channels.json"), "w") as f:
        json.dump(out, f, indent=2)
    md = ["| " + " | ".join(h for h, _ in cols) + " |", "|" + "---|" * len(cols)]
    md += ["| " + " | ".join(c.replace("&pm;", "+-") for c in r) + " |" for r in rows]
    with open(os.path.join(outdir, "drs_pe_scan_channels.md"), "w") as f:
        f.write(f"# Photoelectron scan, plateau {PE_FIT_XMIN:.0f} to {PE_FIT_XMAX:.0f} mm\n\n" + "\n".join(md) + "\n")
    print("\n".join(md))
    esc = lambda t: t.replace("<", "&lt;").replace(">", "&gt;")
    th = "".join(f"<th style='padding:2px 8px;border-bottom:1px solid #999'>{esc(h)}</th>" for h, _ in cols)
    trs = "".join("<tr>" + "".join(f"<td style='padding:2px 8px;text-align:right'>{c}</td>" for c in r) + "</tr>"
                  for r in rows)
    return (f"<p>Per run, from the in-time vs out-of-time analysis (DRS_InTime pages), with the random "
            f"coincidences inside the in-time window subtracted using the sideband far from the peak: light per "
            f"hit = (N in &lt;E&gt; in - N random &lt;E&gt; sideband) / (N in - N random), p.e. per hit = that / "
            f"1 p.e. integral of the fibre, hit fraction = (N in - N random) / N. Plateau = mean over runs with {PE_FIT_XMIN:.0f} &le; X &le; {PE_FIT_XMAX:.0f} mm, "
            f"error = standard error of that mean. The scintillating fibres saturate in part of their hits, so "
            f"their p.e. are lower limits. Attenuation: fit of light = A<sub>E</sub> exp(x / lambda) with one "
            f"normalisation per beam energy and one common lambda, also fitted per energy where the runs span "
            f"{MIN_FIT_SPAN_MM:.0f} mm; lambda &gt; 0 means "
            f"more light towards +X, the readout end, and a slope consistent with zero shows up as a very large "
            f"or negative lambda. Runs whose in-time fraction is below {EDGE_FRACTION:.0%} of the fibre's plateau "
            f"value (the edge runs) are drawn but not fitted. Errors are scaled by sqrt(chi2/ndf) where the "
            f"run-to-run scatter exceeds the point errors.</p>"
            f"<table style='border-collapse:collapse;font-size:13px;margin:10px 0'><thead><tr>{th}</tr></thead>"
            f"<tbody>{trs}</tbody></table>")


def run_pe_scan(run_x, lumi, xmin, xmax, use_jsroot=False):
    data = _collect_intime(run_x)
    if not data:
        print("No intime_summary.json found for any scan run; run check_drs_mcp.py --intime first.")
        return
    pm = _make_pm("DRS_PE_Scan", use_jsroot)
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    names, plateau, expfits = _plot_pe_channels(data, pm, lumi, xmin, xmax)
    table_html = _print_pe_table(plateau, expfits, pm.get_output_dir())

    # plateau summaries per fibre, one marker per fibre type
    def _vals(key):
        return {l: (v[key][0], v[key][1]) for l, v in plateau.items() if key in v}
    summaries = [
        plot_fiber_summary(_vals("mean_energy_hit"), "Light per hit [ADC] (plateau)", 10.0, 200000.0,
                           "drs_energy_hit_summary", pm, lumi, logy=True),
        plot_fiber_summary(_vals("pe_per_hit"), "p.e. per hit (plateau)", 0.05, 200.0,
                           "drs_pe_summary", pm, lumi, logy=True),
        plot_fiber_summary(_vals("hit_fraction"), "Hit fraction (plateau)", 0.0, 1.0,
                           "drs_intime_frac_summary", pm, lumi),
        # attenuation as 1/lambda, which stays linear and can be shown with its
        # error even when the slope is consistent with zero
        plot_fiber_summary({l: (1000.0 * f["mean_energy_hit"]["slope_per_mm"], 1000.0 * f["mean_energy_hit"]["slope_err_per_mm"])
                            for l, f in expfits.items() if "mean_energy_hit" in f},
                           "1 / #lambda [m^{-1}]   (light = A e^{x/#lambda})", -1.0, 5.0,
                           "drs_attenuation_summary", pm, lumi, note="from the light per hit, randoms subtracted"),
    ]
    pm.reset_plots()
    for name in summaries:
        if name:
            pm.add_plot(name)
    pm.add_newline()
    for name in names:
        pm.add_plot(name)
    pm.generate_html(
        "DRS_PE_Scan_channels.html",
        plots_per_row=4,
        title="DRS photoelectron scan — selected channels",
        intro_text=table_html,
    )


def _parse_args():
    parser = argparse.ArgumentParser(
        description="DRS MPV / Energy position scan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mode", choices=["time", "energy", "pe", "both", "all"], default="both",
        help="Which scan to run: time (CFD MPV), energy (SumMean), pe (photoelectrons "
             "from the in-time analysis, needs --channels), both (time + energy), or all",
    )
    parser.add_argument(
        "--jsroot", action="store_true",
        help="Embed interactive JSROOT canvases in HTML output",
    )
    parser.add_argument(
        "--channels", default=None, metavar="FILE",
        help="JSON file mapping labels to channel names "
             "(e.g. channels/data/testingfibers.json). "
             "When given, only the listed channels are processed.",
    )
    return parser.parse_args()


def main():
    args = _parse_args()
    setup_root(batch_mode=True, load_functions=(args.mode != "energy"))

    with open("data/Runlist.json") as f:
        runlist = json.load(f)

    _X_OFFSET_FROM_RUN = 2012
    _X_OFFSET          = -482

    def _parse_x(run):
        raw = runlist[str(run)]["Table position X"]
        x   = float(raw.replace("mm", "").strip())
        if run >= _X_OFFSET_FROM_RUN:
            x += _X_OFFSET
        return x

    run_x = {r: _parse_x(r) for r in SCAN_RUNS if str(r) in runlist}
    if not run_x:
        print("No table positions found for the scan runs in data/Runlist.json.")
        return

    # beam type and energy per run, keyed by table position for the plots
    def _beam(run):
        e = runlist[str(run)]
        try:
            energy = int(float(str(e.get("beam energy", "")).replace("GeV", "").strip()))
        except ValueError:
            energy = 0
        return str(e.get("beam type", "")).strip(), energy
    processed = [r for r in sorted(run_x) if os.path.exists(f"results/root/Run{r}/drs_stats.root")]
    for r in processed:
        x = run_x[r]
        if x in X_BEAM and X_BEAM[x] != _beam(r):
            print(f"Warning: runs at X = {x} mm have different beams; markers follow the last one")
        X_BEAM[x] = _beam(r)

    drsboards     = build_drs_boards(run_number=REFERENCE_RUN)
    channel_order = _build_channel_order(drsboards)

    ch_label_map = {}  # channel_name -> label (populated when --channels is used)
    if args.channels:
        with open(args.channels) as f:
            ch_map = json.load(f)
        ch_label_map = {v: k for k, v in ch_map.items()}
        allowed = set(ch_map.values())
        channel_order = [e for e in channel_order if e[0] in allowed]
        print(f"Channel filter: {len(allowed)} channels from {args.channels}")

    xs_all = list(run_x.values())
    xmin, xmax = min(xs_all) - 20, max(xs_all) + 20

    btype, benergy = getRunInfo(REFERENCE_RUN)
    energies = sorted({e for _, e in X_BEAM.values() if e})
    lumi = (f"{_BTYPES.get(btype.lower(), btype.lower())}, {benergy} GeV, 90^{{#circ}}" if len(energies) <= 1
            else f"{_BTYPES.get(btype.lower(), btype.lower())} {'/'.join(str(e) for e in energies)} GeV, 90^{{#circ}}")

    single_page = args.channels is not None
    if args.mode in ("time", "both", "all"):
        run_time_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                      use_jsroot=args.jsroot, single_page=single_page,
                      ch_label_map=ch_label_map if ch_label_map else None)
    if args.mode in ("energy", "both", "all"):
        run_energy_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                        use_jsroot=args.jsroot, single_page=single_page,
                        ch_label_map=ch_label_map if ch_label_map else None)
    if args.mode in ("pe", "all"):
        run_pe_scan(run_x, lumi, xmin, xmax, use_jsroot=args.jsroot)

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
