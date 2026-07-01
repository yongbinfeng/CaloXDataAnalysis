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
import os
import json
import ROOT

from channels.channel_map import build_drs_boards
from core.plot_manager import PlotManager
from plotting.my_function import DrawHistos
from plotting.calox_plot_helper import create_pave_text
from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var, get_hist_mpv
auto_timer("Total Execution Time")

SCAN_RUNS = list(range(1998, 2012))
REFERENCE_RUN = 1998

PLOTDIR = "results/plots/PositionScan"
HTMLDIR = "results/html/PositionScan"



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
    """Read MPV of fine-binned TS_cfd_mcp histogram for each run/channel."""
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
                hist = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
                if hist and hist.GetEntries() >= 5:
                    mpv, mpv_err = get_hist_mpv(hist)
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

        gr = ROOT.TGraphErrors(len(pts))
        gr.SetName(f"gr_mpv_{ch}")
        gr.SetTitle("")
        for i, (x, y, ye) in enumerate(pts):
            gr.SetPoint(i, x, y)
            gr.SetPointError(i, 0.0, ye)
        gr.SetLineWidth(2)

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
            [gr], [],
            xmin, xmax, "Table position X [mm]",
            ymid - yhalf, ymid + yhalf, "MPV [TS]",
            plot_name,
            outdir=pm.get_output_dir(),
            drawoptions="EPL",
            markerstyles=[20],
            mycolors=[ROOT.kBlue + 1],
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

        gr = ROOT.TGraph(len(pts))
        gr.SetName(f"gr_summean_{ch}")
        gr.SetTitle("")
        for i, (x, mean) in enumerate(pts):
            gr.SetPoint(i, x, mean)
        gr.SetLineWidth(2)

        # Constant (pol0) fit restricted to [-220, 50] mm
        fit = ROOT.TF1(f"fit_sum_{ch}", "pol0", FIT_XMIN, FIT_XMAX)
        fit.SetLineColor(ROOT.kRed + 1)
        fit.SetLineWidth(2)
        fit.SetLineStyle(2)
        gr.Fit(fit, "QR")
        p0 = fit.GetParameter(0)
        e0 = fit.GetParError(0)
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
            [gr], [],
            xmin, xmax, "Table position X [mm]",
            ymid - yhalf, ymid + yhalf, "Mean DRS Sum [ADC]",
            plot_name,
            outdir=pm.get_output_dir(),
            drawoptions="EPL",
            markerstyles=[20],
            mycolors=[ROOT.kBlue + 1],
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

_SPEED_PALETTE = [
    ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kOrange + 1,
    ROOT.kViolet + 2, ROOT.kCyan + 2, ROOT.kMagenta + 1, ROOT.kTeal + 3,
    ROOT.kAzure + 7, ROOT.kSpring + 6,
]
# ROOT marker styles: filled circle, square, triangle-up, diamond, triangle-down,
# open circle, open square, open triangle-up, star, open diamond
_SPEED_MARKERS = [20, 21, 22, 33, 23, 24, 25, 26, 29, 27]


def _fiber_base_type(label):
    """Strip trailing _N suffix: 'Quartz400_1' → 'Quartz400', 'Sapphire' → 'Sapphire'."""
    import re
    return re.sub(r'_\d+$', '', label)


def _plot_speed_summary(ch_label_map, p1_map, e1_map, pm, lumi):
    """TGraphErrors of speed (cm/ns) per fiber, one color per base fiber type.

    Returns the plot name (without extension) so the caller can add it to pm.
    """
    # Compute speed and uncertainty for each labeled channel
    speed_data = {}   # label -> (speed, e_speed)
    for ch, label in ch_label_map.items():
        p1 = p1_map.get(ch)
        e1 = e1_map.get(ch)
        if not p1:
            continue
        s  = 0.5 / abs(p1)
        es = s * abs(e1) / abs(p1)
        speed_data[label] = (s, es)
    if not speed_data:
        return None

    # Sort by (base type, full label) so same-type channels are adjacent
    sorted_labels = sorted(speed_data, key=lambda l: (_fiber_base_type(l), l))
    types_ordered = list(dict.fromkeys(_fiber_base_type(l) for l in sorted_labels))
    type_color    = {t: _SPEED_PALETTE[i % len(_SPEED_PALETTE)]
                     for i, t in enumerate(types_ordered)}
    type_marker   = {t: _SPEED_MARKERS[i % len(_SPEED_MARKERS)]
                     for i, t in enumerate(types_ordered)}

    # Group points by base type
    pts_by_type = {}
    for i, label in enumerate(sorted_labels):
        ft = _fiber_base_type(label)
        s, es = speed_data[label]
        pts_by_type.setdefault(ft, []).append((i + 1, s, es))

    n = len(sorted_labels)
    ymin, ymax = 12.0, 24.0

    # Canvas
    c = ROOT.TCanvas("c_speed_summary", "", 1400, 650)
    c.SetBottomMargin(0.28)
    c.SetLeftMargin(0.09)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)

    # Frame histogram — only purpose is to carry x-axis bin labels
    frame = ROOT.TH1F("h_speed_frame", "", n, 0.5, n + 0.5)
    for i, label in enumerate(sorted_labels):
        frame.GetXaxis().SetBinLabel(i + 1, label)
    frame.GetXaxis().LabelsOption("v")
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetXaxis().SetTickLength(0.0)
    frame.GetYaxis().SetTitle("Speed  v = 1/|p_{1}| [cm/ns]")
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.048)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.SetStats(0)
    frame.Draw("AXIS")

    # One TGraphErrors per base type
    graphs = []
    leg = ROOT.TLegend(0.10, 0.84, 0.90, 0.90)
    leg.SetNColumns(min(len(types_ordered), 8))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.036)

    for ft in types_ordered:
        pts = pts_by_type[ft]
        gr  = ROOT.TGraphErrors(len(pts))
        col = type_color[ft]
        mkr = type_marker[ft]
        for j, (x, s, es) in enumerate(pts):
            gr.SetPoint(j, x, s)
            gr.SetPointError(j, 0.0, es)
        gr.SetMarkerStyle(mkr)
        gr.SetMarkerSize(1.3)
        gr.SetMarkerColor(col)
        gr.SetLineColor(col)
        gr.SetLineWidth(2)
        gr.Draw("P SAME")
        leg.AddEntry(gr, ft, "p")
        graphs.append(gr)

    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.038)
    lat.SetTextFont(42)
    lat.SetTextAlign(31)   # right-align
    lat.DrawLatex(0.96, 0.935, lumi)

    plot_name = "drs_cfd_speed_summary"
    outdir = pm.get_output_dir()

    if pm.use_jsroot:
        pm._canvas_jsons[plot_name] = ROOT.TBufferJSON.ToJSON(c).Data()
    else:
        c.SaveAs(os.path.join(outdir, f"{plot_name}.png"))

    # Keep objects alive until after SaveAs
    c._keep = [frame, leg, lat] + graphs

    return plot_name


# ---------------------------------------------------------------------------
# Energy p0 summary plot
# ---------------------------------------------------------------------------

def _plot_p0_summary_energy(ch_label_map, p0_map, e0_map, pm, lumi):
    """TGraphErrors of mean DRS Sum p0 per fiber, one color+marker per base type.

    Returns the plot name so the caller can add it to pm.
    """
    p0_data = {}
    for ch, label in ch_label_map.items():
        p0 = p0_map.get(ch)
        e0 = e0_map.get(ch)
        if p0 is None:
            continue
        p0_data[label] = (p0, e0 if e0 is not None else 0.0)
    if not p0_data:
        return None

    sorted_labels = sorted(p0_data, key=lambda l: (_fiber_base_type(l), l))
    types_ordered = list(dict.fromkeys(_fiber_base_type(l) for l in sorted_labels))
    type_color  = {t: _SPEED_PALETTE[i % len(_SPEED_PALETTE)]  for i, t in enumerate(types_ordered)}
    type_marker = {t: _SPEED_MARKERS[i % len(_SPEED_MARKERS)]  for i, t in enumerate(types_ordered)}

    pts_by_type = {}
    for i, label in enumerate(sorted_labels):
        ft = _fiber_base_type(label)
        p0, e0 = p0_data[label]
        pts_by_type.setdefault(ft, []).append((i + 1, p0, e0))

    n = len(sorted_labels)
    ymin, ymax = 500.0, 1700.0

    c = ROOT.TCanvas("c_p0_energy_summary", "", 1400, 650)
    c.SetBottomMargin(0.28)
    c.SetLeftMargin(0.09)
    c.SetRightMargin(0.04)
    c.SetTopMargin(0.08)

    frame = ROOT.TH1F("h_p0_energy_frame", "", n, 0.5, n + 0.5)
    for i, label in enumerate(sorted_labels):
        frame.GetXaxis().SetBinLabel(i + 1, label)
    frame.GetXaxis().LabelsOption("v")
    frame.GetXaxis().SetLabelSize(0.042)
    frame.GetXaxis().SetTickLength(0.0)
    frame.GetYaxis().SetTitle("Mean DRS Sum  p_{0} [ADC]")
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetTitleSize(0.048)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    frame.SetStats(0)
    frame.Draw("AXIS")

    graphs = []
    leg = ROOT.TLegend(0.10, 0.84, 0.90, 0.90)
    leg.SetNColumns(min(len(types_ordered), 8))
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.036)

    for ft in types_ordered:
        pts = pts_by_type[ft]
        gr  = ROOT.TGraphErrors(len(pts))
        col = type_color[ft]
        mkr = type_marker[ft]
        for j, (x, p0, e0) in enumerate(pts):
            gr.SetPoint(j, x, p0)
            gr.SetPointError(j, 0.0, e0)
        gr.SetMarkerStyle(mkr)
        gr.SetMarkerSize(1.3)
        gr.SetMarkerColor(col)
        gr.SetLineColor(col)
        gr.SetLineWidth(2)
        gr.Draw("P SAME")
        leg.AddEntry(gr, ft, "p")
        graphs.append(gr)

    leg.Draw()

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.038)
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    lat.DrawLatex(0.96, 0.935, lumi)

    plot_name = "drs_summean_p0_summary"
    outdir = pm.get_output_dir()

    if pm.use_jsroot:
        pm._canvas_jsons[plot_name] = ROOT.TBufferJSON.ToJSON(c).Data()
    else:
        c.SaveAs(os.path.join(outdir, f"{plot_name}.png"))

    c._keep = [frame, leg, lat] + graphs
    return plot_name


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

        # Speed summary — drawn first so it appears at the top of the page
        speed_plot = None
        if ch_label_map and p1_map:
            speed_plot = _plot_speed_summary(ch_label_map, p1_map, e1_map, pm, lumi)

        pm.reset_plots()
        if speed_plot:
            pm.add_plot(speed_plot)
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

def _parse_args():
    parser = argparse.ArgumentParser(
        description="DRS MPV / Energy position scan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mode", choices=["time", "energy", "both"], default="both",
        help="Which scan to run: time (CFD MPV), energy (SumMean), or both",
    )
    parser.add_argument(
        "--jsroot", action="store_true",
        help="Embed interactive JSROOT canvases in HTML output",
    )
    parser.add_argument(
        "--channels", default=None, metavar="FILE",
        help="JSON file mapping labels to channel names "
             "(e.g. data/channel_maps/testingfibers.json). "
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

    _btypes = {"pion": "#pi^{+}", "pions": "#pi^{+}", "pi+": "#pi^{+}",
               "positron": "e^{+}", "positrons": "e^{+}", "e+": "e^{+}",
               "mu+": "#mu^{+}"}
    btype, benergy = getRunInfo(REFERENCE_RUN)
    lumi = f"{_btypes.get(btype.lower(), btype.lower())}, {benergy} GeV, 90^{{#circ}}"

    single_page = args.channels is not None
    if args.mode in ("time", "both"):
        run_time_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                      use_jsroot=args.jsroot, single_page=single_page,
                      ch_label_map=ch_label_map if ch_label_map else None)
    if args.mode in ("energy", "both"):
        run_energy_scan(drsboards, run_x, channel_order, lumi, xmin, xmax,
                        use_jsroot=args.jsroot, single_page=single_page,
                        ch_label_map=ch_label_map if ch_label_map else None)

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
