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
from plotting.calox_plot_helper import create_pave_text, BoardPlotHelper
from utils.visualization import visualizeDRSBoards
from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_channel_var, get_hist_mpv
from variables.drs import subtract_type_mpv

auto_timer("Total Execution Time")

SCAN_RUNS = list(range(1897, 1932))
REFERENCE_RUN = 1897

PLOTDIR = "results/plots/PositionScan"
HTMLDIR = "results/html/PositionScan"

FIT_SHARED_SLOPE = True
SAVE_T0_JSON = True

Y1500_EVAL  = 1500
Y1500_SHIFT = 300

_TYPE_ORDER  = ["CerQuartz", "CerPlastic", "Sci"]
_TYPE_COLORS = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2]


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

def _compute_shared_slopes(channel_order, mpv_data):
    """Pre-compute one weighted slope per channel type (shared-slope model)."""
    shared_fit = {}
    for var_name in dict.fromkeys(v for _, _, _, _, v, _, _, _ in channel_order):
        chs_var = [(ch, sorted(mpv_data.get(ch, [])))
                   for ch, _, _, _, v, _, _, _ in channel_order if v == var_name]
        SXX, SXY = 0.0, 0.0
        for ch, pts in chs_var:
            if len(pts) < 2:
                continue
            w_c    = [1.0 / ye**2 for _, _, ye in pts]
            W_c    = sum(w_c)
            xbar_w = sum(wi * x for wi, (x, _, _) in zip(w_c, pts)) / W_c
            ybar_w = sum(wi * y for wi, (_, y, _) in zip(w_c, pts)) / W_c
            for wi, (x, y, _) in zip(w_c, pts):
                SXX += wi * (x - xbar_w) ** 2
                SXY += wi * (x - xbar_w) * (y - ybar_w)
        if SXX == 0:
            continue
        shared_fit[var_name] = (SXY / SXX, 1.0 / SXX ** 0.5, SXX)
    return shared_fit


def _fit_channel(pts, var, shared_fit):
    """Weighted least-squares linear fit. Returns (p0, p1, e0, e1)."""
    w = [1.0 / ye**2 for _, _, ye in pts]
    if FIT_SHARED_SLOPE and var in shared_fit:
        p1, e1, SXX_cat = shared_fit[var]
        W_c    = sum(w)
        xbar_w = sum(wi * x for wi, (x, _, _) in zip(w, pts)) / W_c
        ybar_w = sum(wi * y for wi, (_, y, _) in zip(w, pts)) / W_c
        p0 = ybar_w - p1 * xbar_w
        e0 = (1.0 / W_c + xbar_w ** 2 / SXX_cat) ** 0.5
    else:
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

def _plot_time_channels(channel_order, mpv_data, shared_fit, pm, lumi, xmin, xmax):
    """Draw one MPV-vs-X TGraph per channel. Returns (p0_map, p1_map, plots_by_type)."""
    plots_by_type = {}
    p0_map, p1_map = {}, {}

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

        p0, p1, e0, e1 = _fit_channel(pts, var, shared_fit)

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

        pave_bgc = create_pave_text(0.18, 0.77, 0.60, 0.87)
        pave_bgc.AddText(f"B: {b}, G: {g}, C: {c}")
        pave_bgc.AddText(f"Tower: ({tx}, {ty})")

        pave_fit = create_pave_text(0.18, 0.18, 0.60, 0.28)
        pave_fit.AddText(f"p0: {p0:.1f} #pm {e0:.1f} TS")
        pave_fit.AddText(f"p1: {p1:.3f} #pm {e1:.3f} TS/mm")

        plot_name = f"drs_cfd_mpv_scan_{ch}_{var}"
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
            extra_text=var,
            extraToDraw=[pave_bgc, pave_fit],
        )
        plots_by_type.setdefault(var, []).append(plot_name)
        p0_map[ch] = p0
        p1_map[ch] = p1

    return p0_map, p1_map, plots_by_type


# ---------------------------------------------------------------------------
# Time scan — summary pages
# ---------------------------------------------------------------------------

def _plot_fit_params_page(drsboards, p0_map, p1_map, channel_order, lumi, use_jsroot=False):
    """Board maps of p0, |p1|×100, 1/|p1| [cm/ns], plus speed 1D histogram."""
    p0_shifted = {ch: v - 400          for ch, v in p0_map.items()}
    p1_abs_map = {ch: abs(v) * 100     for ch, v in p1_map.items()}
    p1_inv_map = {ch: 0.5 / abs(v)     for ch, v in p1_map.items() if v != 0}

    pm_map = _make_pm("DRS_CFD_MPV_Scan", use_jsroot)
    helper = BoardPlotHelper(pm_map, run_number=REFERENCE_RUN)

    cer_hists, sci_hists = visualizeDRSBoards(drsboards, valuemaps=p0_shifted, suffix="MPV_scan_p0")
    helper.plot_cer_sci_pair(cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p0",
                             zmin=None, zmax=None, nTextDigits=1, lumitext=lumi, usePDF=False)

    cer_hists, sci_hists = visualizeDRSBoards(drsboards, valuemaps=p1_abs_map, suffix="MPV_scan_p1_abs")
    helper.plot_cer_sci_pair(cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p1_abs",
                             zmin=None, zmax=None, nTextDigits=2, lumitext=lumi, usePDF=False)

    p1i_sorted = sorted(p1_inv_map.values())
    n_inv = len(p1i_sorted)
    p5  = p1i_sorted[max(0, int(0.10 * n_inv))]
    p95 = p1i_sorted[min(n_inv - 1, int(0.90 * n_inv))]
    cer_hists, sci_hists = visualizeDRSBoards(drsboards, valuemaps=p1_inv_map, suffix="MPV_scan_p1_inv")
    helper.plot_cer_sci_pair(cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p1_inv",
                             zmin=p5, zmax=p95, nTextDigits=1, lumitext=lumi, usePDF=False)

    speed_by_type = {var: [] for var in _TYPE_ORDER}
    for ch, _b, _g, _c, var, _tx, _ty, _is6mm in channel_order:
        if ch in p1_inv_map:
            speed_by_type[var].append(p1_inv_map[ch])

    xlo_h, xhi_h = 12.0, 22.0
    hists_speed, speed_labels, speed_cols = [], [], []
    for tname, col in zip(_TYPE_ORDER, _TYPE_COLORS):
        vals = speed_by_type.get(tname, [])
        if not vals:
            continue
        h = ROOT.TH1F(f"h_speed_{tname}", "", 20, xlo_h, xhi_h)
        for v in vals:
            h.Fill(v)
        hists_speed.append(h)
        speed_labels.append(tname)
        speed_cols.append(col)

    if hists_speed:
        DrawHistos(
            hists_speed, speed_labels,
            xlo_h, xhi_h, "Speed [cm/ns]",
            None, None, "Channels",
            "drs_cfd_speed_distribution",
            outdir=pm_map.get_output_dir(),
            dology=False, usePDF=False,
            mycolors=speed_cols,
            drawoptions=["C"] * len(hists_speed),
            lumitext=lumi,
        )
        pm_map.add_plot("drs_cfd_speed_distribution")

    run_list = ", ".join(str(r) for r in SCAN_RUNS)
    pm_map.generate_html(
        "DRS_CFD_MPV_Scan_FitParams.html",
        plots_per_row=2,
        title="DRS CFD MPV Scan — Fit Parameters",
        intro_text=(
            f"**Linear fit parameters from the MPV position scan** "
            f"(runs {run_list}, {lumi}).\n"
            "For each channel the MPV of the fine-binned MCP-corrected CFD "
            "time slice was fitted as a linear function of table position X: "
            "MPV(X) = p0 + p1 * X.\n\n"
            "* **p0 − 400 [TS]**: intercept shifted by 400 for readability.\n"
            "* **|p1| × 100 [TS/mm × 100]**: absolute slope scaled by 100.\n"
            "* **1/|p1| [cm/ns]**: effective speed of light in the medium "
            "(1 TS = 200 ps).\n\n"
            "Channels with no valid data across all runs are not filled."
        ),
    )


def _plot_y_at_x_page(drsboards, p0_map, p1_map, channel_order, lumi,
                      x_eval=Y1500_EVAL, shift=Y1500_SHIFT, use_jsroot=False):
    """Board maps + 1D distribution of MPV evaluated at x_eval mm."""
    y_map = {ch: p0_map[ch] + p1_map[ch] * x_eval - shift for ch in p0_map}

    pm_y = _make_pm("DRS_CFD_MPV_Scan", use_jsroot)
    helper_y = BoardPlotHelper(pm_y, run_number=REFERENCE_RUN)

    cer_hists, sci_hists = visualizeDRSBoards(
        drsboards, valuemaps=y_map, suffix=f"MPV_scan_y{x_eval}")
    helper_y.plot_cer_sci_pair(
        cer_hists, sci_hists, f"DRS_CFD_MPV_Scan_y{x_eval}",
        zmin=None, zmax=None, nTextDigits=1, lumitext=lumi, usePDF=False)

    _size_order = ["3mm", "6mm"]
    _size_type_colors = {
        ("3mm", "CerQuartz"):  ROOT.kBlue + 1,
        ("3mm", "CerPlastic"): ROOT.kRed + 1,
        ("3mm", "Sci"):        ROOT.kGreen + 2,
        ("6mm", "CerQuartz"):  ROOT.kAzure + 7,
        ("6mm", "CerPlastic"): ROOT.kOrange + 1,
        ("6mm", "Sci"):        ROOT.kTeal + 3,
    }

    y_by_size_type = {}
    y_3mm_cerq_b5 = []
    for ch, b, _g, _c, var, _tx, _ty, is6mm in channel_order:
        if ch in y_map:
            size = "6mm" if is6mm else "3mm"
            y_by_size_type.setdefault((size, var), []).append(y_map[ch])
            if not is6mm and var == "CerQuartz" and b == 5:
                y_3mm_cerq_b5.append(y_map[ch])

    xlo_y, xhi_y = 73, 105
    hists_y, y_labels, y_cols = [], [], []
    group_mpv = {}
    for size in _size_order:
        for tname in _TYPE_ORDER:
            vals = y_by_size_type.get((size, tname), [])
            if not vals:
                continue
            h = ROOT.TH1F(f"h_y{x_eval}_{size}_{tname}", "",
                          int((xhi_y - xlo_y) * 2), xlo_y, xhi_y)
            for v in vals:
                h.Fill(v)
            mpv, _ = get_hist_mpv(h, window_ts=3.0)
            group_mpv[f"{size}_{tname}"] = mpv + shift
            hists_y.append(h)
            y_labels.append(f"{size} {tname} (MPV = {mpv:.1f})")
            y_cols.append(_size_type_colors[(size, tname)])

    if y_3mm_cerq_b5:
        h_b5 = ROOT.TH1F(f"h_y{x_eval}_3mm_CerQuartz_B5", "",
                         int((xhi_y - xlo_y) * 2), xlo_y, xhi_y)
        for v in y_3mm_cerq_b5:
            h_b5.Fill(v)
        mpv_b5, _ = get_hist_mpv(h_b5, window_ts=3.0)
        group_mpv["3mm_CerQuartz_B5"] = mpv_b5 + shift
        hists_y.append(h_b5)
        y_labels.append(f"3mm CerQuartz B5 (MPV = {mpv_b5:.1f})")
        y_cols.append(ROOT.kViolet + 1)

    if hists_y:
        DrawHistos(
            hists_y, y_labels,
            xlo_y, xhi_y, f"MPV at X = {x_eval} mm - {shift} [TS]",
            0, 25, "# Channels",
            f"drs_cfd_mpv_at_{x_eval}mm",
            outdir=pm_y.get_output_dir(),
            dology=False, usePDF=False,
            mycolors=y_cols,
            drawoptions=["HIST"] * len(hists_y),
            lumitext=lumi,
            addOverflow=True,
            addUnderflow=True,
            legendPos=[0.20, 0.60, 0.60, 0.90],
        )
        pm_y.add_plot(f"drs_cfd_mpv_at_{x_eval}mm")

    _cer_merged_colors = {
        ("3mm", "Cer"): ROOT.kBlue + 1,
        ("3mm", "Sci"): ROOT.kGreen + 2,
        ("6mm", "Cer"): ROOT.kAzure + 7,
        ("6mm", "Sci"): ROOT.kTeal + 3,
    }
    y_by_size_cer = {}
    for ch, b, _g, _c, var, _tx, _ty, is6mm in channel_order:
        if ch not in y_map:
            continue
        size = "6mm" if is6mm else "3mm"
        cer_key = "Cer" if var in ("CerQuartz", "CerPlastic") else "Sci"
        y_by_size_cer.setdefault((size, cer_key), []).append(y_map[ch])

    _cer_order = [("3mm", "Cer"), ("3mm", "Sci"), ("6mm", "Cer"), ("6mm", "Sci")]
    hists_cer, cer_labels, cer_cols = [], [], []
    group_mpv_cer = {}
    for size, cer_key in _cer_order:
        vals = y_by_size_cer.get((size, cer_key), [])
        if not vals:
            continue
        h = ROOT.TH1F(f"h_y{x_eval}_{size}_{cer_key}", "",
                      int((xhi_y - xlo_y) * 2), xlo_y, xhi_y)
        for v in vals:
            h.Fill(v)
        mpv, _ = get_hist_mpv(h, window_ts=3.0)
        group_mpv_cer[f"{size}_{cer_key}"] = mpv + shift
        hists_cer.append(h)
        cer_labels.append(f"{size} {cer_key} (MPV = {mpv:.1f})")
        cer_cols.append(_cer_merged_colors[(size, cer_key)])

    if y_3mm_cerq_b5:
        h_b5_cer = ROOT.TH1F(f"h_y{x_eval}_3mm_CerQuartz_B5_cer", "",
                             int((xhi_y - xlo_y) * 2), xlo_y, xhi_y)
        for v in y_3mm_cerq_b5:
            h_b5_cer.Fill(v)
        mpv_b5_cer, _ = get_hist_mpv(h_b5_cer, window_ts=3.0)
        group_mpv_cer["3mm_CerQuartz_B5"] = mpv_b5_cer + shift
        hists_cer.append(h_b5_cer)
        cer_labels.append(f"3mm CerQuartz B5 (MPV = {mpv_b5_cer:.1f})")
        cer_cols.append(ROOT.kViolet + 1)

    if group_mpv_cer:
        out_path_cer = f"data/drs/drs_cfd_mpv_at_{x_eval}mm_by_group.json"
        os.makedirs(os.path.dirname(out_path_cer), exist_ok=True)
        with open(out_path_cer, "w") as f:
            json.dump(group_mpv_cer, f, indent=2)
        print(f"Saved group MPVs to {out_path_cer}")

    if hists_cer:
        DrawHistos(
            hists_cer, cer_labels,
            xlo_y, xhi_y, f"MPV at X = {x_eval} mm - {shift} [TS]",
            0, 25, "# Channels",
            f"drs_cfd_mpv_at_{x_eval}mm_by_group",
            outdir=pm_y.get_output_dir(),
            dology=False, usePDF=False,
            mycolors=cer_cols,
            drawoptions=["HIST"] * len(hists_cer),
            lumitext=lumi,
            addOverflow=True,
            addUnderflow=True,
            legendPos=[0.20, 0.60, 0.60, 0.90],
        )
        pm_y.add_plot(f"drs_cfd_mpv_at_{x_eval}mm_by_group")

    pm_y.generate_html(
        f"DRS_CFD_MPV_Scan_y{x_eval}.html",
        plots_per_row=2,
        title=f"DRS CFD MPV at X = {x_eval} mm",
        intro_text=(
            f"Per-channel MPV evaluated at X = {x_eval} mm, shifted by −{shift} TS. "
            f"Computed from the linear fit: MPV({x_eval}) − {shift} = p0 + p1 × {x_eval} − {shift}."
        ),
    )


def _save_t0_json(p0_map, p1_map, x_eval=Y1500_EVAL, shift=Y1500_SHIFT):
    """Save raw MPV(x_eval) per channel to data/drs/drs_cfd_t0.json."""
    t0_map = {ch: (p0_map[ch] + p1_map[ch] * x_eval - shift) + shift for ch in p0_map}
    out_path = "data/drs/drs_cfd_t0.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(t0_map, f, indent=2)
    print(f"Saved {len(t0_map)} t0 values to {out_path}")


def _plot_p0_corrected2_page(drsboards, p0_map, p1_map, channel_order, lumi, use_jsroot=False):
    """Board maps of (p0 − group_MPV) − 1500×p1: CerQuartz, CerPlastic, Cer, Sci."""
    base_corr_map = subtract_type_mpv(drsboards, p0_map)
    if base_corr_map is None:
        print("Skipping p0-corrected2 page: MPV JSON not found")
        return

    ch_meta = {
        ch: (var, is6mm)
        for ch, _b, _g, _c, var, _tx, _ty, is6mm in channel_order
    }

    _cer_types  = {"CerQuartz", "CerPlastic"}
    _type_board = {"CerQuartz": "cer", "CerPlastic": "cer", "Sci": "sci"}

    def _corrected2(ch):
        val = base_corr_map[ch] - 1500.0 * abs(p1_map.get(ch, 0.0)) + 3.5
        return val * 0.2

    def _zrange(vals):
        vs = sorted(vals)
        n = len(vs)
        return vs[max(0, int(0.12 * n))], vs[min(n - 1, int(0.95 * n))]

    pm2 = _make_pm("DRS_CFD_MPV_Scan", use_jsroot)
    helper = BoardPlotHelper(pm2, run_number=REFERENCE_RUN)

    for tname in _TYPE_ORDER:
        filtered = {
            ch: _corrected2(ch)
            for ch, (var, _is6mm) in ch_meta.items()
            if ch in base_corr_map and var == tname
        }
        if not filtered:
            continue
        zmin_val, zmax_val = _zrange(list(filtered.values()))
        if tname in _cer_types:
            for ch, (var, _is6mm) in ch_meta.items():
                if var in _cer_types - {tname} and ch in p0_map:
                    filtered[ch] = 0.0
        cer_hists, sci_hists = visualizeDRSBoards(
            drsboards, valuemaps=filtered, suffix=f"MPV_scan_p0c2_{tname}")
        hists = cer_hists if _type_board[tname] == "cer" else sci_hists
        helper.plot_board_map(
            hists, f"DRS_CFD_MPV_Scan_p0c2_{tname}",
            extra_text=tname, zmin=zmin_val, zmax=zmax_val, nTextDigits=1,
            zlabel="Time [ns]", lumitext=lumi, usePDF=False)

    cer_combined = {
        ch: _corrected2(ch)
        for ch, (var, is6mm) in ch_meta.items()
        if ch in base_corr_map and var in _cer_types
        and not (is6mm and var == "CerPlastic")
    }
    if cer_combined:
        zmin_val, zmax_val = _zrange(list(cer_combined.values()))
        for ch, (var, is6mm) in ch_meta.items():
            if is6mm and var == "CerPlastic":
                cer_combined[ch] = 0.0
        cer_hists, _ = visualizeDRSBoards(
            drsboards, valuemaps=cer_combined, suffix="MPV_scan_p0c2_Cer")
        _aw = 1.2
        _ah = 0.55
        _tip_x, _tip_y = -6.5, -1.0
        beam_shaft = ROOT.TLine(-9.2, _tip_y, _tip_x, _tip_y)
        beam_top   = ROOT.TLine(_tip_x, _tip_y, _tip_x - _aw, _tip_y + _ah)
        beam_bot   = ROOT.TLine(_tip_x, _tip_y, _tip_x - _aw, _tip_y - _ah)
        for obj in (beam_shaft, beam_top, beam_bot):
            obj.SetLineWidth(2)
            obj.SetLineColor(ROOT.kBlack)
        beam_label = ROOT.TLatex(-7.85, -0.05, "e^{+}, 40 GeV")
        beam_label.SetTextAlign(22)
        beam_label.SetTextSize(0.030)
        from plotting.calox_plot_helper import BoardPlotHelper as _BPH
        helper_cer = _BPH(pm2, xrange=(-10, 10), yrange=(-11, 11),
                          W_ref=1100, H_ref=1000)
        helper_cer.plot_board_map(
            cer_hists, "DRS_CFD_MPV_Scan_p0c2_Cer",
            extra_text="Cer  ", zmin=zmin_val, zmax=zmax_val,
            nTextDigits=2, zlabel="Time [ns]", lumitext=lumi, usePDF=True,
            extraToDraw=[beam_shaft, beam_top, beam_bot, beam_label])
        pm2.add_pdf_link("DRS_CFD_MPV_Scan_p0c2_Cer")

    pm2.generate_html(
        "DRS_CFD_MPV_Scan_p0_corrected2.html",
        plots_per_row=2,
        title="DRS CFD MPV Scan — Channel timing residuals [ns]",
        intro_text=(
            f"Per-channel timing residual (in ns) defined as "
            f"[(p0 − group_MPV) − {Y1500_EVAL}×|p1|] × 0.2 ns/TS, "
            f"where p0 is the per-channel CFD intercept, group_MPV is the "
            f"type-average MPV at {Y1500_EVAL} mm from "
            f"data/drs/drs_cfd_mpv_at_{Y1500_EVAL}mm_by_type.json, "
            f"and p1 is the fitted slope (TS/mm). "
            "Color range set to 10–90% quantile of active channels. "
            "Values shown to 0.1 ns precision."
        ),
    )


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


def _plot_energy_channels(channel_order, data, pm, xmin, xmax):
    """Draw one SumMean-vs-position TGraph per channel. Returns plots grouped by type."""
    plots_by_type = {}
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

        ys = [m for _, m in pts]
        ymid = (max(ys) + min(ys)) / 2
        yhalf = max(50.0, (max(ys) - min(ys)) / 2 * 1.4 + 20)

        pave = create_pave_text(0.18, 0.77, 0.60, 0.87)
        pave.AddText(f"B: {b}, G: {g}, C: {c}")
        pave.AddText(f"Tower: ({tx}, {ty})")

        plot_name = f"drs_summean_scan_{ch}_{var}"
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
            extra_text=var,
            extraToDraw=[pave],
        )
        plots_by_type.setdefault(var, []).append(plot_name)
    return plots_by_type


# ---------------------------------------------------------------------------
# Top-level scan runners
# ---------------------------------------------------------------------------

def run_time_scan(drsboards, run_x, channel_order, lumi, xmin, xmax, use_jsroot=False):
    mpv_data   = _collect_mpv_data(drsboards, run_x)
    shared_fit = _compute_shared_slopes(channel_order, mpv_data) if FIT_SHARED_SLOPE else {}

    pm = _make_pm("DRS_CFD_MPV_Scan", use_jsroot)
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    p0_map, p1_map, plots_by_type = _plot_time_channels(
        channel_order, mpv_data, shared_fit, pm, lumi, xmin, xmax)

    for type_name, plot_names in plots_by_type.items():
        pm.reset_plots()
        for name in plot_names:
            pm.add_plot(name)
        pm.generate_html(
            f"DRS_CFD_MPV_Scan_{type_name}.html",
            plots_per_row=4,
            title=f"DRS CFD MPV Scan — {type_name}",
        )

    if p0_map and p1_map:
        _plot_y_at_x_page(drsboards, p0_map, p1_map, channel_order, lumi, use_jsroot=use_jsroot)
        if SAVE_T0_JSON:
            _save_t0_json(p0_map, p1_map)
        _plot_fit_params_page(drsboards, p0_map, p1_map, channel_order, lumi, use_jsroot=use_jsroot)
        _plot_p0_corrected2_page(drsboards, p0_map, p1_map, channel_order, lumi, use_jsroot=use_jsroot)


def run_energy_scan(drsboards, run_x, channel_order, lumi, xmin, xmax, use_jsroot=False):
    data = _collect_summean(run_x)
    if not any(data.values()):
        print("No drs_energy_stats.json data found for any scan run.")
        return

    pm = _make_pm("DRS_SumMean_Scan", use_jsroot)
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    plots_by_type = _plot_energy_channels(channel_order, data, pm, xmin, xmax)

    run_list = ", ".join(str(r) for r in sorted(run_x))
    for type_name, plot_names in plots_by_type.items():
        pm.reset_plots()
        for name in plot_names:
            pm.add_plot(name)
        pm.generate_html(
            f"DRS_SumMean_Scan_{type_name}.html",
            plots_per_row=4,
            title=f"DRS SumMean Scan — {type_name}",
            intro_text=(
                "Per-channel mean DRS CFD energy (DRS_SumMean) vs table "
                f"position X [mm] (runs {run_list})."
            ),
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
    return parser.parse_args()


def main():
    args = _parse_args()
    setup_root(batch_mode=True, load_functions=(args.mode != "energy"))

    with open("data/Runlist.json") as f:
        runlist = json.load(f)

    def _parse_x(run):
        raw = runlist[str(run)]["Table position X"]
        return float(raw.replace("mm", "").strip())

    run_x = {r: _parse_x(r) for r in SCAN_RUNS if str(r) in runlist}
    if not run_x:
        print("No table positions found for the scan runs in data/Runlist.json.")
        return

    drsboards     = build_drs_boards(run_number=REFERENCE_RUN)
    channel_order = _build_channel_order(drsboards)

    xs_all = list(run_x.values())
    xmin, xmax = min(xs_all) - 20, max(xs_all) + 20

    _btypes = {"pion": "#pi^{+}", "pions": "#pi^{+}", "pi+": "#pi^{+}",
               "positron": "e^{+}", "positrons": "e^{+}", "e+": "e^{+}",
               "mu+": "#mu^{+}"}
    btype, benergy = getRunInfo(REFERENCE_RUN)
    lumi = f"{_btypes.get(btype.lower(), btype.lower())}, {benergy} GeV, 90^{{#circ}}"

    if args.mode in ("time", "both"):
        run_time_scan(drsboards, run_x, channel_order, lumi, xmin, xmax, use_jsroot=args.jsroot)
    if args.mode in ("energy", "both"):
        run_energy_scan(drsboards, run_x, channel_order, lumi, xmin, xmax, use_jsroot=args.jsroot)

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
