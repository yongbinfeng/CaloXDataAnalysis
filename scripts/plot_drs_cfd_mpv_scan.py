"""
DRS CFD MPV position scan.

For each DRS channel, plots the MPV of the fine-binned TS_cfd_mcp histogram
versus table position X across a fixed set of scan runs, one graph per channel.
A linear fit is overlaid and its parameters displayed on the plot.

Usage
-----
  python scripts/plot_drs_cfd_mpv_scan.py
"""

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

auto_timer("Total Execution Time")
setup_root(batch_mode=True, load_functions=True)

SCAN_RUNS = [1501, 1507, 1511, 1513, 15130]
REFERENCE_RUN = 1501   # used only for DRS board map

PLOTDIR = "results/plots/PositionScan"
HTMLDIR = "results/html/PositionScan"

# True  → one shared slope per category (CerQuartz/CerPlastic/Sci);
#          each channel keeps its own intercept.
# False → independent linear fit per channel (original behaviour).
FIT_SHARED_SLOPE = True


def main():
    with open("data/Runlist.json") as f:
        runlist = json.load(f)

    def _parse_x(run):
        raw = runlist[str(run)]["Table position X"]
        return float(raw.replace("mm", "").strip())

    run_x = {r: _parse_x(r) for r in SCAN_RUNS if str(r) in runlist}
    drsboards = build_drs_boards(run=REFERENCE_RUN)

    # Collect (x, mpv) per channel, sorted by table position
    mpv_data = {}   # ch -> [(x, mpv), ...]
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

    xs_all = list(run_x.values())
    xmin = min(xs_all) - 20
    xmax = max(xs_all) + 20

    _btypes = {"pion": "#pi^{+}", "pions": "#pi^{+}", "pi+": "#pi^{+}",
               "positron": "e^{+}", "positrons": "e^{+}", "e+": "e^{+}",
               "mu+": "#mu^{+}"}
    btype, benergy = getRunInfo(REFERENCE_RUN)
    lumi = f"{_btypes.get(btype.lower(), btype.lower())}, {benergy} GeV, 90^{{#circ}}"

    # Build ordered list of channels with B/G/C info and type label
    channel_order = []
    for _, board in drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            var = get_channel_var(chan)           # "CerQuartz", "CerPlastic", "Sci"
            channel_order.append((
                chan.get_channel_name(blsub=False),
                board.board_no, chan.group_no, chan.channel_no,
                var,
            ))

    pm = PlotManager(
        rootdir=f"results/root/Run{REFERENCE_RUN}",
        plotdir=PLOTDIR,
        htmldir=HTMLDIR,
    )
    pm.set_output_dir("DRS_CFD_MPV_Scan")
    os.makedirs(pm.get_output_dir(), exist_ok=True)

    # Pre-compute one shared slope per category when flag is set.
    # Model: y_{c,i} = p0_c + p1_cat * x_{c,i}  (within-group OLS)
    shared_fit = {}   # var -> (p1, e1, SXX_total, s2)
    if FIT_SHARED_SLOPE:
        for var_name in dict.fromkeys(v for _, _, _, _, v in channel_order):
            chs_var = [(ch, sorted(mpv_data.get(ch, [])))
                       for ch, _, _, _, v in channel_order if v == var_name]
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
            p1_shared = SXY / SXX
            shared_fit[var_name] = (p1_shared, 1.0 / SXX ** 0.5, SXX)

    plots_by_type = {}   # var -> [plot_name, ...]
    p0_map, p1_map = {}, {}

    for ch, b, g, c, var in channel_order:
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

        # Weighted least squares (w_i = 1/sigma_y_i^2)
        n = len(pts)
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
        fit = ROOT.TF1(f"fit_{ch}", "pol1", xmin, xmax)
        fit.SetParameters(p0, p1)
        fit.SetParError(0, e0)
        fit.SetParError(1, e1)
        fit.SetLineColor(ROOT.kRed + 1)
        fit.SetLineWidth(2)
        fit.SetLineStyle(2)   # dashed
        gr.GetListOfFunctions().Add(fit)

        ys = [y for _, y, _ in pts]
        ymid = (max(ys) + min(ys)) / 2
        yhalf = max(30.0, (max(ys) - min(ys)) / 2 * 1.5 + 10)

        pave_bgc = create_pave_text(0.15, 0.80, 0.60, 0.88)
        pave_bgc.AddText(f"B: {b}, G: {g}, C: {c}")

        pave_fit = create_pave_text(0.15, 0.13, 0.65, 0.28)
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

    for type_name, plot_names in plots_by_type.items():
        pm.reset_plots()
        for name in plot_names:
            pm.add_plot(name)
        pm.generate_html(
            f"DRS_CFD_MPV_Scan_{type_name}.html",
            plots_per_row=4,
            title=f"DRS CFD MPV Scan — {type_name}",
        )
    # --- Board maps: p0-400, |p1|×100, and 1/|p1| in cm/ns ---
    if p0_map and p1_map:
        p0_shifted = {ch: v - 400
                      for ch, v in p0_map.items()}
        # ×100 to move decimal point; unit becomes (TS/mm)×100
        p1_abs_map = {ch: abs(v) * 100
                      for ch, v in p1_map.items()}
        # 1/|p1| [mm/TS] → cm/ns: ×(0.1 cm/mm)/(0.2 ns/TS) = ×0.5
        p1_inv_map = {ch: 0.5 / abs(v)
                      for ch, v in p1_map.items() if v != 0}

        pm_map = PlotManager(
            rootdir=f"results/root/Run{REFERENCE_RUN}",
            plotdir=PLOTDIR,
            htmldir=HTMLDIR,
        )
        pm_map.set_output_dir("DRS_CFD_MPV_Scan")
        helper = BoardPlotHelper(pm_map)

        cer_hists, sci_hists = visualizeDRSBoards(
            drsboards, valuemaps=p0_shifted, suffix="MPV_scan_p0")
        helper.plot_cer_sci_pair(
            cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p0",
            zmin=None, zmax=None,
            nTextDigits=1, lumitext=lumi, usePDF=False)

        cer_hists, sci_hists = visualizeDRSBoards(
            drsboards, valuemaps=p1_abs_map, suffix="MPV_scan_p1_abs")
        helper.plot_cer_sci_pair(
            cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p1_abs",
            zmin=None, zmax=None,
            nTextDigits=2, lumitext=lumi, usePDF=False)

        p1i_sorted = sorted(p1_inv_map.values())
        n_inv = len(p1i_sorted)
        p5  = p1i_sorted[max(0, int(0.10 * n_inv))]
        p95 = p1i_sorted[min(n_inv - 1, int(0.90 * n_inv))]

        cer_hists, sci_hists = visualizeDRSBoards(
            drsboards, valuemaps=p1_inv_map, suffix="MPV_scan_p1_inv")
        helper.plot_cer_sci_pair(
            cer_hists, sci_hists, "DRS_CFD_MPV_Scan_p1_inv",
            zmin=p5, zmax=p95,
            nTextDigits=1, lumitext=lumi, usePDF=False)

        # Speed distribution 1D histograms by channel type
        speed_by_type = {}
        for ch, _b, _g, _c, var in channel_order:
            if ch in p1_inv_map:
                speed_by_type.setdefault(var, []).append(p1_inv_map[ch])

        xlo_h = 12.0
        xhi_h = 22.0
        type_order = ["CerQuartz", "CerPlastic", "Sci"]
        speed_colors = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2]
        hists_speed, speed_labels, speed_cols_used = [], [], []
        for tname, col in zip(type_order, speed_colors):
            vals = speed_by_type.get(tname, [])
            if not vals:
                continue
            h = ROOT.TH1F(f"h_speed_{tname}", "", 20, xlo_h, xhi_h)
            for v in vals:
                h.Fill(v)
            hists_speed.append(h)
            speed_labels.append(tname)
            speed_cols_used.append(col)

        if hists_speed:
            speed_plot_name = "drs_cfd_speed_distribution"
            DrawHistos(
                hists_speed, speed_labels,
                xlo_h, xhi_h, "Speed [cm/ns]",
                None, None, "Channels",
                speed_plot_name,
                outdir=pm_map.get_output_dir(),
                dology=False,
                usePDF=False,
                mycolors=speed_cols_used,
                drawoptions=["C"] * len(hists_speed),
                lumitext=lumi,
            )
            pm_map.add_plot(speed_plot_name)

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
                "* **p0 − 400 [TS]**: intercept shifted by 400 for readability. "
                "Represents the expected MPV offset at X = 0 mm.\n"
                "* **|p1| × 100 [TS/mm × 100]**: absolute slope scaled by 100 "
                "for readability — magnitude of the timing change per mm.\n"
                "* **1/|p1| [cm/ns]**: inverted absolute slope converted from "
                "mm/TS to cm/ns (1 TS = 200 ps). Interpretable as an effective "
                "speed of light in the medium.\n\n"
                "Channels with no valid data across all runs are not filled."
            ),
        )

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
