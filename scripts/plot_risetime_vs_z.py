#!/usr/bin/env python3
"""TB2026 leading-edge rise time versus scan (Z) position, per channel.

Same run set, event selection, run positions and robust core-MPV statistic as
scripts/plot_velocityfits.py --reference-mpv-logic, but the per-run quantity is
the leading-edge rise time (CFD 90% - CFD 20%, variables.drs.add_rise_time)
instead of the MCP-corrected TOA.

For every active calorimeter DRS channel:
  1. fill one rise-time histogram per run (MCP-clean events, valid rise only);
  2. take the core MPV (utils.utils.get_hist_mpv) as that run's rise time,
     with its standard error as the y uncertainty;
  3. plot rise time versus the run scan coordinate and fit a straight line.

Outputs (no PNG/PDF):
  risetime_vs_z.html        per-channel rise time vs position + family overlays
  risetime_hists.html       per-channel overlay of the per-run rise-time histograms
  risetime_vs_z.root        histograms and canvases
  risetime_vs_z_summary.json  per-channel/per-run points and linear fits
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ROOT

from core.analysis_manager import CaloXAnalysisManager
from utils.data_loader import getRunInfo
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import get_hist_mpv
from variables.drs import add_rise_time

from plot_velocityfits import (
    ROOT_COLORS,
    ChannelInfo,
    active_timing_channels,
    draw_header,
    family_key,
    family_label,
    load_run_positions,
    make_manager_args,
    safe_name,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Per-channel rise time (CFD 90% - CFD 20%) versus scan position."
    )
    parser.add_argument("--runs", nargs="+", type=int, required=True)
    parser.add_argument("--json-file", default="data/datafiles.json")
    parser.add_argument("--first-event", type=int, default=0)
    parser.add_argument("--last-event", type=int, default=-1)
    parser.add_argument("--particle", default=None)
    parser.add_argument("--no-mcp-clean", action="store_true")
    parser.add_argument("--channels", nargs="+", default=None, metavar="DRS_BRANCH")
    parser.add_argument("--skip-sci", action="store_true")
    parser.add_argument("--n-threads", type=int, default=10)
    parser.add_argument("--outdir", default="results/RiseTimevsZpos")

    parser.add_argument("--runlist", default="data/Runlist.json")
    parser.add_argument("--position-key", default="Table position X")
    parser.add_argument("--positions-json", default=None)
    parser.add_argument("--positions", nargs="+", default=None, metavar="RUN:MM")
    parser.add_argument("--position-label", default=None)

    # Rise-time definition; anchor/search/baseline match the PE-calibration
    # steps (plot_pe_vs_event_quantity.py) so the same peak is used.
    parser.add_argument("--rise-frac-lo", type=float, default=0.20)
    parser.add_argument("--rise-frac-hi", type=float, default=0.90)
    parser.add_argument("--anchor-offset", type=int, default=-330)
    parser.add_argument("--search-half", type=int, default=60)
    parser.add_argument("--baseline-pre-lo", type=int, default=80)
    parser.add_argument("--baseline-pre-hi", type=int, default=20)
    parser.add_argument("--min-amplitude", type=float, default=0.0)

    parser.add_argument("--rise-bins", type=int, default=600)
    parser.add_argument("--rise-range", nargs=2, type=float, default=(0.0, 6.0),
                        metavar=("MIN_NS", "MAX_NS"))
    parser.add_argument(
        "--mpv-window", type=float, default=0.3,
        help="Half-width [ns] of the core window in get_hist_mpv (default: 0.3; "
        "the TOA scans use 1.2 ns, too wide for a ~1 ns rise time).",
    )
    parser.add_argument("--fit-min-events", type=int, default=25)

    args = parser.parse_args()
    if args.rise_range[0] >= args.rise_range[1]:
        parser.error("--rise-range requires MIN_NS < MAX_NS")
    if args.position_label is None:
        args.position_label = args.position_key
    # load_run_positions() reads this flag from the velocity script's namespace.
    args.no_velocity_fits = False
    return args


def weighted_line_fit(xs: List[float], ys: List[float], es: List[float]):
    """Weighted straight line; returns (intercept, slope, int_err, slope_err, chi2/ndf)."""
    if len(xs) < 2 or len(set(xs)) < 2:
        return None
    w = [1.0 / max(e, 1.0e-6) ** 2 for e in es]
    sw = sum(w)
    swx = sum(wi * x for wi, x in zip(w, xs))
    swy = sum(wi * y for wi, y in zip(w, ys))
    swxx = sum(wi * x * x for wi, x in zip(w, xs))
    swxy = sum(wi * x * y for wi, x, y in zip(w, xs, ys))
    denom = sw * swxx - swx ** 2
    if not math.isfinite(denom) or denom <= 0.0:
        return None
    slope = (sw * swxy - swx * swy) / denom
    intercept = (swy - slope * swx) / sw
    chi2 = sum(wi * (y - intercept - slope * x) ** 2 for wi, x, y in zip(w, xs, ys))
    ndf = len(xs) - 2
    return (intercept, slope, math.sqrt(swxx / denom), math.sqrt(sw / denom),
            chi2 / ndf if ndf > 0 else None)


def make_graph(points: List[Tuple[float, float, float]], name: str):
    graph = ROOT.TGraphErrors(len(points))
    graph.SetName(name)
    graph.SetTitle("")
    for i, (x, y, e) in enumerate(points):
        graph.SetPoint(i, x, y)
        graph.SetPointError(i, 0.0, max(e, 1.0e-6))
    return graph


def style_axes(graph, xtitle: str) -> None:
    graph.GetXaxis().SetTitle(xtitle)
    graph.GetYaxis().SetTitle("Rise time, CFD 90% #minus CFD 20% [ns]")
    for axis in (graph.GetXaxis(), graph.GetYaxis()):
        axis.SetTitleSize(0.045)
        axis.SetLabelSize(0.04)
    graph.GetYaxis().SetTitleOffset(1.25)


def new_canvas(name: str):
    canvas = ROOT.TCanvas(name, "", 1100, 800)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.12)
    canvas.SetBottomMargin(0.12)
    canvas.SetTicks(1, 1)
    return canvas


def channel_canvas(branch, meta, points, fit, args, selection_label):
    xs = [p[0] for p in points]
    graph = make_graph(points, "gr_risetime_{}".format(safe_name(branch)))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.15)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)
    canvas = new_canvas("risetime_{}".format(safe_name(branch)))
    graph.Draw("AP")
    style_axes(graph, "{} [mm]".format(args.position_label))
    keep: List[Any] = [graph]
    if fit is not None:
        intercept, slope, _, slope_err, chi2_ndf = fit
        xpad = max(1.0, 0.08 * (max(xs) - min(xs)))
        func = ROOT.TF1("fit_risetime_{}".format(safe_name(branch)), "pol1",
                        min(xs) - xpad, max(xs) + xpad)
        func.SetParameters(intercept, slope)
        func.SetLineColor(ROOT.kRed + 1)
        func.SetLineWidth(3)
        func.SetLineStyle(2)
        func.Draw("SAME")
        info = ROOT.TPaveText(0.15, 0.72, 0.55, 0.83, "NDC")
        info.SetBorderSize(0)
        info.SetFillStyle(0)
        info.SetTextFont(42)
        info.SetTextSize(0.026)
        info.SetTextAlign(12)
        info.AddText("d(rise)/dz = ({:.5f} #pm {:.5f}) ns/mm".format(slope, slope_err))
        if chi2_ndf is not None:
            info.AddText("#chi^{{2}}/ndf = {:.2f}".format(chi2_ndf))
        info.Draw()
        keep += [func, info]
    draw_header(canvas, "Rise time vs position: {}".format(meta.short_label),
                selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keep
    return canvas


def family_canvas(family, branches, channel_points, channel_info, args, selection_label):
    canvas = new_canvas("risetime_family_{}".format(safe_name(family)))
    all_x = [p[0] for b in branches for p in channel_points[b]]
    all_y = [p[1] for b in branches for p in channel_points[b]]
    ymax = max(all_y)
    ymin = min(all_y)
    pad = max(0.05, 0.15 * (ymax - ymin))
    frame = ROOT.TH1F("frame_{}".format(safe_name(family)), "", 10,
                      min(all_x) - 5.0, max(all_x) + 5.0)
    frame.SetMinimum(max(0.0, ymin - pad))
    frame.SetMaximum(ymax + 2.0 * pad)
    frame.SetStats(0)
    frame.GetXaxis().SetTitle("{} [mm]".format(args.position_label))
    frame.GetYaxis().SetTitle("Rise time, CFD 90% #minus CFD 20% [ns]")
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetTitleSize(0.045)
        axis.SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw("AXIS")
    keep: List[Any] = [frame]
    legend = ROOT.TLegend(0.14, 0.64, 0.94, 0.86)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.016)
    keep.append(legend)
    for i, branch in enumerate(branches):
        graph = make_graph(channel_points[branch], "gr_fam_{}".format(safe_name(branch)))
        color = ROOT_COLORS[i % len(ROOT_COLORS)]
        graph.SetMarkerStyle(20 + (i // len(ROOT_COLORS)) % 5)
        graph.SetMarkerSize(0.9)
        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.Draw("PL SAME")
        legend.AddEntry(graph, channel_info[branch].short_label, "pl")
        keep.append(graph)
    legend.Draw()
    draw_header(canvas, "Rise time vs position: {}".format(family_label(family)),
                selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keep
    return canvas


def hist_overlay_canvas(branch, meta, hists, args, selection_label):
    entries = [(run, h) for run, h in sorted(hists.items()) if h.GetEntries() > 0]
    if not entries:
        return None
    canvas = new_canvas("risehist_{}".format(safe_name(branch)))
    keep: List[Any] = []
    ymax = 0.0
    scaled = []
    for run, h in entries:
        c = h.Clone("{}_norm".format(h.GetName()))
        c.SetDirectory(0)
        if c.Integral() > 0:
            c.Scale(1.0 / c.Integral())
        scaled.append((run, c))
        ymax = max(ymax, c.GetMaximum())
    legend = ROOT.TLegend(0.60, 0.40, 0.94, 0.86)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.016)
    keep.append(legend)
    for i, (run, c) in enumerate(scaled):
        color = ROOT_COLORS[i % len(ROOT_COLORS)]
        c.SetLineColor(color)
        c.SetLineWidth(2)
        c.SetStats(0)
        c.GetXaxis().SetTitle("Rise time, CFD 90% #minus CFD 20% [ns]")
        c.GetYaxis().SetTitle("Normalized events")
        c.GetXaxis().SetTitleSize(0.045)
        c.GetYaxis().SetTitleSize(0.045)
        c.GetYaxis().SetTitleOffset(1.25)
        c.SetMinimum(0.0)
        c.SetMaximum(1.3 * ymax)
        c.Draw("HIST" if i == 0 else "HIST SAME")
        legend.AddEntry(c, "Run {}".format(run), "l")
        keep.append(c)
    legend.Draw()
    draw_header(canvas, "Rise-time distributions: {}".format(meta.short_label),
                selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keep
    return canvas


def to_json(canvas) -> str:
    return str(ROOT.TBufferJSON.ToJSON(canvas).Data())


def main() -> None:
    auto_timer("Total Execution Time")
    args = parse_args()
    setup_root(n_threads=args.n_threads, batch_mode=True, load_functions=True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    positions = load_run_positions(args)
    if len(set(positions.values())) < 2:
        raise SystemExit("[FATAL] Need at least two distinct run positions "
                         "(check --runlist / --position-key).")

    root_path = outdir / "risetime_vs_z.root"
    root_file = ROOT.TFile(str(root_path), "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Cannot create {}".format(root_path))

    selection_label = "MCP-clean events" if not args.no_mcp_clean else "All events"
    print("\n" + "=" * 78)
    print("CaloX TB2026 rise time (CFD {:.0f}%-{:.0f}%) vs {}".format(
        100 * args.rise_frac_lo, 100 * args.rise_frac_hi, args.position_key))
    print("Runs   : {}".format(args.runs))
    print("Output : {}".format(outdir))
    print("=" * 78)

    requested = set(args.channels) if args.channels else None
    channel_info: Dict[str, ChannelInfo] = {}
    rise_hists: Dict[str, Dict[int, Any]] = {}
    run_stats: Dict[str, Dict[int, Dict[str, Any]]] = {}

    for run in args.runs:
        beam_type, beam_energy = getRunInfo(run)
        print("\n[RUN {}] beam={}, energy={} GeV".format(run, beam_type, beam_energy))
        manager = CaloXAnalysisManager(make_manager_args(args, run)).prepare(do_fers=False)
        manager.apply_mcp_selection(flag_only=True)
        if not args.no_mcp_clean:
            manager.apply_beam_pid_selection(flag_only=False, particle="mcp_clean")
        if args.particle:
            manager.apply_beam_pid_selection(flag_only=False, particle=args.particle)

        info = active_timing_channels(manager, requested)
        if args.skip_sci:
            info = {b: m for b, m in info.items() if m.material != "Sci"}
        if not info:
            print("[WARN] Run {}: no channels; skipped.".format(run))
            continue

        rdf = add_rise_time(
            manager.rdf, sorted(info), anchor_offset=args.anchor_offset,
            search_half=args.search_half, frac_lo=args.rise_frac_lo,
            frac_hi=args.rise_frac_hi, baseline_pre_lo=args.baseline_pre_lo,
            baseline_pre_hi=args.baseline_pre_hi, min_amplitude=args.min_amplitude,
        )
        columns = {str(c) for c in rdf.GetColumnNames()}

        results: Dict[str, Any] = {}
        for branch in sorted(info):
            value_col = "{}_risetime_adaptive".format(branch)
            valid_col = "{}_risetime_valid_adaptive".format(branch)
            if value_col not in columns or valid_col not in columns:
                continue
            model = ROOT.RDF.TH1DModel(
                "rise_run{}_{}".format(run, safe_name(branch)),
                "Run {};Rise time [ns];Events".format(run),
                args.rise_bins, args.rise_range[0], args.rise_range[1],
            )
            results[branch] = rdf.Filter(
                "{} && std::isfinite({})".format(valid_col, value_col)
            ).Histo1D(model, value_col)
        if not results:
            print("[WARN] Run {}: no rise-time columns; skipped.".format(run))
            continue
        ROOT.RDF.RunGraphs(list(results.values()))

        run_dir = root_file.mkdir("run{}".format(run))
        run_dir.cd()
        for branch, result in results.items():
            channel_info.setdefault(branch, info[branch])
            hist = result.GetValue().Clone(safe_name(branch))
            hist.SetDirectory(0)
            hist.SetTitle("Run {}: {};Rise time [ns];Events".format(run, info[branch].label))
            hist.Write()
            rise_hists.setdefault(branch, {})[run] = hist

            n_entries = int(round(hist.Integral(1, hist.GetNbinsX())))
            stats: Dict[str, Any] = {"n_entries": n_entries, "mpv_ns": None,
                                     "mpv_err_ns": None, "mean_ns": None, "ok": False}
            if n_entries >= args.fit_min_events:
                mpv, mpv_err = get_hist_mpv(hist, window_ts=args.mpv_window)
                if (math.isfinite(mpv) and math.isfinite(mpv_err)
                        and 0.0 < mpv_err < 999.0):
                    stats.update(mpv_ns=float(mpv), mpv_err_ns=float(mpv_err),
                                 mean_ns=float(hist.GetMean()), ok=True)
            run_stats.setdefault(branch, {})[run] = stats
            if stats["ok"]:
                print("[RISE] Run {} {}: core MPV={:.4f} +/- {:.4f} ns, N={}".format(
                    run, info[branch].short_label, stats["mpv_ns"],
                    stats["mpv_err_ns"], n_entries))
            else:
                print("[WARN] Run {} {}: rise-time MPV unavailable (N={}).".format(
                    run, info[branch].short_label, n_entries))

    if not rise_hists:
        root_file.Close()
        raise RuntimeError("No rise-time histograms produced.")

    ordered = sorted(rise_hists)
    canvas_dir = root_file.mkdir("canvases")
    z_keys: List[str] = []
    z_jsons: Dict[str, str] = {}
    h_keys: List[str] = []
    h_jsons: Dict[str, str] = {}
    channel_points: Dict[str, List[Tuple[float, float, float]]] = {}
    fits: Dict[str, Any] = {}

    for branch in ordered:
        pts = sorted(
            (positions[run], s["mpv_ns"], s["mpv_err_ns"])
            for run, s in run_stats.get(branch, {}).items()
            if s["ok"] and run in positions
        )
        channel_points[branch] = pts
        if len(pts) < 2:
            continue
        fit = weighted_line_fit([p[0] for p in pts], [p[1] for p in pts],
                                [p[2] for p in pts])
        fits[branch] = fit
        canvas = channel_canvas(branch, channel_info[branch], pts, fit, args,
                                selection_label)
        name = "risetime_{}".format(safe_name(branch))
        canvas_dir.cd()
        canvas.Write(name)
        z_keys.append(name)
        z_jsons[name] = to_json(canvas)

    families = sorted({family_key(channel_info[b], "material_size")
                       for b in ordered if len(channel_points.get(b, [])) >= 2})
    family_keys: List[str] = []
    for family in families:
        members = [b for b in ordered
                   if len(channel_points.get(b, [])) >= 2
                   and family_key(channel_info[b], "material_size") == family]
        canvas = family_canvas(family, members, channel_points, channel_info, args,
                               selection_label)
        name = "risetime_family_{}".format(safe_name(family))
        canvas_dir.cd()
        canvas.Write(name)
        family_keys.append(name)
        z_jsons[name] = to_json(canvas)
    z_keys = family_keys + ["NEWLINE"] + z_keys if family_keys else z_keys

    for branch in ordered:
        canvas = hist_overlay_canvas(branch, channel_info[branch], rise_hists[branch],
                                     args, selection_label)
        if canvas is None:
            continue
        name = "risehist_{}".format(safe_name(branch))
        canvas_dir.cd()
        canvas.Write(name)
        h_keys.append(name)
        h_jsons[name] = to_json(canvas)
    root_file.Close()

    position_text = ", ".join("run {} = {:g} mm".format(r, positions[r])
                              for r in args.runs if r in positions)
    common = [
        "**Runs and positions:** {}".format(position_text),
        "**Selection:** {}{}".format(selection_label,
                                     " + {}".format(args.particle) if args.particle else ""),
        "**Rise time:** CFD {:.0f}% to CFD {:.0f}% of the leading edge "
        "(`variables.drs.add_rise_time`, anchor {} TS, search +/-{} TS), "
        "valid pulses only.".format(100 * args.rise_frac_lo, 100 * args.rise_frac_hi,
                                    args.anchor_offset, args.search_half),
    ]
    if z_keys:
        generate_jsroot_html(
            z_keys, z_jsons, plots_per_row=2,
            output_html=str(outdir / "risetime_vs_z.html"),
            title="CaloX TB2026 rise time vs position",
            intro_text="\n".join(common + [
                "* Y value per run: robust core MPV of the rise-time histogram "
                "(+/-{:g} ns window), error = standard error in the window.".format(
                    args.mpv_window),
                "* First row: all channels of a (material, fiber size) family overlaid; "
                "then one plot per channel with a weighted linear fit.",
            ]),
        )
    if h_keys:
        generate_jsroot_html(
            h_keys, h_jsons, plots_per_row=2,
            output_html=str(outdir / "risetime_hists.html"),
            title="CaloX TB2026 rise-time distributions",
            intro_text="\n".join(common + [
                "* Unit-area rise-time histograms, one curve per run.",
            ]),
        )

    summary = {
        "runs": args.runs,
        "position_key": args.position_key,
        "run_positions_mm": {str(r): positions.get(r) for r in args.runs},
        "rise_fractions": [args.rise_frac_lo, args.rise_frac_hi],
        "mpv_window_ns": args.mpv_window,
        "channels": {
            b: {
                "label": channel_info[b].label,
                "short_label": channel_info[b].short_label,
                "material": channel_info[b].material,
                "size": channel_info[b].size,
                "per_run": {str(r): s for r, s in sorted(run_stats.get(b, {}).items())},
                "linear_fit": (
                    {"intercept_ns": fits[b][0], "slope_ns_per_mm": fits[b][1],
                     "intercept_err_ns": fits[b][2], "slope_err_ns_per_mm": fits[b][3],
                     "chi2_ndf": fits[b][4]}
                    if fits.get(b) is not None else None
                ),
            }
            for b in ordered
        },
    }
    (outdir / "risetime_vs_z_summary.json").write_text(json.dumps(summary, indent=2),
                                                       encoding="utf-8")
    print("\n[DONE] {}".format(outdir))


if __name__ == "__main__":
    main()
