#!/usr/bin/env python3
"""Simultaneous (shared-slope) rise-time vs position fits per fibre family.

Post-processes results/RiseTimevsZpos/risetime_vs_z_summary.json (written by
scripts/plot_risetime_vs_z.py) - no event loop is re-run.

Fit model, same as scripts/plot_velocityfits.py --reference-mpv-logic
--slope-family material_size: within each (material, fibre size) family one
common slope d(rise)/dz is fitted to all channels simultaneously, while every
channel keeps its own intercept (its own absolute rise time).

Families: CerQuartz 3 mm, CerPlastic 3 mm, CerQuartz 6 mm, CerPlastic 6 mm.

Outputs in <outdir> (default results/RiseTimevsZpos/shared_fit):
  risetime_shared_fit_families.html  one raw + one offset-removed plot per family
  risetime_shared_fit_channels.html  per-channel plots with the shared-slope line
  risetime_shared_fit.root           canvases
  risetime_shared_fit_summary.json   shared slopes, per-channel intercepts, chi2
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import ROOT

from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root

from plot_risetime_vs_z import make_graph, new_canvas, style_axes, to_json
from plot_velocityfits import (
    ChannelInfo,
    VelocityPoint,
    compute_shared_family_slopes,
    draw_header,
    family_key,
    family_label,
    safe_name,
)

YTITLE = "Rise time, CFD 90% #minus CFD 20% [ns]"
FAMILY_ORDER = ["CerQuartz 3mm", "CerPlastic 3mm", "CerQuartz 6mm", "CerPlastic 6mm"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--summary",
                        default="results/RiseTimevsZpos/risetime_vs_z_summary.json")
    parser.add_argument("--outdir", default="results/RiseTimevsZpos/shared_fit")
    parser.add_argument("--position-label", default="Table position X")
    parser.add_argument("--selection-label", default="MCP-clean events")
    parser.add_argument("--min-points", type=int, default=3,
                        help="Channels with fewer valid runs are left out (default: 3).")
    return parser.parse_args()


def load(summary_path: str):
    d = json.load(open(summary_path, encoding="utf-8"))
    positions = {int(r): p for r, p in d["run_positions_mm"].items() if p is not None}
    info: Dict[str, ChannelInfo] = {}
    points: Dict[str, List[VelocityPoint]] = {}
    for branch, ch in d["channels"].items():
        info[branch] = ChannelInfo(branch, ch["label"], ch["short_label"],
                                   ch["material"], ch["size"])
        pts = [
            VelocityPoint(int(r), positions[int(r)], s["mpv_ns"], s["mpv_err_ns"], "CoreMPV")
            for r, s in ch["per_run"].items()
            if s["ok"] and int(r) in positions
        ]
        points[branch] = sorted(pts, key=lambda p: p.position_mm)
    return d, info, points


def fmt_slope(slope: float, err: float) -> str:
    """Slopes are O(1e-6) ns/mm, so quote them in units of 1e-6 ns/mm."""
    return "({:+.2f} #pm {:.2f}) #times 10^{{-6}} ns/mm".format(slope * 1e6, err * 1e6)


def intercept_for(points: List[VelocityPoint], slope: float) -> Tuple[float, float, float]:
    """Channel intercept with the family slope fixed: (b, chi2, weight_sum)."""
    w = [1.0 / max(p.point_error_ns, 1.0e-6) ** 2 for p in points]
    sw = sum(w)
    xbar = sum(wi * p.position_mm for wi, p in zip(w, points)) / sw
    ybar = sum(wi * p.mean_ns for wi, p in zip(w, points)) / sw
    b = ybar - slope * xbar
    chi2 = sum(wi * (p.mean_ns - b - slope * p.position_mm) ** 2 for wi, p in zip(w, points))
    return b, chi2, sw


def line(name: str, x0: float, x1: float, b: float, slope: float, color, width, style):
    f = ROOT.TF1(name, "pol1", x0, x1)
    f.SetParameters(b, slope)
    f.SetLineColor(color)
    f.SetLineWidth(width)
    f.SetLineStyle(style)
    return f


def family_canvas(family, branches, points, shared, intercepts, collapsed, args):
    slope, slope_err, _ = shared
    tag = "collapsed" if collapsed else "raw"
    canvas = new_canvas("risetime_shared_{}_{}".format(safe_name(family), tag))
    mean_b = sum(intercepts[b][0] for b in branches) / len(branches)

    def pts_of(b):
        shift = (mean_b - intercepts[b][0]) if collapsed else 0.0
        return [(p.position_mm, p.mean_ns + shift, p.point_error_ns) for p in points[b]]

    all_x = [x for b in branches for x, _, _ in pts_of(b)]
    all_y = [y for b in branches for _, y, _ in pts_of(b)]
    x0, x1 = min(all_x) - 0.05 * (max(all_x) - min(all_x)), max(all_x) + 0.05 * (max(all_x) - min(all_x))
    ylo, yhi = min(all_y), max(all_y)
    pad = max(0.03, 0.12 * (yhi - ylo))
    frame = ROOT.TH1F("frame_{}_{}".format(safe_name(family), tag), "", 10, x0, x1)
    frame.SetMinimum(ylo - pad)
    frame.SetMaximum(yhi + 2.2 * pad)
    frame.SetStats(0)
    frame.GetXaxis().SetTitle("{} [mm]".format(args.position_label))
    frame.GetYaxis().SetTitle(YTITLE if not collapsed else YTITLE + " (offsets removed)")
    for axis in (frame.GetXaxis(), frame.GetYaxis()):
        axis.SetTitleSize(0.045)
        axis.SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(1.25)
    frame.Draw("AXIS")
    keep: List[Any] = [frame]

    for b in branches:
        graph = make_graph(pts_of(b), "gr_shared_{}_{}".format(safe_name(b), tag))
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.55)
        graph.SetMarkerColorAlpha(ROOT.kBlue + 1, 0.55)
        graph.SetLineColorAlpha(ROOT.kBlue + 1, 0.35)
        graph.Draw("P SAME")
        keep.append(graph)
        if not collapsed:
            f = line("ln_{}_{}".format(safe_name(b), tag), x0, x1,
                     intercepts[b][0], slope, ROOT.kAzure + 1, 1, 1)
            f.SetLineColorAlpha(ROOT.kAzure + 1, 0.4)
            f.Draw("SAME")
            keep.append(f)

    shared_line = line("shared_{}_{}".format(safe_name(family), tag), x0, x1,
                       mean_b, slope, ROOT.kRed + 1, 4, 1)
    shared_line.Draw("SAME")
    keep.append(shared_line)

    n_pts = sum(len(points[b]) for b in branches)
    chi2 = sum(intercepts[b][1] for b in branches)
    ndf = n_pts - len(branches) - 1
    info = ROOT.TPaveText(0.15, 0.69, 0.60, 0.83, "NDC")
    info.SetBorderSize(0)
    info.SetFillStyle(0)
    info.SetTextFont(42)
    info.SetTextSize(0.027)
    info.SetTextAlign(12)
    info.AddText("{}: {} channels, {} points".format(family_label(family), len(branches), n_pts))
    info.AddText("shared d(rise)/dz = " + fmt_slope(slope, slope_err))
    if ndf > 0:
        info.AddText("#chi^{{2}}/ndf = {:.1f}/{} = {:.2f}".format(chi2, ndf, chi2 / ndf))
    info.Draw()
    keep.append(info)

    draw_header(canvas, "Rise time vs position, simultaneous fit: {}{}".format(
        family_label(family), " (channel offsets removed)" if collapsed else ""),
        args.selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = keep
    return canvas, chi2, ndf


def channel_canvas(branch, meta, pts, slope, slope_err, b, chi2, args):
    graph = make_graph([(p.position_mm, p.mean_ns, p.point_error_ns) for p in pts],
                       "gr_ch_{}".format(safe_name(branch)))
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.15)
    graph.SetMarkerColor(ROOT.kBlue + 1)
    graph.SetLineColor(ROOT.kBlue + 1)
    graph.SetLineWidth(2)
    canvas = new_canvas("risetime_shared_ch_{}".format(safe_name(branch)))
    graph.Draw("AP")
    style_axes(graph, "{} [mm]".format(args.position_label))
    xs = [p.position_mm for p in pts]
    pad = max(1.0, 0.08 * (max(xs) - min(xs)))
    f = line("fit_ch_{}".format(safe_name(branch)), min(xs) - pad, max(xs) + pad,
             b, slope, ROOT.kRed + 1, 3, 2)
    f.Draw("SAME")
    info = ROOT.TPaveText(0.15, 0.72, 0.58, 0.83, "NDC")
    info.SetBorderSize(0)
    info.SetFillStyle(0)
    info.SetTextFont(42)
    info.SetTextSize(0.026)
    info.SetTextAlign(12)
    info.AddText("shared slope: {}, {}".format(meta.material, meta.size))
    info.AddText("d(rise)/dz = " + fmt_slope(slope, slope_err))
    if len(pts) > 2:
        info.AddText("channel #chi^{{2}}/{} = {:.2f}".format(len(pts) - 1, chi2 / (len(pts) - 1)))
    info.Draw()
    draw_header(canvas, "Rise time vs position: {}".format(meta.short_label),
                args.selection_label)
    canvas.RedrawAxis()
    canvas.Update()
    canvas._keepalive = [graph, f, info]
    return canvas


def main() -> None:
    args = parse_args()
    setup_root(n_threads=1, batch_mode=True, load_functions=False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d, info, points = load(args.summary)
    points = {b: p for b, p in points.items() if len(p) >= args.min_points}
    shared = compute_shared_family_slopes(points, info, "material_size")

    root_file = ROOT.TFile(str(outdir / "risetime_shared_fit.root"), "RECREATE")
    canvas_dir = root_file.mkdir("canvases")
    fam_keys: List[str] = []
    ch_keys: List[str] = []
    jsons: Dict[str, str] = {}
    summary: Dict[str, Any] = {}

    families = [f for f in FAMILY_ORDER if f in shared] + \
               sorted(f for f in shared if f not in FAMILY_ORDER)
    for family in families:
        branches = sorted(b for b in points if family_key(info[b], "material_size") == family)
        slope, slope_err, _ = shared[family]
        intercepts = {b: intercept_for(points[b], slope) for b in branches}
        for collapsed in (False, True):
            canvas, chi2, ndf = family_canvas(family, branches, points, shared[family],
                                              intercepts, collapsed, args)
            name = canvas.GetName()
            canvas_dir.cd()
            canvas.Write(name)
            fam_keys.append(name)
            jsons[name] = to_json(canvas)
        for b in branches:
            bi, chi2_b, _ = intercepts[b]
            canvas = channel_canvas(b, info[b], points[b], slope, slope_err, bi, chi2_b, args)
            name = canvas.GetName()
            canvas_dir.cd()
            canvas.Write(name)
            ch_keys.append(name)
            jsons[name] = to_json(canvas)
        summary[family] = {
            "label": family_label(family),
            "slope_ns_per_mm": slope,
            "slope_err_ns_per_mm": slope_err,
            "n_channels": len(branches),
            "n_points": sum(len(points[b]) for b in branches),
            "chi2": chi2,
            "ndf": ndf,
            "chi2_ndf": chi2 / ndf if ndf > 0 else None,
            "channels": {
                b: {"short_label": info[b].short_label, "intercept_ns": intercepts[b][0],
                    "chi2": intercepts[b][1], "n_points": len(points[b])}
                for b in branches
            },
        }
        print("[SHARED] {:<16s} slope = {:+.2f} +/- {:.2f} e-6 ns/mm  ({} ch, chi2/ndf = {:.2f})".format(
            family_label(family), slope * 1e6, slope_err * 1e6, len(branches),
            chi2 / ndf if ndf > 0 else float("nan")))
    root_file.Close()

    positions = {int(r): p for r, p in d["run_positions_mm"].items() if p is not None}
    common = [
        "**Runs:** {}-{} ({} positions, {:g} to {:g} mm)".format(
            min(positions), max(positions), len(positions),
            min(positions.values()), max(positions.values())),
        "**Rise time:** CFD 90% - CFD 20%; per-run value = core MPV (+/-{:g} ns window).".format(
            d.get("mpv_window_ns", 0.3)),
        "**Fit:** one common slope per (material, fibre size) family fitted to all its "
        "channels simultaneously; each channel keeps its own intercept. "
        "Channels with fewer than {} valid runs are excluded.".format(args.min_points),
    ]
    generate_jsroot_html(
        fam_keys, {k: jsons[k] for k in fam_keys}, plots_per_row=2,
        output_html=str(outdir / "risetime_shared_fit_families.html"),
        title="CaloX TB2026 rise time vs position, simultaneous fits",
        intro_text="\n".join(common + [
            "* Per family: left = raw points with each channel's line (shared slope, own "
            "intercept) and the red family line at the mean intercept; right = same with "
            "every channel's offset removed so all points collapse onto the red line.",
        ]),
    )
    generate_jsroot_html(
        ch_keys, {k: jsons[k] for k in ch_keys}, plots_per_row=2,
        output_html=str(outdir / "risetime_shared_fit_channels.html"),
        title="CaloX TB2026 rise time vs position, per channel (shared slope)",
        intro_text="\n".join(common + ["* Red dashed line: the family slope with this "
                                       "channel's fitted intercept."]),
    )
    (outdir / "risetime_shared_fit_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print("[DONE] {}".format(outdir))


if __name__ == "__main__":
    main()
