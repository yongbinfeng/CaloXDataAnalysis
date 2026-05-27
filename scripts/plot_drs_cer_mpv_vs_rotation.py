"""
DRS Cherenkov corrected CFD MPV vs rotation angle.

For each iX slice in IX_TARGETS, collects all DRS Cherenkov channels at that
i_tower_x (optionally filtered by iY), reads the raw MPV from the
TS_cfd_mcp_finebins histogram, applies the type-average MPV correction,
combines all channels per run (weighted average, uncorrelated), and draws one
TGraphErrors per iX slice on a single canvas vs Table Degree.

HTML layout: combined plots (all iY) first, then one pair per iY slice.
All plots for the same particle type share a common y-axis range.

Usage
-----
  cd /path/to/CaloXDataAnalysis
  python scripts/plot_drs_cer_mpv_vs_rotation.py
"""

import os
import json
import ROOT

from channels.channel_map import build_drs_boards
from plotting.my_function import DrawHistos
from utils.html_generator import generate_jsroot_html
from utils.root_setup import setup_root
from utils.utils import get_hist_mpv
from variables.drs import subtract_type_mpv

setup_root(batch_mode=True, load_functions=True)

RUNS_EPLUS  = [1485, 1486, 1489, 1490, 1493, 1496, 1497]
RUNS_MUON   = [1483, 1484, 1487, 1488, 1491, 1494, 1495, 1498]
REFERENCE_RUN = 1485
IX_TARGETS  = [-3.5, -2.5, 2.5, 3.5]
IY_TARGETS  = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
IX_CELL_CM  = 1.2
IY_CELL_CM  = 1.6
TS_TO_NS    = 0.2
PLOTDIR     = "results/plots/RotationScan"
HTMLDIR     = "results/html/RotationScan"

# indexed by position in IX_TARGETS — constant across all plots
COLORS       = [ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kOrange + 1, ROOT.kRed + 1]
MARKERS      = [20, 21, 22, 23]
MARKER_SIZES = [1.4, 1.6, 1.8, 2.0]
LINE_WIDTH   = 3


def _combine(mpvs, errs):
    weights  = [1.0 / e**2 for e in errs]
    W        = sum(weights)
    return sum(w * m for w, m in zip(weights, mpvs)) / W, 1.0 / W**0.5


def _collect_for_channels(channels, runs, runlist, drsboards):
    """Return (degrees, comb_mpvs, comb_errs) for a channel list over runs."""
    degrees, comb_mpvs, comb_errs = [], [], []
    for run in runs:
        entry = runlist.get(str(run))
        if entry is None:
            continue
        degree    = float(entry["Table Degree"])
        root_path = f"results/root/Run{run}/drs_stats.root"
        if not os.path.exists(root_path):
            continue
        infile = ROOT.TFile(root_path, "READ")
        raw_mpv_map, raw_err_map = {}, {}
        for chan in channels:
            ch = chan.get_channel_name(blsub=False)
            hist = infile.Get(f"hist_{ch}_TS_cfd_mcp_finebins")
            if not hist or hist.GetEntries() < 5:
                continue
            raw_mpv, raw_err = get_hist_mpv(hist)
            raw_mpv_map[ch] = raw_mpv
            raw_err_map[ch] = raw_err
        infile.Close()
        corr_map = subtract_type_mpv(drsboards, raw_mpv_map)
        if not corr_map:
            continue
        run_mpvs, run_errs = [], []
        for ch, corr_ts in corr_map.items():
            err_ts = raw_err_map.get(ch, 0)
            if err_ts <= 0:
                continue
            run_mpvs.append(corr_ts * TS_TO_NS)
            run_errs.append(err_ts * TS_TO_NS)
        if not run_mpvs:
            continue
        comb, comb_e = _combine(run_mpvs, run_errs)
        degrees.append(degree)
        comb_mpvs.append(comb)
        comb_errs.append(comb_e)
    return degrees, comb_mpvs, comb_errs


def _collect_all(runs, drsboards, runlist, iy_target=None, iy_exclude=()):
    """Return {ix_idx: (degrees, mpvs, errs)} for all IX_TARGETS."""
    result = {}
    for ix_idx, ix_target in enumerate(IX_TARGETS):
        channels = [
            chan
            for _, board in drsboards.items()
            for chan in board
            if not chan.is_reference
            and chan.isCer
            and chan.i_tower_x == ix_target
            and (iy_target is None or chan.i_tower_y == iy_target)
            and chan.i_tower_y not in iy_exclude
        ]
        if not channels:
            continue
        degrees, mpvs, errs = _collect_for_channels(channels, runs, runlist, drsboards)
        if degrees:
            result[ix_idx] = (degrees, mpvs, errs)
    return result


def _coord_label(val):
    sign = "m" if val < 0 else "p"
    return f"{sign}{abs(val):.1f}".replace(".", "p")


def _akima_smooth_graph(ix_idx, degrees, mpvs, n_points=300):
    """Akima spline - local slope weighting avoids cubic-spline overshoot near sparse data."""
    pairs = sorted(zip(degrees, mpvs))
    xd = [p[0] for p in pairs]
    yd = [p[1] for p in pairs]
    n = len(xd)

    s = [(yd[i+1] - yd[i]) / (xd[i+1] - xd[i]) for i in range(n-1)]
    if n >= 3:
        s_ext = [3*s[0]-2*s[1], 2*s[0]-s[1]] + s + [2*s[-1]-s[-2], 3*s[-1]-2*s[-2]]
    else:
        s_ext = [s[0], s[0]] + s + [s[0], s[0]]

    t = []
    for i in range(n):
        sm2, sm1, sp0, sp1 = s_ext[i], s_ext[i+1], s_ext[i+2], s_ext[i+3]
        w1 = abs(sp1 - sp0)
        w2 = abs(sm1 - sm2)
        t.append((w1*sm1 + w2*sp0) / (w1 + w2) if w1 + w2 > 1e-14 else (sm1 + sp0) / 2)

    gr = ROOT.TGraph(n_points)
    xmin, xmax = xd[0], xd[-1]
    for k in range(n_points):
        x = xmin + k * (xmax - xmin) / (n_points - 1)
        i = next((j for j in range(n-1) if x <= xd[j+1]), n-2)
        h = xd[i+1] - xd[i]
        dx = x - xd[i]
        dy = yd[i+1] - yd[i]
        c = (3*dy/h - 2*t[i] - t[i+1]) / h
        d = (-2*dy/h + t[i] + t[i+1]) / h**2
        gr.SetPoint(k, x, yd[i] + t[i]*dx + c*dx**2 + d*dx**3)

    gr.SetLineColor(COLORS[ix_idx])
    gr.SetLineWidth(LINE_WIDTH)
    gr.SetName(f"smooth_{ix_idx}")
    return gr


def plot_from_data(data_by_ix, beam_label, lumi_text, output_name,
                   canvas_jsons, ymin, ymax, iy_target=None):
    """Draw one canvas from pre-collected data. Color/marker fixed by IX_TARGETS index."""
    if not data_by_ix:
        return

    all_degrees = [d for ix_idx, (degs, _, _) in data_by_ix.items() for d in degs]

    graphs, labels, colors, markers, smooth_graphs = [], [], [], [], []
    for ix_idx, (degrees, mpvs, errs) in sorted(data_by_ix.items()):
        ix_target = IX_TARGETS[ix_idx]
        gr = ROOT.TGraphErrors(len(degrees))
        gr.SetName(f"gr_{output_name}_iX{ix_target:+.1f}")
        gr.SetTitle("")
        for i, (deg, mpv, err) in enumerate(zip(degrees, mpvs, errs)):
            gr.SetPoint(i, deg, mpv)
            gr.SetPointError(i, 0.0, err)
        gr.SetLineWidth(LINE_WIDTH)
        gr.SetMarkerSize(MARKER_SIZES[ix_idx])
        graphs.append(gr)
        labels.append(f"X = {IX_TARGETS[ix_idx] * IX_CELL_CM:.1f} cm")
        colors.append(COLORS[ix_idx])
        markers.append(MARKERS[ix_idx])
        if len(degrees) >= 2:
            smooth_graphs.append(_akima_smooth_graph(ix_idx, degrees, mpvs))

    extra = ("#it{{Cer, Y = {:.1f} cm}}".format(iy_target * IY_CELL_CM)
             if iy_target is not None else "#it{Cer}")

    DrawHistos(
        graphs, labels,
        min(all_degrees) - 5, max(all_degrees) + 5, "Rotation Angle [degrees]",
        ymin, ymax, "Time [ns]",
        output_name,
        outdir=PLOTDIR,
        drawoptions="EP",
        markerstyles=markers,
        legendoptions=["LPE"] * len(graphs),
        legendPos=[0.33, 0.64, 0.58, 0.88],
        mycolors=colors,
        dology=False,
        usePNG=False,
        lumitext=lumi_text,
        extra_text=extra,
        W_ref=900, H_ref=600,
        canvas_json_store=canvas_jsons,
        extraToDraw=[(gr, "L") for gr in smooth_graphs],
    )


def _global_yrange(all_data_list):
    """Compute a padded y range from a list of data dicts."""
    mpvs = [m for data in all_data_list for _, (_, mvs, _) in data.items() for m in mvs]
    if not mpvs:
        return 0.0, 1.0
    pad  = max(0.05, (max(mpvs) - min(mpvs)) * 0.3 + 0.1)
    ymid = (max(mpvs) + min(mpvs)) / 2
    return ymid - pad, ymid + pad


def main():
    with open("data/Runlist.json") as f:
        runlist = json.load(f)

    drsboards = build_drs_boards(run=REFERENCE_RUN)
    os.makedirs(PLOTDIR, exist_ok=True)
    os.makedirs(HTMLDIR, exist_ok=True)

    # ── pre-collect all data ──────────────────────────────────────────────
    print("Collecting e+ data...")
    eplus_data = {None: _collect_all(RUNS_EPLUS, drsboards, runlist, iy_exclude=(-3.5, 3.5))}
    for iy in IY_TARGETS:
        eplus_data[iy] = _collect_all(RUNS_EPLUS, drsboards, runlist, iy_target=iy)

    print("\nCollecting muon data...")
    muon_data = {None: _collect_all(RUNS_MUON, drsboards, runlist, iy_exclude=(-3.5, 3.5))}
    for iy in IY_TARGETS:
        muon_data[iy] = _collect_all(RUNS_MUON, drsboards, runlist, iy_target=iy)

    # ── shared y ranges ───────────────────────────────────────────────────
    ymin_ep, ymax_ep = _global_yrange(list(eplus_data.values()))
    ymin_mu, ymax_mu = _global_yrange(list(muon_data.values()))
    print(f"\ne+ y range:  [{ymin_ep:.3f}, {ymax_ep:.3f}] ns")
    print(f"mu y range:  [{ymin_mu:.3f}, {ymax_mu:.3f}] ns")

    # ── plot ──────────────────────────────────────────────────────────────
    canvas_jsons = {}

    # combined
    plot_from_data(eplus_data[None], "e^{+}", "e^{+}, 80 GeV",
                   "drs_cer_corr_mpv_vs_rotation_eplus",
                   canvas_jsons, ymin_ep, ymax_ep)
    plot_from_data(muon_data[None], "#mu^{+}", "#mu^{+}, 80 GeV",
                   "drs_cer_corr_mpv_vs_rotation_muon",
                   canvas_jsons, ymin_mu, ymax_mu)

    # per-iY
    for iy in IY_TARGETS:
        tag = _coord_label(iy)
        plot_from_data(eplus_data[iy], "e^{+}", "e^{+}, 80 GeV",
                       f"drs_cer_corr_mpv_vs_rotation_eplus_iY{tag}",
                       canvas_jsons, ymin_ep, ymax_ep, iy_target=iy)
        plot_from_data(muon_data[iy], "#mu^{+}", "#mu^{+}, 80 GeV",
                       f"drs_cer_corr_mpv_vs_rotation_muon_iY{tag}",
                       canvas_jsons, ymin_mu, ymax_mu, iy_target=iy)

    # ── single HTML ───────────────────────────────────────────────────────
    html_path = os.path.join(HTMLDIR, "drs_cer_corr_mpv_vs_rotation.html")
    generate_jsroot_html(
        list(canvas_jsons.keys()),
        canvas_jsons,
        plots_per_row=2,
        output_html=html_path,
        title="DRS Cer Corrected MPV vs Rotation",
        intro_text=(
            "Combined corrected CFD MPV [ns] per iX slice, weighted average "
            "over Cherenkov channels at each iX (and iY) per run.\n"
            f"e^{{+}} runs: {', '.join(str(r) for r in RUNS_EPLUS)}.\n"
            f"#mu^{{+}} runs: {', '.join(str(r) for r in RUNS_MUON)}."
        ),
    )
    print(f"\nSaved JSROOT HTML to {html_path}")


if __name__ == "__main__":
    main()
