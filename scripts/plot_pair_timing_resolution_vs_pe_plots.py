#!/usr/bin/env python3
"""Plots for scripts/plot_pair_timing_resolution_vs_pe.py (no ROOT / event loop).

Reads <indir>/pair_timing_resolution_vs_pe_summary.json and
pair_timing_resolution_vs_pe_pooled_hists.json and writes, under <indir>/plots/:

  res_vs_pos_<mats>_pe<max>[_linfit][_log].png/.pdf
      mats    quartzplastic | withsci
      pe<max> p.e. curves 1..5 or 1..6
      linfit  straight-line fit per curve (weighted least squares), different palette
      log     log-scaled y axis
  dt_hists_<material>.png/.pdf
      pooled dt histograms with the Gaussian fit, rows = p.e., columns = scan
      positions (the "sigma_t histograms at different positions")
  linear_fits.json   slope/intercept of every straight-line fit
  index.html         gallery of all PNGs
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import (  # noqa: E402
    AutoMinorLocator, FuncFormatter, LogFormatter, LogLocator, ScalarFormatter,
)

from utils.html_generator import generate_html  # noqa: E402

LEGEND_NAME = {"CerQuartz": "quartz", "CerPlastic": "plastic", "Sci": "sci"}
PALETTE = {
    "default": {"CerQuartz": "Oranges_r", "CerPlastic": "Blues_r", "Sci": "Greens_r"},
    "linfit": {"CerQuartz": "Reds_r", "CerPlastic": "Purples_r", "Sci": "YlGn_r"},
}
MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
MATERIAL_SETS = {
    "quartzplastic": ["CerQuartz", "CerPlastic"],
    "withsci": ["CerQuartz", "CerPlastic", "Sci"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--indir", required=True,
                        help="Directory with the summary/pooled-hists JSON (e.g. .../edges).")
    parser.add_argument("--extra-indir", nargs="+", default=[],
                        help="Further result directories (e.g. a Sci-only run) merged in.")
    parser.add_argument("--outdir", default=None, help="Default: <indir>/plots")
    parser.add_argument("--n-gallery-positions", type=int, default=6)
    parser.add_argument("--gallery-runs", nargs="+", type=int, default=None,
                        help="Runs for the dt-histogram gallery (default: evenly spaced).")
    return parser.parse_args()


def colors_for(material: str, n_curves: int, palette: str):
    cmap = plt.get_cmap(PALETTE[palette][material])
    return [cmap(0.05 + 0.6 * i / max(1, n_curves - 1)) for i in range(n_curves)]


def style_axes(ax, log: bool) -> None:
    ax.grid(which="major", alpha=0.35)
    ax.grid(which="minor", alpha=0.15, linewidth=0.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if log:
        ax.yaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=12))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(LogFormatter(labelOnlyBase=False, minor_thresholds=(2, 0.4)))
        ax.tick_params(axis="y", which="minor", labelsize=9)
    else:
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", direction="out", top=False, right=False)
    ax.tick_params(which="major", length=7)
    ax.tick_params(which="minor", length=3.5)
    ax.spines[["top", "right"]].set_visible(False)


def weighted_linear_fit(x, y, err):
    x, y, err = (np.asarray(v, float) for v in (x, y, err))
    w = 1.0 / np.clip(err, 1e-3, None)
    coef, cov = np.polyfit(x, y, 1, w=w, cov="unscaled")
    resid = (y - np.polyval(coef, x)) / np.clip(err, 1e-3, None)
    ndf = len(x) - 2
    return {
        "slope_ps_per_mm": float(coef[0]), "slope_err": float(np.sqrt(cov[0, 0])),
        "intercept_ps": float(coef[1]), "intercept_err": float(np.sqrt(cov[1, 1])),
        "chi2": float(np.sum(resid ** 2)), "ndf": int(ndf),
    }


def make_resolution_plot(combined, meta, materials, pe_max, linfit, log, outbase, fit_log):
    palette = "linfit" if linfit else "default"
    pe_values = [n for n in meta["pe_values"] if n <= pe_max]
    fig, ax = plt.subplots(figsize=(14, 8))
    ymin, ymax = np.inf, 0.0
    for material in materials:
        colors = colors_for(material, len(pe_values), palette)
        for i, n in enumerate(pe_values):
            pts = sorted(
                (c for c in combined if c["material"] == material and c["pe"] == n
                 and c["position_mm"] is not None),
                key=lambda c: c["position_mm"])
            if not pts:
                continue
            x = [p["position_mm"] for p in pts]
            y = [p["resolution_ps"] for p in pts]
            e = [p["resolution_err_ps"] for p in pts]
            ymin, ymax = min(ymin, min(y)), max(ymax, max(y))
            label = "{}  {} p.e.".format(LEGEND_NAME[material], n)
            ax.errorbar(
                x, y, yerr=e, color=colors[i], marker=MARKERS[i % len(MARKERS)],
                markersize=6, capsize=3, linewidth=0 if linfit else 1.5,
                elinewidth=1.0, alpha=0.75 if linfit else 1.0, label=label)
            if linfit and len(x) >= 3:
                fit = weighted_linear_fit(x, y, e)
                xx = np.array([min(x), max(x)])
                ax.plot(xx, fit["slope_ps_per_mm"] * xx + fit["intercept_ps"],
                        color=colors[i], linestyle="--", linewidth=2.0)
                fit_log.setdefault(outbase.name, {})["{} {} p.e.".format(
                    LEGEND_NAME[material], n)] = fit

    ax.set_xlabel(meta["position_label"] + " (mm)", fontsize=16)
    ax.set_ylabel("Time resolution (ps)", fontsize=16)
    ax.set_title("Time Resolution - position scan with {} at 90°  [PE: {}{}]".format(
        meta.get("beam") or "", meta["pe_definition"], ", linear fits" if linfit else ""),
        fontsize=15)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.8 * ymin, (6.0 if len(materials) > 2 else 2.5) * ymax)
    else:
        ax.set_ylim(0, 1.4 * ymax)
    style_axes(ax, log)
    ax.legend(ncol=len(materials), fontsize=12, frameon=False, loc="upper right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("{}.{}".format(outbase, ext), dpi=150)
    plt.close(fig)


def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def make_dt_gallery(material, combined, hists, meta, runs, outbase):
    pe_values = meta["pe_values"]
    lo, hi = hists["dt_range"]
    nbins = hists["dt_bins"]
    edges = np.linspace(lo, hi, nbins + 1)
    pos_of = {c["run"]: c["position_mm"] for c in combined}
    fig, axes = plt.subplots(len(pe_values), len(runs),
                             figsize=(3.4 * len(runs), 2.5 * len(pe_values)), squeeze=False)

    for r, n in enumerate(pe_values):
        fits = {c["run"]: c for c in combined if c["material"] == material and c["pe"] == n}
        sigmas = [fits[run]["sigma_ns"] for run in runs if run in fits]
        half = max(0.4, 6.0 * max(sigmas)) if sigmas else 1.5
        half = min(half, hi)
        for c, run in enumerate(runs):
            ax = axes[r][c]
            counts = hists["counts"].get("{}|{}|{}".format(run, material, n))
            fit = fits.get(run)
            if counts is None:
                ax.text(0.5, 0.5, "no events", transform=ax.transAxes, ha="center")
            else:
                factor = fit["rebin"] if fit else 1
                counts = np.asarray(counts, float).reshape(-1, factor).sum(axis=1)
                ax.stairs(counts, edges[::factor], color="tab:blue", linewidth=1.3)
                if fit:
                    xx = np.linspace(-half, half, 400)
                    ax.plot(xx, gaussian(xx, fit["amplitude"], fit["mean_ns"], fit["sigma_ns"]) + fit["background"],
                            color="tab:red", linewidth=1.6)
                    ax.text(0.03, 0.95,
                            "N={}\n$\\sigma_{{dt}}$={:.0f} ps\nres={:.0f} ps".format(
                                fit["n_entries"], 1000 * fit["sigma_ns"], fit["resolution_ps"]),
                            transform=ax.transAxes, va="top", fontsize=8)
                else:
                    ax.text(0.5, 0.5, "no fit\n(N={})".format(int(counts.sum())),
                            transform=ax.transAxes, ha="center", fontsize=9)
            ax.set_xlim(-half, half)
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(which="both", direction="out")
            ax.spines[["top", "right"]].set_visible(False)
            if r == 0:
                ax.set_title("{} = {:.0f} mm (run {})".format(
                    meta["position_label"].replace("Table position ", "").replace("Beam ", ""),
                    pos_of.get(run, float("nan")), run), fontsize=9)
            if c == 0:
                ax.set_ylabel("{} p.e.\nEvents".format(n), fontsize=10)
            if r == len(pe_values) - 1:
                ax.set_xlabel("$\\Delta t$ [ns]", fontsize=10)
    fig.suptitle("Pooled adjacent-pair dt, {} [PE: {}]".format(
        LEGEND_NAME[material], meta["pe_definition"]), fontsize=13)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig("{}.{}".format(outbase, ext), dpi=130)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    with (indir / "pair_timing_resolution_vs_pe_summary.json").open() as handle:
        meta = json.load(handle)
    with (indir / "pair_timing_resolution_vs_pe_pooled_hists.json").open() as handle:
        hists = json.load(handle)
    for extra in args.extra_indir:
        with (Path(extra) / "pair_timing_resolution_vs_pe_summary.json").open() as handle:
            extra_meta = json.load(handle)
        with (Path(extra) / "pair_timing_resolution_vs_pe_pooled_hists.json").open() as handle:
            extra_hists = json.load(handle)
        for key in ("pe_definition", "pair_pe_mode", "pe_values"):
            if extra_meta[key] != meta[key]:
                raise SystemExit("{}: {} differs from {}".format(extra, key, args.indir))
        meta["combined"] += extra_meta["combined"]
        meta["materials"] += [m for m in extra_meta["materials"] if m not in meta["materials"]]
        hists["counts"].update(extra_hists["counts"])
    combined: List[Dict[str, Any]] = meta["combined"]
    meta.setdefault("position_label", meta.get("position_key", "Position"))
    available = set(meta["materials"])

    png_files: List[str] = []
    fit_log: Dict[str, Any] = {}
    for mats_name, materials in MATERIAL_SETS.items():
        if not set(materials) <= available:
            continue
        for pe_max in (5, 6):
            if pe_max > max(meta["pe_values"]):
                continue
            for linfit in (False, True):
                for log in (False, True):
                    name = "res_vs_pos_{}_pe{}{}{}".format(
                        mats_name, pe_max, "_linfit" if linfit else "", "_log" if log else "")
                    make_resolution_plot(combined, meta, materials, pe_max, linfit, log,
                                         outdir / name, fit_log)
                    png_files.append(name + ".png")
                    print("[WROTE] {}".format(outdir / name))

    runs_with_fits = sorted({c["run"] for c in combined if c["position_mm"] is not None},
                            key=lambda r: next(c["position_mm"] for c in combined if c["run"] == r))
    if args.gallery_runs:
        gallery_runs = args.gallery_runs
    else:
        idx = np.unique(np.linspace(0, len(runs_with_fits) - 1, args.n_gallery_positions).round().astype(int))
        gallery_runs = [runs_with_fits[i] for i in idx]
    for material in meta["materials"]:
        name = "dt_hists_{}".format(material)
        make_dt_gallery(material, combined, hists, meta, gallery_runs, outdir / name)
        png_files.append(name + ".png")
        print("[WROTE] {}".format(outdir / name))

    with (outdir / "linear_fits.json").open("w") as handle:
        json.dump(fit_log, handle, indent=2)

    generate_html(
        png_files, str(outdir), plots_per_row=2, output_html=str(outdir / "index.html"),
        title="Pair timing resolution vs. position [PE: {}]".format(meta["pe_definition"]),
        intro_text=(
            "**PE definition:** {}; **pair mode:** {}; pooled over adjacent pairs of each material.\n"
            "**res_vs_pos_*:** `withsci` includes Sci, `pe6` adds the 6 p.e. curve, `linfit` adds "
            "weighted straight-line fits (dashed, alternate palette), `log` is log-y.\n"
            "**dt_hists_*:** pooled dt histograms with Gaussian fit at {} scan positions "
            "(runs {}).".format(meta["pe_definition"], meta["pair_pe_mode"], len(gallery_runs),
                                gallery_runs)),
    )
    print("[WROTE] {}".format(outdir / "index.html"))


if __name__ == "__main__":
    main()
