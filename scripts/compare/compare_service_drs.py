"""Compare service-DRS detector distributions (PSD, TTUMuonVeto, Cer, …) between runs.

Reads pre-computed histograms from drs_services.root produced by check_service_drs.py.
For each PID detector, overlays the energy, peak-value and integral/peak distributions
from all requested runs on the same canvas (normalised to unit area).

Default: example run pair — edit DEFAULT_RUNS or pass --runs on the command line.

Usage
-----
  python scripts/compare_service_drs.py
  python scripts/compare_service_drs.py --runs 1422 1439 --labels "e+_30GeV" "pi+_30GeV"
  python scripts/compare_service_drs.py --runs 1422 1439 --detectors PSD TTUMuonVeto Cer474
"""

import argparse
import os
import sys
import ROOT

from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from configs.plot_config import get_service_drs_processed_info_ranges
from configs.selection_config import get_service_drs_cut

setup_root(batch_mode=True, load_functions=False)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_RUNS = [1422, 1439]

ALL_DETECTORS = [
    "HoleVeto", "PSD", "TTUMuonVeto",
    "Cer474", "Cer519", "Cer537",
    "KT1", "KT2", "T3", "T4",
]

_DET_LABEL = {
    "TTUMuonVeto": "Muon Counter",
}

_EM_BTYPES = {"e+", "e-", "positron", "positrons", "electron", "electrons"}

_WARM_COLORS = [ROOT.kRed+1, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kOrange-3]
_COLD_COLORS = [ROOT.kBlue+1, ROOT.kAzure+7, ROOT.kGreen+2, ROOT.kTeal+1]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _is_em(run):
    try:
        btype, _ = getRunInfo(run)
        return btype.lower() in _EM_BTYPES
    except Exception:
        return False


def _build_run_colors(runs):
    warm_idx, cold_idx = 0, 0
    result = []
    for r in runs:
        if _is_em(r):
            result.append(_WARM_COLORS[warm_idx % len(_WARM_COLORS)])
            warm_idx += 1
        else:
            result.append(_COLD_COLORS[cold_idx % len(_COLD_COLORS)])
            cold_idx += 1
    return result


_BTYPE_LATEX = {
    "e+": "e^{+}", "e-": "e^{-}",
    "positron": "e^{+}", "positrons": "e^{+}",
    "electron": "e^{-}", "electrons": "e^{-}",
    "pi+": "#pi^{+}", "pi-": "#pi^{-}",
    "pion": "#pi^{+}", "pions": "#pi^{+}",
    "proton": "p", "protons": "p",
    "mu+": "#mu^{+}", "mu-": "#mu^{-}",
}


def _build_labels(runs, user_labels):
    if user_labels:
        if len(user_labels) != len(runs):
            print("Error: --labels count must match --runs count")
            sys.exit(1)
        return [l.replace("_", " ") for l in user_labels]
    labels = []
    for r in runs:
        try:
            btype, benergy = getRunInfo(r)
            btype_str = _BTYPE_LATEX.get(btype.lower(), btype)
            labels.append(f"{btype_str} {benergy} GeV")
        except Exception:
            labels.append(f"Run{r}")
    return labels


def _load_hist(run, hist_name, fname="drs_services.root"):
    path = f"results/root/Run{run}/{fname}"
    f = ROOT.TFile(path, "READ")
    if f.IsZombie():
        print(f"  Warning: could not open {path}")
        return None
    h = f.Get(hist_name)
    if h:
        h.SetDirectory(0)
    else:
        print(f"  Warning: histogram '{hist_name}' not found in {path}")
    f.Close()
    return h


def _collect(runs, labels, colors, hist_name):
    """Return (hists, labels, colors) for runs that have the histogram."""
    hs, ls, cs = [], [], []
    for r, lbl, col in zip(runs, labels, colors):
        h = _load_hist(r, hist_name)
        if h:
            hs.append(h)
            ls.append(lbl)
            cs.append(col)
    return hs, ls, cs


def _normalize(hists):
    """Return clones scaled to unit area (overflow + underflow included)."""
    result = []
    for h in hists:
        c = h.Clone()
        c.SetDirectory(0)
        integral = c.Integral(0, c.GetNbinsX() + 1)
        if integral > 1e-10:
            c.Scale(1.0 / integral)
        result.append(c)
    return result


def _norm_style(colors, dology=False):
    return PlotStyle(
        dology=dology,
        drawoptions="HIST",
        mycolors=list(colors),
        addOverflow=True,
        addUnderflow=True,
        legendoptions="L",
    )


def _cut_line(x, dology=False):
    ylo = 1e-5 if dology else 0.0
    yhi = 1.0
    line = ROOT.TLine(x, ylo, x, yhi)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.SetLineStyle(ROOT.kDashed)
    return line


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def run_comparison(runs, labels, colors, detectors, use_jsroot=False):
    out_tag = "_vs_".join(f"Run{r}" for r in runs)

    os.makedirs("results/plots/Compare", exist_ok=True)
    os.makedirs("results/html/Compare", exist_ok=True)

    pm = PlotManager(
        f"results/root/Run{runs[0]}",
        "results/plots/Compare",
        "results/html/Compare",
        run_number=runs[0],
        use_jsroot=use_jsroot,
    )
    pm.set_output_dir(out_tag)

    for det in detectors:
        en_xmin, en_xmax = get_service_drs_processed_info_ranges(det, "sum")
        pv_xmin, pv_xmax = get_service_drs_processed_info_ranges(det, "peak_value")
        _, _, _, _, value_cut, cut_method = get_service_drs_cut(det)

        # ── energy / integral ────────────────────────────────────────────────
        hs, ls, cs = _collect(runs, labels, colors, f"{det}_energy")
        if hs:
            extra = [_cut_line(value_cut, dology=True)] if cut_method == "Sum" else []
            pm.plot_1d(
                _normalize(hs), f"{det}_energy",
                "Integral [ADC]", (en_xmin, en_xmax),
                ylabel="Frac. of Events", yrange=(1e-5, None),
                legends=ls, style=_norm_style(cs, dology=True),
                extra_text=_DET_LABEL.get(det, det), includeRunNumber=False,
                extraToDraw=extra or None,
                legendPos=[0.53, 0.65, 0.82, 0.88],
                usePDF=True,
            )

        # ── peak value ───────────────────────────────────────────────────────
        hs, ls, cs = _collect(runs, labels, colors, f"{det}_peak_value")
        if hs:
            extra = [_cut_line(value_cut)] if cut_method == "PeakValue" else []
            pm.plot_1d(
                _normalize(hs), f"{det}_peak_value",
                "Peak Value [ADC]", (pv_xmin, pv_xmax),
                ylabel="Frac. of Events",
                legends=ls, style=_norm_style(cs),
                extra_text=_DET_LABEL.get(det, det), includeRunNumber=False,
                extraToDraw=extra or None,
                legendPos=[0.53, 0.65, 0.82, 0.88],
            )

        # ── integral / peak ratio ────────────────────────────────────────────
        hs, ls, cs = _collect(runs, labels, colors, f"{det}_integral_to_peak")
        if hs:
            pm.plot_1d(
                _normalize(hs), f"{det}_integral_to_peak",
                "Integral / Peak [TS]", (-5, 20),
                ylabel="Frac. of Events",
                legends=ls, style=_norm_style(cs),
                extra_text=_DET_LABEL.get(det, det), includeRunNumber=False,
                legendPos=[0.53, 0.65, 0.82, 0.88],
            )

        # ── waveform profile (mean ADC vs TS) ────────────────────────────────
        wf_hists = []
        wf_ls, wf_cs = [], []
        for r, lbl, col in zip(runs, labels, colors):
            prof = _load_hist(r, f"{det}_ADC_vs_TS_prof")
            if prof:
                proj = prof.ProjectionX(f"proj_{det}_{r}")
                proj.SetDirectory(0)
                wf_hists.append(proj)
                wf_ls.append(lbl)
                wf_cs.append(col)
        if len(wf_hists) >= 1:
            pm.plot_1d(
                wf_hists, f"{det}_waveform_prof",
                "Time Slice", (0, 1024),
                ylabel="Mean ADC",
                legends=wf_ls,
                style=PlotStyle(
                    dology=False, drawoptions="HIST",
                    mycolors=wf_cs, legendoptions="L",
                    addOverflow=True, addUnderflow=True,
                ),
                extra_text=_DET_LABEL.get(det, det), includeRunNumber=False,
            )

        pm.add_newline()

    output_html = pm.generate_html(
        f"{out_tag}_service_drs.html",
        plots_per_row=4,
        title=f"Service DRS Comparison — {out_tag}",
        intro_text=(
            "Energy, peak value and integral/peak distributions for service DRS "
            "detectors, normalised to unit area. Waveform profiles (mean ADC vs TS) "
            "are shown un-normalised.\n\n"
            "Runs compared: " + ", ".join(f"Run{r} ({lbl})" for r, lbl in zip(runs, labels))
        ),
    )
    print(f"Output: {output_html}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Compare service-DRS detector distributions between runs"
    )
    parser.add_argument(
        "--runs", nargs="+", type=int, default=DEFAULT_RUNS,
        help=f"Run numbers to compare (default: {DEFAULT_RUNS})",
    )
    parser.add_argument(
        "--labels", nargs="+", type=str,
        help="Legend labels for each run (use _ for spaces)",
    )
    parser.add_argument(
        "--detectors", nargs="+", type=str, default=ALL_DETECTORS,
        choices=ALL_DETECTORS, metavar="DET",
        help=f"Detectors to include (default: all). Choices: {ALL_DETECTORS}",
    )
    parser.add_argument(
        "--jsroot", dest="jsroot", action="store_true",
        help="Save interactive JSROOT plots instead of PNG images",
    )
    parser.set_defaults(jsroot=False)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    labels = _build_labels(args.runs, args.labels)
    colors = _build_run_colors(args.runs)

    print(f"\nComparing runs: {args.runs}")
    for r, lbl in zip(args.runs, labels):
        print(f"  Run{r}: {lbl}")
    print(f"Detectors: {args.detectors}\n")

    run_comparison(args.runs, labels, colors, args.detectors, use_jsroot=args.jsroot)
