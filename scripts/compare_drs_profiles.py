"""Compare DRS pulse profiles and CFD timing distributions between runs.

Reads pre-computed combined histograms (drs_prof_combined.root,
drs_prof_corr_combined.root, drs_finebins_combined.root,
drs_finebins_corr_combined.root) and overlays them run-by-run.

Default: run 1422 (e+, 30 GeV) vs run 1439 (π+, 30 GeV).

Usage:
  python scripts/compare_drs_profiles.py
  python scripts/compare_drs_profiles.py --runs 1422 1439 --labels "e+" "pion"
"""

import argparse
import os
import sys
import ROOT

from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle

setup_root(batch_mode=True, load_functions=False)

_SIZES = ("3mm", "6mm")
_VARS = ("CerQuartz", "CerPlastic", "Sci")
_VAR_ORDER = {"CerQuartz": 0, "CerPlastic": 1, "Sci": 2}

_RUN_COLORS = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+2, ROOT.kMagenta+1]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_style(n, dology=False):
    return PlotStyle(
        dology=dology,
        drawoptions="HIST",
        mycolors=_RUN_COLORS[:n],
        addOverflow=False,
        addUnderflow=False,
        legendPos=[0.50, 0.80 - n * 0.05, 0.90, 0.90],
        legendoptions="L",
    )


def _normalize(hists):
    """Return clones each scaled to unit integral (area = 1)."""
    result = []
    for h in hists:
        c = h.Clone()
        c.SetDirectory(0)
        integral = c.Integral(0, c.GetNbinsX() + 1)
        if abs(integral) > 1e-10:
            c.Scale(1.0 / integral)
        result.append(c)
    return result


def _load_hist(run_number, fname, hist_name):
    path = f"results/root/Run{run_number}/{fname}"
    f = ROOT.TFile(path, "READ")
    if f.IsZombie():
        print(f"Warning: could not open {path}")
        return None
    h = f.Get(hist_name)
    if h:
        h.SetDirectory(0)
    f.Close()
    return h


def _load_prof(run, size_tag, var):
    return _load_hist(run, "drs_prof_corr_combined.root",
                      f"combo_prof_corr_{size_tag}_{var}_mcp")


def _load_finebins(run, size_tag, var):
    return _load_hist(run, "drs_finebins_corr_combined.root",
                      f"combo_finebins_corr_{size_tag}_{var}")


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def run_comparison(runs, labels, use_jsroot=False):
    n = len(runs)
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

    # ── Waveform profiles: corr (ns) ─────────────────────────────────────────
    for size_tag in _SIZES:
        for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
            hists = []
            for r in runs:
                h = _load_prof(r, size_tag, var)
                if h:
                    hists.append(h)
            if len(hists) < 2:
                continue
            xmin = min(h.GetXaxis().GetXmin() for h in hists)
            xmax = max(h.GetXaxis().GetXmax() for h in hists)
            pm.plot_1d(
                _normalize(hists), f"Prof_corr_{size_tag}_{var}",
                "Time [ns]", (xmin, xmax),
                "Mean DRS Output (a.u.)", (None, None),
                legends=labels[:len(hists)],
                style=_make_style(len(hists), dology=False),
                extra_text=f"{size_tag} {var}",
                includeRunNumber=False,
            )

    pm.add_newline()

    # ── CFD finebins: corr (ns) ──────────────────────────────────────────────
    for size_tag in _SIZES:
        for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
            hists = []
            for r in runs:
                h = _load_finebins(r, size_tag, var)
                if h:
                    hists.append(h)
            if len(hists) < 2:
                continue
            xmin = min(h.GetXaxis().GetXmin() for h in hists)
            xmax = max(h.GetXaxis().GetXmax() for h in hists)
            pm.plot_1d(
                _normalize(hists), f"FineBins_corr_{size_tag}_{var}",
                "Time [ns]", (xmin, xmax),
                "A.U.", (1e-4, None),
                legends=labels[:len(hists)],
                style=_make_style(len(hists), dology=True),
                extra_text=f"{size_tag} {var}",
                includeRunNumber=False,
            )

    output_html = pm.generate_html(f"{out_tag}.html", plots_per_row=6)
    print(f"Output: {output_html}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

_BTYPE_LATEX = {
    "e+": "e^{+}",
    "e-": "e^{-}",
    "pi+": "#pi^{+}",
    "pi-": "#pi^{-}",
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
            btype = _BTYPE_LATEX.get(btype, btype)
            labels.append(f"{btype} {benergy}GeV")
        except Exception:
            labels.append(f"Run{r}")
    return labels


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Compare DRS pulse profiles and CFD timing between runs"
    )
    parser.add_argument("--runs", nargs="+", type=int, default=[1422, 1439],
                        help="Run numbers to compare (default: 1422 1439)")
    parser.add_argument("--labels", nargs="+", type=str,
                        help="Legend labels for each run (use _ for spaces)")
    parser.add_argument("--no-jsroot", dest="jsroot", action="store_false",
                        help="Save PNG images instead of interactive JSROOT plots")
    parser.set_defaults(jsroot=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    labels = _build_labels(args.runs, args.labels)
    run_comparison(args.runs, labels, use_jsroot=args.jsroot)
