"""Compare DRS pulse profiles and CFD timing distributions between runs.

Reads pre-computed combined histograms:
  drs_prof_combined.root       — mean waveform vs raw TS (per size+var)
  drs_prof_corr_combined.root  — mean waveform vs time [ns] (per size+var)
  drs_finebins_combined.root   — CFD TS fine-bin distributions
  drs_finebins_corr_combined.root — CFD time [ns] fine-bin distributions

Per-channel Profile1D histograms live in drs_profiles.root (not read here).

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

_WARM_COLORS = [ROOT.kRed+1, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kOrange-3]
_COLD_COLORS = [ROOT.kBlue+1, ROOT.kAzure+7, ROOT.kGreen+2, ROOT.kTeal+1]
_EM_BTYPES = {"e+", "e-"}

# Color for 3mm+6mm overlay: (is_em, size_tag) -> color
_SIZE_TYPE_COLORS = {
    (True,  "3mm"): ROOT.kRed+1,
    (True,  "6mm"): ROOT.kOrange+1,
    (False, "3mm"): ROOT.kBlue+1,
    (False, "6mm"): ROOT.kAzure+7,
}

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
    colors = []
    for r in runs:
        if _is_em(r):
            colors.append(_WARM_COLORS[warm_idx % len(_WARM_COLORS)])
            warm_idx += 1
        else:
            colors.append(_COLD_COLORS[cold_idx % len(_COLD_COLORS)])
            cold_idx += 1
    return colors


def _make_style(colors, dology=False):
    return PlotStyle(
        dology=dology,
        drawoptions="HIST",
        mycolors=list(colors),
        addOverflow=False,
        addUnderflow=False,
        legendoptions="L",
    )


def _align_peak(hists):
    """Return clones each scaled so the maximum bin content equals 1."""
    result = []
    for h in hists:
        c = h.Clone()
        c.SetDirectory(0)
        peak = c.GetMaximum()
        if abs(peak) > 1e-10:
            c.Scale(1.0 / peak)
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


def _load_ts_prof(run, size_tag, var):
    return _load_hist(run, "drs_prof_combined.root",
                      f"combo_prof_{size_tag}_{var}_mcp")


def _load_prof(run, size_tag, var):
    return _load_hist(run, "drs_prof_corr_combined.root",
                      f"combo_prof_corr_{size_tag}_{var}_mcp")


def _load_finebins(run, size_tag, var):
    return _load_hist(run, "drs_finebins_corr_combined.root",
                      f"combo_finebins_corr_{size_tag}_{var}")


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def run_comparison(runs, labels, colors, use_jsroot=False):
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

    def _collect(loader, size_tag, var):
        """Return (hists, labels, colors) for runs that have the histogram."""
        hs, ls, cs = [], [], []
        for r, lbl, col in zip(runs, labels, colors):
            h = loader(r, size_tag, var)
            if h:
                hs.append(h)
                ls.append(lbl)
                cs.append(col)
        return hs, ls, cs

    # ── Waveform profiles: raw TS ────────────────────────────────────────────
    for size_tag in _SIZES:
        for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
            hs, ls, cs = _collect(_load_ts_prof, size_tag, var)
            if len(hs) < 2:
                continue
            pm.plot_1d(
                _align_peak(hs), f"Prof_ts_{size_tag}_{var}",
                "TS", (0, 1024),
                "Mean DRS Output (a.u.)", (None, None),
                legends=ls, style=_make_style(cs, dology=False),
                extra_text=f"{size_tag} {var}", includeRunNumber=False,
            )

    pm.add_newline()

    # ── Waveform profiles: corr (ns) ─────────────────────────────────────────
    for size_tag in _SIZES:
        for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
            hs, ls, cs = _collect(_load_prof, size_tag, var)
            if len(hs) < 2:
                continue
            xmin = min(h.GetXaxis().GetXmin() for h in hs)
            xmax = max(h.GetXaxis().GetXmax() for h in hs)
            pm.plot_1d(
                _align_peak(hs), f"Prof_corr_{size_tag}_{var}",
                "Time [ns]", (xmin, xmax),
                "Mean DRS Output (a.u.)", (None, None),
                legends=ls, style=_make_style(cs, dology=False),
                extra_text=f"{size_tag} {var}", includeRunNumber=False,
            )

    pm.add_newline()

    # ── CFD finebins: corr (ns) ──────────────────────────────────────────────
    for size_tag in _SIZES:
        for var in sorted(_VARS, key=lambda v: _VAR_ORDER[v]):
            hs, ls, cs = _collect(_load_finebins, size_tag, var)
            if len(hs) < 2:
                continue
            pm.plot_1d(
                _align_peak(hs), f"FineBins_corr_{size_tag}_{var}",
                "Time [ns]", (5, 15),
                "A.U.", (1e-4, None),
                legends=ls, style=_make_style(cs, dology=True),
                extra_text=f"{size_tag} {var}", includeRunNumber=False,
            )

    pm.add_newline()

    # ── CFD finebins: 3mm+6mm overlay per Cer type ───────────────────────────
    for var in ("CerQuartz", "CerPlastic"):
        hists4, legends4, cols4 = [], [], []
        for size_tag in _SIZES:
            for i, r in enumerate(runs):
                h = _load_finebins(r, size_tag, var)
                if h:
                    hists4.append(h)
                    legends4.append(f"{size_tag} {labels[i]}")
                    cols4.append(_SIZE_TYPE_COLORS[(_is_em(r), size_tag)])
        if len(hists4) < 2:
            continue
        pm.plot_1d(
            _align_peak(hists4), f"FineBins_corr_3mm6mm_{var}",
            "Time [ns]", (5, 15),
            "A.U.", (1e-4, None),
            legends=legends4,
            style=_make_style(cols4, dology=True),
            extra_text=var, includeRunNumber=False,
        )

    pm.add_newline()

    # ── CerQuartz + CerPlastic vs time [ns], per size ────────────────────────
    # Color = run particle type (warm/cold); line style = Cer type (solid/dashed)
    _CER_LINESTYLE = {"CerQuartz": 1, "CerPlastic": 2}
    for size_tag in _SIZES:
        hs, ls, cs, styles = [], [], [], []
        for var in ("CerQuartz", "CerPlastic"):
            for r, lbl, col in zip(runs, labels, colors):
                h = _load_prof(r, size_tag, var)
                if h:
                    hs.append(h)
                    ls.append(f"{var} {lbl}")
                    cs.append(col)
                    styles.append(_CER_LINESTYLE[var])
        if len(hs) < 2:
            continue
        xmin = min(h.GetXaxis().GetXmin() for h in hs)
        xmax = max(h.GetXaxis().GetXmax() for h in hs)
        style_cer = PlotStyle(
            drawoptions="HIST",
            mycolors=cs,
            linestyles=styles,
            addOverflow=False,
            addUnderflow=False,
            legendoptions="L",
        )
        pm.plot_1d(
            _align_peak(hs), f"Prof_corr_{size_tag}_CerTypes",
            "Time [ns]", (xmin, xmax),
            "Mean DRS Output (a.u.)", (None, None),
            legends=ls, style=style_cer,
            extra_text=size_tag, includeRunNumber=False,
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
    colors = _build_run_colors(args.runs)
    run_comparison(args.runs, labels, colors, use_jsroot=args.jsroot)
