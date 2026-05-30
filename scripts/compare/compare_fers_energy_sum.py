"""Compare FERS total energy sum distributions (Cer and Sci) between runs.

Reads pre-computed histograms from fers_energy_sum.root produced by make_fers.py.
For each requested gain, produces two plots: one overlaying the Cer total energy
sum across all runs, one for Sci.

Usage
-----
  python scripts/compare_fers_energy_sum.py
  python scripts/compare_fers_energy_sum.py --runs 1422 1439 1445
  python scripts/compare_fers_energy_sum.py --runs 1422 1439 --gains HG LG Mix
  python scripts/compare_fers_energy_sum.py --runs 1422 1439 --labels "e+_60GeV" "pi+_60GeV"
"""

import argparse
import os
import sys
import ROOT

from utils.data_loader import getRunInfo
from utils.root_setup import setup_root
from utils.utils import get_hist_sigma_eff as _sigma_eff
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from plotting.my_function import DrawHistos

setup_root(batch_mode=True, load_functions=False)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULT_RUNS = [1422, 1439]
DEFAULT_GAINS = ["Mix"]

_GAIN_CALIBS = {"HG": False, "LG": False, "Mix": True}

_CAT_LABEL = {"Cer": "Cer", "Sci": "Sci"}  # font 42 in CMS_lumi renders both correctly in PDF

_EM_BTYPES = {"e+", "e-", "positron", "positrons", "electron", "electrons"}

_WARM_COLORS = [ROOT.kRed+1, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kOrange-3]
_COLD_COLORS = [ROOT.kBlue+1, ROOT.kAzure+7, ROOT.kGreen+2, ROOT.kTeal+1]

_BEAM_TO_PARTICLE = {
    "pion": "pion", "pions": "pion", "pi+": "pion", "pi-": "pion",
    "positron": "electron", "positrons": "electron", "e+": "electron", "e-": "electron",
    "proton": "proton", "protons": "proton",
    "mu+": "muon", "mu-": "muon", "muon": "muon", "muons": "muon",
}

_PARTICLE_LABEL = {
    "pion":     "#pi",
    "proton":   "p",
    "electron": "e^{+}",
    "muon":     "#mu^{+}",
}

# ---------------------------------------------------------------------------
# Histogram name helper  (mirrors FERSBoards.get_energy_sum_name)
# ---------------------------------------------------------------------------

_PION_BTYPES = {"pi+", "pi-", "pion", "pions"}


def _is_pion(run):
    try:
        btype, _ = getRunInfo(run)
        return btype.lower() in _PION_BTYPES
    except Exception:
        return False


def _hist_name(gain, is_cer):
    calib = _GAIN_CALIBS[gain]
    cat = "Cer" if is_cer else "Sci"
    name = f"FERS_energy{gain}_{cat}"
    if gain != "Mix":      # pdsub flag only for HG / LG
        name += "_pdsub"
    if calib:
        name += "_calib"
    name += "_sum"
    return f"hist_{name}"


def _xlabel(gain):
    return "Energy Sum [GeV]" if gain == "Mix" else "Energy Sum [ADC]"

# ---------------------------------------------------------------------------
# Helpers shared with compare_service_drs pattern
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


def _default_subdir(run):
    """Return the particle subdir for a run based on its beam type."""
    try:
        btype, _ = getRunInfo(run)
        return _BEAM_TO_PARTICLE.get(btype.lower())
    except Exception:
        return None


def _build_labels(runs, subdirs, user_labels=None):
    if user_labels:
        if len(user_labels) != len(runs):
            print("Error: --labels count must match --runs count")
            sys.exit(1)
        return [l.replace("_", " ") for l in user_labels]
    labels = []
    for r, sd in zip(runs, subdirs):
        try:
            _, benergy = getRunInfo(r)
            particle_str = _PARTICLE_LABEL.get(sd, sd or "")
            labels.append(f"{particle_str} {benergy} GeV".strip())
        except Exception:
            labels.append(f"Run{r}")
    return labels


def _load_hist(run, hist_name, fname="fers_energy_sum.root", subdir=None):
    base = f"results/root/Run{run}/{subdir}" if subdir else f"results/root/Run{run}"
    path = f"{base}/{fname}"
    f = ROOT.TFile(path, "READ")
    if f.IsZombie():
        print(f"  Warning: could not open {path}")
        return None
    h = f.Get(hist_name)
    if h:
        h.SetDirectory(0)
    else:
        print(f"  Warning: '{hist_name}' not found in {path}")
    f.Close()
    return h


def _collect(runs, labels, colors, hist_name, subdirs=None):
    hs, ls, cs = [], [], []
    subdirs = subdirs or [None] * len(runs)
    for r, lbl, col, sd in zip(runs, labels, colors, subdirs):
        h = _load_hist(r, hist_name, subdir=sd)
        if h:
            hs.append(h)
            ls.append(lbl)
            cs.append(col)
    return hs, ls, cs


def _normalize(hists):
    """Clones scaled to unit area (overflow + underflow included)."""
    result = []
    for h in hists:
        c = h.Clone()
        c.SetDirectory(0)
        integral = c.Integral(0, c.GetNbinsX() + 1)
        if integral > 1e-10:
            c.Scale(1.0 / integral)
        result.append(c)
    return result


def _style(colors):
    return PlotStyle(
        dology=False,
        drawoptions="HIST",
        mycolors=list(colors),
        addOverflow=True,
        addUnderflow=True,
        legendoptions="L",
    )

# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def run_comparison(runs, labels, colors, gains, use_jsroot=False, include_dr=True, subdirs=None):
    subdirs = subdirs or [None] * len(runs)
    out_tag = "_vs_".join(
        f"Run{r}_{sd}" if sd else f"Run{r}"
        for r, sd in zip(runs, subdirs or [None]*len(runs))
    )

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

    for gain in gains:
        # ── Overlay distributions ────────────────────────────────────────────
        for cat, is_cer in [("Cer", True), ("Sci", False)]:
            hname = _hist_name(gain, is_cer)
            hs, ls, cs = _collect(runs, labels, colors, hname, subdirs=subdirs)
            if not hs:
                continue

            xmin = min(h.GetXaxis().GetXmin() for h in hs)
            xmax = max(h.GetXaxis().GetXmax() for h in hs)
            normed = _normalize(hs)
            ymax = max(h.GetMaximum() for h in normed) * 1.8

            pm.plot_1d(
                normed,
                f"FERS_ESum_{gain}_{cat}",
                _xlabel(gain), (xmin, xmax),
                ylabel="Frac. of Events",
                yrange=(0, ymax),
                legends=ls,
                style=_style(cs),
                extra_text=_CAT_LABEL.get(cat, cat),
                includeRunNumber=False,
                legendPos=[0.45, 0.65, 0.88, 0.88],
                usePDF=True,
            )

        # ── DR method 1 overlay (pion runs only, Mix gain only) ──────────────
        if include_dr and gain == "Mix":
            dr_hname = _hist_name("Mix", True).replace("Cer", "DR")
            dr_hs, dr_ls, dr_cs = [], [], []
            for r, lbl, col, sd in zip(runs, labels, colors, subdirs):
                if not _is_pion(r):
                    continue
                h = _load_hist(r, dr_hname, fname="fers_dr.root", subdir=sd)
                if h:
                    dr_hs.append(h)
                    dr_ls.append(lbl)
                    dr_cs.append(col)
            if dr_hs:
                dr_normed = _normalize(dr_hs)
                xmin = min(h.GetXaxis().GetXmin() for h in dr_normed)
                xmax = max(h.GetXaxis().GetXmax() for h in dr_normed)
                ymax = max(h.GetMaximum() for h in dr_normed) * 1.8
                pm.plot_1d(
                    dr_normed,
                    f"FERS_ESum_{gain}_DR",
                    _xlabel(gain), (xmin, xmax),
                    ylabel="Frac. of Events",
                    yrange=(0, ymax),
                    legends=dr_ls,
                    style=_style(dr_cs),
                    extra_text="DR",
                    includeRunNumber=False,
                    legendPos=[0.45, 0.65, 0.88, 0.88],
                    usePDF=True,
                )

        # ── Response and resolution vs beam energy (grouped by particle) ────────
        dr_hname = _hist_name("Mix", True).replace("Cer", "DR")
        outdir = pm.get_output_dir()


        def _resp_graph(es, means):
            gr = ROOT.TGraph(len(es))
            for i, (e, m) in enumerate(zip(es, means)):
                gr.SetPoint(i, e, m / e)
            gr.SetLineWidth(2); gr.SetMarkerSize(1.4)
            return gr

        def _reso_graph(es, rmss, means):
            gr = ROOT.TGraph(len(es))
            for i, (e, rv, m) in enumerate(zip(es, rmss, means)):
                gr.SetPoint(i, e, rv / m if m != 0 else 0)
            gr.SetLineWidth(2); gr.SetMarkerSize(1.4)
            return gr

        # collect per particle
        from collections import defaultdict
        raw_by_p = defaultdict(lambda: {'e': [], 'mc': [], 'rc': [], 'ms': [], 'rs': [],
                                        'e_dr': [], 'm_dr': [], 'r_dr': []})
        for r, sd in zip(runs, subdirs):
            key = sd or "inclusive"
            _, benergy = getRunInfo(r)
            hc  = _load_hist(r, _hist_name(gain, True),  subdir=sd)
            hs_ = _load_hist(r, _hist_name(gain, False), subdir=sd)
            if hc is None or hs_ is None:
                continue
            s = raw_by_p[key]
            mc, rc = hc.GetMean(),  _sigma_eff(hc)
            ms, rs = hs_.GetMean(), _sigma_eff(hs_)
            s['e'].append(benergy); s['mc'].append(mc)
            s['rc'].append(rc);     s['ms'].append(ms); s['rs'].append(rs)
            if include_dr and gain == "Mix" and _is_pion(r):
                hdr = _load_hist(r, dr_hname, fname="fers_dr.root", subdir=sd)
                if hdr:
                    s['e_dr'].append(benergy)
                    s['m_dr'].append(hdr.GetMean())
                    s['r_dr'].append(_sigma_eff(hdr))

        # markers: filled for pion, open for others; colours: Cer=red, Sci=blue, DR=green
        _P_MARKERS = {"pion": (20, 21, 22), "proton": (24, 25, 26),
                      "electron": (23, 32, 27), "inclusive": (20, 21, 22)}
        _DEFAULT_MARKERS = (20, 21, 22)

        all_e = [e for s in raw_by_p.values() for e in s['e']]
        if all_e:
            has_pion = any(k == "pion" for k in raw_by_p)
            resp_lo = 0.4 if has_pion else 0.6
            resp_hi = 1.4 if has_pion else 1.2
            reso_hi = 0.5 if has_pion else 0.3
            emin = min(all_e) - 5
            emax = max(all_e) + 5

            ref_line = ROOT.TLine(emin, 1.0, emax, 1.0)
            ref_line.SetLineColor(ROOT.kBlack)
            ref_line.SetLineWidth(2)
            ref_line.SetLineStyle(ROOT.kDashed)

            gr_list_resp, gr_list_reso = [], []
            lbl_rr, col_rr, mk_rr = [], [], []

            for particle_key, s in raw_by_p.items():
                if not s['e']:
                    continue
                plbl = _PARTICLE_LABEL.get(particle_key, particle_key)
                mk_cer, mk_sci, mk_dr = _P_MARKERS.get(particle_key, _DEFAULT_MARKERS)

                gr_list_resp += [_resp_graph(s['e'], s['mc']), _resp_graph(s['e'], s['ms'])]
                gr_list_reso += [_reso_graph(s['e'], s['rc'], s['mc']),
                                 _reso_graph(s['e'], s['rs'], s['ms'])]
                lbl_rr += [f"{plbl} Cer", f"{plbl} Sci"]
                col_rr += [ROOT.kRed+1, ROOT.kBlue+1]
                mk_rr  += [mk_cer, mk_sci]

                if s['e_dr']:
                    gr_list_resp.append(_resp_graph(s['e_dr'], s['m_dr']))
                    gr_list_reso.append(_reso_graph(s['e_dr'], s['r_dr'], s['m_dr']))
                    lbl_rr.append(f"{plbl} DR")
                    col_rr.append(ROOT.kGreen+2)
                    mk_rr.append(mk_dr)

            for plot_name, gr_list, ylabel, ymin, ymax, extra in [
                (f"FERS_Response_{gain}",   gr_list_resp, "Mean / E_{beam}", resp_lo, resp_hi, [ref_line]),
                (f"FERS_Resolution_{gain}", gr_list_reso, "#sigma_{eff} / Mean",      0.0,    reso_hi, []),
            ]:
                DrawHistos(
                    gr_list, lbl_rr,
                    emin, emax, "Beam Energy [GeV]",
                    ymin, ymax, ylabel,
                    plot_name, outdir=outdir,
                    drawoptions="LP", mycolors=col_rr, markerstyles=mk_rr,
                    legendoptions=["LP"] * len(gr_list),
                    dology=False, usePDF=True,
                    extraToDraw=extra or None,
                )
                pm.add_plot(plot_name, with_pdf=True)

        # ── Overlay distributions with Gaussian fits ─────────────────────────
        outdir = pm.get_output_dir()
        fit_by_p = defaultdict(lambda: {
            'e_cer': [], 'm_cer': [], 's_cer': [],
            'e_sci': [], 'm_sci': [], 's_sci': [],
            'e_dr':  [], 'm_dr':  [], 's_dr':  [],
        })

        for cat, is_cer, f_e_key, f_m_key, f_s_key in [
            ("Cer", True,  "e_cer", "m_cer", "s_cer"),
            ("Sci", False, "e_sci", "m_sci", "s_sci"),
        ]:
            hname = _hist_name(gain, is_cer)
            norm_hists, gaussians, ls_fit, cs_fit = [], [], [], []

            for r, lbl, col, sd in zip(runs, labels, colors, subdirs):
                h = _load_hist(r, hname, subdir=sd)
                if h is None:
                    continue
                _, benergy = getRunInfo(r)

                c = h.Clone(f"hfit_{gain}_{cat}_{r}")
                c.SetDirectory(0)
                integral = c.Integral(0, c.GetNbinsX() + 1)
                if integral > 1e-10:
                    c.Scale(1.0 / integral)

                # Fit around expected signal peak
                if _GAIN_CALIBS[gain]:
                    fit_lo, fit_hi = benergy * 0.5, benergy * 1.5
                else:
                    fit_lo = c.GetXaxis().GetXmin()
                    fit_hi = c.GetXaxis().GetXmax()
                lo_bin = c.FindBin(fit_lo)
                hi_bin = c.FindBin(fit_hi)
                if c.Integral(lo_bin, hi_bin) < 1e-10:
                    norm_hists.append(c); ls_fit.append(lbl); cs_fit.append(col)
                    continue
                try:
                    fr = c.Fit("gaus", "SQN", "", fit_lo, fit_hi)
                    if fr.IsValid():
                        amp, mean, sigma = fr.Parameter(0), fr.Parameter(1), fr.Parameter(2)
                        gf = ROOT.TF1(f"gf_{gain}_{cat}_{r}",
                                      "gaus", mean - 3*sigma, mean + 3*sigma)
                        gf.SetParameters(amp, mean, sigma)
                        gf.SetLineColor(col)
                        gf.SetLineWidth(2)
                        gaussians.append(gf)
                        pkey = sd or "inclusive"
                        fit_by_p[pkey][f_e_key].append(benergy)
                        fit_by_p[pkey][f_m_key].append(mean)
                        fit_by_p[pkey][f_s_key].append(sigma)
                except Exception:
                    pass

                norm_hists.append(c)
                ls_fit.append(lbl)
                cs_fit.append(col)

            if norm_hists:
                xmin = min(h.GetXaxis().GetXmin() for h in norm_hists)
                xmax = max(h.GetXaxis().GetXmax() for h in norm_hists)
                ymax = max(h.GetMaximum() for h in norm_hists) * 1.8
                pm.plot_1d(
                    norm_hists,
                    f"FERS_ESum_Fit_{gain}_{cat}",
                    _xlabel(gain), (xmin, xmax),
                    ylabel="Frac. of Events",
                    yrange=(0, ymax),
                    legends=ls_fit,
                    style=_style(cs_fit),
                    extra_text=_CAT_LABEL.get(cat, cat),
                    includeRunNumber=False,
                    legendPos=[0.45, 0.65, 0.88, 0.88],
                    usePDF=True,
                    extraToDraw=gaussians or None,
                )

        # ── Fit DR for pion runs ──────────────────────────────────────────────
        if include_dr and gain == "Mix":
            for r, lbl, col, sd in zip(runs, labels, colors, subdirs):
                if not _is_pion(r):
                    continue
                h = _load_hist(r, dr_hname, fname="fers_dr.root", subdir=sd)
                if h is None:
                    continue
                _, benergy = getRunInfo(r)
                c = h.Clone(f"hfit_{gain}_DR_{r}")
                c.SetDirectory(0)
                integral = c.Integral(0, c.GetNbinsX() + 1)
                if integral > 1e-10:
                    c.Scale(1.0 / integral)
                lo_bin = c.FindBin(benergy * 0.5)
                hi_bin = c.FindBin(benergy * 1.5)
                if c.Integral(lo_bin, hi_bin) < 1e-10:
                    continue
                try:
                    fr = c.Fit("gaus", "SQN", "", benergy * 0.5, benergy * 1.5)
                    if fr.IsValid():
                        pkey = sd or "inclusive"
                        fit_by_p[pkey]['e_dr'].append(benergy)
                        fit_by_p[pkey]['m_dr'].append(fr.Parameter(1))
                        fit_by_p[pkey]['s_dr'].append(fr.Parameter(2))
                except Exception:
                    pass

        # ── Response/resolution from fitted mean and sigma (per particle) ─────
        def _fit_resp_graph(es, ms):
            gr = ROOT.TGraph(len(es))
            for i, (e, m) in enumerate(zip(es, ms)):
                gr.SetPoint(i, e, m / e)
            gr.SetLineWidth(2); gr.SetMarkerSize(1.4)
            return gr

        def _fit_reso_graph(es, ss, ms):
            gr = ROOT.TGraph(len(es))
            for i, (e, s, m) in enumerate(zip(es, ss, ms)):
                gr.SetPoint(i, e, s / m if m != 0 else 0)
            gr.SetLineWidth(2); gr.SetMarkerSize(1.4)
            return gr

        all_e_fit = [e for fs in fit_by_p.values() for e in fs['e_cer'] + fs['e_sci']]
        if all_e_fit:
            has_pion_fit = any(k == "pion" for k in fit_by_p)
            resp_lo_f = 0.4 if has_pion_fit else 0.6
            resp_hi_f = 1.4 if has_pion_fit else 1.2
            reso_hi_f = 0.5 if has_pion_fit else 0.3
            emin_f = min(all_e_fit) - 5
            emax_f = max(all_e_fit) + 5

            ref_fit = ROOT.TLine(emin_f, 1.0, emax_f, 1.0)
            ref_fit.SetLineColor(ROOT.kBlack)
            ref_fit.SetLineWidth(2)
            ref_fit.SetLineStyle(ROOT.kDashed)

            fg_resp, fg_reso, fl, fc, fm_ = [], [], [], [], []
            for particle_key, fs in fit_by_p.items():
                if not fs['e_cer'] or not fs['e_sci']:
                    continue
                plbl = _PARTICLE_LABEL.get(particle_key, particle_key)
                mk_cer, mk_sci, mk_dr = _P_MARKERS.get(particle_key, _DEFAULT_MARKERS)

                fg_resp += [_fit_resp_graph(fs['e_cer'], fs['m_cer']),
                            _fit_resp_graph(fs['e_sci'], fs['m_sci'])]
                fg_reso += [_fit_reso_graph(fs['e_cer'], fs['s_cer'], fs['m_cer']),
                            _fit_reso_graph(fs['e_sci'], fs['s_sci'], fs['m_sci'])]
                fl += [f"{plbl} Cer", f"{plbl} Sci"]
                fc += [ROOT.kRed+1, ROOT.kBlue+1]
                fm_ += [mk_cer, mk_sci]

                if fs['e_dr']:
                    fg_resp.append(_fit_resp_graph(fs['e_dr'], fs['m_dr']))
                    fg_reso.append(_fit_reso_graph(fs['e_dr'], fs['s_dr'], fs['m_dr']))
                    fl.append(f"{plbl} DR"); fc.append(ROOT.kGreen+2); fm_.append(mk_dr)

            for plot_name, gr_list, ylabel, ymin, ymax, extra in [
                (f"FERS_Fit_Response_{gain}",
                 fg_resp, "Fitted #mu / E_{beam}", resp_lo_f, resp_hi_f, [ref_fit]),
                (f"FERS_Fit_Resolution_{gain}",
                 fg_reso, "Fitted #sigma / #mu", 0.0, reso_hi_f, []),
            ]:
                DrawHistos(
                    gr_list, fl,
                    emin_f, emax_f, "Beam Energy [GeV]",
                    ymin, ymax, ylabel,
                    plot_name, outdir=outdir,
                    drawoptions="LP", mycolors=fc, markerstyles=fm_,
                    legendoptions=["LP"] * len(gr_list),
                    dology=False, usePDF=True,
                    extraToDraw=extra or None,
                )
                pm.add_plot(plot_name, with_pdf=True)

        pm.add_newline()

    output_html = pm.generate_html(
        f"{out_tag}_fers_esum.html",
        plots_per_row=2,
        title=f"FERS Energy Sum Comparison — {out_tag}",
        intro_text=(
            "Total FERS energy sum (Cer and Sci) normalised to unit area, "
            "overlaid across runs. Gains shown: " + ", ".join(gains) + ".\n\n"
            "Runs: " + ", ".join(f"Run{r} ({lbl})" for r, lbl in zip(runs, labels))
        ),
    )
    print(f"Output: {output_html}")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    parser = argparse.ArgumentParser(
        description="Compare FERS total energy sum distributions between runs"
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
        "--gains", nargs="+", type=str, default=DEFAULT_GAINS,
        choices=["HG", "LG", "Mix"],
        help=f"Gains to include (default: all three)",
    )
    parser.add_argument(
        "--particles", nargs="+", type=str, default=None, metavar="TYPE",
        help="Particle subdirectory per run (e.g. pion proton). "
             "One value applies to all runs; one per run for mixed comparisons. "
             "Default: inclusive (no subdirectory).",
    )
    parser.add_argument(
        "--dr", dest="include_dr", action="store_true",
        help="Include dual-readout (DR) histograms for pion runs",
    )
    parser.add_argument(
        "--jsroot", dest="jsroot", action="store_true",
        help="Save interactive JSROOT plots instead of PNG+PDF",
    )
    parser.set_defaults(include_dr=False, jsroot=False)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    colors = _build_run_colors(args.runs)

    # Resolve --particles to a per-run subdir list;
    # default: beam particle type for each run
    if args.particles is None:
        subdirs = [_default_subdir(r) for r in args.runs]
    elif len(args.particles) == 1:
        subdirs = args.particles * len(args.runs)
    else:
        if len(args.particles) != len(args.runs):
            print("Error: --particles count must be 1 or match --runs count")
            sys.exit(1)
        subdirs = args.particles

    labels = _build_labels(args.runs, subdirs, user_labels=args.labels)

    print(f"\nComparing runs: {args.runs}")
    for r, lbl, sd in zip(args.runs, labels, subdirs):
        print(f"  Run{r}: {lbl}" + (f"  [{sd}]" if sd else ""))
    print(f"Gains: {args.gains}\n")

    run_comparison(args.runs, labels, colors, args.gains,
                   use_jsroot=args.jsroot, include_dr=args.include_dr,
                   subdirs=subdirs)
