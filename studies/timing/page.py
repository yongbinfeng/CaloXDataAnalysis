"""Page of the timing-resolution study: breakdown of the timing resolution into the MCP reference, the DRS group / board alignment,
the fibre dispersion and the rest; the electronics measurements (electronics.json) and the dispersion fits (dispersion_<dataset>.json).
Breakdown per class (fibre type, p.e. bin), median over the channels of the class, in quadrature, at the far end of the scan:
  total     measured width at the near end of the scan (fit at the largest X) plus the dispersion over the scan
  MCP       one MCP with its digitization: the 2025 same-chip pair (both MCPs ~125 ADC) / sqrt 2, for all periods (in 2026 the upstream MCP
            is too weak to separate the reference MCP from it)
  DRS       alignment of the channel's group to the MCP's group: 0 in the MCP's group, otherwise the width of the two copies of the same
            MCP signal on two boards after the alignment (interpolated reference), which is what the alignment adds to real pulses
  dispersion  sqrt(|per-photon slope| x scan length / N), per-photon slope of the channel (all p.e. bins), N of the bin: what the
            dispersion adds between the near and the far end of the scan
  other     the rest (SiPM single-photon jitter, CFD and noise, shower and beam-spot spread, photostatistics of the light, and the
            dispersion up to the near end of the scan)
Usage (repository root): python3 studies/timing/page.py      Page: results/html/Timing/Timing_resolution.html"""
import json, os, sys, types
import numpy as np
from config import DATA, PLOTS, HTML, DISP_DATASETS
sys.path.insert(0, ".")
import ROOT
from plotting.my_function import DrawHistos
from plotting.calox_plot_helper import plot_fiber_summary, fiber_base_type
from utils.html_generator import generate_html

ROOT.gROOT.SetBatch(True)
OUTD = PLOTS; OUT = f"{HTML}/Timing_resolution.html"
K = 1e3                                  # ns^2/m -> ps^2/mm
TLAB = {"Sci": "scintillating", "CerQuartz": "Cherenkov quartz", "CerPlastic": "Cherenkov plastic", "Cer": "Cherenkov"}
COL = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 1]
MIN_COPY_AMP = 100                       # [ADC] median amplitude of both copies for the alignment term
PERIOD = {"2025": "2025", "2026A": "2026A", "2026B": "2026B", "testfibres": "2026B"}   # electronics period of each dataset
draw = dict(outdir=OUTD, usePDF=False)


def graph(xs, ys, es=None):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, float(x), float(y)); g.SetPointError(i, 0, float(es[i]) if es is not None else 0)
    return g


def hist(name, edges, counts, norm=True):
    h = ROOT.TH1D(name, "", len(edges) - 1, np.asarray(edges, float))
    for i, c in enumerate(counts): h.SetBinContent(i + 1, float(c))
    if norm and h.Integral() > 0: h.Scale(1 / h.Integral())
    return h


TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
details = lambda title, body: f'<details style="margin:10px 0"><summary><b>{title}</b></summary>{body}</details>'
pm = lambda v, e: ((f"{K * v:+.1f} &plusmn; {K * e:.1f}" if K * e < 5 else f"{K * v:+.0f} &plusmn; {K * e:.0f}") if v is not None and np.isfinite(v) else "&ndash;")
ps = lambda v: f"{v:.0f}" if v is not None and np.isfinite(v) else "&ndash;"
gname = lambda ch: ch.rsplit("_Channel", 1)[0]
bname = lambda ch: gname(ch).rsplit("_Group", 1)[0]


def electronics_terms(E):
    """per period: MCP term, per-group alignment, board-to-board alignment [ps] (medians over the runs)"""
    out = {}
    for p, v in E.items():
        rr = v["runs"]
        same = np.median([r["mcp_same_chip"]["sigma"] for r in rr]) / np.sqrt(2)
        grp = np.median([r["groups"]["ref_fw"]["per_group"] for r in rr]); grp_int = np.median([r["groups"]["ref_int"]["per_group"] for r in rr])
        # pairs of copies of the same MCP on two boards, both with large pulses (the weak 2026A upstream MCP measures its own noise)
        cop = [c for r in rr for k, c in r["copies"].items() if min(r["amp"][m] for m in k.split("-")) > MIN_COPY_AMP]
        board = np.median([c["ref_fw"]["core"] for c in cop]) if cop else np.nan
        board_int = np.median([c["ref_int"]["core"] for c in cop]) if cop else np.nan
        bl = sorted({b for c in cop for b in c["boards"]})
        out[p] = dict(mcp_pair=same * np.sqrt(2), group=grp, group_int=grp_int, board=board, board_int=board_int, boards=bl, n_pairs=len(cop))
    for p in out: out[p]["mcp"] = out["2025"]["mcp_pair"] / np.sqrt(2) if "2025" in out else np.nan
    return out


def drs_term(ch, mcp_ref, et):
    if gname(ch) == gname(mcp_ref): return 0.0, "same group"
    return et["board"], "same board" if bname(ch) == bname(mcp_ref) else "other board"


def breakdown(D, ET):
    """classes (dataset, fibre type, p.e. bin) -> medians of the terms over the channels"""
    rows = []
    for n, T in D.items():
        et = ET[PERIOD[n]]; X = np.array(T["x_m"]); L = 1000 * (X.max() - X.min()); G = sorted(set(T["groups"])); gi = [G.index(g) for g in T["groups"]]
        byclass = {}
        for ch, c in T["channels"].items():
            if c["type"] == "Sci": continue
            cls = fiber_base_type(c["label"]) if n == "testfibres" else TLAB.get(c["type"], c["type"])
            drs, rel = drs_term(ch, T["mcp_ref"], et)
            for k, b in enumerate(c["bins"]):
                if "slope" not in b or k not in (0, 2): continue
                s = np.array(b["sigma"], float); ok = np.where(np.isfinite(s))[0]
                if not len(ok): continue
                i = ok[np.argmax(X[ok])]; key = str(gi[i])
                if key not in b["icpt"]: continue
                near2 = (b["icpt"][key] + b["slope"] * X[i]) * 1e6               # ps^2
                # from the channel's per-photon slope (all p.e. bins together) / N of the bin: less noisy than the bin's own slope
                pp, epp = c.get("per_photon", 0.0), c.get("per_photon_err", 0.0)
                d2 = lambda v: max(-v, 0) * 1e3 * L / b["n_mean"]                 # ps^2/mm x mm
                byclass.setdefault((cls, k), []).append(dict(near2=near2, disp2=d2(pp), disp2_lo=d2(pp + epp), disp2_hi=d2(pp - epp), drs=drs, rel=rel,
                                                             n_mean=b["n_mean"]))
        for (cls, k), v in sorted(byclass.items()):
            med = lambda key: float(np.median([x[key] for x in v]))
            near, disp, drs = np.sqrt(max(med("near2"), 0)), np.sqrt(med("disp2")), med("drs")
            total = np.sqrt(near ** 2 + disp ** 2); other2 = near ** 2 - et["mcp"] ** 2 - drs ** 2
            rows.append(dict(dataset=n, cls=cls, bin=k, pe="1 p.e." if k == 0 else "3-4 p.e.", n_ch=len(v), n_mean=med("n_mean"), L=L,
                             disp_lo=np.sqrt(med("disp2_lo")), disp_hi=np.sqrt(med("disp2_hi")),
                             total=total, near=near, mcp=et["mcp"], drs=drs, disp=disp, other=np.sqrt(max(other2, 0)),
                             rels=", ".join(f"{r} {sum(x['rel'] == r for x in v)}" for r in ("same group", "same board", "other board") if any(x["rel"] == r for x in v))))
    return rows


def sigma_vs_x(T, ch, name):
    c = T["channels"][ch]; X = np.array(T["x_m"]) * 1000; gr = T["groups"]; G = sorted(set(gr)); hs, ls, cols, ex = [], [], [], []
    for k, b in enumerate(c["bins"]):
        s, e = np.array(b["sigma"], float), np.array(b["sigma_err"], float); ok = np.isfinite(s)
        if ok.sum() < 3: continue
        lo, hi = b["pe"]; hs.append(graph(X[ok], 1000 * s[ok], 1000 * e[ok]))
        ls.append(f"{lo + 0.5:.0f}" + (f"-{hi - 0.5:.0f}" if hi - lo > 1.5 and hi < 1e8 else ("+" if hi > 1e8 else "")) + " p.e."); cols.append(COL[k % len(COL)])
        if "slope" in b:
            for gk, g in enumerate(G):                # intercepts keyed by the index of the stretch (sorted names)
                m = ok & (np.array(gr) == g)
                if not m.any() or str(gk) not in b["icpt"]: continue
                f = ROOT.TF1(f"f_{name}_{k}_{gk}", f"1000*sqrt(TMath::Max({b['icpt'][str(gk)]} + {b['slope']}*x/1000., 0.))", X[m].min() - 10, X[m].max() + 10)
                f.SetLineColor(COL[k % len(COL)]); f.SetLineWidth(2); ex.append((f, "SAME"))
    if not hs: return None
    allv = np.concatenate([np.array(b["sigma"], float) for b in c["bins"]]); allv = allv[np.isfinite(allv)] * 1000
    lab = c["label"].replace("DRS_", "").replace("Brg", "B").replace("_Board", " Bd").replace("_Group", " G").replace("_Channel", " C")
    DrawHistos(hs, ls, X.min() - 20, X.max() + 20, "X [mm] (readout #rightarrow)", 0, 1.6 * np.percentile(allv, 98), "Width of t_{CFD} - t_{MCP} [ps]", name,
               dology=False, drawoptions="P", legendoptions="P", mycolors=cols, markerstyles=[20, 21, 22, 23, 33][:len(hs)], legendPos=[0.55, 0.70, 0.92, 0.89],
               extraToDraw=ex, lumitext=f"{T['dataset']}: {lab}", **draw)
    return name


def main():
    os.makedirs(OUTD, exist_ok=True)
    E = json.load(open(f"{DATA}/electronics.json")); ET = electronics_terms(E)
    D = {n: json.load(open(f"{DATA}/dispersion_{n}.json")) for n in DISP_DATASETS if os.path.exists(f"{DATA}/dispersion_{n}.json")}
    B = breakdown(D, ET); sections = []

    # ---- breakdown plot: every term in ps (they add in quadrature to the total), side by side per class ----
    plots = []
    comps = [("total", "total (far end of the scan)", ROOT.kBlack), ("mcp", "MCP reference", ROOT.kGray + 1), ("drs", "DRS alignment", ROOT.kAzure + 1),
             ("disp", "dispersion over the scan", ROOT.kRed + 1), ("other", "other (SiPM, CFD, shower)", ROOT.kGreen - 6)]
    wbar = 0.8 / len(comps)
    for n in D:
        rows = [r for r in B if r["dataset"] == n]
        if not rows: continue
        hs = []
        for k, (key, lab, col) in enumerate(comps):
            h = ROOT.TH1F(f"bd_{n}_{key}", "", len(rows), 0, len(rows))
            for i, r in enumerate(rows):
                h.SetBinContent(i + 1, r[key]); h.GetXaxis().SetBinLabel(i + 1, f"{r['cls']} {r['pe']}")
            h.SetFillColor(col); h.SetLineColor(col); h.SetBarWidth(wbar); h.SetBarOffset(0.1 + k * wbar); h.SetStats(0); hs.append(h)
        ymax = 1.15 * max(r["total"] for r in rows)
        c = ROOT.TCanvas(f"c_bd_{n}", "", 1200, 650); c.SetBottomMargin(0.30); c.SetLeftMargin(0.08); c.SetRightMargin(0.25)
        hs[0].SetMaximum(ymax); hs[0].SetMinimum(0); hs[0].GetYaxis().SetTitle("width [ps]"); hs[0].GetXaxis().LabelsOption("v"); hs[0].GetXaxis().SetLabelSize(0.04)
        hs[0].Draw("BAR")
        for h in hs[1:]: h.Draw("BAR SAME")
        leg = ROOT.TLegend(0.76, 0.55, 0.99, 0.90); leg.SetBorderSize(0); leg.SetTextSize(0.03)
        for h, (_, lab, _) in zip(hs, comps): leg.AddEntry(h, lab, "f")
        leg.Draw(); lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextAlign(11); lat.SetTextSize(0.035)
        lat.DrawLatex(0.08, 0.95, f"{n}: timing resolution and its parts (they add in quadrature)")
        c.SaveAs(f"{OUTD}/breakdown_{n}.png"); plots.append(f"breakdown_{n}")
    sections.append(("Breakdown of the timing resolution", plots))

    # ---- electronics ----
    plots = []
    e10 = np.arange(-600, 601, 10)
    hs, ls = [], []
    for p, v in E.items():
        r = v["runs"][0]; hs.append(hist(f"mcp_{p}", e10, r["mcp_same_chip"]["hist"])); ls.append(f"{p} (run {r['run']}): {r['mcp_same_chip']['sigma']:.0f} ps")
    DrawHistos(hs, ls, -600, 600, "MCP_{US1} - MCP_{DS1} [ps]", 1e-4, 1, "Fraction of events", "mcp_same_chip", dology=True, drawoptions="HIST",
               mycolors=COL[1:4], legendPos=[0.18, 0.74, 0.62, 0.88], lumitext="two MCPs, same DRS chip", **draw)
    plots.append("mcp_same_chip")
    hs, ls = [], []
    for p, v in E.items():
        r = v["runs"][0]; c_ = r["copies"].get("MCP_DS_0-MCP_DS_1")
        if c_: hs.append(hist(f"cop_{p}", e10, c_["hist"])); ls.append(f"{p}: {c_['ref_fw']['core']:.0f} ps")
    DrawHistos(hs, ls, -600, 600, "MCP_{DS0} - MCP_{DS1} (aligned) [ps]", 1e-4, 1, "Fraction of events", "mcp_copies", dology=True, drawoptions="HIST",
               mycolors=COL[1:4], legendPos=[0.18, 0.74, 0.62, 0.88], lumitext="same MCP on two DRS boards", **draw)
    plots.append("mcp_copies")
    labs, wi, wf, used = [], [], [], []
    for p, v in E.items():
        r = v["runs"][0]
        for k, c in r["copies"].items():
            labs.append(f"{p} {k.split('-')[0][4:6]} " + "/".join(b.replace("DRS_", "").replace("Brg", "Bridge").replace("_Board", " B").replace("Board", "B") for b in c["boards"]))
            wi.append(c["ref_int"]["core"]); wf.append(c["ref_fw"]["core"]); used.append(min(r["amp"][m] for m in k.split("-")) > MIN_COPY_AMP)
    hi_ = ROOT.TH1F("cop_int", "", len(labs), 0, len(labs)); hf_ = ROOT.TH1F("cop_fw", "", len(labs), 0, len(labs))
    hx_ = ROOT.TH1F("cop_fw_x", "", len(labs), 0, len(labs))          # pairs left out of the alignment term (weak pulses)
    for i, (l, a1, a2, u) in enumerate(zip(labs, wi, wf, used)):
        hi_.SetBinContent(i + 1, a1); (hf_ if u else hx_).SetBinContent(i + 1, a2); hi_.GetXaxis().SetBinLabel(i + 1, l)
    c_ = ROOT.TCanvas("c_cop", "", 1200, 650); c_.SetBottomMargin(0.36); c_.SetLeftMargin(0.11); c_.SetRightMargin(0.03); c_.SetTopMargin(0.08)
    for h, lc, fc in ((hi_, ROOT.kGray + 2, ROOT.kGray), (hf_, ROOT.kBlue + 1, ROOT.kAzure + 1), (hx_, ROOT.kAzure - 4, ROOT.kAzure - 9)):
        h.SetLineColor(lc); h.SetFillColor(fc); h.SetBarWidth(0.8); h.SetBarOffset(0.1); h.SetStats(0)
    hi_.SetMaximum(160); hi_.SetMinimum(0); hi_.GetYaxis().SetTitle("width of the copy difference [ps]"); hi_.GetYaxis().SetTitleSize(0.045)
    hi_.GetYaxis().SetTitleOffset(1.1); hi_.GetYaxis().SetLabelSize(0.04); hi_.GetXaxis().LabelsOption("v"); hi_.GetXaxis().SetLabelSize(0.035)
    hi_.Draw("BAR"); hf_.Draw("BAR SAME"); hx_.Draw("BAR SAME")
    lg = ROOT.TLegend(0.13, 0.76, 0.75, 0.90); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetNColumns(2); lg.SetTextSize(0.032)
    lg.AddEntry(hi_, "integer reference (former)", "f"); lg.AddEntry(hf_, "interpolated reference", "f")
    lg.AddEntry(hx_, f"interpolated, weak MCP (< {MIN_COPY_AMP} ADC), not used", "f"); lg.Draw()
    tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextSize(0.04); tl.DrawLatex(0.11, 0.94, "same MCP signal on two DRS boards, after the group alignment")
    c_.SaveAs(f"{OUTD}/copy_pairs.png"); plots.append("copy_pairs")
    r = E["2026B"]["runs"][0] if "2026B" in E else next(iter(E.values()))["runs"][0]; ex = r["groups"]["example"]; e5 = np.arange(-500, 501, 10)
    hs = [hist(f"grp_{v}", e5, ex[v]) for v in ("ref_int", "ref_fw", "ref_fit")]
    DrawHistos(hs, [f"integer sample ({r['groups']['ref_int']['pair']:.0f} ps)", f"interpolated ({r['groups']['ref_fw']['pair']:.0f} ps)", f"line fit ({r['groups']['ref_fit']['pair']:.0f} ps)"],
               -500, 500, "Reference time, group A - group B [ps]", 1e-4, 2, "Fraction of events", "groups", dology=True, drawoptions="HIST",
               mycolors=COL[1:4], legendPos=[0.18, 0.74, 0.66, 0.88], lumitext=f"run {r['run']}, two groups of one board", **draw)
    plots.append("groups")
    cr = next((x for v in E.values() for x in v["runs"] if "cells" in x), None)
    if cr:
        prof = cr["cells"]["profile"]; cc = [i for i, v in enumerate(prof) if v is not None]
        DrawHistos([graph(cc, [prof[i] for i in cc])], [f"{cr['cells']['raw']:.0f} ps; per-cell correction: {cr['cells']['corrected']:.0f} ps"], 0, 1024, "DRS cell of the edge (group A)",
                   -300, 300, "median group A - group B [ps]", "cells", dology=False, drawoptions="P", legendoptions="P", mycolors=[ROOT.kBlack],
                   legendPos=[0.18, 0.80, 0.80, 0.88], lumitext=f"run {cr['run']}: DRS cell pattern", **draw)
        plots.append("cells")
    sections.append(("Electronics: MCP and DRS alignment", plots))

    # ---- dispersion ----
    plots = []
    for n, T in D.items():
        hs, ls, cols = [], [], []
        for i, (t, v) in enumerate(T["types"].items()):
            if t == "Sci": continue
            q = [(np.mean([c["bins"][k].get("n_mean", np.nan) for c in T["channels"].values() if c["type"] == t and "n_mean" in c["bins"][k]]), b)
                 for k, b in enumerate(v["per_bin"]) if b]
            q = [(x, b) for x, b in q if np.isfinite(x)]
            if not q: continue
            hs.append(graph([x for x, _ in q], [K * b["value"] for _, b in q], [K * b["err"] * np.sqrt(max(b["chi2"], 1)) if np.isfinite(b["chi2"]) else K * b["err"] for _, b in q]))
            ls.append(TLAB.get(t, t)); cols.append(COL[(i + 1) % len(COL)])
        if not hs or n == "testfibres": continue
        z = ROOT.TLine(0, 0, 12, 0); z.SetLineStyle(2)
        DrawHistos(hs, ls, 0, 12, "Mean p.e. of the bin", -400, 200, "N #times d#sigma^{2}/dX [ps^{2}/mm]", f"per_photon_{n}", dology=False, drawoptions="P",
                   legendoptions="P", mycolors=cols, markerstyles=[20, 21, 22][:len(hs)], legendPos=[0.50, 0.75, 0.92, 0.89], extraToDraw=[z],
                   lumitext=f"{n}: per-photon slope by p.e. bin", **draw)
        plots.append(f"per_photon_{n}")
    if "testfibres" in D:
        T = D["testfibres"]; pmo = types.SimpleNamespace(use_jsroot=False, get_output_dir=lambda: OUTD)
        vals = {c["label"]: (K * c["per_photon"], K * c["per_photon_err"] * np.sqrt(max(c["per_photon_chi2"], 1) if np.isfinite(c["per_photon_chi2"]) else 1))
                for c in T["channels"].values() if "per_photon" in c and c["type"] != "Sci"}
        plots.append(plot_fiber_summary(vals, "N #times d#sigma^{2}/dX [ps^{2}/mm]", -500, 200, "per_photon_testfibres", pmo, "test fibres, per-photon slope"))
        byl = {c["label"]: ch for ch, c in T["channels"].items()}
        for lab in ("Multiclad_0", "Quartz700_0", "Quartz_0"):
            if lab in byl: plots.append(sigma_vs_x(T, byl[lab], f"sigma_x_{lab}"))
    for n in ("2026A",):
        if n not in D: continue
        T = D[n]
        for t in ("CerPlastic", "CerQuartz"):
            v = [(ch, c) for ch, c in T["channels"].items() if c["type"] == t and "per_photon" in c]
            if v: plots.append(sigma_vs_x(T, min(v, key=lambda kv: kv[1]["per_photon_err"])[0], f"sigma_x_{n}_{t}"))
    for n, T in D.items():
        c = T["control"]; X = np.array(T["x_m"]) * 1000; s = np.array(c["sigma"], float) * 1000; e = np.array(c["sigma_err"], float) * 1000; ok = np.isfinite(s)
        if ok.sum() < 2 or n == "testfibres": continue
        DrawHistos([graph(X[ok], s[ok], e[ok])], [f"slope {pm(c.get('slope'), c.get('slope_err'))} ps^{{2}}/mm".replace("&plusmn;", "#pm")], X.min() - 20, X.max() + 20,
                   "X [mm]", 0, 1.8 * np.nanmax(s), "Width of MCP_{US1} - MCP_{DS1} [ps]", f"control_{n}", dology=False, drawoptions="P", legendoptions="P",
                   mycolors=[ROOT.kBlack], legendPos=[0.20, 0.80, 0.80, 0.88], lumitext=f"{n}: control (fixed MCPs)", **draw)
        plots.append(f"control_{n}")
    sections.append(("Fibre dispersion: width vs position at fixed p.e.", [p for p in plots if p]))

    # ---- tables ----
    bd_tab = table(["Dataset", "Fibres", "p.e. bin", "channels (group relation to the MCP)", "total far end", "near end", "MCP", "DRS alignment",
                    "dispersion over the scan (1&sigma; range)", "other"],
                   [[r["dataset"], r["cls"], r["pe"], f"{r['n_ch']} ({r['rels']})", f"<b>{ps(r['total'])}</b>", ps(r["near"]), ps(r["mcp"]), ps(r["drs"]),
                     f"{ps(r['disp'])} ({ps(r['disp_lo'])}&ndash;{ps(r['disp_hi'])})", ps(r["other"])] for r in B])
    el_tab = table(["Period", "MCP amplitude US / DS [ADC]", "two MCPs, same chip", "MCP term used",
                    "same MCP on two boards (alignment term, median of the pairs): integer &rarr; interpolated", "boards",
                    "reference edges of two groups, per group: integer &rarr; interpolated (line fit)"],
                   [[p, " / ".join(f"{E[p]['runs'][0]['amp'].get(m, np.nan):.0f}" for m in ("MCP_US_1", "MCP_DS_1")), ps(t["mcp_pair"]), f"<b>{ps(t['mcp'])}</b>",
                     f"{ps(t['board_int'])} &rarr; <b>{ps(t['board'])}</b> ({t['n_pairs']} pairs)", ", ".join(b.replace("DRS_", "") for b in t["boards"]) if t["boards"] else "&ndash;",
                     f"{ps(t['group_int'])} &rarr; {ps(t['group'])} ({ps(np.median([r['groups']['ref_fit']['per_group'] for r in E[p]['runs']]))})"]
                    for p, t in ET.items()])
    pb_rows = []
    for p, v in E.items():
        for r in v["runs"]:
            for vv, lab in (("ref_int", "integer"), ("ref_fw", "interpolated")):
                if vv in r.get("per_board", {}):
                    pb_rows.append([p, r["run"], lab, *[f"{b.replace('DRS_', '')}: {x:.0f}" for b, x in r["per_board"][vv].items()]])
            pb_rows.append([p, r["run"], "same-chip MCP pair", *[f"{b.replace('DRS_', '')}: {x['sigma']:.0f}" for b, x in r.get("same_chip_boards", {}).items()]])
    pb_tab = table(["Period", "run", "", "per board [ps]", "", "", ""], pb_rows)
    ch_rows = []
    for p, v in E.items():
        for r in v["runs"]:
            for rel in ("same group", "same board", "other board"):
                if rel in r["channels"]:
                    c_ = r["channels"][rel]; ch_rows.append([p, r["run"], rel, c_["n"], ps(c_["int"]), ps(c_["interp"]), ps(c_["gain"])])
    ch_tab = table(["Period", "run", "channel group vs the MCP", "channels", "width, integer reference", "interpolated", "difference in quadrature"], ch_rows)
    rows = []
    for n, T in D.items():
        c = T["control"]
        for t, v in T["types"].items():
            rows.append([n, TLAB.get(t, t), v["n"], f"<b>{pm(v['per_photon'], v['err'])}</b>", f"{v['chi2']:.1f}",
                         *[pm(b["value"], b["err"]) if b else "&ndash;" for b in v["per_bin"]], pm(c.get("slope"), c.get("slope_err"))])
    bins = [f"{lo + 0.5:.0f}" + (f"&ndash;{hi - 0.5:.0f}" if hi - lo > 1.5 and hi < 1e8 else ("+" if hi > 1e8 else "")) + " p.e." for lo, hi in next(iter(D.values()))["pe_bins"]]
    disp_tab = table(["Dataset", "Fibres", "channels", "per-photon slope [ps<sup>2</sup>/mm]", "&chi;<sup>2</sup>/ndf between channels", *bins, "MCP control"], rows)
    tf_tab = ""
    if "testfibres" in D:
        byt = {}
        for c in D["testfibres"]["channels"].values():
            if "per_photon" in c and c["type"] != "Sci": byt.setdefault(fiber_base_type(c["label"]), []).append(c)
        rows = []
        for t, cs in sorted(byt.items()):
            w = np.array([1 / c["per_photon_err"] ** 2 for c in cs]); v = np.array([c["per_photon"] for c in cs]); m = np.sum(w * v) / w.sum()
            chi2 = np.sum(w * (v - m) ** 2) / (len(v) - 1) if len(v) > 1 else 1.0; e = np.sqrt(max(chi2, 1) / w.sum())
            rows.append([t, len(cs), f"<b>{pm(m, e)}</b>", " / ".join(f"{K * c['per_photon']:+.0f}" for c in cs),
                         *[f"{np.sqrt(abs(m) * 1e6 / (2 * dd)):.0f}" if m < 0 else "&ndash;" for dd in (1.0, 2.0)]])
        tf_tab = table(["Fibre", "bundles", "per-photon slope [ps<sup>2</sup>/mm]", "per bundle", "k [ps/m] if d&#772; = 1 m", "k if d&#772; = 2 m"], rows)

    et26 = ET.get("2026B", next(iter(ET.values())))
    intro = f"""
<p>Breakdown of the DRS timing resolution (t<sub>CFD</sub> of a channel &minus; t<sub>CFD</sub> of the reference MCP_DS_1) into its parts, from
dedicated measurements; all terms add in quadrature. Study: <code>studies/timing</code>. Related: <a href="../ZScan/ZScan_overview.html">z scans</a>.</p>
<ul style="margin:5px 0;padding-left:20px">
<li><b>MCP reference</b>: two MCPs digitized on the same DRS chip (no group alignment) differ by {ps(ET['2025']['mcp_pair']) if '2025' in ET else '&ndash;'} ps
in 2025, where both give ~125 ADC: about {ps(et26['mcp'])} ps per MCP including its digitization. This value is used for all periods: in 2026 the
upstream MCP gives about half the amplitude of the reference (~80 against ~160 ADC) and dominates the pair ({ps(et26['mcp_pair'])} ps), so the reference MCP cannot be separated from it.</li>
<li><b>DRS alignment</b>: each channel is aligned to the MCP through the reference signal (Channel8) of its own group and of the MCP's group.
The two digitized copies of the same MCP signal on two boards measure what this alignment adds to real pulses, without any MCP jitter:
{ps(et26['board_int'])} ps with the former integer reference time, <b>{ps(et26['board'])} ps</b> with the interpolated crossing now used by the framework
(the integer sample added 0.2&nbsp;ns/&radic;12 = 58&nbsp;ps per group). Channels in the MCP's own group have no alignment term. In 2025 each MCP was digitized on all four boards: all 12 pairs of copies give
the same width ({ps(ET['2025']['board_int']) if '2025' in ET else '&ndash;'} &rarr; {ps(ET['2025']['board']) if '2025' in ET else '&ndash;'} ps), i.e. about the same term for every board, so
the alignment does not depend on which boards are combined; in 2026 it is larger. All copies are in Group3, so the MCPs do not separate
two groups of one board from two boards; channels of the same type in another group of the MCP's board and on other boards have similar
widths. The rest is
mostly the non-uniform DRS4 cell widths: the group-to-group difference of the reference edges follows the cell of the edge (plot
&ldquo;cells&rdquo;), and a per-cell time calibration could remove most of it. (The reference edges of two groups compared directly scatter more,
{ps(et26['group'])} ps per group, because a pulse and the reference of its own group share part of the chip's cell-timing error, which cancels
in a channel time.)</li>
<li><b>Fibre dispersion</b>: at fixed number of photoelectrons N the width grows with the distance to the readout, with
N&nbsp;&times;&nbsp;d&sigma;<sup>2</sup>/dX the same in every p.e. bin (dispersion acting on each photon). The term in the breakdown is what it adds between the
near and the far end of the scan, &radic;(|slope|&nbsp;&times;&nbsp;scan length / N); the dispersion up to the near end is in &ldquo;other&rdquo;.</li>
<li><b>Other</b>: what remains of the measured width at the near end of the scan after the MCP and the DRS alignment: SiPM single-photon
jitter, CFD and noise, the spread of the shower and of the beam spot along the fibres, the photostatistics of the light, the dispersion
up to the near end.</li>
</ul>
<h3 style="margin:14px 0 4px 0">Breakdown [ps], median over the channels of each class</h3>
{bd_tab}
<p>Test fibres: Multiclad, Singleclad, Quartz700 disperse, Quartz and Quartz400 do not; at one photoelectron the dispersion over the 84&nbsp;cm scan is
comparable to everything else, at 3&ndash;4 p.e. it is diluted by 1/N. Calorimeter channels: the electronics (MCP and alignment) is a sizeable part
of the width of the multi-p.e. hits, the dispersion over their 33&nbsp;cm scans a small one.</p>
<h3 style="margin:14px 0 4px 0">Electronics [ps]</h3>
{el_tab}
{details("Per board: alignment term from the copy pairs (sigma<sup>2</sup>(i &minus; j) = a<sub>i</sub><sup>2</sup> + a<sub>j</sub><sup>2</sup>, 2025: four copies of each MCP) and the same-chip MCP pair on every board", pb_tab)}
{details("Channel widths with the integer and the interpolated reference, by group relation to the MCP", ch_tab)}
<h3 style="margin:14px 0 4px 0">Fibre dispersion: per-photon slope</h3>
{disp_tab}
{tf_tab}
<p>Fits per channel and p.e. bin: &sigma;<sup>2</sup> = a<sub>stretch</sub> + s&middot;X; errors scaled by &radic;(&chi;<sup>2</sup>/ndf) where above 1, type averages again by
the scatter between channels. Hits within &plusmn;3&nbsp;ns (Cherenkov) / &plusmn;8&nbsp;ns (scintillating) of the arrival peak, random hits of the same p.e.
bin from the sidebands removed. Scintillating channels are listed for completeness only (decay-time dominated, no 1/N scaling). To turn
the slope into k (ps per m of fibre) the mean distance d&#772; from the beam spot to the SiPMs is needed: k = &radic;(|slope|/2d&#772;).</p>
"""
    pl = sum(([f"{p}.png" for p in ps_] + ["NEWLINE"] for _, ps_ in sections if ps_), [])[:-1]
    generate_html(pl, OUTD, plots_per_row=4, output_html=OUT, title="Timing resolution: breakdown", intro_text=intro)
    html = open(OUT).read(); parts = html.split('    <div class="grid">\n')
    html = parts[0] + "".join(f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{t}</h3>\n    <div class="grid">\n' + part
                              for (t, _), part in zip([s for s in sections if s[1]], parts[1:]))
    open(OUT, "w").write(html)
    print(OUT)


if __name__ == "__main__":
    main()
