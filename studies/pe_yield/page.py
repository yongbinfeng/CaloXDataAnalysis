"""Page of the p.e. yield study (gain.json, yields.json): gain from the single-p.e. structure, p.e. per GeV per fibre type vs energy,
event selection and the runs left out.
Usage (repository root): python3 studies/pe_yield/page.py      Page: results/html/PEYield/PEYield_2025.html"""
import json, os, sys
sys.path.insert(0, ".")
import numpy as np
import ROOT
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html
from config import DATA, PLOTS, HTML, PED_RUN, GAIN_RUN, GAIN_BOARDS, MAP_RUN, SEL_MIN_FRACTION, FAR

ROOT.gROOT.SetBatch(True)
OUTD = f"{PLOTS}/2025"; os.makedirs(OUTD, exist_ok=True); OUT = f"{HTML}/PEYield_2025.html"
G = json.load(open(f"{DATA}/gain.json")); Y = json.load(open(f"{DATA}/yields.json"))
GAIN, SPREAD = G["summary"]["gain"], G["summary"]["spread"]
LAB = {"Sci": "scintillating", "CerQ": "quartz", "CerP": "plastic", "Cer": "Cherenkov (quartz + plastic)"}
COL = {"Sci": ROOT.kGreen + 2, "CerQ": ROOT.kBlue + 1, "CerP": ROOT.kRed + 1, "Cer": ROOT.kBlack}
MK = {"Sci": 20, "CerQ": 21, "CerP": 22, "Cer": 24}
HEADLINE = [y for y in Y["runs"] if y["energy"] >= 20]     # 10 GeV: channel offsets matter (plastic sum pulled down)
draw = dict(outdir=OUTD, usePDF=False)


def hist(name, edges, counts):
    h = ROOT.TH1D(name, "", len(edges) - 1, np.asarray(edges, float))
    for i, c in enumerate(counts): h.SetBinContent(i + 1, float(c))
    return h


def graph(xs, ys, es=None):
    g = ROOT.TGraphErrors(len(xs))
    for i, (x, y) in enumerate(zip(xs, ys)):
        g.SetPoint(i, float(x), float(y)); g.SetPointError(i, 0, float(es[i]) if es is not None else 0)
    return g


TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
details = lambda title, body: f'<details style="margin:10px 0"><summary><b>{title}</b></summary>{body}</details>'

sections = []

# ---- gain ----
P = np.load(f"{DATA}/raw/raw_{PED_RUN}.npz"); meta_raw = P["meta"].astype(int)
ref = [c for c in G["channels"] if c["board"] in GAIN_BOARDS and c["good"]]
other = [c for c in G["channels"] if c["board"] not in GAIN_BOARDS]
ex_ref = [min((c for c in ref if c["type"] == t), key=lambda c: abs(c["light"] - 150)) for t in ("Sci", "CerP")]
ex_other = [max((c for c in other if c["board"] == b and c["chi2ndf"] < 1.5), key=lambda c: c["light"]) for b in (5, 13)]
PB = {int(b): v for b, v in G["per_board"].items()}
plots = []
# pedestal-run spectra: a reference-board channel vs channels of other boards
hs, ls = [], []
for c in ex_ref[:1] + ex_other:
    i = int(np.where((meta_raw[:, 0] == c["board"]) & (meta_raw[:, 1] == c["channel"]))[0][0]); v = P["hg"][:, i].astype(float)
    e = np.arange(-80, 401, 2.0); h = hist(f"ped_{i}", e, np.histogram(v - np.median(v), e)[0]); h.Scale(1 / h.Integral())
    hs.append(h); ls.append(f"board {c['board']} ch {c['channel']} ({LAB[c['type']]})")
DrawHistos(hs, ls, -80, 400, "HG - median [ADC]", 1e-6, 1, "Fraction of events", "ped_spectra", dology=True, drawoptions="HIST",
           mycolors=[ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1], legendPos=[0.40, 0.72, 0.92, 0.88], lumitext=f"pedestal run {PED_RUN}", **draw)
plots.append("ped_spectra")
# example fits and chi2 profiles
edges = np.array(G["edges"])
for c in ex_ref + ex_other:
    nm = f"fit_B{c['board']}_{c['channel']}"
    h = hist(f"h_{nm}", edges, c["hist"]); m = hist(f"m_{nm}", edges, c["model"])
    txt = ROOT.TLatex(0.22, 0.80, f"g = {c['gain']:.1f} #pm {c['err']:.1f} ADC/p.e., #chi^{{2}}/ndf {c['chi2ndf']:.2f}"); txt.SetNDC(); txt.SetTextSize(0.035)
    DrawHistos([h, m], ["data", "fit"], edges[0], edges[-1], "HG - pedestal-run median [ADC]", 0.5, 60 * max(c["hist"]), "Events / ADC", nm,
               dology=True, drawoptions=["HIST", "HIST"], mycolors=[ROOT.kBlack, ROOT.kRed + 1], legendPos=[0.70, 0.74, 0.92, 0.88],
               extraToDraw=[txt], lumitext=f"run {GAIN_RUN}, board {c['board']} ch {c['channel']} ({LAB[c['type']]})", **draw)
    plots.append(nm)
bl = sorted(PB); cols = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 1, ROOT.kCyan + 2, ROOT.kGray + 2]
gs = [graph(G["gains"], np.clip((np.array(PB[b]["profile_sum"]) - min(PB[b]["profile_sum"])) / PB[b]["n"], 0, 30)) for b in bl]
DrawHistos(gs, [f"board {b} ({PB[b]['n']} ch)" for b in bl], G["gains"][0], G["gains"][-1], "Gain [ADC/p.e.]", 0, 20,
           "#Delta#chi^{2} summed / channels", "profiles", dology=False, drawoptions="L", legendoptions="L", mycolors=cols[:len(bl)],
           legendPos=[0.20, 0.60, 0.50, 0.90], lumitext=f"run {GAIN_RUN}, summed channel profiles", **draw)
plots.append("profiles")
gb = ROOT.TGraphErrors(len(bl))
for i, b in enumerate(bl): gb.SetPoint(i, b, PB[b]["gain"]); gb.SetPointError(i, 0, PB[b]["spread"])
other_g = Y["unfitted_boards"]; lo = ROOT.TLine(1.5, other_g, 15.5, other_g); lo.SetLineStyle(2); lo.SetLineColor(ROOT.kGray + 2)
DrawHistos([gb], ["channel median, bar: spread"], 1.5, 15.5, "FERS board", 40, 100, "Gain [ADC/p.e.]", "gain_boards", dology=False,
           drawoptions="P", legendoptions="P", mycolors=[ROOT.kBlack], legendPos=[0.20, 0.80, 0.75, 0.88], extraToDraw=[lo],
           lumitext=f"run {GAIN_RUN}; dashed: boards without a fit", **draw)
plots.append("gain_boards")
e = np.arange(60, 95.1, 1.0)
hs = [hist(f"g_{t}", e, np.histogram([c["gain"] for c in ref if c["type"] == t], e)[0]) for t in ("Sci", "CerP", "CerQ")]
line = ROOT.TLine(GAIN, 0, GAIN, max(h.GetMaximum() for h in hs) * 1.3); line.SetLineStyle(2)
DrawHistos(hs, [f"{LAB[t]} ({int(h.Integral())})" for t, h in zip(("Sci", "CerP", "CerQ"), hs)], e[0], e[-1], "Gain [ADC/p.e.]", 0,
           max(h.GetMaximum() for h in hs) * 1.5, "Channels", "gain_hist", dology=False, drawoptions="HIST",
           mycolors=[COL["Sci"], COL["CerP"], COL["CerQ"]], legendPos=[0.20, 0.72, 0.55, 0.88], extraToDraw=[line],
           lumitext=f"boards {', '.join(map(str, GAIN_BOARDS))}, good fits", **draw)
plots.append("gain_hist")
sections.append(("Gain from the single-p.e. peaks", plots))

# ---- yields vs energy ----
plots = []
E = [y["energy"] for y in Y["runs"]]
for name, ts, ymax, lp in (("pe_sci", ("Sci",), 80, [0.20, 0.80, 0.60, 0.88]), ("pe_cer", ("CerQ", "CerP", "Cer"), 18, [0.20, 0.72, 0.70, 0.88])):
    gs = [graph(E, [y["types"][t]["pe_per_gev"] for y in Y["runs"]], [y["types"][t].get("stat", 0) for y in Y["runs"]]) for t in ts]
    DrawHistos(gs, [LAB[t] for t in ts], 0, 110, "Beam energy [GeV]", 0, ymax, "p.e. per GeV", name, dology=False, drawoptions="P",
               legendoptions="P", mycolors=[COL[t] for t in ts], markerstyles=[MK[t] for t in ts], legendPos=lp,
               lumitext="2025 e^{+}, gain per board", **draw)
    plots.append(name)
bins = np.array(Y["sce_bins"]); hs, ls = [], []
for y in Y["runs"][-3:] + Y["bad_runs"][:3]:
    h = hist(f"sce_{y['run']}", bins, y["sce_hist"]); h.Scale(1 / max(h.Integral(), 1)); hs.append(h); ls.append(f"{y['run']}, {y['energy']:.0f} GeV")
cut = ROOT.TLine(SEL_MIN_FRACTION, 0, SEL_MIN_FRACTION, 0.06); cut.SetLineStyle(2)
DrawHistos(hs, ls, -0.2, 3.0, "Calibrated (S+C) / E", 0, 0.08, "Fraction of events", "sce", dology=False, drawoptions="HIST",
           mycolors=[ROOT.kGreen + 2, ROOT.kBlue + 1, ROOT.kBlack, ROOT.kRed + 1, ROOT.kOrange + 1, ROOT.kMagenta + 1],
           linestyles=[1, 1, 1, 2, 2, 2], legendPos=[0.55, 0.60, 0.92, 0.88], extraToDraw=[cut], lumitext="positron PID, hodoscope 1-3", **draw)
plots.append("sce")
# light map of the map run: p.e. per event per channel, by tower
mp = next(y for y in Y["runs"] if y["run"] == MAP_RUN)["map"]; mm = np.array(mp["meta"]); pe = np.array(mp["pe_per_event"])
for name, sel, zt in (("map_sci", mm[:, 2] == 0, "scintillating"), ("map_cer", mm[:, 2] == 1, "Cherenkov")):
    h2 = ROOT.TH2D(f"h_{name}", "", 27, -13.5, 13.5, 19, -9.5, 9.5)
    for (ix, iy), v in zip(mm[sel][:, 4:6], pe[sel]):
        if np.isfinite(v): h2.Fill(ix, iy, max(v, 1e-3))
    DrawHistos([h2], [], -13.5, 13.5, "Tower X", -9.5, 9.5, "Tower Y", name, dology=False, doth2=True, drawoptions=["COLZ"], dologz=True,
               zmin=0.1, zmax=max(1.0, h2.GetMaximum()), zlabel="p.e. per event", lumitext=f"run {MAP_RUN}, {zt}", W_ref=800, H_ref=600, **draw)
    plots.append(name)
sections.append(("p.e. per GeV vs energy, event selection, light map", plots))

# ---- energy resolution (resolution.py) ----
RS = json.load(open(f"{DATA}/resolution.json")) if os.path.exists(f"{DATA}/resolution.json") else None
res_html = ""
if RS:
    plots = []; RR = RS["runs"]; En = [r["energy"] for r in RR]; FT = RS["fits"]; NZ = RS["noise_gev"]; PS = RS["photostat"]
    SLAB = {"S": "S", "C": "C", "SC": "(S+C)/2"}; SCOL = {"S": ROOT.kGreen + 2, "C": ROOT.kBlue + 1, "SC": ROOT.kBlack}; SMK = {"S": 20, "C": 21, "SC": 24}
    def rcurve(name, f, col, style=1, xlo=8, xhi=110):
        g = ROOT.TF1(name, f"100*sqrt({f['a']}**2/x + {max(f['b2'], 0)} + ({f['c']}/x)**2)", xlo, xhi); g.SetLineColor(col); g.SetLineStyle(style); g.SetLineWidth(2); return g
    # resolution vs E, main configuration, with the fits and the photostatistics
    hs, ls, ex = [], [], []
    for k in ("S", "C", "SC"):
        hs.append(graph(En, [100 * r["spot_reg"]["signals"][k]["res_impact"] for r in RR], [100 * r["spot_reg"]["signals"][k]["err_impact"] for r in RR]))
        f = FT[f"spot_reg_{k}_impact"]; ls.append(f"{SLAB[k]}: {100 * f['a']:.0f}%/#sqrt{{E}} #oplus {100 * f['b']:.1f}% #oplus {f['c']:.2f}/E")
        ex.append((rcurve(f"rf_{k}", f, SCOL[k]), "SAME"))
        ph = ROOT.TF1(f"ph_{k}", f"100*{PS[k]}/sqrt(x)", 8, 110); ph.SetLineColor(SCOL[k]); ph.SetLineStyle(2); ex.append((ph, "SAME"))
    DrawHistos(hs, ls, 0, 110, "Beam energy [GeV]", 0, 30, "#sigma / E [%]", "res_vs_e", dology=False, drawoptions="P", legendoptions="P",
               mycolors=[SCOL[k] for k in ("S", "C", "SC")], markerstyles=[SMK[k] for k in ("S", "C", "SC")], legendPos=[0.32, 0.68, 0.92, 0.88],
               extraToDraw=ex, lumitext="2025 e^{+}; dashed: photostatistics", **draw)
    plots.append("res_vs_e")
    # module vs region sum (noise)
    hs, ls, cols, mks, ex = [], [], [], [], []
    for k in ("S", "C"):
        for sel_, lab, mk in (("spot", "whole module", 24 if k == "C" else 20), ("spot_reg", f"{RS['region_towers']} towers", 21 if k == "C" else 22)):
            hs.append(graph(En, [100 * r[sel_]["signals"][k]["res"] for r in RR])); ls.append(f"{k}, {lab}"); cols.append(SCOL[k]); mks.append(mk)
        for sel_, st in (("", 2), ("_reg", 3)):
            nzf = ROOT.TF1(f"nz_{k}{sel_}", f"100*{NZ[k + sel_]}/x", 8, 110); nzf.SetLineColor(SCOL[k]); nzf.SetLineStyle(st); ex.append((nzf, "SAME"))
    DrawHistos(hs, ls, 0, 110, "Beam energy [GeV]", 0, 40, "#sigma / E [%]", "res_noise", dology=False, drawoptions="P", legendoptions="P", mycolors=cols,
               markerstyles=mks, legendPos=[0.45, 0.66, 0.92, 0.88], extraToDraw=ex, lumitext="no position correction; lines: noise / E", **draw)
    plots.append("res_noise")
    # response vs impact point and centre of gravity
    pr = next(r for r in RR if "profiles" in r)["profiles"]; pb = np.array(RS["profile_bins"]); pc = 0.5 * (pb[1:] + pb[:-1])
    for ax, xt in (("impact", "impact point - median [cm]"), ("cog", "shower centre - median [towers]")):
        hs, ls, cols, mks = [], [], [], []
        for k in ("S", "C"):
            for d_, mk in (("x", 20), ("y", 24)):
                v = pr[f"{k}_{ax}_{d_}"]; ok = [i for i, x in enumerate(v) if x is not None]
                hs.append(graph(pc[ok], [v[i] for i in ok])); ls.append(f"{k} vs {d_}"); cols.append(SCOL[k]); mks.append(mk)
        DrawHistos(hs, ls, pb[0], pb[-1], xt, 0.85, 1.15 if ax == "impact" else 1.35, "response / median", f"profile_{ax}", dology=False, drawoptions="LP", legendoptions="LP",
                   mycolors=cols, markerstyles=mks, legendPos=[0.20, 0.72, 0.60, 0.88], lumitext=f"run {MAP_RUN}, whole beam spot", **draw)
        plots.append(f"profile_{ax}")
    # position corrections, whole beam spot
    hs, ls, cols, mks = [], [], [], []
    for k in ("S", "C"):
        for t, lab, mk in (("res", "no correction", 24), ("res_cog", "centre of gravity", 20), ("res_impact", "impact point", 21)):
            hs.append(graph(En, [100 * r["all_reg"]["signals"][k][t] for r in RR])); ls.append(f"{k}, {lab}"); cols.append(SCOL[k]); mks.append(mk)
    DrawHistos(hs, ls, 0, 110, "Beam energy [GeV]", 0, 30, "#sigma / E [%]", "res_corrections", dology=False, drawoptions="P", legendoptions="P", mycolors=cols,
               markerstyles=mks, legendPos=[0.42, 0.62, 0.92, 0.88], lumitext=f"whole beam spot, {RS['region_towers']}-tower sums", **draw)
    plots.append("res_corrections")
    sections.append(("Energy resolution (positrons)", plots))
    rows = [[r["run"], f"{r['energy']:.0f}", r["spot_reg"]["n"], *[f"{100 * r['spot_reg']['signals'][k]['res_impact']:.2f}" for k in ("S", "C", "SC")],
             f"{100 * r['spot_reg']['s_minus_c']:.2f}", f"{r['spot_reg']['corr_sc']:+.2f}",
             *[f"{100 * r['spot']['signals'][k]['res']:.2f}" for k in ("S", "C")], *[f"{100 * r['all_reg']['signals'][k][t]:.2f}" for k in ("S", "C") for t in ("res", "res_cog")]]
            for r in RR]
    res_tab = table(["Run", "E [GeV]", "events", "S", "C", "(S+C)/2", "S &minus; C", "corr(S, C)", "S, module", "C, module",
                     "S, spot", "S, spot + CoG", "C, spot", "C, spot + CoG"], rows)
    frows = []
    for k in ("S", "C", "SC"):
        f = FT[f"spot_reg_{k}_impact"]
        frows.append([SLAB[k], f"{100 * f['a']:.1f} &plusmn; {100 * f['a_err']:.1f}", f"{100 * f['b']:.1f} (&lt; {100 * f['b_upper']:.1f})", f"{f['c']:.2f}",
                      f"{f['chi2ndf']:.1f}", f"{100 * PS[k]:.1f}", f"{100 * (PS[k] / f['a']) ** 2:.0f}%"])
    fsc = FT["spot_reg_s_minus_c"]
    frows.append(["S &minus; C", f"{100 * fsc['a']:.1f} &plusmn; {100 * fsc['a_err']:.1f}", f"{100 * fsc['b']:.1f}", "&ndash;", f"{fsc['chi2ndf']:.1f}",
                  f"{100 * np.sqrt(PS['S'] ** 2 + PS['C'] ** 2):.1f}", f"{100 * (PS['S'] ** 2 + PS['C'] ** 2) / fsc['a'] ** 2:.0f}%"])
    fit_tab = table(["Signal", "a [%/&radic;GeV]", "b [%]", "c [GeV] (pedestal run, fixed)", "&chi;<sup>2</sup>/ndf", "photostatistics [%/&radic;GeV]",
                     "photostatistics share of a<sup>2</sup>"], frows)
    res_html = f"""
<h3 style="margin:14px 0 4px 0">Energy resolution</h3>
<p>Positrons of the same runs, beam spot within &plusmn;{RS['beam_spot_cm']:.0f} cm of the run median, sums over the towers within {RS['region_towers']} towers
of the shower, corrected for the impact point (quadratic surface fitted on the even events, applied to the odd ones); width of the gaussian
core. Fits over {RS['fit_emin']}&ndash;100 GeV: &sigma;/E = a/&radic;E &oplus; b &oplus; c/E, with the noise c measured in the pedestal run
(widths of the calibrated sums: whole module S {NZ['S']:.2f} / C {NZ['C']:.2f} GeV, {RS['region_towers']} towers S {NZ['S_reg']:.2f} / C {NZ['C_reg']:.2f} GeV; the channels are correlated,
without correlations the whole module would give S {NZ['S_uncorrelated']:.2f} / C {NZ['C_uncorrelated']:.2f} GeV). Photostatistics: 1/&radic;(p.e. per GeV).</p>
{fit_tab}
<ul style="margin:5px 0;padding-left:20px">
<li><b>Photostatistics</b> is a minor part of the stochastic term: about a sixth of it (in a<sup>2</sup>) for S and a third for C; the rest is sampling.</li>
<li><b>Noise</b>: summed over the whole module the pedestal and dark counts of the ~450 Cherenkov channels give {NZ['C']:.1f} GeV, 25% at 10 GeV; summing
only near the shower halves it without losing signal.</li>
<li><b>S &minus; C</b>: the fluctuations common to S and C (correlation 0.3&ndash;0.8) cancel; what remains, beyond the photostatistics of both, is the
independent sampling of the two sets of fibres.</li>
<li><b>Position corrections</b>: on the whole beam spot the response depends on the impact point, mostly for C (quartz and plastic towers); the
centre of gravity of the shower corrects it as well as the hodoscope impact point (C at 100 GeV: {100 * RR[-1]['all_reg']['signals']['C']['res']:.1f}% &rarr;
{100 * RR[-1]['all_reg']['signals']['C']['res_cog']:.1f}%).</li>
</ul>
{details("Resolution per run [%]", res_tab)}
"""

# ---- tables ----
f1 = lambda v: f"{v:.1f}"; f2 = lambda v: f"{v:.2f}"
rows = []
for y in Y["runs"]:
    t = y["types"]
    rows.append([y["run"], f"{y['energy']:.0f}", y["n_sel"], f1(t["Sci"]["pe_per_gev"]), f2(t["CerQ"]["pe_per_gev"]), f2(t["CerP"]["pe_per_gev"]),
                 f2(t["Cer"]["pe_per_gev"]), f1(t["Sci"]["pe_per_gev"] / t["Cer"]["pe_per_gev"]),
                 " / ".join(f"{t[k]['share_ref']:.2f}" for k in ("Sci", "CerQ", "CerP")), " / ".join(f"{t[k]['far']:+.2f}" for k in ("Sci", "CerQ", "CerP")),
                 f"{t['Sci']['pe_per_gev_common']:.1f} / {t['CerQ']['pe_per_gev_common']:.2f} / {t['CerP']['pe_per_gev_common']:.2f}"])
avg = {t: np.mean([y["types"][t]["pe_per_gev"] for y in HEADLINE]) for t in ("Sci", "CerQ", "CerP", "Cer")}
rms = {t: np.std([y["types"][t]["pe_per_gev"] for y in HEADLINE]) for t in avg}
rows.append([f"<b>mean {min(y['energy'] for y in HEADLINE):.0f}&ndash;{max(y['energy'] for y in HEADLINE):.0f} GeV</b>", "", "",
             *[f"<b>{avg[t]:.{1 if t == 'Sci' else 2}f}</b> &plusmn; {rms[t]:.{1 if t == 'Sci' else 2}f}" for t in ("Sci", "CerQ", "CerP", "Cer")],
             f"<b>{avg['Sci'] / avg['Cer']:.1f}</b>", "", "", ""])
yield_tab = table(["Run", "E [GeV]", "events", "scintillating", "quartz", "plastic", "Cherenkov", "S/C",
                   f"on boards {'+'.join(map(str, GAIN_BOARDS))} (S/Q/P)", f"beyond {FAR} towers (S/Q/P)", f"one gain {GAIN:.0f} for all boards (S/Q/P)"], rows)
gt = G["summary"]["per_type"]
gain_tab = table(["Fibre type", "channels (good fits)", "median gain [ADC/p.e.]"], [[LAB[t], gt[t]["n"], f1(gt[t]["median"])] for t in gt])
board_tab = table(["Board", "peaks", "channels fitted", "gain [ADC/p.e.]", "&plusmn; (median)", "channel spread", "summed-profile minimum"],
                  [[b, "resolved" if PB[b]["resolved"] else "smeared", PB[b]["n"], f1(PB[b]["gain"]), f1(PB[b]["err"]), f1(PB[b]["spread"]),
                    f1(PB[b]["profile_min"])] for b in bl]
                  + [[", ".join(str(b) for b in range(2, 16) if b not in PB), "", "too few lit channels", f"{other_g:.1f}", "", "", "median of the smeared boards"]])
bad_tab = table(["Run", "E [GeV]", "PID events", "median (S+C)/E"], [[y["run"], f"{y['energy']:.0f}", y["n_pid"], f2(y["sce_median"])] for y in Y["bad_runs"]]
                + [[f"<i>{y['run']}</i>", f"{y['energy']:.0f}", y["n_pid"], f2(y["sce_median"])] for y in Y["runs"]])

intro = f"""
<p>Photoelectrons per GeV of beam energy in the 2025 FERS readout (6 mm SiPMs), for positrons at position #1 (table X &asymp; 0,
Y = 1032 mm, 3.3&deg;): the Mix ADC of each channel (pedestal-subtracted HG; LG scaled to HG where HG saturates) divided by the gain,
summed over the channels of a fibre type and averaged over the events. No correction for leakage. Study: <code>studies/pe_yield</code>.</p>
<h3 style="margin:14px 0 4px 0">p.e. per GeV</h3>
{yield_tab}
<p>Flat from 20 to 100 GeV: about <b>{avg['Sci']:.0f} p.e./GeV scintillating</b>, <b>{avg['CerQ']:.1f} quartz</b> and <b>{avg['CerP']:.1f} plastic</b>
(<b>{avg['Cer']:.0f} Cherenkov</b>), S/C &asymp; {avg['Sci'] / avg['Cer']:.1f}. At 10 GeV the small offsets of the many channels away from
the shower matter (plastic: boards {'+'.join(map(str, GAIN_BOARDS))} alone give more than the total). Uncertainties beyond the statistics:</p>
<ul style="margin:5px 0;padding-left:20px">
<li><b>Gain</b> per FERS board (table below): {GAIN:.0f} ADC/p.e. on boards {', '.join(map(str, GAIN_BOARDS))} where the single-p.e. peaks are resolved
(channel spread &plusmn;{SPREAD:.0f}), 53&ndash;61 on the other fitted boards (shallow profiles, &plusmn;2&ndash;3 per board), {Y['unfitted_boards']:.0f} assumed for
the boards with too few lit channels at {next(y['energy'] for y in Y['runs'] if y['run'] == GAIN_RUN):.0f} GeV. With one gain ({GAIN:.0f}) for all boards
(last column) the scintillating and quartz yields are 4% and 12% lower: boards {'+'.join(map(str, GAIN_BOARDS))} carry most of the scintillating and plastic
light but only about 2/3 of the quartz light (column &ldquo;on boards&rdquo;).</li>
<li><b>Crosstalk</b>: ADC / gain counts fired SiPM cells, including optical crosstalk; primary photoelectrons are fewer, by the crosstalk
probability of the SiPMs (typically 5&ndash;25%).</li>
<li><b>Offsets</b>: the light beyond {FAR} towers from the shower centre is 1&ndash;3% above 10 GeV, so channel offsets of a few ADC matter little.</li>
</ul>
<h3 style="margin:14px 0 4px 0">Gain</h3>
<p>Fit of the HG spectrum of each channel with a little light in run {GAIN_RUN} ({next(y['energy'] for y in Y['runs'] if y['run'] == GAIN_RUN):.0f} GeV):
the channel's own spectrum in pedestal run {PED_RUN} (pedestal + dark counts) is the zero-p.e. shape, shifted by n &times; gain with free
non-negative weights per n; the gain is scanned and the &chi;<sup>2</sup> profiled over a common offset and a per-p.e. smearing. All 2025
FERS channels are 6 mm SiPMs. On boards {', '.join(map(str, GAIN_BOARDS))} (about 5&times; lower SiPM current than the others at the same temperature)
the single-p.e. peaks are resolved and every channel gives a sharp minimum (good fit: profile width below 2 ADC, &chi;<sup>2</sup>/ndf below
1.5). On the other boards the dark counts, which the FERS hold samples at random times, smear the peaks: single channels give shallow
profiles, so the gain is taken per board as the median of the channel minima (error from their spread); the minimum of the summed
profiles agrees. Those boards come out about 25% lower in gain than boards {'+'.join(map(str, GAIN_BOARDS))}, consistent with a different SiPM batch.</p>
{board_tab}
{details("Boards " + "+".join(map(str, GAIN_BOARDS)) + " by fibre type", gain_tab)}
{details("Runs left out: the energy scan of the night Sep 25&ndash;26", "<p>In runs 1409&ndash;1415 the calibrated (S+C)/E is 0.3&ndash;0.5 with a tail to zero, against &asymp; 2 (S and C each calibrated to E) in all other runs, and the central channels sit at the pedestal in most events: the FERS signal is not right in these runs (sampling off the pulse peak?). The logbook does not say. Runs used in italics.</p>" + bad_tab)}
{res_html}
{details("Selection", f"<p>Positron PID (beam Cherenkov counters, muon veto), 1&ndash;3 hodoscope fibres per axis, calibrated (S+C)/E above {SEL_MIN_FRACTION} (removes the few events with a partial shower).</p>")}
"""
pl = sum(([f"{n}.png" for n in ns] + ["NEWLINE"] for _, ns in sections), [])[:-1]
generate_html(pl, OUTD, plots_per_row=4, output_html=OUT, title="p.e. per GeV, 2025 FERS (positrons)", intro_text=intro)
html = open(OUT).read(); parts = html.split('    <div class="grid">\n')
html = parts[0] + "".join(f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{t}</h3>\n    <div class="grid">\n' + part for (t, _), part in zip(sections, parts[1:]))
open(OUT, "w").write(html)
print(OUT)
