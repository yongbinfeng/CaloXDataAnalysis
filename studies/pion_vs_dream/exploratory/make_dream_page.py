"""HTML page: CaloX tb2025 pion response and resolution vs the original DREAM paper."""
import json
import os
import sys

REPO = REPO_DIR
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
SCR = DATA_DIR
sys.path.insert(0, REPO)
os.chdir(REPO)
from utils.html_generator import generate_html

PLOTDIR = "results/plots/FERS_PionResponse"
OUT = "results/html/FERS_PionResponse/PionResponse_vs_DREAM.html"

resp = json.load(open(f"{PLOTDIR}/pion_response.json"))
res = json.load(open(f"{SCR}/calox_res_holeveto.json"))
d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json"))
d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
dS = dict(zip(d11["energy_GeV"], d11["S_over_E"]))
dQ = dict(zip(d11["energy_GeV"], d11["QS_over_E"]))
dres = {e: (s, c, q) for e, s, c, q in zip(d32["energy_GeV"], d32["S_res_pct"], d32["C_res_pct"], d32["QS_res_pct"])}
NO_TARGET = ["1441", "1439", "1438", "1437", "1429", "1432", "1433", "1434"]

TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"
TD = "padding:3px 10px;text-align:right"


def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


resp_rows = []
for run in NO_TARGET:
    r = resp[run]; E = int(r["E"])
    resp_rows.append([run, E, f"{r['N']}", f"{r['S'][0]:.3f}", f"{r['C'][0]:.3f}", f"{r['C'][0] / r['S'][0]:.3f}",
                      f"{r['DR'][0]:.3f}", f"{r['QS'][0]:.3f}",
                      f"{dS[E]:.3f}" if E in dS else "&ndash;",
                      "0.640" if E == 100 else "&ndash;", "0.78" if E == 100 else "&ndash;",
                      f"{dQ[E]:.3f}" if E in dQ else "&ndash;"])
resp_tab = table(["Run", "E [GeV]", "pions", "CaloX S/E", "CaloX C/E", "CaloX C/S", "CaloX DR/E", "CaloX Q/S-corr./E",
                  "DREAM S/E", "DREAM C/E", "DREAM C/S", "DREAM Q/S-corr./E"], resp_rows)

res_rows = []
for run in NO_TARGET:
    r = res[run]; E = int(r["E"])
    dd = dres.get(E)
    dcells = [f"{v:.1f}" for v in dd] if dd else ["&ndash;"] * 3
    res_rows.append([run, E, f"{r['N']}", f"{r['S_rms']:.1f}", f"{r['C_rms']:.1f}", f"{r['DR_rms']:.1f}",
                     f"{r['QS_rms']:.1f}", *dcells])
res_tab = table(["Run", "E [GeV]", "events", "CaloX S [%]", "CaloX C [%]", "CaloX DR [%]", "CaloX Q/S [%]",
                 "DREAM S [%]", "DREAM C [%]", "DREAM Q/S [%]"], res_rows)

# ---- sections 6-9: attenuation, depth reference, electron scale, corrected pion response
att = json.load(open(f"{SCR}/att90_Y1064.json"))
fin = json.load(open(f"{SCR}/pion_final.json"))
t12 = json.load(open(f"{SCR}/test12.json"))
TYPE_LAB = {"Sci": "Scintillation", "CerQuartz": "Cherenkov, quartz", "CerPlastic": "Cherenkov, plastic"}
att_rows = [[TYPE_LAB[t], f"{att['FERS'][t]['lam']:.2f} &plusmn; {att['FERS'][t]['elam']:.2f}",
             " ".join(f"{x:+.1f}" for x in att['FERS'][t]['res_pct']), f"{att['FERS'][t]['r1501_pct']:+.1f}",
             f"{att['DRS'][t]['lam']:.2f} &plusmn; {att['DRS'][t]['elam']:.2f}", f"{att['DRS'][t]['r1501_pct']:+.1f}"]
            for t in ("Sci", "CerQuartz", "CerPlastic")]
att_tab = table(["Fibres", "FERS &lambda; [m]", "FERS residuals [%] (1513, 1507, 1511, 15130)", "FERS 1501 off fit [%]",
                 "DRS &lambda; [m]", "DRS 1501 off fit [%]"], att_rows)
cent_tab = table(["Run", "Table X [mm]", "Angle", "shower centre ix [towers]", "iy [towers]", "DRS arrival-time X [mm]"],
                 [["1513", "&minus;54.5", "93.3&deg;", "&minus;2.45", "&minus;0.44", "&minus;55.5"],
                  ["1501", "&minus;168", "90&deg;", "<b>&minus;3.22</b>", "&minus;0.56", "&minus;168 (ref.)"],
                  ["1507", "&minus;218", "90&deg;", "&minus;2.43", "&minus;0.60", "&minus;218.1"],
                  ["1511", "&minus;268", "93.3&deg;", "&minus;2.45", "&minus;0.65", "&minus;262.3"],
                  ["15130", "&minus;400.3", "93.3&deg;", "&minus;2.43", "&minus;0.85", "&minus;388.2"]])
ref_tab = table(["Positron run", "E [GeV]", "Period", "median displacement dX [cm]", "dY [cm]"],
                [["1423", "20", "before the pion scan", "&minus;0.54", "&minus;0.09"], ["1416", "30", "before", "&minus;0.70", "&minus;0.12"],
                 ["1409", "80", "before", "&minus;0.58", "&minus;0.04"], ["1411", "120", "before", "&minus;0.53", "&minus;0.03"],
                 ["1514", "40", "after the 13.3&ndash;93.3&deg; rotations", "+0.04", "&minus;0.01"],
                 ["1526", "60", "after", "&minus;0.04", "&minus;0.03"], ["1527", "100", "after (map defined here)", "&minus;0.01", "&minus;0.01"]])
esc_tab = table(["E [GeV]", "Run", "e<sup>+</sup> S/E", "e<sup>+</sup> C/E", "C/S"],
                [[E, POSRUN, f"{v[0]:.3f}", f"{v[1]:.3f}", f"{v[1] / v[0]:.3f}"]
                 for (E, v), POSRUN in zip(sorted(((int(k), v) for k, v in fin["esc"].items())), (1423, 1416, 1514, 1526, 1409, 1527, 1411))])
pos = t12["pos"]
off_tab = table(["Run", "Table (X, Y) [mm]", "shower centre (x, y) [cm]", "S/E", "rel. to centre", "C/E", "rel. to centre"],
                [[r, f"({v['tx']:+.0f}, {v['ty']:.0f})", f"({v['ex']:+.1f}, {v['ey']:+.1f})", f"{v['S']:.3f}", f"{v['S'] / pos['1409']['S']:.3f}",
                  f"{v['C']:.3f}", f"{v['C'] / pos['1409']['C']:.3f}"]
                 for r, v in sorted(pos.items(), key=lambda kv: (kv[1]["ty"], kv[1]["tx"])) if r != "1399"])
def fin_rows(yd):
    rows = []
    for E, v in sorted(((int(k), v) for k, v in fin["rows"][yd].items())):
        rows.append([E, v["run"], v["N"], f"{v['fq']:.2f}", f"{v['S']:.3f} &plusmn; {v['eS']:.3f}", f"{dS[E]:.3f}" if E in dS else "&ndash;",
                     f"{v['C']:.3f} &plusmn; {v['eC']:.3f}", "0.640" if E == 100 else "&ndash;", f"{v['CS']:.3f}", "0.78" if E == 100 else "&ndash;",
                     f"{v['Q']:.3f} &plusmn; {v['eQ']:.3f}", f"{dQ[E]:.3f}" if E in dQ else "&ndash;"])
    return table(["E [GeV]", "Run", "pions", "quartz share of C", "CaloX S/e", "DREAM S/E", "CaloX C/e", "DREAM C/E", "CaloX C/S",
                  "DREAM C/S", "CaloX Q/S-corr./e", "DREAM Q/S-corr./E"], rows)
fin_tab0, fin_tab1 = fin_rows("0"), fin_rows("1")
fres = json.load(open(f"{SCR}/pion_final_res.json"))
def pair(v, k, i):
    return f"{v['final_dY0'][k][i]:.1f} / {v['final_dY1'][k][i]:.1f}"
fres_rows = []
for E, v in sorted(((int(k), v) for k, v in fres.items())):
    dd = dres.get(E)
    fres_rows.append([E, pair(v, "S", 0), pair(v, "C", 0), pair(v, "DR", 0), pair(v, "QS", 0), pair(v, "QS", 1),
                      f"{v['insample']['QS'][0]:.1f} / {v['insample']['QS'][1]:.1f}", f"{v['raw']['QS'][0]:.1f} / {v['raw']['QS'][1]:.1f}",
                      *([f"{dd[0]:.1f}", f"{dd[1]:.1f}", f"{dd[2]:.1f}"] if dd else ["&ndash;"] * 3)])
fres_tab = table(["E [GeV]", "S rms", "C rms", "DR rms", "Q/S rms", "Q/S core", "Q/S in-sample (rms / core)", "Q/S uncorr. (rms / core)",
                  "DREAM S rms", "DREAM C rms", "DREAM Q/S"], fres_rows)
L = fin["lam"]

intro = f"""
<p>CaloX pions (SPS 2025, FERS readout) compared with the original DREAM paper on the same module:
N. Akchurin et al., NIM A 537 (2005) 537 (<a href="http://www.phys.ttu.edu/~dream/resources/publications/537_537.pdf">pdf</a>).
&pi;<sup>+</sup> runs 1429&ndash;1441 (20&ndash;160 GeV, no target), table X = 0, Y = 1032 mm, 3.3&deg;. Pion selection: muon counter and
PSD vetoed, three Cherenkov counters fired. DREAM values are read off their Figs. 11, 32a and 35.</p>

<div style="border:1px solid #888;background:#f4f4f4;padding:6px 12px;margin:10px 0">
<b>Main points</b>
<ul style="margin:4px 0;padding-left:20px">
<li>Raw S looks almost compensating, so dual readout does not help (section 1).</li>
<li>Cause: deeper, hadron-rich showers lose less light to attenuation. Corrected event by event with the hodoscope,
&chi; &asymp; 0.3 (e/h<sub>S</sub> &asymp; 1.35) as in DREAM, and Q/S improves from 16.6% to 12.6% (rms) at 100 GeV (section 2).</li>
<li>Attenuation measured on 90&deg; positron runs: &lambda; = {L['Sci'][0]:.1f} m (scintillating), {L['CerQuartz'][0]:.1f} m (quartz),
{L['CerPlastic'][0]:.1f} m (plastic) (section 3).</li>
<li>Scaled to positron depth with these &lambda; and put on our own electron scale, the 100 GeV response is
S/e {fin['rows']['1']['100']['S']:.2f} (DREAM 0.83), C/e {fin['rows']['1']['100']['C']:.2f} (0.64), Q/S/e {fin['rows']['1']['100']['Q']:.2f} (1.01):
about 6&ndash;7% above DREAM, 4&ndash;6% with the hottest-channel &lambda;. Q/S core resolution 7.9% (8.2%) (section 6).</li>
<li>Open: the remaining ~6% response excess; S 9&ndash;16% above DREAM at 20&ndash;40 GeV; S and C resolutions 3&ndash;8 points worse
than DREAM.</li>
</ul></div>

<h3 style="margin:16px 0 4px 0">1. Raw response and resolution</h3>
{resp_tab}
{res_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>Response per beam energy (DREAM Fig. 11) and &sigma;<sub>rms</sub>/mean; resolution sample with the hole veto and S &gt; 0.3 E.
DR = (S &minus; 0.29 C)/0.71; Q/S = DREAM's method (Section 5.3).</li>
<li>S is high and nearly flat in C/S (e/h<sub>S</sub> &asymp; 1.0 from DREAM's Fig. 27 check); DR and Q/S give no improvement over S.</li>
</ul>

<h3 style="margin:16px 0 4px 0">2. Shower depth and the hodoscope correction</h3>
<ul style="margin:4px 0;padding-left:20px">
<li>With the 3.3&deg; tilt, a deeper shower is displaced sideways: displacement = shower centre &minus; hodoscope impact point
(map from 100 GeV positrons, 1.3&ndash;1.4 mm residual). S/E rises with it.</li>
<li>Here the slopes of ln S and ln C vs displacement (and impact position) are <i>fitted in the pion data</i>; each run keeps its
own average depth. The fitted slope mixes attenuation with physics (deeper showers are hadron-rich: less S, much less C), so it
is right for the resolution but not for the absolute response. Section 6 replaces it with the measured attenuation.</li>
<li>&chi; becomes 0.22&ndash;0.36 (three methods) against DREAM's 0.25&ndash;0.29, and &lt;S&gt; vs f<sub>em</sub> follows e/h<sub>S</sub> = 1.3.</li>
</ul>

<h3 style="margin:16px 0 4px 0">3. Attenuation (90&deg; positron runs, Y &asymp; 1064 mm)</h3>
{att_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>Runs 1513, 1501, 1507, 1511, 15130 (X = &minus;54.5 to &minus;400 mm); fit of ln(type-summed FERS signal) vs X. Run 1501 is excluded:
its shower lands in different towers (all X = &minus;168 mm runs do).</li>
<li>The scintillating fibres attenuate twice as fast as the Cherenkov fibres. The hottest channels alone give 2.5 / 3.3 / 5.2 m
(used as a systematic). DRS reads few towers and is less reliable here.</li>
<li>Consistent with the depth slope in the pion data at fixed C/S (0.059 vs 0.056 per cm of displacement).</li>
</ul>

<h3 style="margin:16px 0 4px 0">4. Positron-depth reference</h3>
{ref_tab}
<p>The hodoscope zero moved by 0.56 cm after the table rotations; the pions use the positron runs taken just before them.</p>

<h3 style="margin:16px 0 4px 0">5. Electron scale</h3>
{esc_tab}
{off_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>Our positron response is 0.98&ndash;1.04; the pion signals in section 6 are divided by it.</li>
<li>Off-centre positrons respond 2&ndash;10% lower; weighted with the wider pion showers this lowers the pion electron scale by
only 0.5%.</li>
</ul>

<h3 style="margin:16px 0 4px 0">6. Final pion response and resolution</h3>
<p>One chain per event: displacement (dX, dY) &rarr; depth z relative to the positrons (section 4) &rarr; factor exp(z/&lambda;) removed,
with &lambda; from section 3 for S, quartz and plastic separately; then divided by the positron response (section 5). Errors from
the &lambda; uncertainties. Depth from dX alone: z = dX / tan 3.3&deg;. The pion dY/dX is 0.71 at all energies and dX, dY are
correlated event by event, i.e. a vertical tilt of about 2.3&deg; as well; projecting on that direction,
z = (dX + 0.71 dY) / (1.5 tan 3.3&deg;), gives the same mean depth and a slightly better per-event estimate.</p>
<p><b>Depth from dX only:</b></p>
{fin_tab0}
<p><b>Depth from dX and dY (projected):</b></p>
{fin_tab1}
<p><b>Resolution</b> [%], "dX only / dX and dY":</p>
{fres_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>At 80&ndash;100 GeV S is 6&ndash;8% above DREAM, C 6% (100 GeV), Q/S 7&ndash;8%; C/S 0.76&ndash;0.77 (0.78). The hottest-channel
&lambda; values lower all responses by 1.5&ndash;2.5%.</li>
<li>Resolution: the dX + dY depth is 0.1&ndash;0.2 points better than dX alone; the fitted slopes of section 2 are about 1 point
better still in Q/S (they also absorb effects other than attenuation). C alone gets slightly worse than uncorrected, as the
f<sub>em</sub> fluctuation is no longer masked.</li>
</ul>

<h3 style="margin:16px 0 4px 0">7. Open questions</h3>
<ul style="margin:4px 0;padding-left:20px">
<li><b>Remaining ~6% excess</b>: &lambda; extraction (summed vs hottest channels), the 2.3&deg; vertical tilt (not logged; check with the
muon runs 1443&ndash;1446), readout differences (FERS shaping vs DREAM's 120 ns charge gate).</li>
<li><b>20&ndash;40 GeV</b>: S 9&ndash;16% above DREAM (pion purity, beam energy, depth reference at low energy).</li>
<li><b>S and C resolutions</b> (3&ndash;8 points worse, before and after correction): beam spot, 3.3&deg; vs 2&deg; orientation, noise summed over ~900 SiPM channels, &pi;<sup>+</sup> vs &pi;<sup>&minus;</sup>.</li>
</ul>

<h3 style="margin:16px 0 4px 0">Plots (rows)</h3>
<ol style="margin:4px 0;padding-left:20px">
<li>Raw response (with and without leakage correction) and resolution vs DREAM; raw distributions per energy (two rows).</li>
<li>Hodoscope map check and S vs displacement; response and resolution before/after the section-2 correction; distributions (two rows).</li>
<li>&chi;: S vs C at 100 GeV before/after, &lt;S&gt; vs f<sub>em</sub>, &chi; vs energy.</li>
<li>Attenuation fits (FERS, DRS); final response; final resolution (rms, core); final distributions per energy (two rows).
Circles / dashed lines: depth from dX; crosses / solid lines: from dX and dY.</li>
</ol>
"""

plots = ["pion_response_vs_DREAM.png", "pion_response_vs_DREAM_leakcorr.png", "pion_resolution_vs_DREAM.png", "NEWLINE"]
plots += [f"pion_dist_E{e}.png" for e in (20, 30, 40, 60, 80, 100, 120, 160)]
plots += ["NEWLINE", "hodo_map_check_X.png", "hodo_profile_S_dX.png", "hodo_profile_S_dY.png",
          "NEWLINE", "hodo_response.png", "hodo_resolution_rms.png", "hodo_resolution_core.png",
          "NEWLINE"] + [f"hodo_dist_E{e}.png" for e in (20, 30, 40, 60, 80, 100, 120, 160)]
plots += ["NEWLINE", "chi_SvsC_E100_impact.png", "chi_SvsC_E100_edepth.png", "chi_S_vs_fem_E100.png", "chi_vs_energy.png"]
plots += ["NEWLINE", "att90_Y1064_FERS.png", "att90_Y1064_DRS.png", "pion_response_final.png", "pion_resolution_final_rms.png",
          "pion_resolution_final_core.png"]
plots += ["NEWLINE"] + [f"pion_dist_final_E{e}.png" for e in (20, 30, 40, 60, 80, 100, 120, 160)]
os.makedirs(os.path.dirname(OUT), exist_ok=True)
generate_html(plots, PLOTDIR, plots_per_row=4, output_html=OUT,
              title="CaloX pion response and resolution vs DREAM (NIM A 537, 2005)", intro_text=intro)
print(OUT)
