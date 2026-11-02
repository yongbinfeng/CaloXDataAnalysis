"""HTML page: CaloX tb2025 pions vs DREAM (NIM A 537, 2005), final version. Built from final_analysis.json."""
import os as _os, sys as _sys
_here = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:0] = [_here, _os.path.dirname(_here)]
from paths import STUDY_DIR, DATA_DIR, REPO_DIR, code
import json, os, sys
REPO = REPO_DIR; SCR = DATA_DIR
sys.path.insert(0, REPO); os.chdir(REPO)
from utils.html_generator import generate_html
F = json.load(open(f"{SCR}/final_analysis.json")); R = F["runs"]; L = F["lam"]
d11 = json.load(open(f"{STUDY_DIR}/dream_fig11.json")); d32 = json.load(open(f"{STUDY_DIR}/dream_fig32.json"))
dS = {int(e): v for e, v in zip(d11["energy_GeV"], d11["S_over_E"])}; dQ = {int(e): v for e, v in zip(d11["energy_GeV"], d11["QS_over_E"])}
dres = {int(e): (s, c, q) for e, s, c, q in zip(d32["energy_GeV"], d32["S_res_pct"], d32["C_res_pct"], d32["QS_res_pct"])}
TH = "padding:3px 10px;border-bottom:1px solid #999;text-align:right"; TD = "padding:3px 10px;text-align:right"
def table(cols, rows):
    head = "".join(f"<th style='{TH}'>{c}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td style='{TD}'>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f"<table style='border-collapse:collapse;font-size:13px;margin:8px 0 14px 0'><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"
D = "&ndash;"; Es = sorted(int(e) for e in R); r100 = R["100"]
lam_tab = table(["Fibres", "&lambda; [m]"], [["Scintillating", f"{L['Sci'][0]:.2f} &plusmn; {L['Sci'][1]:.2f}"],
                                            ["Cherenkov, quartz", f"{L['CerQuartz'][0]:.1f} &plusmn; {L['CerQuartz'][1]:.1f}"],
                                            ["Cherenkov, plastic", f"{L['CerPlastic'][0]:.1f} &plusmn; {L['CerPlastic'][1]:.1f}"]])
resp_tab = table(["E [GeV]", "Run", "Pions", "S/e", "DREAM S/E", "C/e", "DREAM C/E", "C/S", "DREAM C/S", "Q/S/e", "DREAM Q/S/E"],
    [[e, R[str(e)]["run"], R[str(e)]["N"], f"{R[str(e)]['resp']['S']:.3f} &plusmn; {R[str(e)]['eresp']['S']:.3f}", f"{dS[e]:.3f}" if e in dS else D,
      f"{R[str(e)]['resp']['C']:.3f} &plusmn; {R[str(e)]['eresp']['C']:.3f}", "0.640" if e == 100 else D, f"{R[str(e)]['CS']:.3f}", "0.78" if e == 100 else D,
      f"{R[str(e)]['resp']['QS']:.3f} &plusmn; {R[str(e)]['eresp']['Q']:.3f}", f"{dQ[e]:.3f}" if e in dQ else D] for e in Es])
res_tab = table(["E [GeV]", "S", "C", "DR", "Q/S", "Q/S core", "Q/S uncorrected", "DREAM S", "DREAM C", "DREAM Q/S"],
    [[e] + [f"{R[str(e)]['res'][k][0]:.1f}" for k in ("S", "C", "DR", "QS")] + [f"{R[str(e)]['res']['QS'][1]:.1f}",
      f"{R[str(e)]['res_raw']['QS'][0]:.1f}"] + ([f"{x:.1f}" for x in dres[e]] if e in dres else [D] * 3) for e in Es])
unc_tab = table(["E [GeV]", "Pions", "S/E", "DREAM S/E", "C/E", "DREAM C/E", "DR/E", "DREAM DR (Q/S)/E", "S rms", "DREAM S",
                 "C rms", "DREAM C", "DR rms", "DREAM DR (Q/S)"],
    [[e, R[str(e)]["N_unc"], f"{R[str(e)]['resp_unc']['S']:.3f}", f"{dS[e]:.3f}" if e in dS else D, f"{R[str(e)]['resp_unc']['C']:.3f}",
      "0.640" if e == 100 else D, f"{R[str(e)]['resp_unc']['DR']:.3f}", f"{dQ[e]:.3f}" if e in dQ else D]
     + sum(([f"{R[str(e)]['res_unc'][k][0]:.1f}", f"{dres[e][i]:.1f}" if e in dres else D] for i, k in ((0, "S"), (1, "C"), (2, "DR"))), [])
     for e in Es])
RF = json.load(open(f"{SCR}/res_fit.json"))
def pm(v, i): return f"{v[i]:.0f} &plusmn; {v[i + 1]:.0f}" if i == 0 else f"{v[i]:.1f} &plusmn; {v[i + 1]:.1f}"
fit_tab = table(["", "CaloX a [%]", "CaloX b [%]", "DREAM a [%]", "DREAM b [%]", "CaloX a with b = DREAM&#39;s (&chi;<sup>2</sup>/ndf)", "CaloX a &oplus; b", "DREAM a &oplus; b"],
    [[lab, pm(RF[k]["calox"]["lin"], 0), pm(RF[k]["calox"]["lin"], 2), pm(RF[k]["dream"]["lin"], 0), pm(RF[k]["dream"]["lin"], 2), f"{RF[k]['calox']['fixb'][0]:.0f} &plusmn; {RF[k]['calox']['fixb'][1]:.0f} ({RF[k]['calox']['fixb'][3]:.0f})",
      f"{RF[k]['calox']['quad'][0]:.0f}% &oplus; {RF[k]['calox']['quad'][2]:.1f}%", f"{RF[k]['dream']['quad'][0]:.0f}% &oplus; {RF[k]['dream']['quad'][2]:.1f}%"]
     for k, lab in (("S", "S"), ("C", "C"), ("QS", "Q/S (DR)"))])
fit_unc_tab = table(["", "CaloX a [%]", "CaloX b [%]", "DREAM a [%]", "DREAM b [%]", "CaloX a with b = DREAM&#39;s (&chi;<sup>2</sup>/ndf)", "CaloX a &oplus; b", "DREAM a &oplus; b"],
    [[lab, pm(RF["unc"][k]["calox"]["lin"], 0), pm(RF["unc"][k]["calox"]["lin"], 2), pm(RF["unc"][k]["dream"]["lin"], 0), pm(RF["unc"][k]["dream"]["lin"], 2), f"{RF['unc'][k]['calox']['fixb'][0]:.0f} &plusmn; {RF['unc'][k]['calox']['fixb'][1]:.0f} ({RF['unc'][k]['calox']['fixb'][3]:.0f})",
      f"{RF['unc'][k]['calox']['quad'][0]:.0f}% &oplus; {RF['unc'][k]['calox']['quad'][2]:.1f}%", f"{RF['unc'][k]['dream']['quad'][0]:.0f}% &oplus; {RF['unc'][k]['dream']['quad'][2]:.1f}%"]
     for k, lab in (("S", "S"), ("C", "C"), ("DR", "DR (DREAM: Q/S)"))])
chi_tab = table(["E [GeV]", "&chi; from &lt;S&gt; vs f<sub>em</sub> (DREAM Fig. 27)", "e/h<sub>S</sub>", "&chi; for the best DR resolution"],
    [[e, f"{R[str(e)]['chi']['fig27']:.2f}", f"{R[str(e)]['chi']['eh_S']:.2f}", f"{R[str(e)]['chi']['best']:.2f}"] for e in Es])
intro = f"""
<p>CaloX pions (SPS 2025, FERS readout) compared with the original DREAM paper on the same module: N. Akchurin et al., NIM A 537
(2005) 537 (<a href="http://www.phys.ttu.edu/~dream/resources/publications/537_537.pdf">pdf</a>). &pi;<sup>+</sup> runs 1429&ndash;1441,
20&ndash;160 GeV, no target. DREAM values are read off their Figs. 11, 32a and 35.</p>

<div style="border:1px solid #888;background:#f4f4f4;padding:6px 12px;margin:10px 0"><b>Main points</b>
<ul style="margin:4px 0;padding-left:20px">
<li>Without depth information the dual-readout combinations hardly improve on S (100 GeV, all events: Q/S 16.6% vs S 16.9%):
deeper, hadron-rich showers lose less light to the attenuation in the (aged) fibres, which hides the f<sub>em</sub> dependence of S.</li>
<li>Correcting each event for its shower depth (and cutting the deepest 10%, as DREAM did) restores DREAM's behaviour:
&chi; = 0.21&ndash;0.30 above 40 GeV (DREAM 0.25&ndash;0.29); Q/S at 100 GeV {r100['res']['QS'][0]:.1f}% rms, {r100['res']['QS'][1]:.1f}% Gaussian core.</li>
<li>Response at 80&ndash;160 GeV within 5% of DREAM; at 100 GeV S/e {r100['resp']['S']:.3f} (0.830), C/e {r100['resp']['C']:.3f} (0.640),
C/S {r100['CS']:.2f} (0.78), Q/S/e {r100['resp']['QS']:.3f} (1.009).</li>
<li>Resolutions (rms) remain worse than DREAM's because of the stochastic term: S {RF['S']['calox']['lin'][0]:.0f}% vs {RF['S']['dream']['lin'][0]:.0f}%,
C {RF['C']['calox']['lin'][0]:.0f}% vs {RF['C']['dream']['lin'][0]:.0f}%, Q/S {RF['QS']['calox']['lin'][0]:.0f}% vs {RF['QS']['dream']['lin'][0]:.0f}% per &radic;E; the constant terms are
similar ({RF['S']['calox']['lin'][2]:.1f} / {RF['C']['calox']['lin'][2]:.1f} / {RF['QS']['calox']['lin'][2]:.1f}% vs {RF['S']['dream']['lin'][2]:.1f} / {RF['C']['dream']['lin'][2]:.1f} / {RF['QS']['dream']['lin'][2]:.1f}%).
Also open: S 6&ndash;15% above DREAM at 20&ndash;40 GeV.</li>
</ul></div>

<h3 style="margin:16px 0 4px 0">Method</h3>
<ol style="margin:4px 0;padding-left:20px">
<li><b>Selection:</b> muon counter and PSD vetoed, three Cherenkov counters fired; 1&ndash;3 hit fibres per hodoscope axis; S &gt; 0.3 E;
light produced less than {F['depth_cut_cm']:.0f} cm beyond the positron depth (DREAM: centre of gravity in the first metre, their positrons at
30 cm), which removes 7&ndash;10% and improves the rms resolutions by up to 0.7 points (Q/S at 100 GeV: 12.6% &rarr; 11.9%).</li>
<li><b>Depth:</b> the module is tilted, so a deeper shower sits further sideways. Displacement = shower centre (energy-weighted S) minus
impact point (fibre hodoscope, mapped with 100 GeV positrons), relative to the positron runs taken just before the pion scan.
Depth z = displacement along the tilt direction / tan {F['tilt_deg']}&deg;. The direction is measured (dY/dX = {F['r_tilt']:.2f}); the angle
is not known precisely, {F['tilt_deg']}&deg; gives the best Q/S resolution (2&ndash;4&deg; changes the response by up to &plusmn;4%). Without any correction, the pion S/E rises by 4&ndash;8% per cm of displacement in both X and Y: deeper showers
lose less light (plots: "Attenuation").</li>
<li><b>Attenuation:</b> each event's S, quartz-C and plastic-C are scaled by exp(&minus;z/&lambda;), with &lambda; measured on 40 GeV
positrons entering the side of the module (90&deg;) at five positions along the fibres (plots: "Attenuation"). The
impact-position dependence is fitted per run.</li>
<li><b>Electron scale:</b> divided by our positron response at the same energy (S/E and C/E 0.96&ndash;1.06). Errors: from &lambda;.</li>
</ol>
<p>All numbers and plots below use the {F['tilt_deg']}&deg; tilt and the depth cut.</p>
{lam_tab}

<h3 style="margin:16px 0 4px 0">Response</h3>
{resp_tab}
<p>On the calibrated energy scale instead (signal / beam energy, without dividing by our positron response; second plot), at 100 GeV:
S/E {r100['resp_abs']['S']:.3f}, C/E {r100['resp_abs']['C']:.3f}, Q/S/E {r100['resp_abs']['QS']:.3f}, i.e. 3&ndash;4% higher from 60 GeV on,
where our positron S response is 1.03&ndash;1.04.</p>
<h3 style="margin:16px 0 4px 0">Resolution [%]</h3>
<p>&sigma;<sub>rms</sub>/mean (Q/S also Gaussian core); DR = (S &minus; 0.29 C)/0.71, Q/S = DREAM's method (Section 5.3), their DR.
DREAM's S and C are rms; for Q/S the paper does not say (probably rms too), so it is drawn on both plots. "Q/S uncorrected": same
events without the depth correction.</p>
{res_tab}
<p><b>Fits</b> &sigma;/E = a/&radic;E + b (DREAM's linear sum; straight lines vs 1/&radic;E, plots: "Resolution", free b and b fixed to DREAM's), CaloX at
30&ndash;160 GeV (20 GeV left out), DREAM fitted to their Fig. 32a points (20&ndash;250 GeV) the same way; last two columns: quadratic sum.
Errors scaled with &radic;(&chi;<sup>2</sup>/ndf); CaloX &chi;<sup>2</sup>/ndf is 9&ndash;16, i.e. run-to-run scatter beyond statistics.</p>
{fit_tab}
<p>The constant terms are similar to DREAM's (or smaller); the stochastic terms are 1.5&ndash;2.3 times larger. With b fixed to DREAM's
value the stochastic terms are 1.5 (S), 1.6 (C) and 2.0 (Q/S) times DREAM's.</p>
<h3 style="margin:16px 0 4px 0">&chi;</h3>
{chi_tab}

<h3 style="margin:16px 0 4px 0">Without corrections</h3>
<p>All selected events (no depth correction, no impact-position correction, no depth cut), signal / beam energy; DR = (S &minus; 0.29 C)/0.71.
Compared with DREAM's S, C and their DR (the Q/S method).</p>
{unc_tab}
<p>Fits as above (30&ndash;160 GeV; plots: "Without corrections"):</p>
{fit_unc_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>S is 0.87&ndash;1.0 of the beam energy, 17&ndash;25% above DREAM, and DR overshoots to 1.05&ndash;1.11 from 60 GeV on.</li>
<li>DR does not improve on S (100 GeV: 17.2% vs 16.9% rms), where DREAM's DR improves S from 12.2% to 8.2%.</li>
<li>The correction mainly lowers the constant terms (S 10.4 &rarr; 6.8%, DR 10.5 &rarr; 5.6%); the stochastic terms stay 1.3&ndash;2.3 times DREAM's.</li>
</ul>

<h3 style="margin:16px 0 4px 0">Open questions</h3>
<ul style="margin:4px 0;padding-left:20px">
<li>20&ndash;40 GeV response: pion purity, beam energy, or the depth estimate at low energy.</li>
<li>Larger stochastic terms: photostatistics, noise summed over ~900 SiPM channels, sampling; beam spot, &pi;<sup>+</sup> vs &pi;<sup>&minus;</sup>.</li>
<li>The tilt could be measured directly from the DRS arrival time (light from deeper showers arrives earlier).</li>
</ul>

"""
GROUPS = [("Response (corrected): relative to the positron response, and on the energy scale", ["final_response.png", "final_response_abs.png"]),
          ("Resolution (corrected): rms (with DREAM S, C, DR), Gaussian core (with DREAM Q/S), fits vs 1/&radic;E (free b, b = DREAM's)",
           ["final_resolution_rms.png", "final_resolution_core.png", "final_resolution_fit.png", "final_resolution_fit_fixb.png"]),
          ("Without corrections: S, C, DR response, rms, Gaussian core, fits vs 1/&radic;E (free b, b = DREAM's)",
           ["final_uncorr_response.png", "final_uncorr_resolution_rms.png", "final_uncorr_resolution_core.png",
            "final_uncorr_resolution_fit.png", "final_uncorr_resolution_fit_fixb.png"]),
          ("Attenuation: uncorrected pion S vs shower displacement in X and Y; 90&deg; positron runs at Y &asymp; 1064 mm (FERS)",
           ["final_uncorr_S_vs_dX.png", "final_uncorr_S_vs_dY.png", "att90_Y1064_FERS.png"]),
          ("Distributions per energy in GeV (calibrated energy scale): corrected (solid), uncorrected (dashed); legends mean / E, rms / mean", [f"final_dist_E{e}.png" for e in Es])]
plots = sum((g[1] + ["NEWLINE"] for g in GROUPS), [])[:-1]
OUT = "results/html/FERS_PionResponse/PionResponse_vs_DREAM.html"
generate_html(plots, "results/plots/FERS_PionResponse", plots_per_row=4, output_html=OUT,
              title="CaloX pion response and resolution vs DREAM (NIM A 537, 2005)", intro_text=intro)
# a heading above each group of plots (one grid per group)
html = open(OUT).read(); parts = html.split('    <div class="grid">\n')
assert len(parts) == len(GROUPS) + 1, (len(parts), len(GROUPS))
html = parts[0] + "".join(f'    <h3 style="margin:22px 0 6px 0;font-size:16px">{t}</h3>\n    <div class="grid">\n' + part for (t, _), part in zip(GROUPS, parts[1:]))
open(OUT, "w").write(html)
print(OUT)
