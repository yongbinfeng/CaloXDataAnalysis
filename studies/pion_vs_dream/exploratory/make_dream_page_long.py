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
<p>Comparison of the CaloX test-beam pion data (SPS, 2025) with the original DREAM paper,
N. Akchurin et al., <i>Hadron and jet detection with a dual-readout calorimeter</i>, NIM A 537 (2005) 537
(<a href="http://www.phys.ttu.edu/~dream/resources/publications/537_537.pdf">pdf</a>).
The calorimeter module is the same; the readout changed from photomultipliers to SiPMs with FERS electronics.
The question behind it: why the dual-readout (DR) correction does not improve the CaloX pion resolution over the
scintillation signal (S) alone, while it did in DREAM.</p>

<p><b>Data.</b> &pi;<sup>+</sup> energy scan without the interaction target, all at table position X = 0 mm,
Y = 1032 mm, 3.3&deg;: runs 1441 (20 GeV), 1439 (30), 1438 (40), 1437 (60), 1429 (80), 1432 (100), 1433 (120),
1434 (160), first 25 000 events each. Pion selection (tb2025): muon counter and preshower (PSD) vetoed, all three
Cherenkov counters fired. Signals are the calibrated Mix-gain FERS sums on the electron scale (constants from a fit of
the positron calibration runs to the beam energy). DREAM values are read off their figures (Fig. 11 for the response,
Fig. 32a for the resolutions) or taken from the text. In the plots DREAM is drawn with light filled markers
(open markers are drawn in wrong colours by this ROOT version).</p>

<div style="border:1px solid #888;background:#f4f4f4;padding:6px 12px;margin:10px 0">
<p style="margin:4px 0"><b>Current status (updated 27 Oct 2026; sections 6&ndash;9 and the last row of plots).</b></p>
<ul style="margin:4px 0;padding-left:20px">
<li><b>Why DR did not help:</b> deeper (hadron-rich) showers produce their light closer to the SiPMs and lose less to
attenuation, which hides the f<sub>em</sub> dependence of S. Corrected event by event with the hodoscope displacement,
&chi; &asymp; 0.3 (e/h<sub>S</sub> &asymp; 1.35) and the Q/S resolution improves as in DREAM (sections 4, 4b).</li>
<li><b>Attenuation, measured</b> on the 2025 90&deg; positron runs at Y &asymp; 1064 mm (five table positions, 35 cm along
the fibres, FERS): &lambda; = {L['Sci'][0]:.1f} &plusmn; {L['Sci'][1]:.1f} m (scintillating), {L['CerQuartz'][0]:.1f} &plusmn;
{L['CerQuartz'][1]:.1f} m (quartz), {L['CerPlastic'][0]:.1f} &plusmn; {L['CerPlastic'][1]:.1f} m (plastic). The depth slope
seen in the pion data at fixed C/S matches the scintillation value (section 6).</li>
<li><b>Positron-depth reference:</b> the hodoscope zero differs by 0.56 cm between the positron runs taken before the pion
scan and those after the 13.3&ndash;93.3&deg; rotations; the pions use the earlier one (section 7).</li>
<li><b>Response vs DREAM</b> at positron depth with the measured &lambda; (quartz and plastic separately) and on our own
positron scale: at 100 GeV S/e = {fin['rows']['1']['100']['S']:.3f}&ndash;{fin['rows']['0']['100']['S']:.3f} (DREAM 0.830),
C/e = {fin['rows']['1']['100']['C']:.3f}&ndash;{fin['rows']['0']['100']['C']:.3f} (0.640), Q/S/e =
{fin['rows']['1']['100']['Q']:.3f}&ndash;{fin['rows']['0']['100']['Q']:.3f} (1.009), with a further &minus;1.5 to &minus;2.5%
from how the attenuation lengths are extracted (section 6). The range is whether the vertical
displacement dY is also depth (only the 3.3&deg; X tilt is logged). The earlier ~12% excess is mostly the depth boost
(&asymp; +7&ndash;10% in S, +3&ndash;5% in C) plus our electron scale (&asymp; +3%); the off-centre calibration does not
contribute (section 8). Statements in sections 1 and 4b that are superseded are marked.</li>
</ul></div>

<h3 style="margin:16px 0 4px 0">1. Response</h3>
<p>Mean signal per beam energy, raw (no leakage correction), as in DREAM Fig. 11. DREAM quotes the single-pion
Cherenkov response only at 100 GeV (64.0 GeV, text); its response-vs-energy figures (Fig. 11, Fig. 35) show
the scintillation signal only. "Q/S-corr." is the scintillation signal after DREAM's Q/S correction (Section 5.3,
described in Section 2 below); the DREAM values are read off Fig. 35 ("Pions (after Q/S)").</p>
{resp_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>CaloX pions give more signal relative to electrons at every energy; S reaches 0.96&ndash;1.0 of the beam energy
from 80 GeV on, i.e. the scintillation signal behaves as if nearly compensating.</li>
<li>At 100 GeV, C is 12% above DREAM (0.717 vs 0.640). <i>[Superseded: most of this is the attenuation boost of the
deeper pion showers plus our electron scale, see sections 6&ndash;9.]</i></li>
<li>S is up by more than C: C/S is 0.73 vs DREAM's 0.78 at 100 GeV, independent of the overall scale. About 6% of the
excess is specific to the scintillation signal.</li>
<li>DREAM's own e/h check (slope of the leakage-corrected S vs the em fraction, their Fig. 27) gives p<sub>1</sub>/p<sub>0</sub>
= 0.25&ndash;0.30, e/h<sub>S</sub> &asymp; 1.25&ndash;1.3. The same check on CaloX run 1465 (80 GeV, pion-selected)
gives 0.003, e/h<sub>S</sub> &asymp; 1.0: the mean S does not change with C/S.</li>
<li>After the Q/S correction DREAM's pion response is 0.97&ndash;1.04, equal to the electron response within a few
percent, as the method intends. On CaloX the same correction overshoots to 1.10&ndash;1.24, increasing with energy:
it assumes e/h<sub>S</sub> = 1.3 and scales S up for its hadronic part, but the CaloX S is already close to the
electron scale.</li>
<li><b>DR</b> is the dual-readout energy of make_fers.py (getFERSEnergyDR in variables/fers.py),
(S &minus; 0.29 C)/0.71 on the raw S and C. The Q/S-corrected signal is the same formula with &chi; = 0.288 but on the
leakage-corrected S (DREAM's choice), hence about 10% higher. The raw-S DR response is 0.96&ndash;1.10, above 1 from
60 GeV on.</li>
</ul>

<h3 style="margin:16px 0 4px 0">2. Energy resolution</h3>
<p>&sigma;<sub>rms</sub>/mean, as in DREAM. CaloX: pion selection plus hole veto (in place of DREAM's 1 cm radius
beam-spot cut) and S &gt; 0.3 E against non-showering pions (removes 0.1&ndash;4%); the 20 GeV point (104 events)
is left out of the plot. Q/S = DREAM's beam-energy-free method (Section 5.3): S corrected for lateral leakage,
S(1 + f<sub>out</sub>), em fraction from the leakage-corrected C/S with e/h = 1.3 (S) and 5 (C), then Eq. 7.</p>
{res_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>CaloX raw S and C are already 4&ndash;5 points worse than DREAM, before any DR correction.</li>
<li>The raw-S DR energy is 0.3&ndash;0.7 points worse than S alone at every energy; the Q/S-corrected one is about equal
to S.</li>
<li>DREAM's Q/S correction improves its single pions considerably (12.2% &rarr; 8.2% at 100 GeV); on CaloX it gives
nothing (16.6% &rarr; 17.5%), consistent with the nearly compensating S.</li>
<li>DREAM's headline single-pion result (12.3% &rarr; 2.6% at 100 GeV, Eq. 4) uses the nominal beam energy, like the
CaloX "method 2", and is not a measurable resolution.</li>
</ul>

<h3 style="margin:16px 0 4px 0">2b. Energy distributions</h3>
<p>Per energy (second row of plots): S, C and DR (&chi; = 0.29, raw S) in GeV, up to 1.6 times the beam energy,
as event counts; the legends give mean / beam energy and rms / mean, pion selection plus hole veto, <i>without</i> the S &gt; 0.3 E cut (so the low tails are visible; the numbers
in the legends therefore differ slightly from the table above). Overflow is shown in the last bin.</p>

<h3 style="margin:16px 0 4px 0">3. Electron-scale checks</h3>
<ul style="margin:4px 0;padding-left:20px">
<li><b>HG/LG gain switch</b> (switch at 7500 HG counts): the converted LG value agrees with HG to 0&ndash;2% just below
the switch (median 1.009 / 1.001 for S / C channels in 60 GeV positrons, 1.020 / 1.020 in 80 GeV pions). Positrons put
more of S on low gain (37&ndash;67% from 20 to 100 GeV vs 32% for pions), so this shifts the pion scale by &lt; 1%.</li>
<li><b>SiPM saturation</b>: the hottest S tower keeps a constant share of the positron S (0.22&ndash;0.23 from 20 to
100 GeV), so no sign of saturation.</li>
<li><b>Linearity</b>: positron S/E rises from 0.95 to 1.03 and C/E from 0.92 to 1.01 between 20 and 100 GeV; S and C move
together, pointing to a common offset of about 1&ndash;1.5 GeV (pedestals or thresholds), a couple of percent at 80 GeV.</li>
</ul>

<h3 style="margin:16px 0 4px 0">4. Shower position (depth) dependence and the hodoscope correction</h3>
<p>The detector is tilted with respect to the beam (nominally 3.3&deg;), so a shower that develops deeper is also
displaced sideways. The fibre hodoscope in front of the calorimeter (TTU, FERS boards 0/1, 64 fibres per axis, 0.64 mm
pitch; amplitude-weighted centroid of the fibres above threshold) gives the impact point of every particle, so the
displacement between the shower centre (energy-weighted, scintillation) and the impact point measures the sideways
shift independently of where the particle entered. Selection: pion selection, 1&ndash;3 hit fibres per hodoscope axis and
S &gt; 0.3 E, <i>no</i> hole veto. Plots in the third row onwards.</p>
<ul style="margin:4px 0;padding-left:20px">
<li><b>Hodoscope-to-detector map</b> (from 100 GeV positrons, run 1527, whose showers are shallow): shower centre =
A &middot; hodoscope + b with scale 0.94 (X) and 0.975 (Y), no measurable rotation (off-diagonal &le; 0.004), offset
(&minus;0.69, &minus;1.57) cm; residual 1.4 mm (X) and 1.3 mm (Y), against a beam-spot spread of 8 mm and 14 mm. Checked on
60 GeV positrons (run 1526).</li>
<li><b>Pions</b>: the shower centre sits about +6 mm (X) and +7 mm (Y) beside the impact point (positrons: 0 by
construction), and S/E rises steadily with this displacement. A linear dependence of ln S and ln C on the displacement is
fitted per run and divided out, for all events, together with the impact position (linear and quadratic terms), which
takes over the role of the hole veto: at 100 GeV, Q/S goes from 16.6% (rms) / 13.5% (core) to 12.6% / 9.8% with all
4926 events, against 12.0% / 8.9% with the hole veto (1609 events).</li>
<li><b>Response and resolution</b> are shown before (light) and after (dark) the correction for S, C, DR (&chi; = 0.29 on
raw S) and Q/S (DREAM method, leakage-corrected S), with DREAM. DREAM's raw S and C are &sigma;<sub>rms</sub>/mean;
its corrected (Q/S) resolutions are most likely Gaussian-fit widths (the corrected distributions in the paper are quoted
from Gaussian fits), so DREAM Q/S is drawn on the Gaussian-core plot.</li>
</ul>

<h3 style="margin:16px 0 4px 0">4b. &chi; before and after the depth correction</h3>
<p>Last row of plots. Same pion sample. "Impact" = S and C corrected for the impact position only; "+ depth" =
additionally normalised event by event to positron depth (hodoscope displacement = 0), as DREAM normalised its signals to
the depth of its 40 GeV calibration electrons (their Section 4.3).</p>
<ul style="margin:4px 0;padding-left:20px">
<li><b>S vs C plane</b> (100 GeV): without the depth correction the pion cloud is nearly flat in S (best &chi; = 0.05: S looks
compensating); with it, the cloud lines up with the &chi; = 0.29 line through the positron point, best &chi; = 0.31.</li>
<li><b>&lt;S&gt; vs f<sub>em</sub></b> (DREAM's Fig. 27): flat without the depth correction, and following the e/h<sub>S</sub> = 1.3
expectation with it. Deep showers, which are hadron-rich, get the largest attenuation boost; this cancelled the f<sub>em</sub>
dependence of S and made it look compensating.</li>
<li><b>&chi; vs energy</b>, three ways (centroid of the cloud relative to the positron point; slope of DREAM's Fig. 27; the &chi;
that minimises the DR resolution): 0.05&ndash;0.15 without, 0.22&ndash;0.36 with the depth correction (e/h<sub>S</sub> = 1.24&ndash;1.40),
against DREAM's 0.25&ndash;0.29. The Fig. 27 slope and the best-resolution &chi; come from the event-by-event part of the
correction; the centroid also needs the overall scale at positron depth.</li>
<li>Normalised to positron depth, C/S = 0.787 at 100 GeV (DREAM: 0.78) and the DR energy is 0.98&ndash;1.005 of the beam
energy from 60 GeV on; S and C remain ~12% above DREAM's values together, which points to how the electron scale is set.
<i>[Superseded: this normalisation used the in-sample slope (biased low by the f<sub>em</sub>&ndash;depth correlation,
&asymp; 0 for C) and the positron depth of runs taken after the rotations. Section 9 redoes it with the measured
attenuation lengths and the right reference. &chi; and the resolutions are unchanged by this: the correction is nearly a
per-event scale.]</i></li>
</ul>

<h3 style="margin:16px 0 4px 0">6. Attenuation along the fibres (2025 90&deg; positron runs)</h3>
<p>40 GeV positrons entering the side of the module (table at 90&deg; or 93.3&deg;), so the shower sits at one point along the
fibres and moving the table in X moves it along them. Runs at Y &asymp; 1064 mm: 1513 (X = &minus;54.5 mm), 1501 (&minus;168),
1507 (&minus;218), 1511 (&minus;268) and 15130 (&minus;400.3); 1513 and 15130 are the two data files of DAQ run 1513, with their own
runlist entries, analysed with the tb2025 setup. Positron selection as above. Signal = mean over events of the summed
calibrated FERS signal per fibre type (and the DRS integrals of the channels it reads). Fit: ln(signal) linear in X.
First two plots of the last row.</p>
{cent_tab}
<ul style="margin:4px 0;padding-left:20px">
<li><b>Positions:</b> the DRS arrival time moves by &minus;5.4 to &minus;6.2 ns/m (2026 scan), so it gives the position along
the fibres independently of the table readout; it agrees with the table X within 1&ndash;12 mm.</li>
<li><b>Run 1501 is irregular:</b> its shower lands in different towers (ix = &minus;3.22 against &minus;2.43 to &minus;2.45).
All runs at X = &minus;168 mm (1499, 1501, 1502, taken first) show this shift in their rows; the runlist notes an error and
restart at run 1505&ndash;1507. Possibly irregularities of the side surface, where the beam enters. 1501 is left out of the fit.</li>
</ul>
{att_tab}
<ul style="margin:4px 0;padding-left:20px">
<li><b>Method:</b> per run, all calibrated FERS channels of one type (448 scintillation, 182 quartz, 266 plastic) are summed
event by event, the sum is averaged over the positron events, and ln(average) is fitted linearly in table X. In the four fitted
runs the shower lands in the same towers, so the change of the sum is the light-weighted attenuation of the fibres that carry
the shower. Single channels change much more between runs (16&ndash;84% range of the slopes &minus;0.6 to +1.6 /m for the
channels with &ge; 1% of the sum): small lateral shifts move light between neighbouring channels, which the sum cancels.</li>
<li><b>Cross-check with the hottest channels</b> (each &ge; 3% of its type's sum, tightly grouped slopes): quartz 3.3&ndash;3.4 m
(5 channels), plastic 5.2 m (5 channels); scintillation 2.5 m from the channels with &ge; 1%. The shower core attenuates faster
than the full sum suggests, most clearly for quartz (3.3 vs 6.7 m); possibly the outer channels pick up changing lateral
leakage. Used as a systematic in section 9 (both sets of lengths).</li>
<li><b>FERS</b> (whole module): residuals &le; 0.5% (scintillation) and &le; 1.2% (Cherenkov) over 35 cm. The quartz and plastic
Cherenkov fibres attenuate at about half the rate of the scintillating fibres (summed signals). These are the nominal values used
below.</li>
<li><b>DRS</b> reads only a few towers, so it is much more sensitive to where the shower lands (1501 is 12&ndash;14% off) and
gives shorter lengths; in the other Y rows the shower mostly misses the DRS towers (DRS sums 15&ndash;25 times smaller).</li>
<li><b>Two-run comparisons</b> (50 mm apart) used earlier gave 0.7&ndash;2.5 m because most of them involve the shifted
&minus;168 mm runs; they are superseded by this fit.</li>
<li><b>2026 90&deg; DRS scan</b> (runs 1897&ndash;1931, 34 positions over 33 cm, table Y = &minus;2577 mm): type-summed
&lambda; = 6.0 m (Sci), 19 m (plastic), &gt; 30 m (quartz), residuals 0.3&ndash;1.2%. Different table frame and DRS channel map,
possibly another stretch of fibre; not used for the 2025 pions.</li>
<li><b>Cross-check in the pion data:</b> the slope of ln S against the displacement dX at fixed C/S is 0.059 per cm; with the
3.3&deg; tilt (1 cm of dX = 17 cm of depth) the measured &lambda;<sub>S</sub> = 3.07 m predicts 0.056 per cm.</li>
</ul>

<h3 style="margin:16px 0 4px 0">7. The positron-depth reference</h3>
<p>Median displacement (shower centre minus hodoscope-predicted impact, on the map from run 1527) of the central-position
positron runs. All are at table X = 0, Y = 1032 mm, 3.3&deg;.</p>
{ref_tab}
<p>The hodoscope zero differs by 0.56 cm between the two periods; between them the table was rotated through
13.3&ndash;93.3&deg; and set back to 3.3&deg;. The pion scan (runs 1429&ndash;1441) directly follows run 1423, so the pions are
referred to the mean of the earlier runs, ({fin['ref0'][0]:+.2f}, {fin['ref0'][1]:+.2f}) cm. With the later reference, 20 GeV
pion showers would be only 2 cm deeper than the positrons, which is implausible; with the earlier one about 12 cm.</p>

<h3 style="margin:16px 0 4px 0">8. Electron scale (tests 1 and 2)</h3>
<p><b>Test 1, our own positron response</b> (Gaussian-core peak, central position, positron selection plus 1&ndash;3 hodoscope
fibres per axis). The pion signals in section 9 are divided by these at the same energy (interpolated in ln E; 160 GeV uses
the 120 GeV value).</p>
{esc_tab}
<p><b>Test 2, calibration away from the centre</b>: 80 GeV positrons at other table positions (core peak), relative to the
central run 1409.</p>
{off_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>9&ndash;14 cm off the centre the response is 2&ndash;10% <i>lower</i> (S on average 0.96, C 0.95 of the centre).</li>
<li>Pion showers are wider (25% of S beyond 6 cm of the axis vs 8% for positrons, 100 GeV pions vs 80 GeV positrons, per-channel
FERS). Weighting this map with where the energy goes, the electron scale seen by pions is 0.5% <i>lower</i> than for
central positrons (S and C alike): negligible, and of the opposite sign to the excess.</li>
</ul>

<h3 style="margin:16px 0 4px 0">9. Pion response at positron depth with the measured attenuation</h3>
<p>Same pion sample as section 4 (hodoscope selection, S &gt; 0.3 E). Each event is scaled to the positron depth of section 7
with the attenuation lengths of section 6, separately for the scintillation signal, the quartz and the plastic Cherenkov
sums (1 cm of displacement dX = 1/tan 3.3&deg; = 17.3 cm of depth); the impact-position terms are fitted with the depth
terms fixed. The result is divided by our positron response at the same energy (section 8). Errors: from the
attenuation-length uncertainties. Third plot of the last row; the resolutions are the fourth and fifth.</p>
<p><b>dY not depth</b> (only the 3.3&deg; horizontal tilt is logged):</p>
{fin_tab0}
<p><b>dY also depth</b> (the pion dY displacement grows with energy like dX, which hints at a small vertical angle too):</p>
{fin_tab1}
<ul style="margin:4px 0;padding-left:20px">
<li>At 80&ndash;100 GeV: S is 1&ndash;8% above DREAM, C 3&ndash;6% (100 GeV), Q/S 1&ndash;8%; the lower ends are with dY
treated as depth, which is the largest open item.</li>
<li><b>Attenuation-length systematic:</b> with the hottest-channel lengths (2.5 / 3.3 / 5.2 m) instead of the summed ones, S/e, C/e
and Q/S/e drop by a further 1.5&ndash;2.5% at every energy (100 GeV, dY = depth: S/e 0.814, C/e 0.638, Q/S/e 0.993).</li>
<li>C/S is 0.77&ndash;0.79 at 100 GeV (DREAM 0.78).</li>
<li>At 20&ndash;40 GeV S stays 4&ndash;17% above DREAM (the S &gt; 0.3 E cut moves the mean by &le; 1.3%; 20 GeV has 750 events).</li>
</ul>

<p><b>Resolutions with this correction</b> [%], each cell "dY not depth / dY = depth" (S, C, DR = (S &minus; 0.29 C)/0.71 and
Q/S, all on the electron scale): &sigma;<sub>rms</sub>/mean, and for Q/S also the Gaussian-core &sigma;/mean. For comparison: Q/S
with the in-sample correction of section 4 (fitted slopes, each run moved to its own median depth) and without any correction,
same events. Plots: last row (rms with DREAM S and C; core with DREAM Q/S).</p>
{fres_tab}
<ul style="margin:4px 0;padding-left:20px">
<li>The resolutions hardly depend on which depth the events are moved to; what matters is the event-by-event slope. With dY
treated as depth the measured attenuation lengths give the same resolutions as the in-sample fit (100 GeV Q/S 12.6% rms,
7.3% core); with dX only they are 0.5&ndash;1.1 points worse.</li>
<li>C gets slightly <i>worse</i> than uncorrected (100 GeV: 22.7% &rarr; 23.6&ndash;23.9% rms): deep showers are hadron-rich
and give less C, which the attenuation boost had partly compensated. The f<sub>em</sub> fluctuation is real, and the DR and Q/S
combinations use it.</li>
<li>Against DREAM: Q/S core 7.3&ndash;7.9% at 100 GeV vs 8.2%; S (15.0&ndash;15.2% vs 12.2%) and C (23.6&ndash;23.9% vs 18.6%)
remain worse in rms by 3&ndash;5 points (S) and 5&ndash;8 points (C) from 40 to 100 GeV.</li>
</ul>

<p><b>Distributions</b> (two rows of plots after the resolutions, one per energy): S, C, DR and Q/S with this correction, as
signal / (beam energy &times; positron response); solid lines with dY treated as depth, dashed without. The legends give the mean
and rms/mean; overflow and underflow are in the end bins. Same events as the tables (hodoscope selection, S &gt; 0.3 E).</p>

<h3 style="margin:16px 0 4px 0">10. Open questions and next checks</h3>
<ul style="margin:4px 0;padding-left:20px">
<li><b>Vertical displacement dY</b>: is it depth (a vertical beam angle) or lateral shower containment? This sets the remaining
1&ndash;7% difference to DREAM at 100 GeV. The muon runs at the pion position (1443&ndash;1446), with signal along the full
depth, could tell.</li>
<li><b>Low energies</b>: the 4&ndash;17% S excess at 20&ndash;40 GeV (pion purity, beam energy, or the depth reference at low energy).</li>
<li><b>Readout timing</b>: DREAM integrated the photomultiplier charge in a 120 ns gate (signal arriving after ~30 ns);
FERS shapes the signal and samples its peak (shaping time, hold delay). If they weigh the slow scintillation from
nuclear recoils differently, the hadronic part of S changes relative to electrons and e/h<sub>S</sub> moves towards 1.
Test: integrate the DRS waveforms of calorimeter towers over 20&ndash;200 ns windows for pions and positrons.</li>
<li><b>Worse raw resolution</b>: beam-spot size behind the hole veto, 3.3&deg; vs 2&deg;/0.7&deg; orientation, the
electronic noise summed over ~900 SiPM channels for wide hadron showers, and the &pi;<sup>+</sup> beam composition
(DREAM used &pi;<sup>&minus;</sup>).</li>
<li><b>Jets</b>: the runs with the interaction target could reproduce DREAM's "jet" selection, where their Q/S method
worked best, if the interaction-counter signal is available.</li>
</ul>
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
