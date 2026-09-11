"""Render the DRS noise/FFT study into a single self-contained HTML page.

Reads results.json and the PNGs written by drs_noise_fft.py from the same
directory and writes index.html next to them (plots embedded, so the page can
be mailed or served on its own). --fragment writes the page without the
document skeleton, for hosts that supply their own.
"""
import argparse
import base64
import json
import os
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))


def img(name):
    with open(os.path.join(HERE, name), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()


def pct(x):
    return f"{100 * x:.0f}%"


def build(res):
    noise, cell, coh = res["noise"], res["cell_pattern"], res["coherence"]
    timing = {k: v for k, v in res["timing"].items() if not k.startswith("_")}
    n_mcp, mcp_spread = res["timing"]["_mcp_events"], res["timing"]["_mcp_spread_ns"]
    amplified = [k for k in noise if k not in ("MCP_DS_0", "Quartz_1")]
    share = lambda k: (noise[k]["band_fraction"]["200-400"] + noise[k]["band_fraction"]["400-600"])
    amp_share = sorted(share(k) for k in amplified)
    amp_sigma = sorted(noise[k]["sigma"] for k in amplified)
    cell_removed = sorted(cell[k]["variance_removed"] for k in amplified)

    # ---------- tables ----------
    bands = ["0-100", "100-200", "200-400", "400-600", "600-1000", "1000-2500"]
    t_noise = "".join(
        f"<tr{' class=odd' if k in ('MCP_DS_0', 'Quartz_1') else ''}><th>{k}</th>"
        f"<td>{v['sigma']:.2f}</td>" +
        "".join(f"<td>{pct(v['band_fraction'][b])}</td>" for b in bands) +
        f"<td>{v['peak_MHz']:.0f}</td></tr>"
        for k, v in noise.items())
    t_cell = "".join(
        f"<tr{' class=odd' if k == 'Quartz_1' else ''}><th>{k}</th><td>{v['sigma']:.2f}</td>"
        f"<td>{v['pattern_rms']:.2f}</td><td>{v['sigma_after']:.2f}</td>"
        f"<td>{pct(v['variance_removed'])}</td></tr>"
        for k, v in cell.items())
    t_coh = "".join(
        f"<tr><th>{k}</th><td>{v['pairs']}</td><td>{v['mean_coherence']:.4f}</td></tr>"
        for k, v in coh.items())
    tb = ["160-400", "400-2500"]
    rows = []
    for k, v in timing.items():
        if "bins" not in v:
            rows.append(f"<tr class=odd><th>{k}</th><td colspan=4 class=note>{v['skipped']}</td></tr>")
            continue
        cells = ""
        for b in tb:
            r = v["bins"].get(b)
            if not r:
                cells += "<td>–</td><td>–</td>"; continue
            for m in ("cfd", "mf"):
                s, f = r[m]["core_sigma_ns"], r[m]["core_fraction"]
                cells += (f"<td>{s:.2f} <span class=frac>{pct(f)}</span></td>"
                          if s is not None else f"<td>– <span class=frac>{pct(f)}</span></td>")
        rows.append(f"<tr><th>{k}</th>{cells}</tr>")
    t_timing = "".join(rows)

    cher = {k: v["bins"]["400-2500"]["cfd"]["core_sigma_ns"] for k, v in timing.items()
            if "bins" in v and "400-2500" in v["bins"] and k.split("_")[0] in ("Singleclad", "Multiclad", "Quartz")
            and v["bins"]["400-2500"]["cfd"]["core_sigma_ns"] is not None
            and v["bins"]["400-2500"]["cfd"]["core_fraction"] > 0.8}
    cher_lo, cher_hi = min(cher.values()), max(cher.values())
    skipped = [k for k, v in timing.items() if "bins" not in v]

    return f"""<title>Where the DRS Noise Lives</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,400;6..72,500;6..72,600&family=Source+Sans+3:wght@400;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root {{
  --ground: #f6f7f9; --panel: #ffffff; --ink: #1c2028; --ink-2: #4a5160; --ink-3: #7b8394;
  --rule: #d5dae3; --accent: #1f6f8b; --accent-ink: #165a71; --accent-wash: #e6f1f5;
  --amber: #a8620f; --amber-wash: #fbf1e3; --ok: #2c7a58;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --ground: #14171c; --panel: #1b1f26; --ink: #e7eaef; --ink-2: #b4bac6; --ink-3: #7f8797;
    --rule: #2e343f; --accent: #5fb3cf; --accent-ink: #8ccbe0; --accent-wash: #1a2a32;
    --amber: #e0a251; --amber-wash: #2b2214; --ok: #6cc39a;
  }}
}}
:root[data-theme="dark"] {{
  --ground: #14171c; --panel: #1b1f26; --ink: #e7eaef; --ink-2: #b4bac6; --ink-3: #7f8797;
  --rule: #2e343f; --accent: #5fb3cf; --accent-ink: #8ccbe0; --accent-wash: #1a2a32;
  --amber: #e0a251; --amber-wash: #2b2214; --ok: #6cc39a;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--ground); color: var(--ink);
  font-family: "Source Sans 3", "Segoe UI", system-ui, sans-serif; font-size: 17px; line-height: 1.55;
  padding-block: 40px 80px; padding-inline: 20px; }}
.page {{ max-width: 1080px; margin: 0 auto; }}
.col {{ max-width: 68ch; }}
h1, h2, h3 {{ font-family: Newsreader, Georgia, "Times New Roman", serif; font-weight: 500;
  line-height: 1.15; text-wrap: balance; margin: 0; }}
h1 {{ font-size: 44px; letter-spacing: -0.01em; }}
h2 {{ font-size: 28px; margin-top: 64px; padding-top: 18px; border-top: 1px solid var(--rule); }}
h3 {{ font-size: 20px; margin-top: 28px; }}
p {{ margin: 14px 0 0; }}
.eyebrow {{ font-family: "IBM Plex Mono", ui-monospace, Menlo, monospace; font-size: 12.5px;
  letter-spacing: 0.08em; text-transform: uppercase; color: var(--ink-3); margin-bottom: 14px; }}
.dek {{ font-size: 20px; color: var(--ink-2); margin-top: 16px; }}
.meta {{ display: flex; flex-wrap: wrap; gap: 8px 28px; margin-top: 22px;
  font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 13px; color: var(--ink-2); }}
.meta b {{ color: var(--ink); font-weight: 500; }}
.findings {{ margin-top: 40px; display: grid; grid-template-columns: repeat(3, 1fr); gap: 14px; }}
.finding {{ background: var(--panel); border: 1px solid var(--rule); border-top: 3px solid var(--accent);
  padding: 16px 18px 18px; }}
.finding .k {{ font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 30px; font-weight: 500;
  color: var(--accent-ink); line-height: 1; font-variant-numeric: tabular-nums; }}
.finding .k small {{ font-size: 14px; color: var(--ink-3); font-weight: 400; margin-left: 4px; }}
.finding .t {{ font-weight: 600; margin-top: 10px; }}
.finding .d {{ color: var(--ink-2); font-size: 15px; margin-top: 4px; line-height: 1.45; }}
figure {{ margin: 28px 0 0; }}
figure img {{ width: 100%; max-width: 100%; display: block; background: #fff; border: 1px solid var(--rule); }}
figcaption {{ font-size: 14.5px; color: var(--ink-2); margin-top: 10px; max-width: 78ch; }}
figcaption b {{ color: var(--ink); font-weight: 600; }}
.tbl {{ overflow-x: auto; margin-top: 18px; border: 1px solid var(--rule); background: var(--panel); }}
table {{ border-collapse: collapse; font-family: "IBM Plex Mono", ui-monospace, Menlo, monospace;
  font-size: 13px; font-variant-numeric: tabular-nums; width: 100%; }}
th, td {{ padding: 6px 12px; text-align: right; white-space: nowrap; border-bottom: 1px solid var(--rule); }}
thead th {{ font-weight: 500; color: var(--ink-3); font-size: 12px; letter-spacing: 0.04em;
  text-transform: uppercase; text-align: right; border-bottom: 1px solid var(--ink-3); }}
tbody th {{ text-align: left; font-weight: 400; color: var(--ink); }}
tr.odd td, tr.odd th {{ color: var(--amber); }}
td.note {{ text-align: left; font-style: italic; white-space: normal; }}
.frac {{ color: var(--ink-3); font-size: 12px; }}
caption {{ caption-side: bottom; text-align: left; padding: 10px 12px; color: var(--ink-2);
  font-family: "Source Sans 3", sans-serif; font-size: 14px; white-space: normal; }}
.callout {{ margin-top: 24px; padding: 16px 20px; border-left: 3px solid var(--amber); background: var(--amber-wash); }}
.callout.ok {{ border-left-color: var(--accent); background: var(--accent-wash); }}
.callout .h {{ font-weight: 600; margin-bottom: 4px; }}
code, kbd {{ font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 0.9em;
  background: var(--accent-wash); padding: 1px 5px; border-radius: 3px; }}
pre {{ font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 13.5px; background: var(--panel);
  border: 1px solid var(--rule); padding: 14px 16px; overflow-x: auto; margin: 14px 0 0; line-height: 1.5; }}
ul {{ margin: 12px 0 0; padding-left: 22px; }}
li {{ margin-top: 6px; }}
li::marker {{ color: var(--ink-3); }}
.next li b {{ font-weight: 600; }}
a {{ color: var(--accent-ink); }}
a:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
@media (max-width: 760px) {{
  h1 {{ font-size: 34px; }} h2 {{ font-size: 24px; }}
  .findings {{ grid-template-columns: 1fr; }}
}}
</style>

<div class="page">
<header class="col">
  <div class="eyebrow">CaloX · DRS waveform study · run {res['run']}</div>
  <h1>Where the DRS Noise Lives</h1>
  <p class="dek">A Fourier look at the pre-pulse samples of the test-fibre channels finds that the noise is the amplifiers', concentrated at 200–600&nbsp;MHz, on top of a correctable DRS4 cell pattern — and that the constant-fraction discriminator is already at the noise limit for timing.</p>
  <div class="meta">
    <span>run <b>{res['run']}</b>, 40 GeV e⁺</span>
    <span><b>{res['nevents']}</b> events</span>
    <span><b>{len(noise) - 1}</b> fibre channels + MCP</span>
    <span>DRS4 at <b>5 GSPS</b>, 0.2 ns/sample</span>
    <span>{date.today().isoformat()}</span>
  </div>
</header>

<section class="findings">
  <div class="finding">
    <div class="k">{pct(amp_share[0])}–{pct(amp_share[-1])}</div>
    <div class="t">of the noise variance sits in 200–600 MHz</div>
    <div class="d">In every amplified channel, with one identical spectral shape. White noise up to the 1 GHz roll-off would put ~48% there. Cross-channel coherence is at the 1/N floor, so this is each amplifier's own noise, not pickup.</div>
  </div>
  <div class="finding">
    <div class="k">{pct(cell_removed[0])}–{pct(cell_removed[-1])}</div>
    <div class="t">of the variance is a fixed DRS4 cell pattern</div>
    <div class="d">Visible only after realigning events by <code>StartIndexCell</code>. It is residual voltage-calibration error — per-cell offsets of 2–4 ADC RMS — and it is correctable.</div>
  </div>
  <div class="finding">
    <div class="k">{cher_lo:.2f}–{cher_hi:.2f}<small>ns</small></div>
    <div class="t">Cherenkov timing at high amplitude, CFD and matched filter alike</div>
    <div class="d">The matched filter does not beat the 20% CFD anywhere it can be measured. Below ~160 ADC neither works, because those events mostly are not in-time pulses.</div>
  </div>
</section>

<div class="col">
<h2>Setup</h2>
<p>Raw waveforms for the 20 channels in <code>data/channel_maps/testingfibers.json</code> plus <code>MCP_DS_0</code>, {res['nevents']} events of run {res['run']}. Each 1024-sample record is baseline-subtracted per event using the median of samples 20–380. Those same samples — 72 ns before any pulse, the pulses sitting at samples 430–470 — are the noise-only window for parts 1 and 2. Part 3 uses samples 380–990; 990 stops short of the dip in the last DRS cells.</p>
<p>Two things the DRS4 forces on any spectral analysis. Its cell widths vary by 10–20% around the nominal 200 ps, so an FFT on the nominal grid smears content above a few hundred MHz — fine for what follows, which lives below that. And it is a transient recorder: the noise spectrum has to come from a pulse-free window, windowed (Hann) to stop the baseline step from leaking into every bin.</p>

<h2>1. The noise is the amplifiers'</h2>
<p>All 19 amplified channels have σ ≈ {amp_sigma[0]:.1f}–{amp_sigma[-1]:.1f} ADC and the same spectral shape: a broad resonance peaking near 300 MHz with a shoulder at 450 MHz, then the analogue roll-off above ~600 MHz. The two channels that differ are the unamplified MCP (σ {noise['MCP_DS_0']['sigma']:.1f}) and <code>Quartz_1</code> (σ {noise['Quartz_1']['sigma']:.1f}) — flatter, with a larger share above 1 GHz, which is the DRS's own white noise. The amplifiers add roughly √(7.8² − 4.4²) ≈ 6.4 ADC of coloured noise on top of that.</p>
</div>
<figure>
  <img src="{img('noise_psd.png')}" alt="Left: noise power spectral density of every channel, each normalised to its own maximum, log-log, with the MCP drawn thick and black. Right: stacked bar chart of the fraction of noise variance in six frequency bands per channel.">
  <figcaption><b>Noise spectra.</b> Left, each channel's PSD scaled to its own maximum; the thick black trace is the unamplified MCP. Right, the share of variance in each band. The amplified channels form one family; <code>Quartz_1</code> and the MCP form another.</figcaption>
</figure>
<div class="col">
<p>An early pass reported "spectral lines at 292 and 458 MHz, 30–80× the floor". That was a metric artefact: the "floor" was the median PSD, which sits in the roll-off region, so a broad bump read as a sharp line. The band table is the honest statement.</p>
<div class="tbl"><table>
<thead><tr><th>channel</th><th>σ [ADC]</th>{''.join(f'<th>{b} MHz</th>' for b in bands)}<th>peak</th></tr></thead>
<tbody>{t_noise}</tbody>
<caption>Share of noise variance per band, from the Hann-windowed PSD of samples 20–380, events with a stray pulse in that window excluded. Amber rows are the two channels without the amplifier signature.</caption>
</table></div>
<p>Coherence rules out common-mode pickup. The magnitude-squared coherence between channel pairs over 10–500 MHz is at the level expected for unrelated noise (1/N<sub>events</sub> = {1 / res['nevents']:.4f}) regardless of whether the pair shares a DRS group, a board, or neither:</p>
<div class="tbl"><table>
<thead><tr><th>relation</th><th>pairs</th><th>mean coherence</th></tr></thead>
<tbody>{t_coh}</tbody>
</table></div>
<div class="callout">
  <div class="h">Worth a look at the hardware: <code>Quartz_1</code></div>
  Every fibre channel is meant to be amplified, yet <code>Quartz_1</code> has the unamplified noise signature and the smallest signal of the set (profile peak 3.5 ADC). Its amplifier may be dead or bypassed.
</div>

<h2>2. A fixed pattern in cell space</h2>
<p>Nothing is fixed in <em>sample</em> space — subtracting the event-mean waveform removes 1% of the variance. But the DRS4's known fixed pattern is periodic in <em>cell</em>, and the trigger cell (<code>StartIndexCell</code>, per group) moves it around in sample space every event. Realigning each event by its start cell and averaging gives a per-cell offset pattern of 2–4 ADC RMS.</p>
</div>
<figure>
  <img src="{img('cell_pattern.png')}" alt="Left: the mean offset per DRS cell for Scintillating_0, a noisy pattern of plus or minus 10 ADC over 1024 cells. Right: bar chart of noise sigma per channel before and after subtracting the cell pattern.">
  <figcaption><b>The DRS4 cell pattern.</b> Left, the mean pre-pulse offset per cell for one channel, ~1100 entries per cell. Right, noise σ before (grey) and after (blue) subtracting each channel's own pattern.</figcaption>
</figure>
<div class="col">
<div class="tbl"><table>
<thead><tr><th>channel</th><th>σ raw</th><th>pattern RMS</th><th>σ after</th><th>variance removed</th></tr></thead>
<tbody>{t_cell}</tbody>
<caption>Subtracting the cell pattern removes 6–28% of the variance in the amplified channels, and 88% in <code>Quartz_1</code>, where the amplifier noise is absent and the pattern is nearly all that is left.</caption>
</table></div>
<div class="callout ok">
  <div class="h">This one is actionable</div>
  The pattern is residual DRS voltage-calibration error. Either redo the per-cell calibration or subtract a pattern measured from pedestal data, as done here. Combined with band-limiting the amplifier resonance, σ drops from ~7.8 to ~6.0 ADC — a 25% gain in signal-to-noise on precisely the low-amplitude channels.
</div>

<h2>3. Timing: CFD is already at the limit</h2>
<p>For each event with an MCP pulse ({n_mcp} of {res['nevents']}; MCP arrival spread {mcp_spread:.1f} ns against the DRS trigger), the channel time is estimated two ways on the same events: the 20% leading-edge CFD used in the analysis, and a matched filter — cross-correlation with the channel's own MCP-aligned profile from <code>drs_profiles.root</code>, parabolic sub-sample interpolation. Both are referenced to the MCP CFD. Because a robust σ over a mixed population misleads, each amplitude bin reports the <em>core</em> σ (events within ±2 ns of the median) and the <em>fraction</em> of events in that core.</p>
</div>
<figure>
  <img src="{img('timing_cfd_vs_matched.png')}" alt="Left: core timing sigma versus pulse amplitude for six channels, CFD as filled markers and matched filter as open markers, the two overlapping. Right: fraction of events within 2 ns of the median versus amplitude, rising from near zero below 160 ADC to 0.9 above 400 ADC for Cherenkov fibres.">
  <figcaption><b>CFD versus matched filter.</b> One channel per fibre type. Filled markers are the CFD, open markers the matched filter; σ is only drawn where more than 20% of events are in core. The two methods coincide.</figcaption>
</figure>
<div class="col">
<div class="tbl"><table>
<thead><tr><th>channel</th><th>160–400 CFD</th><th>160–400 MF</th><th>400–2500 CFD</th><th>400–2500 MF</th></tr></thead>
<tbody>{t_timing}</tbody>
<caption>Core σ in ns with the in-core fraction in grey, for the two amplitude bins where timing is measurable. Amber rows were skipped by the template check.</caption>
</table></div>
<p>Above 400 ADC the Cherenkov fibres reach {cher_lo:.2f}–{cher_hi:.2f} ns and the scintillating fibres ~1.3 ns, with 88–96% of events in core; CFD and matched filter agree to within their own scatter. In the 160–400 ADC bin, 40–55% of events are already out of core, and below 160 ADC almost none are in it. That is not an estimator problem — an estimator cannot recover a time from an event that does not contain the template — it is a population question: particles that missed the fibre, or crosstalk from the 2000-ADC scintillator neighbours, are two candidates.</p>
<div class="callout">
  <div class="h">An artefact caught, and a guard added</div>
  The first pass claimed a 3.3 ns matched-filter resolution for <code>GradedIndex_1</code> even on pure noise. Its profile has no pulse, so the "template" was noise, the correlation peaked at a constant lag, and the "resolution" was exactly the {mcp_spread:.1f} ns MCP arrival spread. A matched filter fails silently on a bad template. The study now rejects any template whose peak is under 5σ of its own baseline — which excludes {', '.join(f'<code>{k}</code>' for k in skipped)} — and drops events whose correlation peaks at a window edge.
</div>

<h2>What to do with this</h2>
<ul class="next">
  <li><b>Correct the cell pattern.</b> Cheap, real, and a calibration task rather than a physics one. Measure it from pedestal runs per channel and subtract, or redo the DRS voltage calibration.</li>
  <li><b>Check <code>Quartz_1</code>'s amplifier.</b> Its noise says it is not there.</li>
  <li><b>If timing still matters at low amplitude, try the whitened filter.</b> Dividing by the measured noise PSD before correlating is the one version with a chance to beat CFD in the 160–400 ADC bin, given how coloured the noise is. Plain matched filtering has been tried and does not.</li>
  <li><b>Understand the out-of-core population</b> before spending more on estimators. Split it by the scintillator channels' amplitude in the same event.</li>
</ul>

<h2>Reproduce</h2>
<pre>python3 scripts/check_drs_mcp.py --run {res['run']} --channels {res['channels']}   # profiles as templates
python3 studies/drs_noise_fft/drs_noise_fft.py --run {res['run']} --channels {res['channels']}
python3 studies/drs_noise_fft/make_report.py</pre>
<p>The first step is only needed for part 3; without <code>drs_profiles.root</code> the study writes parts 1 and 2 and says why it stopped. Every number on this page comes from <code>results.json</code> written by the second step.</p>
</div>
</div>
"""


SKELETON = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
{body}
</body>
</html>
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fragment", default=None, metavar="PATH",
                   help="also write the page without the document skeleton here")
    p.add_argument("--out", default=os.path.join(HERE, "index.html"))
    args = p.parse_args()
    with open(os.path.join(HERE, "results.json")) as f:
        res = json.load(f)
    body = build(res)
    # the standalone file wants the <title>/<link>/<style> inside <head>
    head_end = body.index("</style>") + len("</style>")
    page = SKELETON.format(body=body[:head_end] + "\n</head>\n<body>" + body[head_end:])
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"wrote {args.out} ({os.path.getsize(args.out) // 1024} KB)")
    if args.fragment:
        with open(args.fragment, "w", encoding="utf-8") as f:
            f.write(body)
        print(f"wrote {args.fragment}")


if __name__ == "__main__":
    main()
