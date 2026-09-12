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


CLEANING_PLOTS = ["summary_noise_sigma", "summary_peak_snr", "summary_peak_ratio",
                  "summary_dt_core_sigma", "summary_ringing_coherence", "example_waveform_Sapphire"]


def harvest_cleaning(run):
    """Copy the cleaning study's summary and key plots from results/ next to
    this script, so the report (and the repo) carry them. Returns the summary
    or None if the cleaning study has not been run."""
    import shutil
    root = os.path.join("results", "root", f"Run{run}")
    src = os.path.join(root, "drs_cleaning_summary.json")
    if not os.path.exists(src):
        return None
    shutil.copy(src, os.path.join(HERE, "cleaning_summary.json"))
    for name in CLEANING_PLOTS:
        fn = os.path.join("results", "plots", f"Run{run}", "DRS_Cleaning", f"{name}.png")
        if os.path.exists(fn):
            shutil.copy(fn, os.path.join(HERE, f"cleaning_{name}.png"))
    with open(src) as f:
        return json.load(f)


DECONV_PLOTS = ["response_h", "Scintillating_0_pulse_fit", "Scintillating_0_unfolded",
                "Sapphire_unfolded", "Singleclad_2_unfolded"]


def harvest_deconvolution(run):
    """Copy the unfolding study's summary and key plots next to this script."""
    import shutil
    src = os.path.join("results", "root", f"Run{run}", "drs_deconvolution.json")
    if not os.path.exists(src):
        return None
    shutil.copy(src, os.path.join(HERE, "deconvolution.json"))
    for name in DECONV_PLOTS:
        fn = os.path.join("results", "plots", f"Run{run}", "DRS_Deconvolution", f"{name}.png")
        if os.path.exists(fn):
            shutil.copy(fn, os.path.join(HERE, f"deconvolution_{name}.png"))
    with open(src) as f:
        return json.load(f)


def deconvolution_section(dc):
    if dc is None:
        return ""
    ch = dc["channels"]
    resp = "; ".join(f"<code>{b.replace('DRS_Brg0_', '')}</code> from {', '.join(r['from'])}, control {r['control']}"
                     for b, r in dc["responses"].items())
    def num(v, err=None, dp=2):
        if v is None or v != v:
            return "–"
        return f"{v:.{dp}f}" + (f" ± {err:.{dp}f}" if err is not None and err == err and err < 10 else "")
    rows = ""
    for lab, r in ch.items():
        f, n = r["fit"], r["nnls"]
        spike = n["fwhm_ns"] <= 0.8
        rows += (f"<tr{' class=odd' if r['kind'] == 'control' else ''}><th>{lab}</th><td>{r['kind']}</td>"
                 f"<td>{num(f['tau_rise_ns'], f['tau_rise_err']) if not spike else '–'}</td>"
                 f"<td>{num(f['tau_decay_ns'], f['tau_decay_err']) if not spike else 'spike'}</td>"
                 f"<td>{f['chi2'] / f['ndf']:.2f}</td><td>{n['fwhm_ns']:.1f}</td>"
                 f"<td>{num(n['tau_decay_ns'])}</td><td>{n['chi2_ndf']:.2f}</td></tr>")
    sc = {k: v for k, v in ch.items() if k.startswith("Scintillating")}
    sa = ch.get("Sapphire")
    controls = [v for v in ch.values() if v["nnls"]["fwhm_ns"] <= 0.8]
    floor = max(v["nnls"]["fwhm_ns"] for v in controls) if controls else float("nan")
    sc_txt = " and ".join(f"<code>{k}</code> {v['fit']['tau_decay_ns']:.2f} ± {v['fit']['tau_decay_err']:.2f} ns"
                          for k, v in sc.items())
    return f"""
<h2>5. Unfolding the light: n(t) from the pulses</h2>
<p>The measured pulse is the photon arrival profile n(t) convolved with the response h(t) of SiPM, amplifier and DRS. The Cherenkov fibres see prompt light, so their mean pulse <em>is</em> h(t); it was taken from them and used to unfold the others, two ways. <b>fit</b>: a parametric n(t) — exponential rise into an exponential decay — convolved with h(t) and fitted to the profile with χ² on the profile errors (plus a 2% shape systematic, the channel-to-channel spread of h). Nothing is divided, so noise is never amplified. <b>NNLS</b>: n(t) as a free histogram of arrival times constrained only to be non-negative, projected-gradient least squares stopped at χ²/ndf = 1 — non-parametric, so it can show what the exponential misses. Both work in a −6…+12 ns window around the peak, because the bright scintillating pulses droop beyond that (−13% of the peak at +40 ns) in a way the 60 ADC Cherenkov pulses do not. Responses: {resp}; one Cherenkov channel per board is held out and unfolded as a control.</p>
</div>
<figure>
  <img src="{img('deconvolution_Scintillating_0_pulse_fit.png')}" alt="Scintillating_0 mean pulse with the fitted model n(t) convolved with h(t) overlaid; the two agree across the pulse.">
  <figcaption><b>Forward fit.</b> <code>Scintillating_0</code>'s mcp_clean profile and n(t) ⊛ h(t) with n(t) an exponential rise and decay: χ²/ndf {sc['Scintillating_0']['fit']['chi2'] / sc['Scintillating_0']['fit']['ndf']:.2f}.</figcaption>
</figure>
<figure>
  <img src="{img('deconvolution_Scintillating_0_unfolded.png')}" alt="Unfolded photon arrival profile for Scintillating_0 from the fit model and from NNLS: a fast rise into a roughly 4 ns exponential decay, the two methods agreeing.">
  <figcaption><b>The light.</b> n(t) for <code>Scintillating_0</code> from the fit (red) and from NNLS (blue): a ~0.5 ns rise into a {sc['Scintillating_0']['fit']['tau_decay_ns']:.1f} ns decay. NNLS also shows a faint bump at +8–10 ns — at the edge of the window, where the droop begins, so not something to read physics into yet.</figcaption>
</figure>
<div class="col">
<div class="tbl"><table>
<thead><tr><th>channel</th><th>kind</th><th>τ<sub>rise</sub> [ns]</th><th>τ<sub>decay</sub> [ns]</th><th>fit χ²/ndf</th><th>NNLS FWHM [ns]</th><th>NNLS tail τ</th><th>NNLS χ²/ndf</th></tr></thead>
<tbody>{rows}</tbody>
<caption>Unfolding results. Controls (amber) are Cherenkov channels held out of h(t); "spike" means the unfolded n(t) is a single peak within the resolution floor. NNLS tail τ is an exponential fit to the unfolded tail and includes the rise, so it sits above the model's τ<sub>decay</sub>.</caption>
</table></div>
<p><b>Scintillating fibres.</b> {sc_txt}, with rise times of ~0.5 ns; NNLS agrees on the shape without assuming it. <b>Sapphire</b> is not purely Cherenkov: its unfolded light is {sa['nnls']['fwhm_ns']:.1f} ns wide against ≤{floor:.1f} ns for the pure-Cherenkov controls, and fits a {sa['fit']['tau_decay_ns']:.2f} ± {sa['fit']['tau_decay_err']:.2f} ns decay — a fast luminescence component. <b>Resolution floor</b>: the controls unfold to spikes of {floor:.1f} ns FWHM, which is what the ~0.3 ns alignment jitter and the 0.2 ns sampling allow.</p>
</div>
<figure>
  <img src="{img('deconvolution_Sapphire_unfolded.png')}" alt="Unfolded light profile of the Sapphire channel: a prompt rise and a roughly 2 ns decay, wider than the Cherenkov controls.">
  <figcaption><b>Sapphire.</b> Prompt, then a ~2 ns tail: a slow component the quartz and clad fibres do not have.</figcaption>
</figure>
<div class="col">
<div class="callout">
  <div class="h">Three things that had to be right, or nothing fits</div>
  <b>The response must come from the same DRS board.</b> With h(t) from Board 2 the Board 0 scintillators sat at χ²/ndf 3 and τ<sub>decay</sub> came out 20% short; the residual inter-board jitter after the reference correction (~0.3 ns) is enough to blur a 1.4 ns edge. <b>The pedestal before the pulse must be removed first.</b> The profiles sit 1–5% of the peak above their far baseline in the 10 ns before the pulse; left in, it sets the onset of h(t) instead of the pulse and the fits fail at χ²/ndf ≈ 50. <b>The window must stop before the droop.</b> Beyond +12 ns the bright pulses are not in the regime the Cherenkov response describes. The first two were found the hard way; the script now enforces all three.
</div>
<p>What these numbers are and are not: the shape of n(t) in relative units — absolute photon numbers would need the single-photoelectron gain per channel and a SiPM saturation correction. The Cherenkov light is taken as prompt; chromatic dispersion in a metre of fibre spreads it by a few hundred ps, so h(t) is very slightly too wide and τ correspondingly slightly under. The +3 ns shoulder in h(t) is present in every Cherenkov channel on both boards and is treated as instrumental; if it is in fact optical, the scintillator τ moves by a fraction of that.</p>
"""


def cleaning_section(cl):
    if cl is None:
        return ""
    ch = cl["channels"]
    amp = [k for k in ch if k != "Quartz_1" and ch[k]["n_selected"] >= 30]
    med = lambda key, v: sorted(ch[k][v][key] for k in amp)[len(amp) // 2]
    sig = [med("noise_sigma", v) for v in ("raw", "cell", "cell+LP")]
    snr = [med("peak_snr", v) for v in ("raw", "cell", "cell+LP")]
    cher = [k for k in amp if k.split("_")[0] in ("Singleclad", "Multiclad", "Quartz", "Sapphire")]
    peak_cost = sorted(ch[k]["cell+LP"]["peak_median"] / ch[k]["raw"]["peak_median"] for k in cher)
    integ = sorted(ch[k]["cell+LP"]["integral_median"] / ch[k]["raw"]["integral_median"] for k in cher)
    gain = sorted(ch[k]["ringing_extrapolation_gain"] for k in amp)[len(amp) // 2]
    with_prof = [k for k in amp if "analysis_profile_peak" in ch[k]]
    hit_frac = ch[amp[0]]["mcp_clean_fraction"]
    ratios = sorted((ch[k]["mean_pulse_peak_mcp_clean"] / ch[k]["mean_pulse_peak_all"], k) for k in amp if ch[k]["mean_pulse_peak_all"] > 0)
    hit_ratio = ratios[len(ratios) // 2][0]
    ex = "Multiclad_1" if "Multiclad_1" in ch else amp[0]
    hit_example = f"<code>{ex}</code> {ch[ex]['mean_pulse_peak_all']:.0f} → {ch[ex]['mean_pulse_peak_mcp_clean']:.0f} ADC"
    matched = [ch[k]["profile_matches"] for k in with_prof]
    prof_kind = "raw, mcp_clean" if matched.count("raw, mcp_clean") > len(matched) / 2 else "raw"
    diff_key = "mcp_clean_vs_profile_max_abs_diff" if prof_kind != "raw" else "raw_vs_profile_max_abs_diff"
    strong = [k for k in with_prof if ch[k]["analysis_profile_peak"] > 10]
    prof_agree = max(ch[k][diff_key] for k in with_prof) if with_prof else float("nan")
    holdout = max(abs(ch[k]["cell"]["noise_sigma"] - ch[k]["cell"]["noise_sigma_in_sample"]) for k in amp)
    prof_sentence = ("was booked with <code>--mcp-clean</code>, and the study's <code>raw</code> curve, computed on those same events, reproduces it"
                     if prof_kind != "raw" else "was booked over every event, and the study's <code>raw</code> curve reproduces it")
    fiqr = sorted(ch[k]["ringing_fit_iqr_MHz"] for k in amp)[len(amp) // 2]
    rows = "".join(
        f"<tr><th>{k}</th><td>{ch[k]['n_selected']}</td>"
        + "".join(f"<td>{ch[k][v]['noise_sigma']:.2f}</td>" for v in ("raw", "cell", "cell+LP"))
        + "".join(f"<td>{ch[k][v]['peak_snr']:.0f}</td>" for v in ("raw", "cell", "cell+LP"))
        + "".join(f"<td>{ch[k][v]['dt_core_sigma_ns']:.2f}</td>" if ch[k][v]['dt_core_sigma_ns'] == ch[k][v]['dt_core_sigma_ns']
                  else "<td>–</td>" for v in ("raw", "cell", "cell+LP"))
        + f"<td>{100 * ch[k]['ringing_extrapolation_gain']:+.0f}%</td></tr>"
        for k in ch if ch[k]["n_selected"] >= 30)
    return f"""
<h2>4. Cleaning: what it buys and what it cannot</h2>
<p>Two corrections were applied to the raw waveforms of the whole run and every pulse feature recomputed on the <em>same</em> events (the analysis's <code>mcp_clean</code> selection; timing and rise time additionally need a raw pulse above {cl['min_snr_raw']:.0f}σ within {cl['in_time_ns']:.0f} ns of where the channel's pulses sit). The cell pattern is derived from the first half of the run; the noise σ quoted below is from the held-out second half, and it agrees with the in-sample half to {holdout:.2f} ADC, so the split changes nothing and is kept only as the check. <code>cell</code> subtracts the per-cell offsets; <code>cell+LP</code> adds a zero-phase low-pass at {cl['lowpass_MHz']:.0f} MHz aimed at the amplifier resonance.</p>
<p>The mean pulse per channel is built exactly as the analysis builds its <code>_VS_ts_mcp</code> profiles — every event, shifted by its group reference channel's 50% crossing and by the MCP's CFD from the same <code>[400, 550)</code> window with the same threshold, the float shift floored as <code>Profile1D</code> bins it — and the analysis's own profile from <code>drs_profiles.root</code> is overlaid as a check. The profile in the file {prof_sentence} to within {prof_agree:.2f} ADC at every bin near the peak, on every channel — same events, same alignment, same estimator, so whatever differs between the other curves and the dashed one is the cleaning and nothing else. Getting there took two corrections to the study: a first version averaged only bright events, which on dim channels are mostly out-of-time pulses, and a second took the MCP CFD over the whole record, where events without an MCP pulse lock onto the last-cell dip and land its shape just after the peak as a spurious undershoot.</p>
<p>The check exposes what the profile is made of. Only {hit_frac:.0%} of events pass <code>mcp_clean</code>; the rest get a noise time from the CFD window and smear their fibre pulses into a plateau of a few ADC. On the <code>mcp_clean</code> events alone — where the three treatments are compared, so that the panel compares like with like — the same estimator gives a peak {hit_ratio:.1f}× higher — {hit_example} — with the same shape, and no undershoot. The shoulder ~3 ns after the peak is present in both, so it is a feature of the pulse, not of the averaging. <code>--mcp-clean</code> is the analysis switch that keeps only those events.</p>
<div class="callout ok">
  <div class="h">Found and fixed on the way: the "MCP-aligned" profiles were not MCP-aligned</div>
  <code>variables/drs.py</code> used to set <code>MCP_REF = "MCP_DS_0"</code> and overwrite it on the next line with <code>"MCP_1"</code>, an MCP that runs after 1828 do not have — so <code>_ts_mcp</code> silently fell back to the group reference for every tb2026 run (and, for the same reason, the Sep-2024 runs). The choice now lives in one place, <code>channels.maps.services.get_mcp_reference()</code>: <code>MCP_DS_1</code> where the run has it, else <code>MCP_1</code>, else none. With real MCP alignment the profiles are sharper (<code>Multiclad_0</code> FWHM 10 → 7 samples) and lower, because pulses in events without an MCP hit no longer line up — use <code>--mcp-clean</code> to select them.
</div>
</div>
<figure>
  <img src="{img('cleaning_example_waveform_Sapphire.png')}" alt="One event of the Sapphire channel: raw, cell-corrected and low-passed waveforms overlaid. A 2 ns wide double pulse of 200 ADC sits on a visible 300 MHz ringing of about 15 ADC that continues before and after the pulse.">
  <figcaption><b>One event, three treatments.</b> The Sapphire Cherenkov pulse is ~2 ns wide. The ~300 MHz ringing is present before and after it with the same amplitude — it is the noise of part 1 seen in the time domain, and its period is about the width of the pulse.</figcaption>
</figure>
<div class="col">
<p>Noise σ falls from {sig[0]:.1f} to {sig[1]:.1f} ADC with the cell correction and to {sig[2]:.1f} with the low-pass (medians over the amplified channels). The cell correction costs nothing. The low-pass costs {100 * (1 - peak_cost[-1]):.0f}–{100 * (1 - peak_cost[0]):.0f}% of the Cherenkov peak amplitude — those pulses have real content at the cut-off — while leaving the integral within {100 * (1 - min(integ)):.0f}% and the rise time unchanged at the 0.2 ns sample granularity. Net, the median peak-to-noise ratio goes {snr[0]:.0f} → {snr[1]:.0f} → {snr[2]:.0f}.</p>
<p>Timing does not move. With amplitude above {cl['min_snr_raw']:.0f}σ and within {cl['in_time_ns']:.0f} ns of where the channel's pulses sit the core σ against the MCP is the same to two decimals across all three treatments; at these amplitudes the CFD is not noise-limited, and where noise would matter the events are mostly not in-time pulses (part 3).</p>
</div>
<figure>
  <img src="{img('cleaning_summary_peak_snr.png')}" alt="Bar chart per channel of median peak amplitude divided by noise sigma, for raw, cell-corrected and low-passed waveforms; the corrected versions sit 15 to 25 percent above raw for every channel.">
  <figcaption><b>Peak-to-noise per channel.</b> The scintillating channels start near 140 and the Cherenkov fibres near 35; both gain 15–25% from <code>cell+LP</code>.</figcaption>
</figure>
<div class="col">
<div class="tbl"><table>
<thead><tr><th>channel</th><th>events</th><th>σ raw</th><th>σ cell</th><th>σ cell+LP</th><th>SNR raw</th><th>SNR cell</th><th>SNR cell+LP</th><th>σ<sub>t</sub> raw</th><th>σ<sub>t</sub> cell</th><th>σ<sub>t</sub> cell+LP</th><th>ringing extrap.</th></tr></thead>
<tbody>{rows}</tbody>
<caption>Pre-pulse noise σ [ADC], median peak / σ, and core timing σ [ns] for each treatment; channels with fewer than 30 selected events omitted. Last column: change in post-pulse variance when the pre-pulse sinusoid fit is extrapolated under it.</caption>
</table></div>
<div class="callout">
  <div class="h">The ringing cannot be subtracted, only filtered</div>
  If the 300 MHz ringing kept its phase across the 205 ns record, one could fit it on the pulse-free samples and remove it from under the pulse without touching the signal. It does not: extrapolating a sinusoid fitted on samples 20–380 into samples 620–990 makes the residual <em>worse</em> (median {100 * gain:+.0f}%), and the best-fit frequency scatters by {fiqr:.0f} MHz from event to event. It is a narrow-band random process with a coherence time of a few cycles. That leaves filtering, and the filter band is the signal band — so the real fix for this noise is the amplifier, not the analysis.
</div>
</div>
<figure>
  <img src="{img('cleaning_summary_ringing_coherence.png')}" alt="Bar chart per channel of the post-pulse variance removed by extrapolating the pre-pulse sinusoid fit; every bar is negative, between minus 5 and minus 25 percent.">
  <figcaption><b>Extrapolation test.</b> Negative everywhere: the fitted sinusoid predicts the later samples worse than assuming nothing.</figcaption>
</figure>
<div class="col">
<p>Per-channel distributions of every feature under the three treatments are in <code>results/html/Run{cl['run']}/DRS/DRS_Cleaning_Channels.html</code>, with the summary bars in <code>DRS_Cleaning.html</code> alongside the other DQM pages.</p>
"""


def build(res, cl=None, dc=None):
    noise, cell, coh = res["noise"], res["cell_pattern"], res["coherence"]
    timing = {k: v for k, v in res["timing"].items() if not k.startswith("_")}
    n_mcp, mcp_spread = res["timing"]["_mcp_events"], res["timing"]["_mcp_spread_ns"]
    mcp_label = res.get("mcp_label", "MCP")
    amplified = [k for k in noise if k not in (mcp_label, "Quartz_1")]
    share = lambda k: (noise[k]["band_fraction"]["200-400"] + noise[k]["band_fraction"]["400-600"])
    amp_share = sorted(share(k) for k in amplified)
    amp_sigma = sorted(noise[k]["sigma"] for k in amplified)
    cell_removed = sorted(cell[k]["variance_removed"] for k in amplified)

    # ---------- tables ----------
    bands = ["0-100", "100-200", "200-400", "400-600", "600-1000", "1000-2500"]
    t_noise = "".join(
        f"<tr{' class=odd' if k in (mcp_label, 'Quartz_1') else ''}><th>{k}</th>"
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
    scint = [v["bins"]["400-2500"]["cfd"]["core_sigma_ns"] for k, v in timing.items()
             if "bins" in v and "400-2500" in v["bins"] and k.startswith("Scintillating")
             and v["bins"]["400-2500"]["cfd"]["core_sigma_ns"] is not None]
    scint_lo, scint_hi = (min(scint), max(scint)) if scint else (float("nan"), float("nan"))
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
<p>Raw waveforms for the 20 channels in <code>data/channel_maps/testingfibers.json</code> plus <code>{mcp_label}</code> — the MCP the analysis references its timing to, chosen per run by <code>get_mcp_reference()</code> — {res['nevents']} events of run {res['run']}. Each 1024-sample record is baseline-subtracted per event using the median of samples 20–380. Those same samples — 72 ns before any pulse, the pulses sitting at samples 430–470 — are the noise-only window for parts 1 and 2. Part 3 uses samples 380–990; 990 stops short of the dip in the last DRS cells.</p>
<p>Two things the DRS4 forces on any spectral analysis. Its cell widths vary by 10–20% around the nominal 200 ps, so an FFT on the nominal grid smears content above a few hundred MHz — fine for what follows, which lives below that. And it is a transient recorder: the noise spectrum has to come from a pulse-free window, windowed (Hann) to stop the baseline step from leaking into every bin.</p>

<h2>1. The noise is the amplifiers'</h2>
<p>All 19 amplified channels have σ ≈ {amp_sigma[0]:.1f}–{amp_sigma[-1]:.1f} ADC and the same spectral shape: a broad resonance peaking near 300 MHz with a shoulder at 450 MHz, then the analogue roll-off above ~600 MHz. The two channels that differ are the reference MCP (σ {noise[mcp_label]['sigma']:.1f}) and the entry labelled <code>Quartz_1</code> (σ {noise['Quartz_1']['sigma']:.1f}) — flatter, with a larger share above 1 GHz, which is the DRS's own white noise. <code>Quartz_1</code> looks like an MCP because it is one: see the box below. The amplifiers add roughly √(7.8² − 4.4²) ≈ 6.4 ADC of coloured noise on top of the DRS floor.</p>
</div>
<figure>
  <img src="{img('noise_psd.png')}" alt="Left: noise power spectral density of every channel, each normalised to its own maximum, log-log, with the MCP drawn thick and black. Right: stacked bar chart of the fraction of noise variance in six frequency bands per channel.">
  <figcaption><b>Noise spectra.</b> Left, each channel's PSD scaled to its own maximum; the thick black trace is the reference MCP. Right, the share of variance in each band. The amplified fibre channels form one family; the entry labelled <code>Quartz_1</code> and the MCP form another — because both are MCPs.</figcaption>
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
  <div class="h">The channel listed as <code>Quartz_1</code> is <code>MCP_US_1</code></div>
  <code>data/channel_maps/testingfibers.json</code> maps <code>Quartz_1</code> to <code>DRS_Brg0_Board0_Group3_Channel6</code>, which the channel map assigns to <code>MCP_US_1</code> for runs from 1994. The data say the channel map is right: its pulses are negative and ~130 ADC where every fibre is positive and ~500, its noise is the MCP's (σ 4.4, no amplifier resonance), the analysis inverts it as an MCP, and its "timing against <code>MCP_DS_1</code>" of 0.15 ns is one MCP against the other. Earlier drafts of this page read those symptoms as a fibre with a dead amplifier; they are the upstream MCP. The <code>Quartz_1</code> entry of the JSON needs correcting — or, if a quartz fibre really was on that channel for some runs, the MCP map does.
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
<caption>Subtracting the cell pattern removes 6–28% of the variance in the amplified channels, and 88% in the <code>Quartz_1</code> entry — an MCP channel with no amplifier noise, where the cell pattern is nearly all there is.</caption>
</table></div>
<div class="callout ok">
  <div class="h">This one is actionable</div>
  The pattern is residual DRS voltage-calibration error. Either redo the per-cell calibration or subtract a pattern measured from pedestal data, as done here. Combined with band-limiting the amplifier resonance, σ drops from ~7.8 to ~6.0 ADC — a 25% gain in signal-to-noise on precisely the low-amplitude channels.
</div>

<h2>3. Timing: CFD is already at the limit</h2>
<p>For each event with an MCP pulse ({n_mcp} of {res['nevents']}; MCP arrival spread {mcp_spread:.1f} ns against the DRS trigger), the channel time is estimated two ways on the same events: the 20% leading-edge CFD used in the analysis, and a matched filter — cross-correlation with the channel's own MCP-aligned profile from <code>drs_profiles.root</code>, parabolic sub-sample interpolation. Both are referenced to the MCP CFD, and every time is first corrected by its DRS group's reference channel exactly as the analysis does. That correction is not optional: without it, channels on a different DRS board from the MCP show ~1.3 ns of extra spread that is board jitter, not detector timing. With it, the resolution is the same whichever MCP is the reference. Because a robust σ over a mixed population misleads, each amplitude bin reports the <em>core</em> σ (events within ±2 ns of the median) and the <em>fraction</em> of events in that core.</p>
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
<p>Above 400 ADC the Cherenkov fibres reach {cher_lo:.2f}–{cher_hi:.2f} ns and the scintillating fibres {scint_lo:.2f}–{scint_hi:.2f} ns, with ~90% of events in core; CFD and matched filter agree to within their own scatter. The one place the matched filter shows something is the 160–400 ADC bin of the scintillating fibres, where it puts about twice as many events in core as the CFD (36–40% against 15–19%) at similar width — it finds moderate pulses the leading-edge CFD mis-times, without improving the ones both find. In the 160–400 ADC bin, 40–55% of events are already out of core, and below 160 ADC almost none are in it. That is not an estimator problem — an estimator cannot recover a time from an event that does not contain the template — it is a population question: particles that missed the fibre, or crosstalk from the 2000-ADC scintillator neighbours, are two candidates.</p>
<div class="callout">
  <div class="h">An artefact caught, and a guard added</div>
  The first pass claimed a 3.3 ns matched-filter resolution for <code>GradedIndex_1</code> even on pure noise. Its profile has no pulse, so the "template" was noise, the correlation peaked at a constant lag, and the "resolution" was exactly the {mcp_spread:.1f} ns MCP arrival spread. A matched filter fails silently on a bad template. The study now rejects any template whose peak is under 5σ of its own baseline — which excludes {', '.join(f'<code>{k}</code>' for k in skipped)} — and drops events whose correlation peaks at a window edge.
</div>

{cleaning_section(cl)}
{deconvolution_section(dc)}
<h2>What to do with this</h2>
<ul class="next">
  <li><b>Correct the cell pattern in the pipeline.</b> It is free — no signal cost — and worth ~1 ADC of σ. Measure per channel from pedestal or pre-pulse samples and subtract by cell, or redo the DRS voltage calibration.</li>
  <li><b>Treat the 300 MHz noise as a hardware problem.</b> It cannot be subtracted and sits on the Cherenkov signal band; a low-pass buys 15–25% in peak SNR at the cost of 7–10% of the peak, and nothing in timing. Look at the amplifier's stability and bandwidth.</li>
  <li><b>Regenerate the timing pages</b> for the tb2026 and Sep-2024 runs now that they really align to an MCP.</li>
  <li><b>Fix the <code>Quartz_1</code> entry in <code>testingfibers.json</code></b> — it names the <code>MCP_US_1</code> channel. If a fibre was meant to be there, find where it actually went.</li>
  <li><b>Use the unfolding on the fibre types that still lack signal</b> (GradedIndex, Quartz400) once they have it; the scintillator decay times and Sapphire's slow component are already measurable at the 2% level.</li>
  <li><b>Understand the out-of-time population</b> before spending more on timing estimators: on the dim channels most bright pulses are 20–80 ns off the MCP. Split it by the scintillator channels' amplitude in the same event.</li>
</ul>

<h2>Reproduce</h2>
<pre>python3 scripts/check_drs_mcp.py --run {res['run']} --channels {res['channels']}   # profiles as templates
python3 studies/drs_noise_fft/drs_noise_fft.py --run {res['run']} --channels {res['channels']}
python3 studies/drs_noise_fft/drs_noise_cleaning.py --run {res['run']} --channels {res['channels']}
python3 studies/drs_noise_fft/drs_pulse_deconvolution.py --run {res['run']} --channels {res['channels']}
python3 studies/drs_noise_fft/make_report.py</pre>
<p>The first step is only needed for part 3; without <code>drs_profiles.root</code> the study writes parts 1 and 2 and says why it stopped. Parts 1–3 come from <code>results.json</code>, part 4 from <code>drs_cleaning_summary.json</code>; the report harvests both.</p>
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
    body = build(res, harvest_cleaning(res["run"]), harvest_deconvolution(res["run"]))
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
