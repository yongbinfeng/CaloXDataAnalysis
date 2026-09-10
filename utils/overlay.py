"""Interactive overlay pages: many channels on one set of axes, pick which.

The page is nothing but y-vs-x line series, so it carries the arrays rather
than a rendered ROOT canvas and is drawn by Plotly in the browser. Building it
needs no RDataFrame and no channel map -- only histograms that already exist
in a ROOT file -- which is why scripts/make_drs_overlay.py can regenerate a
page without redoing the analysis.
"""
import os
import statistics

from utils.html_generator import generate_plotly_picker_html

# Plain CSS colours: the page draws from the raw arrays, so ROOT's colour
# table is not involved. tab10 first, then ten further well-separated shades.
OVERLAY_COLOURS = [
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
    "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f",
    "#393b79", "#843c39", "#5254a3", "#8c6d31", "#a55194",
    "#637939", "#3182bd", "#e6550d", "#31a354", "#756bb1",
]


def build_channel_overlay(entries, output_html, xlabel, ylabel, title,
                          intro_text, filename, y_pad=0.15):
    """Write an overlay page for [(hist, label), ...].

    All histograms must share one uniform x binning; the x values are then
    implied by (first bin centre, bin width) and never stored. Returns the
    path written, or None if there was nothing worth plotting.
    """
    if len(entries) < 2:
        print(f"Overlay needs at least 2 channels, got {len(entries)}; skipping")
        return None

    ref = entries[0][0]
    nbins = ref.GetNbinsX()
    x0, dx = ref.GetBinCenter(1), ref.GetBinWidth(1)

    series = []
    lo, hi = float("inf"), float("-inf")            # raw
    hi_n, norm_minima = float("-inf"), []           # peak-normalised
    for i, (hist, label) in enumerate(entries):
        if hist.GetNbinsX() != nbins:
            print(f"  {label}: {hist.GetNbinsX()} bins, expected {nbins}; skipped")
            continue
        if abs(hist.GetBinWidth(1) - dx) > 1e-9 or abs(hist.GetBinCenter(1) - x0) > 1e-9:
            print(f"  {label}: x binning differs from {entries[0][1]}; skipped")
            continue
        ys, vals, peak = [], [], 0.0
        for b in range(1, nbins + 1):
            v = hist.GetBinContent(b)
            ys.append(round(v, 4))
            vals.append(v)
            peak = max(peak, v)
            if hist.GetBinError(b) > 0 or v != 0:
                lo, hi = min(lo, v), max(hi, v)
        # The page divides by this peak to compare shapes; sending the peak
        # rather than a second array keeps the payload the same size and keeps
        # the normalisation identical to the range computed here.
        usable = abs(peak) > 1e-10
        if usable:
            series_min = float("inf")
            for b, v in enumerate(vals):
                if hist.GetBinError(b + 1) > 0 or v != 0:
                    series_min = min(series_min, v / peak)
                    hi_n = max(hi_n, v / peak)
            if series_min < float("inf"):
                norm_minima.append(series_min)
        else:
            print(f"  {label}: peak is ~0, cannot normalise; shown unscaled")
        series.append({"name": label,
                       "color": OVERLAY_COLOURS[i % len(OVERLAY_COLOURS)],
                       "peak": peak if usable else 0.0,
                       "y": ys})

    if len(series) < 2:
        print("Overlay: fewer than 2 usable channels; skipping")
        return None
    if lo == float("inf"):
        lo, hi = -1.0, 1.0
    # Normalising divides by the peak, so a low-amplitude channel turns into
    # huge excursions and its minimum alone would set the axis. Take the median
    # of the per-channel minima instead: half the curves then fit entirely, and
    # the reader can type a range or auto-scale for the rest.
    if norm_minima:
        lo_n = max(statistics.median(norm_minima), -3.0)
    else:
        lo_n, hi_n = -1.0, 1.0
    if hi_n == float("-inf"):
        hi_n = 1.0
    pad = (hi - lo) * y_pad or 1.0
    pad_n = (hi_n - lo_n) * y_pad or 1.0

    spec = {"x0": x0, "dx": dx, "xlabel": xlabel, "ylabel": ylabel,
            "xrange": [x0 - dx / 2.0, x0 + (nbins - 0.5) * dx],
            "yrange": [lo - pad, hi + pad],
            "yrangeNorm": [lo_n - pad_n, hi_n + pad_n],
            "filename": filename,
            "series": series}
    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
    return generate_plotly_picker_html(spec, output_html=output_html,
                                       title=title, intro_text=intro_text)
