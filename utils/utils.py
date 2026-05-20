import math


def number_to_string(n):
    s = str(n)
    return s.replace('-', 'm').replace('.', 'p')


def string2number(s):
    return float(s.replace('m', '-').replace('p', '.'))


def round_up_to_1eN(x):
    """Round a number up to the nearest 10^N."""
    if x <= 0:
        return 0
    return 10 ** math.ceil(math.log10(x))


def getBranchStats(rdf, branches):
    stats = {
        br: {
            "mean": rdf.Mean(br),
            "min": rdf.Min(br),
            "max": rdf.Max(br)
        } for br in branches
    }
    return stats


def get_hist_mpv(hist, window_ts=6.0):
    """Return (mpv, mpv_err) for a ROOT TH1 using a robust core-window method.

    1. Scans the histogram with a sliding window to find the true signal cluster.
    2. Defines a symmetric window (± window_ts) around the cluster center.
    3. Calculates the weighted mean and standard error within this core window.
    """
    n_bins = hist.GetNbinsX()
    xaxis  = hist.GetXaxis()

    if n_bins < 1 or hist.GetEntries() == 0:
        return 0.0, 999.0

    bin_width = xaxis.GetBinWidth(1)

    # 1. Find the densest ±2-bin cluster 
    # Uses Sum of Squares to prioritize distinct peaks over invisible 1-count noise.
    max_score = -1.0
    best_bin  = 1
    for b in range(1, n_bins + 1):
        local_score = 0.0
        # Replaces hist.Integral with a squared sum over the local window
        for i in range(max(1, b - 2), min(n_bins, b + 2) + 1):
            w = hist.GetBinContent(i)
            local_score += w * w
            
        if local_score > max_score:
            max_score = local_score
            best_bin  = b
            
    anchor_x = xaxis.GetBinCenter(best_bin)

    # 2. Fixed geometric window — avoids FWHM walk failures on jagged/Landau tails
    bin_low  = max(1,      xaxis.FindBin(anchor_x - window_ts))
    bin_high = min(n_bins, xaxis.FindBin(anchor_x + window_ts))

    # 3. Weighted mean and RMS inside the window
    sum_w, sum_wx, sum_wx2 = 0.0, 0.0, 0.0
    max_w = 0.0
    for b in range(bin_low, bin_high + 1):
        w = hist.GetBinContent(b)
        x = xaxis.GetBinCenter(b)
        sum_w   += w
        sum_wx  += w * x
        sum_wx2 += w * x * x
        if w > max_w:
            max_w = w

    if sum_w <= 1.0:
        return anchor_x, 999.0

    mpv      = sum_wx / sum_w
    variance = sum_wx2 / sum_w - mpv ** 2
    rms      = math.sqrt(max(0.0, variance))
    stat_err = rms / math.sqrt(sum_w)

    # 4. Poisson plateau penalty
    plateau_err = 0.0
    if max_w > 0:
        # Explicit override guarantees 1-count hits are caught in ultra-low stats
        if max_w <= 4.0:
            threshold = 0.5  
        else:
            threshold = max_w - math.sqrt(max_w)
            
        candidate_x = [xaxis.GetBinCenter(b)
                       for b in range(bin_low, bin_high + 1)
                       if hist.GetBinContent(b) >= threshold]
                       
        if candidate_x:
            spread      = max(candidate_x) - min(candidate_x) + bin_width
            plateau_err = spread / math.sqrt(12.0)

    # 5. Final uncertainty: largest of stat, plateau, and single-bin floor
    mpv_err = max(stat_err, plateau_err, bin_width / math.sqrt(12.0))

    return mpv, mpv_err


def get_channel_var(channel):
    if channel.is_reference:
        return "Ref"
    elif channel.isCer:
        return "CerQuartz" if channel.isQuartz else "CerPlastic"
    elif channel.isSci:
        return "Sci"
    else:
        return "Unknown"
