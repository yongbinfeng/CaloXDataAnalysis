#include "ROOT/RVec.hxx"
#include <algorithm>
#include <cmath>

ROOT::RVec<int> FillIndices(size_t n)
{
    ROOT::RVec<int> out(n);
    for (size_t i = 0; i < n; ++i)
        out[i] = i;
    return out;
}

float compute_baseline_median(const ROOT::RVec<float> &vec,
                              size_t ibegin = 0,
                              size_t iend = -1)
{
    if (vec.empty())
        return -9999.0f;

    if (iend > vec.size())
        iend = vec.size();

    if (ibegin >= iend)
        return -9999.0f;

    size_t n = iend - ibegin;

    std::vector<float> slice(vec.begin() + ibegin, vec.begin() + iend);

    size_t mid = n / 2;

    //    O(N) partial sort: guarantees slice[mid] is the correct upper-median,
    //    and all elements in [0, mid) are <= slice[mid]
    std::nth_element(slice.begin(), slice.begin() + mid, slice.end());

    return slice[mid];
}

// Standard deviation (RMS about the mean) of vec[ibegin, iend) — used to
// quantify the per-event DRS baseline noise level.
float compute_baseline_rms(const ROOT::RVec<float> &vec,
                           size_t ibegin = 0,
                           size_t iend = -1)
{
    if (vec.empty())
        return -9999.0f;
    if (iend > vec.size())
        iend = vec.size();
    if (ibegin >= iend)
        return -9999.0f;

    size_t n = iend - ibegin;
    double sum = 0.0, sum2 = 0.0;
    for (size_t i = ibegin; i < iend; ++i)
    {
        sum += vec[i];
        sum2 += double(vec[i]) * vec[i];
    }
    double mean = sum / n;
    double var = sum2 / n - mean * mean;
    return var > 0.0 ? std::sqrt(var) : 0.0f;
}

float compute_baseline_average(const ROOT::RVec<float> &waveform, size_t baseline_start = 0, size_t baseline_end = 300)
{
    if (waveform.size() < baseline_end)
        return -9999.0f;

    float baseline_sum = std::accumulate(
        waveform.begin() + baseline_start,
        waveform.begin() + baseline_end,
        0.0f);
    return baseline_sum / (baseline_end - baseline_start);
}

size_t findTrigFireTime(const ROOT::VecOps::RVec<float> &vec, float val_min)
{
    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (vec[i] < val_min / 2.0)
        {
            return i; // return the index of the first value below the threshold
        }
    }
    return -1; // return -1 if no value is below the threshold
}

ROOT::VecOps::RVec<float> clipToZero(const ROOT::VecOps::RVec<float> &vec)
{
    ROOT::VecOps::RVec<float> out;
    for (float v : vec)
    {
        if (fabs(v) < 3.0f)
            v = 0.0f; // clip to zero if below threshold
        out.push_back(v);
    }
    return out;
}

float SumRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j)
{
    j = std::min(j, v.size());
    if (i >= j)
        return 0.0;
    return std::accumulate(v.begin() + i, v.begin() + j, 0.0f);
}

float MaxRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j)
{
    j = std::min(j, v.size());
    if (i >= j)
        return 0.0;
    float maxVal = v[i];
    for (size_t k = i + 1; k < j; ++k)
        if (v[k] > maxVal)
            maxVal = v[k];
    return maxVal;
}

float MinRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j)
{
    j = std::min(j, v.size());
    if (i >= j)
        return 0.0;
    float minVal = v[i];
    for (size_t k = i + 1; k < j; ++k)
        if (v[k] < minVal)
            minVal = v[k];
    return minVal;
}

size_t ArgMinRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j, float threshold = 0.0f)
{
    j = std::min(j, v.size());
    if (i >= j)
        return -1;
    auto minIt = std::min_element(v.begin() + i, v.begin() + j);
    if (*minIt > threshold)
        return -1;
    return std::distance(v.begin(), minIt);
}

size_t ArgMaxRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j, float threshold = 0.0f)
{
    j = std::min(j, v.size());
    if (i >= j)
        return -1;
    auto maxIt = std::max_element(v.begin() + i, v.begin() + j);
    if (*maxIt < threshold)
        return -1;
    return std::distance(v.begin(), maxIt);
}

// --- Hard-coded Inverse Mappings ---
// These arrays are "Inverted":
// The array INDEX is the value you are looking for.
// The array VALUE is the original index/position.

static const int inverse_x[64] = {
    56, 24, 40, 8, 55, 23, 39, 7,
    57, 25, 41, 9, 54, 22, 38, 6,
    58, 26, 42, 10, 53, 21, 37, 5,
    59, 27, 43, 11, 52, 20, 36, 4,
    60, 28, 44, 12, 51, 19, 35, 3,
    61, 29, 45, 13, 50, 18, 34, 2,
    62, 30, 46, 14, 49, 17, 33, 1,
    63, 31, 47, 15, 48, 16, 32, 0};

static const int inverse_y[64] = {
    63, 31, 47, 15, 48, 16, 32, 0,
    62, 30, 46, 14, 49, 17, 33, 1,
    61, 29, 45, 13, 50, 18, 34, 2,
    60, 28, 44, 12, 51, 19, 35, 3,
    59, 27, 43, 11, 52, 20, 36, 4,
    58, 26, 42, 10, 53, 21, 37, 5,
    57, 25, 41, 9, 54, 22, 38, 6,
    56, 24, 40, 8, 55, 23, 39, 7};

// --- Mapping Functions ---

int get_x_index(int val)
{
    if (val < 0 || val > 63)
        return -1;
    return inverse_x[val];
}

int get_y_index(int val)
{
    if (val < 0 || val > 63)
        return -1;
    return inverse_y[val];
}

struct RefHit
{
    float time_slice;
    float peak_value;
    int peak_position;
    bool is_valid;
};

// =================================================================
// Dynamic Leading Edge Discriminator (For Square Reference Pulses)
// =================================================================
RefHit process_dynamic_led(const ROOT::RVec<float> &waveform, bool is_positive = true)
{
    constexpr float kInvalidTime = -9999.0f;
    constexpr float kMinAmplitude = 500.0f;
    constexpr float kLedFraction = 0.50f;

    RefHit hit = {kInvalidTime, 0.0f, -1, false};
    if (waveform.empty())
        return hit;

    const float sign = is_positive ? 1.0f : -1.0f;

    // --- 1. Find Peak (unified via sign) ---
    // Multiply by sign so max_element always finds the extremum we care about
    auto peak_iter = std::max_element(waveform.begin(), waveform.end(),
                                      [sign](float a, float b)
                                      { return sign * a < sign * b; });

    const int peak_idx = std::distance(waveform.begin(), peak_iter);
    const float current_amplitude = std::abs(*peak_iter);

    if (current_amplitude < kMinAmplitude)
        return hit;

    hit.peak_value = current_amplitude;
    hit.peak_position = peak_idx;

    // --- 2. Dynamic Threshold (sign carries polarity) ---
    const float dynamic_threshold = sign * kLedFraction * current_amplitude;

    // --- 3. Walk back from peak while signal is beyond threshold ---
    // Condition: sign * waveform[i] > sign * threshold  (works for both polarities)
    int scan_idx = peak_idx;
    while (scan_idx > 0 && sign * waveform[scan_idx] > sign * dynamic_threshold)
        scan_idx--;

    // --- 4. Anchor to nearest integer slice ---
    hit.time_slice = static_cast<float>(scan_idx + 1);
    hit.is_valid = true;
    return hit;
}

struct CaloHit
{
    float energy;
    float time_slice;
    float peak_value;
    int peak_position;
    float integral_to_peak;
    bool is_valid;
};

// =================================================================
// Constant Fraction Discriminator (CFD) for Physics Pulses
// =================================================================
CaloHit compute_cfd_integral(const ROOT::RVec<float> &waveform, bool is_positive = true, const float min_amplitude = 5.0f,
                             size_t i_start = 0, size_t i_end = static_cast<size_t>(-1),
                             int window_pre = 5, int window_post = 45)
{
    constexpr float kInvalidTime = -9999.0f;
    constexpr float kCfdFraction = 0.20f;

    CaloHit hit = {0.0f, kInvalidTime, 0.0f, -1, 0.0f, false};
    if (waveform.size() < static_cast<size_t>(window_pre + window_post))
        return hit;

    if (i_end > waveform.size())
        i_end = waveform.size();
    if (i_start >= i_end)
        return hit;

    const float sign = is_positive ? 1.0f : -1.0f;

    // --- 1. Find Peak within [i_start, i_end) ---
    auto peak_iter = std::max_element(waveform.begin() + i_start, waveform.begin() + i_end,
                                      [sign](float a, float b)
                                      { return sign * a < sign * b; });

    const int peak_idx = std::distance(waveform.begin(), peak_iter);
    const float peak_amplitude = std::abs(*peak_iter);

    if (peak_amplitude < min_amplitude)
        return hit;

    hit.peak_value = peak_amplitude;
    hit.peak_position = peak_idx;

    // --- 2. CFD Threshold (signed) ---
    const float cfd_thresh = sign * kCfdFraction * peak_amplitude;

    // --- 3. Walk Backward to Crossing ---
    int scan_idx = peak_idx;
    while (scan_idx > static_cast<int>(i_start) && sign * waveform[scan_idx] > sign * cfd_thresh)
        scan_idx--;

    // --- 4. Sub-Slice Linear Interpolation ---
    const int idx_below = scan_idx;
    const int idx_above = scan_idx + 1;
    const float v_below = waveform[idx_below];
    const float v_above = waveform[idx_above];

    hit.time_slice = (v_above != v_below)
                         ? idx_below + (cfd_thresh - v_below) / (v_above - v_below)
                         : static_cast<float>(idx_above);

    // --- 5. Fixed-Window Energy Integration ---
    const int anchor = idx_above;
    const int start_idx = std::max(0, anchor - window_pre);
    const int end_idx = std::min(static_cast<int>(waveform.size()), anchor + window_post);

    const float raw_sum = std::accumulate(waveform.begin() + start_idx,
                                          waveform.begin() + end_idx, 0.0f);
    // Multiply by sign so energy is always a positive magnitude
    hit.energy = sign * raw_sum;
    hit.integral_to_peak = (hit.peak_value > 0.0f) ? hit.energy / hit.peak_value : 0.0f;
    hit.is_valid = true;
    return hit;
}
// =================================================================
// Fixed-window integral with a local pre-pulse baseline
// =================================================================
// Motivation (tb2026 6 mm channels): compute_cfd_integral anchors its window
// on the channel's own 20% CFD crossing and subtracts the global
// median(0,200) baseline. Both assumptions break for the 6 mm fibers:
//
//   * Amplitude. The 6 mm median peak is 16-24 ADC against ~4.3 ADC/sample
//     noise, so the CFD (min_amplitude = 5) runs near threshold and latches
//     onto noise far from the beam window in a large fraction of events. The
//     window then integrates a random stretch of record.
//   * Baseline. The Brg1 (6 mm) boards drift by -2 to -3.8 ADC/sample between
//     the median(0,200) window and the late record. Over a 50-sample window
//     that is a -100 to -190 ADC bias on a ~100-500 ADC signal, and over the
//     full 1024 samples it is the reason every 6 mm full-record integral comes
//     out negative.
//
// This routine instead takes the window position from a caller-supplied
// anchor - use a high-amplitude time reference (the group Channel8 LED
// crossing, or the MCP) rather than the channel itself - and rebaselines on a
// pre-pulse window immediately before the signal, so the drift cancels
// locally. There is no amplitude threshold, so the pedestal (0 photoelectron)
// peak survives; that peak is what makes a photoelectron spectrum readable.
struct FixedHit
{
    float energy;        // baseline-corrected sum over the signal window
    float baseline;      // local pre-pulse level that was subtracted
    float peak_value;    // largest baseline-corrected sample in the window
    int   peak_position; // absolute slice of that sample
    bool  is_valid;
};

FixedHit compute_fixed_window_integral(const ROOT::RVec<float> &waveform,
                                       float anchor,
                                       int window_pre = 4, int window_post = 35,
                                       int baseline_pre_lo = 80,
                                       int baseline_pre_hi = 20)
{
    FixedHit hit = {0.0f, 0.0f, 0.0f, -1, false};
    const int n = static_cast<int>(waveform.size());
    if (n == 0 || !std::isfinite(anchor))
        return hit;

    const int a = static_cast<int>(std::lround(anchor));
    const int sig_lo = a - window_pre;
    const int sig_hi = a + window_post;              // exclusive
    const int bl_lo  = a - baseline_pre_lo;
    const int bl_hi  = a - baseline_pre_hi;          // exclusive
    if (sig_lo < 0 || sig_hi > n || bl_lo < 0 || bl_hi > n || bl_hi <= bl_lo)
        return hit;

    // Mean, not median: the pre-pulse window is signal-free by construction,
    // and the mean is the lower-variance estimator of the level there.
    double bl = 0.0;
    for (int i = bl_lo; i < bl_hi; ++i)
        bl += waveform[i];
    bl /= static_cast<double>(bl_hi - bl_lo);

    double sum = 0.0;
    float best = -1.0e30f;
    int best_idx = sig_lo;
    for (int i = sig_lo; i < sig_hi; ++i)
    {
        const double v = waveform[i] - bl;
        sum += v;
        if (v > best) { best = static_cast<float>(v); best_idx = i; }
    }

    hit.energy = static_cast<float>(sum);
    hit.baseline = static_cast<float>(bl);
    hit.peak_value = best;
    hit.peak_position = best_idx;
    hit.is_valid = true;
    return hit;
}

// =================================================================
// CFD start, baseline-crossing end: a pulse-adaptive integration window
// =================================================================
// The two fixed-length definitions both assume one window length fits every
// channel. Measured on tb2026 run 1917, the pulse width above 10% of peak is
// ~13 slices for 3 mm fibers but ~42-45 for 6 mm, so a 50-slice window is
// generous for one and marginal for the other, and no single choice is right.
//
// Here the window starts at the 20% CFD crossing, as compute_cfd_integral
// does, but ends where the pulse actually returns to baseline: walking forward
// from the peak until the corrected waveform stays at or below stop_level for
// stop_consecutive consecutive slices. The width then adapts per pulse.
//
// stop_level = max(stop_abs, stop_fraction * peak). The absolute floor matters
// because on a small pulse a pure fraction-of-peak threshold sits inside the
// noise and the walk terminates immediately; stop_abs should be set near the
// per-sample noise RMS (~4 ADC for the 6 mm channels, ~10 for 3 mm).
//
// The peak search can be restricted to [search_center +/- search_half]. Use
// that for low-amplitude channels: searching the whole record lets the peak,
// and therefore the whole window, land on a noise excursion.
struct AdaptiveHit
{
    float energy;
    float baseline;
    float time_slice;     // sub-sample CFD crossing
    float peak_value;
    int   peak_position;
    int   start_idx;
    int   end_idx;        // exclusive
    int   width;
    bool  hit_max_width;  // true if the walk was stopped by the cap, not by the pulse
    bool  is_valid;
};

AdaptiveHit compute_cfd_to_baseline_integral(const ROOT::RVec<float> &waveform,
                                             float search_center = -1.0f,
                                             int search_half = 60,
                                             float cfd_fraction = 0.20f,
                                             float stop_fraction = 0.05f,
                                             float stop_abs = 4.0f,
                                             int stop_consecutive = 3,
                                             int max_width = 250,
                                             int baseline_pre_lo = 80,
                                             int baseline_pre_hi = 20,
                                             float min_amplitude = 0.0f,
                                             int smooth_half = 3,
                                             int min_width = 40)
{
    AdaptiveHit hit = {0.0f, 0.0f, -9999.0f, 0.0f, -1, -1, -1, 0, false, false};
    const int n = static_cast<int>(waveform.size());
    if (n < 16)
        return hit;

    // --- 1. Peak search, optionally restricted to a window around an anchor ---
    int lo = 0, hi = n;
    if (search_center >= 0.0f)
    {
        const int c = static_cast<int>(std::lround(search_center));
        lo = std::max(0, c - search_half);
        hi = std::min(n, c + search_half + 1);
        if (lo >= hi)
            return hit;
    }
    int peak_idx = lo;
    float peak_raw = waveform[lo];
    for (int i = lo + 1; i < hi; ++i)
        if (waveform[i] > peak_raw) { peak_raw = waveform[i]; peak_idx = i; }

    // --- 2. Local pre-pulse baseline (same rationale as the fixed-window one) ---
    const int bl_lo = peak_idx - baseline_pre_lo;
    const int bl_hi = peak_idx - baseline_pre_hi;
    if (bl_lo < 0 || bl_hi <= bl_lo || bl_hi > n)
        return hit;
    double bl = 0.0;
    for (int i = bl_lo; i < bl_hi; ++i)
        bl += waveform[i];
    bl /= static_cast<double>(bl_hi - bl_lo);

    const float peak = static_cast<float>(peak_raw - bl);
    if (peak <= 0.0f || peak < min_amplitude)
        return hit;

    // --- 3. Walk back to the CFD crossing ---
    const float cfd_thresh = cfd_fraction * peak;
    int scan = peak_idx;
    while (scan > lo && (waveform[scan] - bl) > cfd_thresh)
        --scan;
    const int idx_below = scan;
    const int idx_above = std::min(scan + 1, n - 1);
    const float v_below = static_cast<float>(waveform[idx_below] - bl);
    const float v_above = static_cast<float>(waveform[idx_above] - bl);
    hit.time_slice = (v_above != v_below)
                         ? idx_below + (cfd_thresh - v_below) / (v_above - v_below)
                         : static_cast<float>(idx_above);

    // --- 4. Walk forward from the peak to the baseline crossing ---
    // The stop test runs on a moving average, not raw samples. On the 6 mm
    // channels the pulse tail (~20 ADC peak) is comparable to the per-sample
    // noise (~4.3 ADC), so an unsmoothed test trips on the first downward
    // fluctuation and the window closes after a few slices instead of the ~45
    // the pulse actually occupies.
    const float stop_level = std::max(stop_abs, stop_fraction * peak);
    auto smoothed = [&](int i) -> float {
        const int a = std::max(0, i - smooth_half);
        const int b = std::min(n - 1, i + smooth_half);
        double acc = 0.0;
        for (int k = a; k <= b; ++k)
            acc += waveform[k] - bl;
        return static_cast<float>(acc / static_cast<double>(b - a + 1));
    };
    int end = peak_idx + 1;
    int run_below = 0;
    while (end < n && (end - idx_above) < max_width)
    {
        if (smoothed(end) <= stop_level)
        {
            if (++run_below >= stop_consecutive)
            {
                // End the window at the first of the consecutive low samples.
                end = end - stop_consecutive + 1;
                break;
            }
        }
        else
        {
            run_below = 0;
        }
        ++end;
    }
    if (end >= n || (end - idx_above) >= max_width)
    {
        hit.hit_max_width = true;
        end = std::min(n, idx_above + max_width);
    }

    const int start = idx_above;
    // Floor the width. Without it the window collapses to a few slices on
    // low-amplitude channels -- the smoothed trace is already at the stop
    // level when the walk begins -- which under-counts small pulses and makes
    // the response non-linear exactly where it needs to stay linear. The floor
    // should be about the pulse length of the fiber type being integrated.
    if (end - start < min_width)
        end = std::min(n, start + min_width);
    if (end <= start)
        return hit;

    double sum = 0.0;
    for (int i = start; i < end; ++i)
        sum += waveform[i] - bl;

    hit.energy = static_cast<float>(sum);
    hit.baseline = static_cast<float>(bl);
    hit.peak_value = peak;
    hit.peak_position = peak_idx;
    hit.start_idx = start;
    hit.end_idx = end;
    hit.width = end - start;
    hit.is_valid = true;
    return hit;
}

// =================================================================
// Leading-edge rise time: CFD frac_lo -> frac_hi crossing (default 20% -> 90%)
// =================================================================
// Same peak search / local baseline as compute_cfd_to_baseline_integral, so
// the two share an anchor and a pedestal treatment. The walk backward from
// the peak passes through frac_hi before frac_lo (it is closer to the peak),
// so both crossings come from one pass: continue the same backward scan used
// for the first (frac_hi) crossing on to the second (frac_lo) one, using the
// same linear-interpolation convention as the CFD crossing above.
struct RiseTimeHit
{
    float time_lo;      // sub-sample TS of the frac_lo crossing
    float time_hi;      // sub-sample TS of the frac_hi crossing
    float rise_time;    // time_hi - time_lo [TS]
    float peak_value;
    int   peak_position;
    float baseline;
    bool  is_valid;
};

RiseTimeHit compute_rise_time(const ROOT::RVec<float> &waveform,
                              float search_center = -1.0f,
                              int search_half = 60,
                              float frac_lo = 0.20f,
                              float frac_hi = 0.90f,
                              int baseline_pre_lo = 80,
                              int baseline_pre_hi = 20,
                              float min_amplitude = 0.0f)
{
    RiseTimeHit hit = {-9999.0f, -9999.0f, -9999.0f, 0.0f, -1, 0.0f, false};
    const int n = static_cast<int>(waveform.size());
    if (n < 16)
        return hit;

    // --- 1. Peak search, optionally restricted to a window around an anchor ---
    int lo = 0, hi = n;
    if (search_center >= 0.0f)
    {
        const int c = static_cast<int>(std::lround(search_center));
        lo = std::max(0, c - search_half);
        hi = std::min(n, c + search_half + 1);
        if (lo >= hi)
            return hit;
    }
    int peak_idx = lo;
    float peak_raw = waveform[lo];
    for (int i = lo + 1; i < hi; ++i)
        if (waveform[i] > peak_raw) { peak_raw = waveform[i]; peak_idx = i; }

    // --- 2. Local pre-pulse baseline ---
    const int bl_lo = peak_idx - baseline_pre_lo;
    const int bl_hi = peak_idx - baseline_pre_hi;
    if (bl_lo < 0 || bl_hi <= bl_lo || bl_hi > n)
        return hit;
    double bl = 0.0;
    for (int i = bl_lo; i < bl_hi; ++i)
        bl += waveform[i];
    bl /= static_cast<double>(bl_hi - bl_lo);

    const float peak = static_cast<float>(peak_raw - bl);
    if (peak <= 0.0f || peak < min_amplitude)
        return hit;

    // --- 3. Walk back from the peak, recording the frac_hi crossing first
    //         (higher up the rising edge, closer to the peak) and the
    //         frac_lo crossing second, from the same backward scan. ---
    const float thresh_hi = frac_hi * peak;
    const float thresh_lo = frac_lo * peak;
    bool found_hi = false;
    int scan = peak_idx;
    while (scan > lo && (waveform[scan] - bl) > thresh_lo)
    {
        if (!found_hi && (waveform[scan] - bl) <= thresh_hi)
        {
            const int idx_below = scan;
            const int idx_above = std::min(scan + 1, n - 1);
            const float v_below = static_cast<float>(waveform[idx_below] - bl);
            const float v_above = static_cast<float>(waveform[idx_above] - bl);
            hit.time_hi = (v_above != v_below)
                              ? idx_below + (thresh_hi - v_below) / (v_above - v_below)
                              : static_cast<float>(idx_above);
            found_hi = true;
        }
        --scan;
    }
    if (!found_hi)
        return hit;
    {
        const int idx_below = scan;
        const int idx_above = std::min(scan + 1, n - 1);
        const float v_below = static_cast<float>(waveform[idx_below] - bl);
        const float v_above = static_cast<float>(waveform[idx_above] - bl);
        hit.time_lo = (v_above != v_below)
                          ? idx_below + (thresh_lo - v_below) / (v_above - v_below)
                          : static_cast<float>(idx_above);
    }

    hit.rise_time = hit.time_hi - hit.time_lo;
    hit.peak_value = peak;
    hit.peak_position = peak_idx;
    hit.baseline = static_cast<float>(bl);
    hit.is_valid = (hit.rise_time > 0.0f);
    return hit;
}
