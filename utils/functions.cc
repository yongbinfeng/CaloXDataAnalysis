#include "ROOT/RVec.hxx"
#include <algorithm>

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

    RefHit hit = {kInvalidTime, false};
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
    bool is_valid;
};

// =================================================================
// Constant Fraction Discriminator (CFD) for Physics Pulses
// =================================================================
CaloHit compute_cfd_integral(const ROOT::RVec<float> &waveform)
{
    // --- Configuration Parameters ---
    constexpr float kInvalidTime = -9999.0f;
    constexpr float kMinAmplitude = 5.0f; // Reject noise fluctuations
    constexpr float kCfdFraction = 0.20f; // 20% CFD threshold
    constexpr int kIntWindowPre = 5;      // Slices to integrate BEFORE the CFD anchor
    constexpr int kIntWindowPost = 45;    // Slices to integrate AFTER the CFD anchor

    CaloHit hit = {0.0f, kInvalidTime, false};

    // Safety check: ensure waveform has enough data to integrate
    if (waveform.size() < static_cast<size_t>(kIntWindowPre + kIntWindowPost))
        return hit;

    // --- Find Peak and Amplitude ---
    // (Assuming baseline is 0.0f / already subtracted upstream)
    auto peak_iter = std::max_element(waveform.begin(), waveform.end());
    int peak_idx = std::distance(waveform.begin(), peak_iter);
    float peak_amplitude = *peak_iter;

    if (peak_amplitude < kMinAmplitude)
        return hit;

    // --- Calculate CFD Threshold ---
    float cfd_thresh = kCfdFraction * peak_amplitude;

    // --- Walk Backward to Find Crossing ---
    int scan_idx = peak_idx;
    while (scan_idx > 0 && waveform[scan_idx] > cfd_thresh)
    {
        scan_idx--;
    }

    // Identify the bounding slices around the 20% threshold
    int idx_below = scan_idx;
    int idx_above = scan_idx + 1;

    float v_below = waveform[idx_below];
    float v_above = waveform[idx_above];

    // --- Sub-Slice Linear Interpolation ---
    // Added safety check to prevent division by zero on flat edges
    if (v_above != v_below)
    {
        float frac_idx = idx_below + (cfd_thresh - v_below) / (v_above - v_below);
        hit.time_slice = frac_idx; // Float for sub-slice precision
    }
    else
    {
        hit.time_slice = static_cast<float>(idx_above);
    }

    // --- Fixed-Window Energy Integration ---
    // Anchor the integration window to the nearest integer slice
    int anchor = idx_above;
    int start_idx = std::max(0, anchor - kIntWindowPre);
    int end_idx = std::min(static_cast<int>(waveform.size()), anchor + kIntWindowPost);

    float raw_sum = std::accumulate(waveform.begin() + start_idx, waveform.begin() + end_idx, 0.0f);

    hit.energy = raw_sum;
    hit.is_valid = true;

    return hit;
}