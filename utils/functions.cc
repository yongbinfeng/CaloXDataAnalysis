#include "ROOT/RVec.hxx"
#include <algorithm>

ROOT::RVec<int> FillIndices(size_t n)
{
    ROOT::RVec<int> out(n);
    for (size_t i = 0; i < n; ++i)
        out[i] = i;
    return out;
}

float compute_median(ROOT::RVec<float> vec)
{
    if (vec.empty())
        return -9999;
    std::sort(vec.begin(), vec.end());
    size_t n = vec.size();
    if (n % 2 == 0)
        return 0.5 * (vec[n / 2 - 1] + vec[n / 2]);
    else
        return vec[n / 2];
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
    56, 24, 40, 8, 55, 23, 39, 7, 57, 25, 41, 9, 54, 22, 38, 6,
    58, 26, 42, 10, 53, 21, 37, 5, 59, 27, 43, 11, 52, 20, 36, 4,
    60, 28, 44, 12, 51, 19, 35, 3, 61, 29, 45, 13, 50, 18, 34, 2,
    62, 30, 46, 14, 49, 17, 33, 1, 63, 31, 47, 15, 48, 16, 32, 0};

static const int inverse_y[64] = {
    0, 32, 16, 48, 15, 47, 31, 63, 1, 33, 17, 49, 14, 46, 30, 62,
    2, 34, 18, 50, 13, 45, 29, 61, 3, 35, 19, 51, 12, 44, 28, 60,
    4, 36, 20, 52, 11, 43, 27, 59, 5, 37, 21, 53, 10, 42, 26, 58,
    6, 38, 22, 54, 9, 41, 25, 57, 7, 39, 23, 55, 8, 40, 24, 56};

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