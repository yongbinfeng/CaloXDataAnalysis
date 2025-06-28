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
    if (i >= v.size() || j > v.size() || i >= j)
        return 0.0;
    return std::accumulate(v.begin() + i, v.begin() + j, 0.0f);
}

float MaxRange(const ROOT::VecOps::RVec<float> &v, size_t i, size_t j)
{
    if (i >= v.size() || j > v.size() || i >= j)
        return 0.0;
    float maxVal = v[i];
    for (size_t k = i + 1; k < j; ++k)
        if (v[k] > maxVal)
            maxVal = v[k];
    return maxVal;
}