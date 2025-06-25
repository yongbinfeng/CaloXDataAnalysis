def number2string(n):
    s = str(n)
    return s.replace('-', 'm').replace('.', 'p')


def string2number(s):
    return float(s.replace('m', '-').replace('p', '.'))


def getDataFile(runNumber):
    runNum = str(runNumber)
    import json
    jsonFile = "data/datafiles.json"
    with open(jsonFile, 'r') as f:
        data = json.load(f)
    if runNum in data:
        return data[runNum]
    else:
        raise ValueError(f"Run number {runNum} not found in datafiles.json")


def getBranchStats(rdf, branches):
    stats = {
        br: {
            "mean": rdf.Mean(br),
            "min": rdf.Min(br),
            "max": rdf.Max(br)
        } for br in branches
    }
    return stats


def processDRSBoards(rdf, DRSBoards):
    import re
    import ROOT
    # Get the list of all branch names
    branches = [str(b) for b in rdf.GetColumnNames()]
    pattern = re.compile(r"DRS.*Group.*Channel.*")
    drs_branches = [b for b in branches if pattern.search(b)]
    stats = getBranchStats(rdf, drs_branches)
    print("DRS branches statistics:")
    for br, res in stats.items():
        print(f"{br}: mean = {res['mean'].GetValue():.4f}, "
              f"min = {res['min'].GetValue():.4f}, "
              f"max = {res['max'].GetValue():.4f}")
        stats[br] = {
            "mean": res['mean'].GetValue(),
            "min": res['min'].GetValue(),
            "max": res['max'].GetValue()
        }

    # get the mean of DRS outputs per channel
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for channel in DRSBoard:
            varname = channel.GetChannelName()
            rdf = rdf.Define(
                f"{varname}_mean",
                f"ROOT::VecOps::Mean({varname})"
            )

    ROOT.gInterpreter.Declare("""
    ROOT::RVec<int> FillIndices(size_t n) {
        ROOT::RVec<int> out(n);
        for (size_t i = 0; i < n; ++i) out[i] = i;
        return out;
    }
    """)
    # Create an array of indices for DRS outputs
    rdf = rdf.Define("TS", "FillIndices(1024)")

    ROOT.gInterpreter.Declare("""
    #include "ROOT/RVec.hxx"
    #include <algorithm>

    float compute_median(ROOT::RVec<float> vec) {
        if (vec.empty()) return -9999;
        std::sort(vec.begin(), vec.end());
        size_t n = vec.size();
        if (n % 2 == 0)
            return 0.5 * (vec[n / 2 - 1] + vec[n / 2]);
        else
            return vec[n / 2];
    }
    """)
    # get the mean of DRS outputs per channel
    for varname in drs_branches:
        rdf = rdf.Define(
            f"{varname}_median",
            f"compute_median({varname})"
        )
        rdf = rdf.Define(
            f"{varname}_subtractMedian",
            f"{varname} - {varname}_median"
        )

    return rdf
