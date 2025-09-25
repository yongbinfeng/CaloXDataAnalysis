def IsScanRun(runNumber):
    import json
    f_scanruns = "data/scanruns.json"
    with open(f_scanruns, 'r') as f:
        temp = json.load(f)
        scanruns = temp["scanruns"]
    return runNumber in scanruns


def getDataFile(runNumber, jsonFile):
    runNum = str(runNumber)
    import json
    with open(jsonFile, 'r') as f:
        data = json.load(f)
    if runNum in data:
        return data[runNum]
    else:
        raise ValueError(f"Run number {runNum} not found in {jsonFile}")


def loadRDF(runNumber, firstEvent=0, lastEvent=-1, jsonFile="data/datafiles.json"):
    import ROOT
    # Open the input ROOT file
    ifile = getDataFile(runNumber, jsonFile)
    infile = ROOT.TFile(ifile, "READ")
    rdf_org = ROOT.RDataFrame("EventTree", infile)
    nevents = rdf_org.Count().GetValue()
    if lastEvent < 0 or lastEvent > nevents:
        lastEvent = nevents
    print(
        f"\033[94mfiltering events from {firstEvent} to {lastEvent} in run {runNumber}, total {nevents} events in the file.\033[0m")
    # Apply the event range filter
    rdf = rdf_org.Filter(f"event_n >= {firstEvent} && event_n < {lastEvent}")

    return rdf, rdf_org


def getRunInfo(runNumber):
    import json
    f_runinfo = "data/RunlistAugust.json"
    with open(f_runinfo, 'r') as f:
        temp = json.load(f)
        runinfo = temp
    runNum = str(runNumber)
    if runNum not in runinfo:
        # default is positrons 80 GeV
        return "positrons", 80
        # raise ValueError(f"Run number {runNum} not found in {f_runinfo}")
    btype = runinfo[runNum]['beam type']
    benergy = int(runinfo[runNum]['beam energy'].replace('GeV', ''))
    # print(f"Run {runNum}: beam type = {btype}, beam energy = {benergy} GeV")
    return btype, benergy
