import json
import os

# Global caches for performance
_SCAN_RUNS_CACHE = None
_RUN_INFO_CACHE = None


def IsScanRun(runNumber):
    global _SCAN_RUNS_CACHE
    if _SCAN_RUNS_CACHE is None:
        f_scanruns = "data/scanruns.json"
        if os.path.exists(f_scanruns):
            with open(f_scanruns, 'r') as f:
                _SCAN_RUNS_CACHE = json.load(f).get("scanruns", [])
    return runNumber in (_SCAN_RUNS_CACHE or [])


def getRunInfo(runNumber):
    global _RUN_INFO_CACHE
    if _RUN_INFO_CACHE is None:
        f_runinfo = "data/RunlistAugust.json"
        if os.path.exists(f_runinfo):
            with open(f_runinfo, 'r') as f:
                _RUN_INFO_CACHE = json.load(f)

    runinfo = _RUN_INFO_CACHE or {}
    runNum = str(runNumber)
    if runNum not in runinfo:
        return "positrons", 80

    info = runinfo[runNum]
    btype = info.get('beam type', 'positrons')
    benergy = int(info.get('beam energy', '80GeV').replace('GeV', ''))
    return btype, benergy
