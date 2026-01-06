import ROOT
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


class CaloXDataLoader:
    """
    Handles data loading and stores TChain/RDataFrame references 
    as instance attributes to prevent PyROOT memory issues.
    """

    def __init__(self, json_file="data/datafiles.json"):
        if not os.path.exists(json_file):
            raise FileNotFoundError(
                f"Configuration file {json_file} not found.")

        with open(json_file, 'r') as f:
            self._data_map = json.load(f)

        self.tchain = None
        self.rdf_org = None

    def load_rdf(self, run_number, first_event=0, last_event=-1):
        """
        Initializes the RDataFrame and handles TChain lifetime.
        Checks if ROOT files exist and if the TChain is loaded successfully.
        """
        run_num_str = str(run_number)
        if run_num_str not in self._data_map:
            raise ValueError(
                f"Run {run_num_str} not found in {self.json_file}.")

        filelist = self._data_map[run_num_str]
        if not isinstance(filelist, list):
            filelist = [filelist]

        self.tchain = ROOT.TChain("EventTree")
        files_added = 0

        for f in filelist:
            # Check if file exists on disk
            if not os.path.exists(f):
                print(
                    f"\033[91mError: ROOT file {f} does not exist. Skipping.\033[0m")
                continue

            # Add to TChain and check return value (returns number of files added)
            status = self.tchain.Add(f)
            if status > 0:
                print(f"Successfully added file {f} to TChain.")
                files_added += 1
            else:
                print(
                    f"\033[91mError: Failed to load ROOT file {f} into TChain.\033[0m")

        # Verify if any files were successfully loaded
        if files_added == 0:
            raise RuntimeError(
                f"Could not load any ROOT files for run {run_number}.")

        self.rdf_org = ROOT.RDataFrame(self.tchain)

        # Check if the TChain/RDF is actually usable
        try:
            total_events = self.rdf_org.Count().GetValue()
        except Exception as e:
            raise RuntimeError(
                f"RDataFrame initialization failed for run {run_number}: {e}")

        if last_event < 0 or last_event > total_events:
            last_event = total_events

        print(f"\033[94mFiltering events {first_event} to {last_event} "
              f"in run {run_number} (Total: {total_events})\033[0m.")

        return self.rdf_org.Filter(f"event_n >= {first_event} && event_n < {last_event}")
