import ROOT
import os
import json
from channels.channel_map import build_fers_boards, build_drs_boards
from variables.fers import (
    vectorizeFERS, subtractFERSPedestal, mixFERSHGLG,
    calibrateFERSChannels, get_fers_energy_sum,
    getFERSEnergyWeightedCenter, getFERSEnergyDR, addFERSPosXY
)
from variables.drs import preProcessDRSBoards, calibrateDRSPeakTS, get_drs_stats
from core.selection_manager import SelectionManager
from utils.data_loader import getRunInfo
from utils.plot_helper import get_run_paths
from utils.timing import register_manager


class CaloXAnalysisManager:
    """
    Unified manager for CaloX data analysis. 
    Handles data loading, hardware mapping, calibration, and selections.
    """

    def __init__(self, args):
        self.args = args
        self.run_number = args.run
        register_manager(self)

        self.beam_type, self.beam_energy = getRunInfo(self.run_number)
        self.fersboards = build_fers_boards(run=self.run_number)
        self.drsboards = build_drs_boards(run=self.run_number)
        self.paths = get_run_paths(self.run_number)

        self.tchain = None
        self.rdf_org = None
        self._data_map = None
        self.rdf = self._load_rdf()

        self.branches = {}  # Registry for branched particle managers
        # Track state to prevent re-definitions
        self._steps_applied = set()

    def _load_rdf(self):
        """
        Initializes the RDataFrame and handles TChain lifetime.
        Checks if ROOT files exist and if the TChain is loaded successfully.
        """
        if self._data_map is None:
            if not os.path.exists(self.args.json_file):
                raise FileNotFoundError(
                    f"Configuration file {self.args.json_file} not found.")

            with open(self.args.json_file, 'r') as f:
                self._data_map = json.load(f)

            run_num_str = str(self.run_number)
            if run_num_str not in self._data_map:
                raise ValueError(
                    f"Run {run_num_str} not found in {self.args.json_file}.")

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
                f"Could not load any ROOT files for run {self.run_number}.")
        self.rdf_org = ROOT.RDataFrame(self.tchain)

        first_event = self.args.first_event
        last_event = self.args.last_event
        if last_event < 0:
            last_event = 1e6

        ROOT.RDF.Experimental.AddProgressBar(self.rdf_org)
        return self.rdf_org.Filter(f"event_n >= {first_event} && event_n < {last_event}")

    def prepare(self, do_drs=True, do_fers=True):
        """Initializes standard baseline subtractions and vectorization."""
        if do_drs and "drs_init" not in self._steps_applied:
            self.rdf = preProcessDRSBoards(
                self.rdf, run_number=self.run_number, drsboards=self.drsboards)
            self._steps_applied.add("drs_init")

        if do_fers and "fers_init" not in self._steps_applied:
            self.rdf = vectorizeFERS(self.rdf, self.fersboards)
            self._steps_applied.add("fers_init")
        return self

    def _get_calibration_paths(self, version="Sep", pedestal_run=None):
        """Helper to centralize calibration file path logic."""
        # Logic moved from calibrate_fers to a central helper
        p_run = pedestal_run if pedestal_run else (
            1425 if self.run_number >= 1350 else 1259)

        return {
            "pedestal": f"data/fers/FERS_pedestals_run{p_run}.json",
            "response": f"data/fers/FERS_response_{version}.json",
            "mixing": f"data/fers/FERS_HG2LG_{version}.json"
        }

    def calibrate_fers(self, pedestal_run=None, version="Sep"):
        """Chains FERS pedestal subtraction, HG/LG mixing, and gain calibrations."""
        if f"fers_calib_{version}" in self._steps_applied:
            return self

        paths = self._get_calibration_paths(version, pedestal_run)

        # The updated functions in variables/fers.py now check for existing columns internally
        self.rdf = subtractFERSPedestal(
            self.rdf, self.fersboards, gain="HG", file_pedestals=paths["pedestal"])
        self.rdf = subtractFERSPedestal(
            self.rdf, self.fersboards, gain="LG", file_pedestals=paths["pedestal"])
        self.rdf = mixFERSHGLG(self.rdf, self.fersboards,
                               file_HG2LG=paths["mixing"])
        self.rdf = calibrateFERSChannels(
            self.rdf, self.fersboards, file_calibrations=paths["response"], gain="Mix")

        self._steps_applied.add(f"fers_calib_{version}")
        return self

    def define_physics_variables(self, gain="Mix", calib=True, pdsub=True):
        """Computes sums, weighted centers, and dual-readout variables for a specific gain."""
        # Spatial positions only need to be added once
        if "fers_pos_added" not in self._steps_applied:
            self.rdf = addFERSPosXY(self.rdf, self.fersboards)
            self._steps_applied.add("fers_pos_added")

        # Define basic sums and centers for the requested gain
        self.rdf = get_fers_energy_sum(
            self.rdf, self.fersboards, gain=gain, pdsub=pdsub, calib=calib)
        self.rdf = getFERSEnergyWeightedCenter(
            self.rdf, self.fersboards, gain=gain, pdsub=pdsub, calib=calib)

        # Dual-Readout (DR) logic typically only applies to the 'Mix' gain
        if gain == "Mix" and "fers_dr_defined" not in self._steps_applied:
            self.rdf = getFERSEnergyDR(
                self.rdf, self.fersboards, energy=self.beam_energy)
            self._steps_applied.add("fers_dr_defined")

        return self

    def apply_hole_veto(self, flag_only=False):
        """
        Applies hole vetoes via SelectionManager.
        """
        if "hole_veto_applied" in self._steps_applied:
            return self

        self.sel_mgr = SelectionManager(self.rdf, self.run_number)
        self.rdf = (self.sel_mgr
                    .apply_hole_veto(flag_only=flag_only)
                    .get_rdf())

        self._steps_applied.add("hole_veto_applied")
        return self

    def get_particle_analysis(self, particle_type, flag_only=False):
        """
        Creates an independent branch. Inherits any global filters 
        (like hole veto) already applied to self.rdf.
        """
        branch_mgr = SelectionManager(self.rdf, self.run_number)
        branch_mgr.apply_particle_selection(particle_type, flag_only)

        # Register for automatic reporting
        self.branches[particle_type.lower()] = branch_mgr
        return branch_mgr.get_rdf()

    def report_all(self):
        """Prints the global cutflow (hole veto) and all particle branches."""
        print(f"\n{'#'*72}\n{'COMBINED ANALYSIS SUMMARY':^72}\n{'#'*72}")

        if hasattr(self, 'sel_mgr'):
            print("\n>>> GLOBAL SELECTIONS (Applied to all branches)")
            self.sel_mgr.print_cutflow()

        for name, mgr in self.branches.items():
            print(f"\n>>> PARTICLE BRANCH: {name.upper()}")
            mgr.print_cutflow()

    def get_rdf(self):
        """Returns the final processed RDataFrame node for booking."""
        return self.rdf
