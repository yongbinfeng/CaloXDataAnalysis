import ROOT
import os
import json
from channels.channel_map import build_fers_boards, build_drs_boards

_BEAM_TO_PARTICLE = {
    "pion":      "pion",
    "pions":     "pion",
    "pi+":       "pion",
    "positron":  "electron",
    "positrons": "electron",
    "e+":        "electron",
}

from variables.fers import (
    vectorizeFERS, subtractFERSPedestal, mixFERSHGLG,
    calibrateFERSChannels, get_fers_energy_sum,
    getFERSEnergyWeightedCenter, getFERSEnergyDR, addFERSPosXY,
    buildTTUHodo
)
from variables.drs import get_psd_energy_deposit, process_drs_data
from core.selection_manager import SelectionManager
from utils.data_loader import getRunInfo
from utils.plot_helper import get_run_paths
from utils.timing import register_manager
from core.hist_book import HistBook


class CaloXAnalysisManager:
    """
    Unified manager for CaloX data analysis. 
    Handles data loading, hardware mapping, calibration, and selections.
    """

    def __init__(self, args, load_data=True):
        self.args = args
        self.run_number = args.run
        register_manager(self)

        self.beam_type, self.beam_energy = getRunInfo(self.run_number)
        self.fersboards = build_fers_boards(run_number=self.run_number)
        self.drsboards = build_drs_boards(run_number=self.run_number)

        self.paths = get_run_paths(self.run_number)
        for path in self.paths.values():
            if not os.path.exists(path):
                os.makedirs(path)
                print(f"Created directory: {path}")

        self.tchain = None
        self.rdf_org = None
        self._data_map = None
        self.do_detailed_plots = True
        self.include_pedestals = True
        self.rdf = self._load_rdf() if load_data else None

        self.hbook = HistBook(self.paths["root"])
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

    def prepare(self, do_drs=True, do_fers=True, do_hodo=True, do_mcp=True,
                linear_baseline=False):
        """Initializes standard baseline subtractions and vectorization.

        linear_baseline=True switches DRS baseline subtraction to the two-anchor
        linear method (removes post-pulse pedestal droop); default is the
        constant median(0,200) baseline.
        """
        if do_drs and "drs_init" not in self._steps_applied:
            self.rdf = process_drs_data(
                self.rdf, self.run_number, drsboards=self.drsboards,
                do_mcp=do_mcp, linear_baseline=linear_baseline)
            if do_mcp:
                self.rdf = get_psd_energy_deposit(self.rdf, self.run_number)
            self._steps_applied.add("drs_init")

        if do_fers and "fers_init" not in self._steps_applied:
            self.rdf = vectorizeFERS(self.rdf, self.fersboards)
            self._steps_applied.add("fers_init")

        if do_hodo and "hodo_init" not in self._steps_applied:
            self.rdf = buildTTUHodo(self.rdf)
            self._steps_applied.add("hodo_init")
        return self

    def _get_calibration_paths(self, version="Sep", pedestal_run=None):
        """Helper to centralize calibration file path logic."""
        # tb2026 (run_number > 1700): dedicated calibration set.
        if self.run_number > 1700:
            pedestal_path = (f"data/fers/FERS_pedestals_run{pedestal_run}.json"
                             if pedestal_run else
                             "data/fers/FERS_pedestals_tb2026.json")
            return {
                "pedestal": pedestal_path,
                "response": "data/fers/FERS_response_tb2026.json",
                "mixing": "data/fers/FERS_HG2LG_tb2026.json",
                "deadchannels": "data/fers/deadchannels_tb2026.json",
            }

        p_run = pedestal_run if pedestal_run else (
            1425 if self.run_number >= 1350 else 1259)
        return {
            "pedestal": f"data/fers/FERS_pedestals_run{p_run}.json",
            "response": f"data/fers/FERS_response_{version}.json",
            "mixing": f"data/fers/FERS_HG2LG_{version}.json",
            "deadchannels": "data/fers/deadchannels.json",
        }

    def calibrate_fers(self, pedestal_run=None, version="Sep"):
        """Chains FERS pedestal subtraction, HG/LG mixing, and gain calibrations."""
        if f"fers_calib_{version}" in self._steps_applied:
            return self

        paths = self._get_calibration_paths(version, pedestal_run)

        print("\033[92mFERS calibration files:\033[0m")
        print(f"  pedestal     : {paths['pedestal']}")
        print(f"  HG->LG mixing: {paths['mixing']}")
        print(f"  calibration  : {paths['response']}")
        print(f"  dead channels: {paths['deadchannels']}")

        # The updated functions in variables/fers.py now check for existing columns internally
        self.rdf = subtractFERSPedestal(
            self.rdf, self.fersboards, gain="HG", file_pedestals=paths["pedestal"])
        self.rdf = subtractFERSPedestal(
            self.rdf, self.fersboards, gain="LG", file_pedestals=paths["pedestal"])
        self.rdf = mixFERSHGLG(self.rdf, self.fersboards,
                               file_HG2LG=paths["mixing"])
        self.rdf = calibrateFERSChannels(
            self.rdf, self.fersboards, file_calibrations=paths["response"], gain="Mix",
            file_deadchannels=paths["deadchannels"])

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

    def _get_or_create_sel_mgr(self):
        """Return the shared SelectionManager, creating it from current rdf if needed."""
        if not hasattr(self, 'sel_mgr'):
            self.sel_mgr = SelectionManager(self.rdf, self.run_number)
        return self.sel_mgr

    def apply_hole_veto(self, flag_only=False):
        """Applies hole veto via SelectionManager."""
        if "hole_veto_applied" in self._steps_applied:
            return self

        self.rdf = (self._get_or_create_sel_mgr()
                    .apply_hole_veto(flag_only=flag_only)
                    .get_rdf())

        self._steps_applied.add("hole_veto_applied")
        return self

    def apply_beam_pid_selection(self, flag_only=False, particle=None):
        """Applies particle ID selection based on beam type.

        Electron beams: electron selection (PSD fired, Cherenkovs fired).
        Pion beams: pion selection (PSD vetoed, Cherenkovs fired).
        No-op for unrecognised beam types (e.g. muon).

        Pass particle= to override the beam-type lookup (e.g. select protons
        from a pion run).
        """
        particle = particle or _BEAM_TO_PARTICLE.get((self.beam_type or "").lower())
        if particle is None:
            return self

        # Per-particle guard so distinct selections stack (e.g. mcp_clean +
        # electron), while re-applying the same particle stays idempotent.
        step_key = f"beam_pid_applied_{particle.lower()}"
        if step_key in self._steps_applied:
            return self

        self.rdf = (self._get_or_create_sel_mgr()
                    .apply_particle_selection(particle, flag_only=flag_only)
                    .get_rdf())

        self._steps_applied.add(step_key)
        return self

    def apply_mcp_selection(self, flag_only=False,
                            min_integral_to_peak=4, min_peak_value=10):
        """Applies MCP clean-pulse selection via SelectionManager with cutflow tracking."""
        if "mcp_selection_applied" in self._steps_applied:
            return self

        self.rdf = (self._get_or_create_sel_mgr()
                    .apply_mcp_selection(flag_only=flag_only,
                                         min_integral_to_peak=min_integral_to_peak,
                                         min_peak_value=min_peak_value)
                    .get_rdf())

        self._steps_applied.add("mcp_selection_applied")
        return self

    def apply_mcp_diff_selection(self, flag_only=False,
                                 min_integral_to_peak=4, min_peak_value=10,
                                 ts_diff_min=2, ts_diff_max=5):
        """Applies MCP clean-pulse + timing-difference selection with cutflow tracking."""
        if "mcp_diff_selection_applied" in self._steps_applied:
            return self

        self.rdf = (self._get_or_create_sel_mgr()
                    .apply_mcp_diff_selection(flag_only=flag_only,
                                              min_integral_to_peak=min_integral_to_peak,
                                              min_peak_value=min_peak_value,
                                              ts_diff_min=ts_diff_min,
                                              ts_diff_max=ts_diff_max)
                    .get_rdf())

        self._steps_applied.add("mcp_diff_selection_applied")
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

    def get_particle_type(self, override=None):
        """Return the effective particle type: override if given, else beam-type default."""
        return override or _BEAM_TO_PARTICLE.get((self.beam_type or "").lower())

    def set_output_subdir(self, subdir):
        """Append a subdirectory to all output paths (root, plots, html)
        and reinitialise the hbook. Useful for separating particle selections."""
        for key in ("root", "plots", "html"):
            self.paths[key] = f"{self.paths[key]}/{subdir}"
        self.hbook = HistBook(self.paths["root"])
        return self

    def report_all(self):
        """Prints the global cutflow (hole veto) and all particle branches."""
        print(f"\n{'#'*72}\n{'COMBINED ANALYSIS SUMMARY':^72}\n{'#'*72}")

        if hasattr(self, 'sel_mgr'):
            print("\n>>> GLOBAL SELECTIONS (Applied to all branches)")
            self.sel_mgr.print_cutflow()

        for name, mgr in self.branches.items():
            print(f"\n>>> PARTICLE BRANCH: {name.upper()}")
            mgr.print_cutflow()

    @property
    def selection_summary(self) -> str:
        """Markdown-formatted selections and cutflow for the HTML intro block."""
        if not hasattr(self, 'sel_mgr'):
            return ""
        lines = []
        if self.sel_mgr.selection_descriptions:
            lines.append("**Selections applied:**")
            lines.extend(f"- {d}" for d in self.sel_mgr.selection_descriptions)
        try:
            data = self.sel_mgr.get_cutflow_dict()
            W = 72
            sep  = '=' * W
            dash = '-' * W
            title = f"Selection Cutflow (Run {self.run_number})"
            header = f"{'Step':<30} | {'Count':>10} | {'Total Eff.':>11} | {'Step Eff.':>11}"
            pre_lines = [
                sep,
                f"{title:^{W}}",
                sep,
                header,
                dash,
                f"{'Initial Events':<30} | {data['initial']:>10,} | {'100.0%':>11} | {'-':>11}",
            ]
            for s in data["steps"]:
                pre_lines.append(
                    f"{s['label']:<30} | {s['count']:>10,} | "
                    f"{s['total_eff']:>10.2f}% | {s['step_eff']:>10.2f}%"
                )
            pre_lines.append(sep)
            content = "&#10;".join(pre_lines)
            lines.append(
                f'<pre style="background:#eef4fb;border:1px solid #a8c8f0;'
                f'border-radius:4px;padding:8px 12px;font-size:12px;'
                f'margin:8px 0;overflow-x:auto">{content}</pre>'
            )
        except Exception:
            pass
        return "\n".join(lines)

    def get_rdf(self):
        """Returns the final processed RDataFrame node for booking."""
        return self.rdf

    def run_sequences(self, sequences):
        """Run hist phase then plot phase for all sequences, using self as context."""
        from core.sequence import run_hist_phase, run_plot_phase
        run_hist_phase(sequences, self)
        run_plot_phase(sequences, self)
        return self
