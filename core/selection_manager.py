from channels.channel_map import get_service_drs_channels
from configs.selection_config import get_service_drs_cut, get_particle_selection


class SelectionManager:
    """
    Manages RDataFrame selections for CaloX data analysis.
    Supports physical filtering and boolean labeling for ML or N-1 studies.
    """

    def __init__(self, rdf, run_number):
        self.run_number = run_number
        self._nodes = [rdf]
        self._defined_columns = set(str(c) for c in rdf.GetColumnNames())

        # Stats tracking
        self.initial_count_proxy = rdf.Count()
        self.stats_proxies = []  # List of (label, result_proxy)
        self._selection_descriptions: list[str] = []

        self._apply_define("passNone", "1.0")

    @property
    def rdf(self):
        """Returns the current head of the RDataFrame graph."""
        return self._nodes[-1]

    def _update_chain(self, new_node):
        """Internal helper to safely advance the RDF graph and prevent memory issues."""
        self._nodes.append(new_node)
        return new_node

    def _apply_define(self, col_name, expression):
        """Defines a column only if it does not already exist."""
        if col_name not in self._defined_columns:
            new_rdf = self.rdf.Define(col_name, expression)
            self._update_chain(new_rdf)
            self._defined_columns.add(col_name)
        return self.rdf

    def _register_stats(self, label, flag_only):
        """Tracks selection statistics regardless of filter/label mode."""
        if flag_only:
            # Sum the boolean column to count 'True' occurrences
            self.stats_proxies.append((label, self.rdf.Sum(label)))
        else:
            # Count the remaining events in the filtered graph
            self.stats_proxies.append((label, self.rdf.Count()))

    def _apply_selection(self, det, flag_only=False, apply_veto=False):
        channel = get_service_drs_channels(self.run_number).get(det)
        if not channel:
            raise ValueError(
                f"Channel for {det} not found for run {self.run_number}.")

        ts_min, ts_max, _, _, val_cut, method = get_service_drs_cut(
            det, self.run_number)

        # 1. Define Calculation
        calc_col = f"{det}_firingVal"
        self._apply_define(
            calc_col, f"{method}Range({channel}_blsub, {ts_min}, {ts_max})")

        # 2. Define Boolean Label (Is the detector firing?)
        firing_label = f"is{det}Fired"
        self._apply_define(firing_label, f"{calc_col} > {val_cut}")

        # 3. Determine Cut Logic (Selection vs Veto)
        # We create a specific 'pass' column to ensure Sum() works in flag_only mode
        pass_label = f"is_{det}_{'vetoed' if apply_veto else 'fired'}"
        cut_str = f"{firing_label} == 0" if apply_veto else f"{firing_label} == 1"
        self._apply_define(pass_label, cut_str)

        # 4. Filter or Label
        if not flag_only:
            self._update_chain(self.rdf.Filter(pass_label, f"{det}_Filter"))

        self._register_stats(pass_label, flag_only)
        return self

    @property
    def selection_descriptions(self) -> list[str]:
        return list(self._selection_descriptions)

    def apply_hole_veto(self, flag_only=False, apply_veto=True):
        suffix = " (flag only)" if flag_only else ""
        self._selection_descriptions.append(f"Hole Veto{suffix}")
        return self._apply_selection("HoleVeto", flag_only, apply_veto)

    def apply_mcp_selection(self, flag_only=False,
                            min_integral_to_peak=4, min_peak_value=10,
                            _describe=True):
        """Apply clean-pulse filter over all MCP channels with cutflow tracking."""
        from channels.channel_map import get_mcp_channels
        if _describe:
            suffix = " (flag only)" if flag_only else ""
            self._selection_descriptions.append(f"MCP clean pulse{suffix}")
        for det in get_mcp_channels(self.run_number):
            pass_label = f"is_{det}_clean"
            cut_expr = (f"{det}_integral_to_peak > {min_integral_to_peak}"
                        f" && {det}_peak_value > {min_peak_value}")
            self._apply_define(pass_label, cut_expr)
            if not flag_only:
                self._update_chain(self.rdf.Filter(pass_label, f"{det}_MCP_filter"))
            self._register_stats(pass_label, flag_only)
        return self

    def apply_mcp_diff_selection(self, flag_only=False,
                                 min_integral_to_peak=4, min_peak_value=10,
                                 ts_diff_min=2, ts_diff_max=5):
        """Apply MCP clean-pulse + timing-difference selection with cutflow tracking."""
        from channels.channel_map import get_mcp_channels
        suffix = " (flag only)" if flag_only else ""
        self._selection_descriptions.append(
            f"MCP clean pulse + timing difference [{ts_diff_min}, {ts_diff_max}] TS{suffix}")
        self.apply_mcp_selection(flag_only=flag_only, _describe=False,
                                 min_integral_to_peak=min_integral_to_peak,
                                 min_peak_value=min_peak_value)
        dets = list(get_mcp_channels(self.run_number).keys())
        det1, det2 = dets[0], dets[4]
        pass_label = f"is_{det1}_{det2}_timing_diff"
        cut_expr = (f"abs({det1}_TS_cfd_ref - {det2}_TS_cfd_ref) > {ts_diff_min}"
                    f" && abs({det1}_TS_cfd_ref - {det2}_TS_cfd_ref) < {ts_diff_max}")
        self._apply_define(pass_label, cut_expr)
        if not flag_only:
            self._update_chain(self.rdf.Filter(pass_label, "MCP_timing_diff_filter"))
        self._register_stats(pass_label, flag_only)
        return self

    def apply_particle_selection(self, particle_type, flag_only=False):
        """Applies a suite of cuts to select a specific particle type."""
        requirements = get_particle_selection(particle_type)
        if not requirements:
            raise ValueError(f"Unknown particle type: {particle_type}")

        suffix = " (flag only)" if flag_only else ""
        self._selection_descriptions.append(f"{particle_type} particle selection{suffix}")

        for detector, should_fire in requirements.items():
            # If should_fire is True, we want 'fired'.
            # If should_fire is False, we want 'vetoed' (apply_veto=True).
            apply_veto = not should_fire
            self._apply_selection(
                detector, flag_only=flag_only, apply_veto=apply_veto)

        return self

    @property
    def all_cutflow_proxies(self):
        """All lazy count/sum proxies for cutflow — include in RunGraphs."""
        return [self.initial_count_proxy] + [proxy for _, proxy in self.stats_proxies]

    def get_rdf(self):
        """Returns the final RDataFrame node."""
        return self.rdf

    def get_cutflow_dict(self):
        """Returns results as a dictionary. Triggers the event loop."""
        initial = self.initial_count_proxy.GetValue()
        stats = {"initial": initial, "steps": []}

        prev_count = initial
        for label, proxy in self.stats_proxies:
            count = int(proxy.GetValue())
            stats["steps"].append({
                "label": label,
                "count": count,
                "total_eff": (count / initial * 100) if initial > 0 else 0,
                "step_eff": (count / prev_count * 100) if prev_count > 0 else 0
            })
            prev_count = count
        return stats

    def print_cutflow(self):
        """UI wrapper for get_cutflow_dict."""
        data = self.get_cutflow_dict()
        header = f"{'Step':<30} | {'Count':>10} | {'Total Eff.':>11} | {'Step Eff.':>11}"
        print(
            f"\n{'='*72}\n{'Selection Cutflow (Run ' + str(self.run_number) + ')':^72}\n{'='*72}")
        print(header + f"\n{'-'*72}")
        print(
            f"{'Initial Events':<30} | {data['initial']:>10} | {'100.0%':>11} | {'-':>11}")

        for s in data["steps"]:
            print(
                f"{s['label']:<30} | {s['count']:>10} | {s['total_eff']:>10.2f}% | {s['step_eff']:>10.2f}%")
        print(f"{'='*72}\n")
