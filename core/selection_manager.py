from channels.channel_map import get_service_drs_channels
from configs.selection_values import get_service_drs_cut


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

    def _register_stats(self, label, condition, label_only):
        """Tracks selection statistics regardless of filter/label mode."""
        if label_only:
            # If we aren't filtering, we count how many 'true' values exist in the label column
            self.stats_proxies.append((label, self.rdf.Sum(label)))
        else:
            # If we are filtering, the count of the resulting node is the statistic
            self.stats_proxies.append((label, self.rdf.Count()))

    def apply_muon_counter_veto(self, label_only=False):
        muon_channel = get_service_drs_channels(
            self.run_number).get("TTUMuonVeto")
        if not muon_channel:
            raise ValueError(
                "Muon counter channel not found for this run. Can not apply muon counter veto.")

        ts_min, ts_max, val_cut, method = get_service_drs_cut("TTUMuonVeto")

        calc_col = "MuonCounterMin"
        self._apply_define(
            calc_col, f"{method}Range({muon_channel}_blsub, {ts_min}, {ts_max})")
        label = "passMuonCounterVeto"
        condition = f"{calc_col} >= {val_cut}"
        self._apply_define(label, condition)

        if not label_only:
            self._update_chain(self.rdf.Filter(condition, f"{label}_Filter"))

        self._register_stats(label, condition, label_only)
        return self

    def apply_psd_selection(self, is_hadron=False, label_only=False):
        psd_chan = get_service_drs_channels(self.run_number).get("PSD")
        if not psd_chan:
            raise ValueError(
                "PSD channel not found for this run. Can not apply PSD selection.")

        ts_min, ts_max, value_cut, method = get_service_drs_cut("PSD")

        calc_col = f"{psd_chan}_sum"
        self._apply_define(
            calc_col, f"{method}Range({psd_chan}_blsub, {ts_min}, {ts_max})")

        label = "passPSDSelection"
        condition = f"{calc_col} < {value_cut}" if not is_hadron else f"{calc_col} >= {value_cut}"
        self._apply_define(label, condition)

        if not label_only:
            self._update_chain(self.rdf.Filter(condition, f"{label}_Filter"))

        self._register_stats(label, condition, label_only)
        return self

    def apply_hole_veto(self, label_only=False):
        chan_holeveto = get_service_drs_channels(
            self.run_number).get("HoleVeto")
        if not chan_holeveto:
            raise ValueError(
                "Hole veto channel not found for this run. Can not apply hole veto.")

        ts_min, ts_max, value_cut, method = get_service_drs_cut("HoleVeto")

        calc_col = f"{chan_holeveto}_peak_value"
        self._apply_define(
            calc_col, f"{method}Range({chan_holeveto}_blsub, {ts_min}, {ts_max})")

        # 2. Boolean Label
        label = "passHoleVeto"
        condition = f"{calc_col} > {value_cut}"
        self._apply_define(label, condition)

        # 3. Action
        if not label_only:
            self._update_chain(self.rdf.Filter(condition, f"{label}_Filter"))

        self._register_stats(label, condition, label_only)
        return self

    def get_rdf(self):
        """Returns the final RDataFrame node."""
        return self.rdf

    def print_cutflow(self):
        """Triggers the event loop and prints results."""
        print(f"\n{'='*50}")
        print(f"{'Selection Cutflow (Run ' + str(self.run_number) + ')':^50}")
        print(f"{'='*50}")

        initial = self.initial_count_proxy.GetValue()
        print(f"{'Initial Events':<35} | {initial:>10}")
        print(f"{'-'*50}")

        for label, proxy in self.stats_proxies:
            # proxy.GetValue() returns either Count or Sum result
            count = int(proxy.GetValue())
            print(f"{label:<35} | {count:>10}")

        print(f"{'='*50}\n")
