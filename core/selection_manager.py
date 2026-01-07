class SelectionManager:
    """
    A class to manage RDataFrame selections for CaloX data analysis.
    Maintains a list of nodes to prevent PyROOT memory issues.
    """

    def __init__(self, rdf, run_number):
        self.run_number = run_number
        # Store nodes in a list to prevent garbage collection
        self._nodes = [rdf]
        # Book the initial count lazily
        self.initial_count_proxy = rdf.Count()
        self.filter_reports = []

    @property
    def rdf(self):
        """Always returns the most recent node in the chain."""
        return self._nodes[-1]

    def _apply_filter(self, filter_string, label, invert=False):
        """
        Internal helper to apply filters and manage memory.
        Added 'invert' parameter to flip the selection logic.
        """
        # If invert is True, wrap the filter string in !( ... )
        final_filter = f"!({filter_string})" if invert else filter_string
        final_label = f"NOT_{label}" if invert else label

        new_rdf = self.rdf.Filter(final_filter, final_label)
        self._nodes.append(new_rdf)

        # Book the count for this specific filter stage lazily
        self.filter_reports.append((final_label, new_rdf.Count()))

    def _column_exists(self, name):
        """Check if a column already exists to prevent re-definition errors."""
        return name in [str(c) for c in self.rdf.GetColumnNames()]

    def veto_muon_counter(self, TSmin=200, TSmax=700, cut=-100, invert=False):
        """
        Vetoes events with muon counter signal.
        Set invert=True to keep ONLY events with a muon signal.
        """
        from channels.channel_map import getDownStreamTTUMuonChannel
        muon_channel = getDownStreamTTUMuonChannel(run=self.run_number)

        if muon_channel is None:
            print("Muon counter channel not found, skipping veto.")
            return self

        col_name = "MuonCounterMin"
        if not self._column_exists(col_name):
            rdf_def = self.rdf.Define(
                col_name, f"MinRange({muon_channel}_blsub, {TSmin}, {TSmax})")
            self._nodes.append(rdf_def)

        # Apply filter with the inversion option
        self._apply_filter(
            f"{col_name} >= {cut}", "MuonVeto", invert=invert)
        return self

    def apply_psd_selection(self, is_hadron=False, invert=False):
        """
        Applies Pre-shower Detector (PSD) selection.
        invert=True will flip the logic regardless of is_hadron setting.
        """
        from channels.channel_map import getPreShowerChannel
        preshower_channel = getPreShowerChannel(self.run_number)

        if preshower_channel is None:
            return self

        col_name = f"{preshower_channel}_sum"
        if not self._column_exists(col_name):
            rdf_def = self.rdf.Define(
                col_name, f"SumRange({preshower_channel}_blsub, 100, 400)")
            self._nodes.append(rdf_def)

        val_cut = -1000.0
        # Standard logic: Electrons (is_hadron=False) pass if sum < cut
        condition = f"({col_name} < {val_cut}) == {0 if is_hadron else 1}"

        label = "PSD_Hadron" if is_hadron else "PSD_Electron"
        self._apply_filter(condition, label, invert=invert)
        return self

    def apply_upstream_veto(self, cut=-1000.0, invert=False):
        """
        Applies upstream veto based on peak value.
        invert=True selects events that FAILED the veto (the 'halo' events).
        """
        from channels.channel_map import getUpstreamVetoChannel
        chan_upveto = getUpstreamVetoChannel(self.run_number)

        if not chan_upveto:
            return self

        col_name = f"{chan_upveto}_peak_value"
        if not self._column_exists(col_name):
            rdf_def = self.rdf.Define(
                col_name, f"ROOT::VecOps::Min({chan_upveto}_blsub)")
            self._nodes.append(rdf_def)

        self._apply_filter(f"{col_name} > {cut}",
                           "UpstreamVeto", invert=invert)
        return self

    def get_rdf(self):
        """Returns the final filtered RDataFrame node."""
        return self.rdf

    def print_cutflow(self):
        """
        Call this ONLY at the end of your script. 
        The first .GetValue() call here will trigger the ONE slow analysis loop.
        """
        print("\n--- Selection Cutflow ---")
        # Accessing the first proxy triggers the loop for all booked results
        initial = self.initial_count_proxy.GetValue()
        print(f"Initial: {initial} events")

        for label, proxy in self.filter_reports:
            # These are now ready instantly without a new loop
            print(f"{label}: {proxy.GetValue()} events")
