class SelectionManager:
    """
    A class to manage RDataFrame selections for CaloX data analysis.
    Maintains a list of nodes to prevent PyROOT memory issues.
    """

    def __init__(self, rdf, run_number):
        self.run_number = run_number
        # Store nodes in a list to prevent garbage collection
        self._nodes = [rdf]
        self.initial_count = rdf.Count().GetValue()
        self.current_count = self.initial_count

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

        # Create the new filtered node
        new_rdf = self.rdf.Filter(final_filter, final_label)

        # Trigger an action to print cut-flow and keep node alive
        new_count = new_rdf.Count().GetValue()
        print(
            f"Filter [{final_label}]: {self.current_count} -> {new_count} events.")

        # Append to node list to maintain reference
        self._nodes.append(new_rdf)
        self.current_count = new_count

    def veto_muon_counter(self, TSmin=400, TSmax=600, cut=-80, invert=False):
        """
        Vetoes events with muon counter signal.
        Set invert=True to keep ONLY events with a muon signal.
        """
        from channels.channel_map import getDownStreamTTUMuonChannel
        muon_channel = getDownStreamTTUMuonChannel(run=self.run_number)

        if muon_channel is None:
            print("Muon counter channel not found, skipping veto.")
            return self

        # Definition remains the same
        rdf_def = self.rdf.Define(
            "MuonCounterMin", f"MinRange({muon_channel}_blsub, {TSmin}, {TSmax})")
        self._nodes.append(rdf_def)

        # Apply filter with the inversion option
        self._apply_filter(
            f"MuonCounterMin >= {cut}", "MuonVeto", invert=invert)
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

        rdf_def = self.rdf.Define(
            f"{preshower_channel}_sum", f"SumRange({preshower_channel}_blsub, 100, 400)")
        self._nodes.append(rdf_def)

        val_cut = -1000.0
        # Standard logic: Electrons (is_hadron=False) pass if sum < cut
        condition = f"({preshower_channel}_sum < {val_cut}) == {0 if is_hadron else 1}"

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

        rdf_def = self.rdf.Define(
            f"{chan_upveto}_peak_value", f"ROOT::VecOps::Min({chan_upveto}_blsub)")
        self._nodes.append(rdf_def)
        self._apply_filter(f"{chan_upveto}_peak_value > {cut}",
                           "UpstreamVeto", invert=invert)
        return self

    def get_rdf(self):
        """Returns the final filtered RDataFrame node."""
        return self.rdf
