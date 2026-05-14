import ROOT
from utils.plot_helper import save_hists_to_file


class HistBook:
    """
    Accumulates histogram groups registered by book_* calls.

    Usage (hist phase):
        ctx.hbook.add("output.root", lazy_hists)
        ROOT.RDF.RunGraphs(ctx.hbook.lazy_objects)
        ctx.hbook.save_all()

    Usage (plot phase):
        infile = ctx.hbook.open_file("output.root")
    """

    def __init__(self, root_dir):
        self._root_dir = root_dir
        self._entries = []  # (rel_path, [lazy_hists], post_save_fn | None)

    def add(self, rel_path, hists, post_save=None):
        """Register lazy objects. rel_path=None skips the ROOT file write."""
        self._entries.append((rel_path, list(hists), post_save))

    @property
    def lazy_objects(self):
        return [h for _, hists, _ in self._entries for h in hists]

    def save_all(self):
        for rel_path, hists, post_save in self._entries:
            if rel_path is not None:
                save_hists_to_file(hists, f"{self._root_dir}/{rel_path}")
            if post_save is not None:
                post_save()

    def open_file(self, rel_path):
        path = f"{self._root_dir}/{rel_path}"
        f = ROOT.TFile(path, "READ")
        if not f or f.IsZombie():
            raise FileNotFoundError(f"ROOT file not found: {path}")
        return f
