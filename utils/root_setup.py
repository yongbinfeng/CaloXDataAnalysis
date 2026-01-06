"""Auto-compile C++ functions if needed."""
from pathlib import Path
import ROOT


def ensure_compiled():
    utils_dir = Path(__file__).parent
    cc_file = utils_dir / "functions.cc"
    so_file = utils_dir / "functions_cc.so"

    if not so_file.exists() or cc_file.stat().st_mtime > so_file.stat().st_mtime:
        print(f"Compiling {cc_file}...")
        ROOT.gROOT.ProcessLine(f'.L {cc_file}+')
    else:
        print(f"Loading pre-compiled {so_file}...")
        ROOT.gSystem.Load(str(so_file))


def setup_root(n_threads: int = 10, batch_mode: bool = True, load_functions: bool = True):
    """Initialize ROOT with standard settings.
    """
    if batch_mode:
        print("Setting ROOT to batch mode.")
        ROOT.gROOT.SetBatch(True)
    if n_threads > 1:
        print(
            f"Enabling ROOT implicit multi-threading with {n_threads} threads.")
        ROOT.ROOT.EnableImplicitMT(n_threads)
    # Prevent histograms from auto-registering in gDirectory
    ROOT.TH1.AddDirectory(False)
    if load_functions:
        ensure_compiled()
