"""Auto-compile C++ functions if needed."""
from pathlib import Path
import ROOT
import os


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


ensure_compiled()
