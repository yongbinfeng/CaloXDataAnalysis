import os
import ROOT
from utils.utils import number2string


def get_run_paths(run_number):
    """Centralizes path management for consistency across all scripts."""
    return {
        "root": f"results/root/Run{run_number}",
        "plots": f"results/plots/Run{run_number}",
        "html": f"results/html/Run{run_number}"
    }


def get_standard_info_pave(fersboard=None, iTowerX=None, iTowerY=None, drsboard=None, extra_lines=None):
    """Generates the standard info box (TPaveText) used in histograms."""
    extra = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
    extra.SetTextAlign(11)
    extra.SetFillColorAlpha(0, 0)
    extra.SetBorderSize(0)
    extra.SetTextFont(42)
    extra.SetTextSize(0.04)

    if fersboard:
        extra.AddText(f"FERS Board: {fersboard.boardNo}")
    if drsboard:
        extra.AddText(f"DRS Board: {drsboard.boardNo}")
    if iTowerX is not None and iTowerY is not None:
        extra.AddText(f"Tower: ({iTowerX}, {iTowerY})")

    if extra_lines:
        for line in extra_lines:
            extra.AddText(line)
    return extra


def save_hists_to_file(hist_list, filename):
    """Wrapper to safely save a list of ROOT histograms to a file."""
    if not hist_list:
        return
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    outfile = ROOT.TFile(filename, "RECREATE")
    for h in hist_list:
        # Handle cases where h is a pointer from RDataFrame
        if hasattr(h, 'GetPtr'):
            h = h.GetPtr()
        h.SetDirectory(outfile)
        h.Write()
    outfile.Close()
