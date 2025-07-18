import ROOT
from utils.html_generator import generate_html
from utils.fitter import runFit

ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing


if __name__ == "__main__":

    infile = "results/root/Run1040/fers_all_channels_1D.root"
    outdir = "results/plots/Run1040/fit_plots"
    ifile = ROOT.TFile(infile)
    hists = ifile.GetListOfKeys()
    hists = [h.GetName() for h in hists if "hist_FERS_Board11" in h.GetName()
             or "hist_FERS_Board3" in h.GetName()]
    # hists = [h for h in hists if "hist_FERS_Board11_Sci_1p5_0p625" in h]

    plots = []
    for hname in hists:
        print(f"Processing histogram: {hname}")
        h = ifile.Get(hname).Clone(f"{hname}_clone")
        runFit(h, outdir, hname)
        plots.append(hname + ".png")

    generate_html(plots, outdir, plots_per_row=4,
                  output_html="results/html/Run1040/sipm_fit_results.html")
