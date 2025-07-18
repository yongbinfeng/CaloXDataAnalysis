import ROOT
import os
from utils.html_generator import generate_html


def runFit(hist, outdir, outname):
    hist = ifile.Get(hname)

    x = ROOT.RooRealVar("x_" + outname, "Pulse Amplitude (a.u.)", 20, 1000)
    data_hist = ROOT.RooDataHist(
        "data_hist_" + outname, "SiPM data", ROOT.RooArgList(x), hist)

    mu0_guess = 145.0
    dpe_guess = 140.0

    sigma_guess = 20.0

    npe_max = 3

    print(
        f"Auto-init: mu0 ≈ {mu0_guess:.2f}, dpe ≈ {dpe_guess:.2f}, sigma ≈ {sigma_guess:.2f}, npe_max = {npe_max}")

    # --------------------------
    # 3. Build Multi-Gaussian Model
    # --------------------------
    mu0 = ROOT.RooRealVar("mu0_" + outname, "Pedestal mean", mu0_guess,
                          mu0_guess-20, mu0_guess+20)
    sigma = ROOT.RooRealVar(
        "sigma_" + outname, "Single p.e. sigma", sigma_guess, 10, 30)
    dpe = ROOT.RooRealVar(
        "dpe_" + outname, "1 p.e. spacing", dpe_guess, 100, 160)

    mus = []
    sigmas = []
    pdfs = []
    coefs = []
    pdf_list = ROOT.RooArgList()
    coef_list = ROOT.RooArgList()

    mu0 = ROOT.RooRealVar("mu0_" + outname, "Pedestal mean", mu0_guess,
                          mu0_guess - 10, mu0_guess + 10)
    sigma0 = ROOT.RooRealVar(
        "sigma0_" + outname, "Pedestal sigma", sigma_guess, 10, 40)
    alpha_cb = ROOT.RooRealVar(
        "alpha_cb_" + outname, "CrystalBall alpha", 1.5, 0.5, 5.0)
    n_cb = ROOT.RooRealVar("n_cb_" + outname, "CrystalBall n", 5.0, 0.01, 50)

    pedestal_cb = ROOT.RooCBShape(
        "pedestal_cb_" + outname, "Pedestal CB", x, mu0, sigma0, alpha_cb, n_cb)

    pdf_list.add(pedestal_cb)

    for i in range(1, npe_max):
        mu_i = ROOT.RooFormulaVar(
            f"mu_{i}_" + outname, f"@0 + {i}*@1", ROOT.RooArgList(mu0, dpe))
        sigma_i = ROOT.RooFormulaVar(
            f"sigma_{i}_" + outname, f"sqrt({i+1})*@0", ROOT.RooArgList(sigma))
        gauss_i = ROOT.RooGaussian(
            f"gauss_{i}_" + outname, f"Gaussian {i}", x, mu_i, sigma_i)
        mus.append(mu_i)
        sigmas.append(sigma_i)
        pdfs.append(gauss_i)
        pdf_list.add(gauss_i)

        coef_i = ROOT.RooRealVar(
            f"coef_{i}_" + outname, f"Coefficient {i}", 1.0/npe_max, 0.0, 1.0)
        coefs.append(coef_i)
        coef_list.add(coef_i)

    total_pdf = ROOT.RooAddPdf(
        "total_pdf_" + outname, "Multi-Gaussian SiPM model", pdf_list, coef_list)

    tau = ROOT.RooRealVar(
        "tau_" + outname, "Background slope", -0.01, -2.0, -1e-3)
    bkg = ROOT.RooExponential("bkg_" + outname, "Background", x, tau)
    frac_bkg = ROOT.RooRealVar(
        "frac_bkg_" + outname, "Background fraction", 0.05, 0.0, 1.0)

    final_pdf = ROOT.RooAddPdf("final_pdf_" + outname, "Peaks + Background",
                               ROOT.RooArgList(total_pdf, bkg),
                               ROOT.RooArgList(frac_bkg))

    # --------------------------
    # 4. Fit
    # --------------------------
    fit_result = final_pdf.fitTo(
        data_hist, ROOT.RooFit.Save())
    fit_result.Print()

    # --------------------------
    # 5. Plot
    # --------------------------
    frame = x.frame()
    frame.SetTitle("")
    data_hist.plotOn(frame)
    final_pdf.plotOn(frame, ROOT.RooFit.Name("final_pdf_" + outname))
    final_pdf.plotOn(frame, ROOT.RooFit.Components("bkg_" + outname),
                     ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed))
    final_pdf.plotOn(frame, ROOT.RooFit.Components("pedestal_cb_" + outname),
                     ROOT.RooFit.LineStyle(ROOT.kDotted), ROOT.RooFit.LineColor(ROOT.kGreen))
    for i in range(1, npe_max):
        final_pdf.plotOn(frame, ROOT.RooFit.Components(f"gauss_{i}_" + outname),
                         ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(6))

    c = ROOT.TCanvas("c_" + outname, "SiPM Fit", 800, 600)
    c.SetLogy()
    frame.Draw()

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outputname = f"{outdir}/{outname}.png"
    c.SaveAs(outputname)


if __name__ == "__main__":

    infile = "results/root/Run1040/fers_all_channels_1D.root"
    outdir = "results/plots/Run1040/fit_plots"
    ifile = ROOT.TFile(infile)
    hists = ifile.GetListOfKeys()
    hists = [h.GetName() for h in hists if "hist_FERS_Board11" in h.GetName()]

    plots = []
    for hname in hists:
        print(f"Processing histogram: {hname}")
        h = ifile.Get(hname).Clone(f"{hname}_clone")
        runFit(h, outdir, hname)
        plots.append(hname + ".png")

    generate_html(plots, outdir, plots_per_row=4,
                  output_html="results/html/Run1040/sipm_fit_results.html")
