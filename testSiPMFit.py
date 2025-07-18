from math import exp, factorial
import ROOT


def estimate_number_of_peaks(hist, threshold_fraction=0.05):
    """
    Estimate the number of visible photoelectron peaks in the histogram.
    """
    npeaks = 0
    maxval = hist.GetMaximum()
    threshold = maxval * threshold_fraction
    for b in range(2, hist.GetNbinsX()-1):
        y0 = hist.GetBinContent(b-1)
        y1 = hist.GetBinContent(b)
        y2 = hist.GetBinContent(b+1)
        if y1 > y0 and y1 > y2 and y1 > threshold:
            npeaks += 1
    return max(1, min(npeaks, 15))  # reasonable limits


def poisson_val(n, lam_val):
    return (lam_val**n) * exp(-lam_val) / factorial(n)


def runFit(h, npe_max=5):
    # RooFit variable and data
    x = ROOT.RooRealVar("x", "Integrated Pulse [a.u.]", 0, 1000)
    dh = ROOT.RooDataHist("dh", "Data", ROOT.RooArgList(x), h)

    pedestal_mean = ROOT.RooRealVar(
        "pedestal_mean", "Pedestal mean", 150.0, 0, 200)
    pedestal_sigma = ROOT.RooRealVar(
        "pedestal_sigma", "Pedestal sigma", 2.0, 0.1, 10.0)
    gain = ROOT.RooRealVar("gain", "Gain", 70.0, 30.0, 150.0)
    sigma = ROOT.RooRealVar("sigma", "PE sigma", 30.0, 10.0, 50.0)
    lam = ROOT.RooRealVar("lambda", "Mean PE", 2.0, 0.1, 10.0)
    xtalk_prob = ROOT.RooRealVar(
        "xtalk_prob", "Crosstalk prob", 0.05, 0.0, 0.4)
    tail_slope = ROOT.RooRealVar(
        "tail_slope", "Tail slope", -1.0 / 5.0, -10.0, -0.01)
    f_tail = ROOT.RooRealVar("f_tail", "Tail fraction", 0.2, 0.0, 0.5)

    # Components for pedestal
    exp_tail = ROOT.RooExponential("exp_tail", "Exp Tail", x, tail_slope)
    pedestal_gaus = ROOT.RooGaussian(
        "pedestal_gaus", "Pedestal Gaussian", x, pedestal_mean, pedestal_sigma)
    pedestal_pdf = ROOT.RooAddPdf("pedestal_pdf", "Pedestal",
                                  ROOT.RooArgList(pedestal_gaus, exp_tail),
                                  ROOT.RooArgList(f_tail))

    lam_eff_init = lam.getVal() / (1.0 - xtalk_prob.getVal())
    raw_weights = [poisson_val(n, lam_eff_init) for n in range(npe_max + 1)]
    norm = sum(raw_weights)
    weights = [w / norm for w in raw_weights]

    # PDFs and coefficients (constant RooRealVars)
    pdfs = [pedestal_pdf]
    means = []
    coeffs = []
    w0 = ROOT.RooRealVar("w_0", "Weight for 0 PE", weights[0])
    w0.setConstant(True)
    coeffs.append(w0)

    for n in range(1, npe_max + 1):
        mean_n = ROOT.RooFormulaVar(
            f"mean_{n}", f"pedestal_mean + {n} * gain", ROOT.RooArgList(pedestal_mean, gain))
        gaus_n = ROOT.RooGaussian(
            f"gaus_{n}", f"{n} PE Gaussian", x, mean_n, sigma)
        means.append(mean_n)
        pdfs.append(gaus_n)

        wn = ROOT.RooRealVar(f"w_{n}", f"Weight for {n} PE", weights[n])
        wn.setConstant(True)
        coeffs.append(wn)

    pdf_list = ROOT.RooArgList()
    coef_list = ROOT.RooArgList()

    for pdf in pdfs:
        pdf_list.add(pdf)

    for coef in coeffs:
        coef_list.add(coef)

    total_pdf = ROOT.RooAddPdf("total_pdf", "total pdf", pdf_list, coef_list)

    # Afterpulse component
    afterpulse_mean = ROOT.RooFormulaVar(
        "afterpulse_mean", "pedestal_mean + 0.5 * gain", ROOT.RooArgList(pedestal_mean, gain))
    afterpulse_sigma = ROOT.RooRealVar(
        "afterpulse_sigma", "Afterpulse sigma", 3.0, 0.5, 10.0)
    afterpulse_frac = ROOT.RooRealVar(
        "afterpulse_frac", "Afterpulse fraction", 0.02, 0.0, 0.3)
    afterpulse_gaus = ROOT.RooGaussian(
        "afterpulse_gaus", "Afterpulse", x, afterpulse_mean, afterpulse_sigma)

    model = ROOT.RooAddPdf("full_model", "Model + Afterpulse",
                           ROOT.RooArgList(total_pdf, afterpulse_gaus),
                           ROOT.RooArgList(afterpulse_frac))

    params = [pedestal_mean, pedestal_sigma, gain, sigma, lam,
              xtalk_prob, f_tail, tail_slope, afterpulse_frac]

    model.fitTo(dh, ROOT.RooFit.Save())

    # Plot
    c = ROOT.TCanvas("c", "Fit", 800, 600)
    c.SetLogy()
    frame = x.frame(ROOT.RooFit.Title(
        "SiPM Multi-PE Fit with Crosstalk and Afterpulse"))
    dh.plotOn(frame)
    model.plotOn(frame)
    model.plotOn(frame, ROOT.RooFit.Components(
        "pedestal_pdf"), ROOT.RooFit.LineStyle(2))
    model.plotOn(frame, ROOT.RooFit.Components(
        "afterpulse_gaus"), ROOT.RooFit.LineStyle(3))
    frame.Draw()
    c.SaveAs("sipm_fit_with_xtalk_afterpulse.png")

    print("\nFit Results:")
    for p in params:
        print(f"{p.GetName():<20}: {p.getVal():.4f} ± {p.getError():.4f}")


if __name__ == "__main__":
    infile = "results/root/Run1033/fers_all_channels_1D.root"
    ifile = ROOT.TFile(infile)
    hname = "hist_FERS_Board11_Sci_1p5_m0p875"
    h = ifile.Get(hname)
    npe_estimate = estimate_number_of_peaks(h, threshold_fraction=0.1)
    print(f"\nEstimated number of visible peaks: {npe_estimate}")
    runFit(h, 3)
