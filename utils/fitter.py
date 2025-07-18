import ROOT
from CMSPLOTS.tdrstyle import setTDRStyle
import os
import math


def runFit(hist, outdir, outname, npe_max=3, is3mm=False):
    # --------------------------
    # 1. Define Observable
    # --------------------------
    x = ROOT.RooRealVar("x_" + outname, "Pulse Amplitude (a.u.)", 0, 1000)
    if is3mm:
        x.setRange("fitrange", 10, 700)
    else:
        x.setRange("fitrange", 20, 500)
    x.setRange("plotrange", 0, 1000)
    data_hist = ROOT.RooDataHist(
        "data_hist_" + outname, "SiPM data", ROOT.RooArgList(x), hist)

    # --------------------------
    # 2. Initial Guesses
    # --------------------------
    if is3mm:
        mu0_guess = 145.0
        dpe_guess = 140.0
        dpe_min = 115
        dpe_max = 170
        sigma_guess = 20.0
    else:
        # Default values for 6mm SiPMs
        mu0_guess = 140.0
        dpe_guess = 60.0
        dpe_min = 51
        dpe_max = 75
        sigma_guess = 15.0

    npe_max = npe_max  # Number of peaks (0 = pedestal)
    mu_pe_guess = 1.0   # Mean photon count
    p_xt_guess = 0.05   # Crosstalk probability

    print(
        f"Auto-init: mu0 ≈ {mu0_guess:.2f}, dpe ≈ {dpe_guess:.2f}, sigma ≈ {sigma_guess:.2f}, npe_max = {npe_max}")

    # --------------------------
    # 3. Build Pedestal (Crystal Ball)
    # --------------------------
    mu_ped = ROOT.RooRealVar(
        "mu_ped_" + outname, "Pedestal mean", mu0_guess, mu0_guess - 10, mu0_guess + 10)
    dpe = ROOT.RooRealVar(
        "dpe_" + outname, "1 p.e. spacing", dpe_guess, dpe_min, dpe_max)
    sigmaL = ROOT.RooRealVar(
        "sigmaL_" + outname, "Single p.e. sigma", sigma_guess, 10, 30)
    sigmaR = ROOT.RooRealVar(
        "sigmaR_" + outname, "Single p.e. sigma right", sigma_guess, 10, 30)

    alphaL = ROOT.RooRealVar("alphaL", "Left tail alpha", 1.5, 0.5, 5.0)
    nL = ROOT.RooRealVar("nL", "Left tail n", 5.0, 0.5, 50.0)
    alphaR = ROOT.RooRealVar("alphaR", "Right tail alpha", 1.5, 0.5, 5.0)
    nR = ROOT.RooRealVar("nR", "Right tail n", 5.0, 0.5, 50.0)

    # --------------------------
    # 4. Build Multi-Peak Crystal Ball Model
    # --------------------------
    pdf_list = ROOT.RooArgList()
    coef_list = ROOT.RooArgList()

    # Keep references to prevent garbage collection
    mu_list = []
    sigma_list = []
    cb_list = []
    coef_refs = []
    alpha_list = []
    n_cb_list = []

    # Poisson-based initial fractions
    poisson_weights = []
    total_weight = 0.0
    for i in range(npe_max):
        weight = (mu_pe_guess ** i / math.factorial(i)) * \
            math.exp(-mu_pe_guess) * ((1 + p_xt_guess) ** i)
        poisson_weights.append(weight)
        total_weight += weight
    poisson_weights = [w / total_weight for w in poisson_weights]

    for i in range(npe_max):
        mu_i = ROOT.RooFormulaVar(
            f"mu_{i}_{outname}", f"@0 + {i}*@1", ROOT.RooArgList(mu_ped, dpe))
        sigma_L_i = ROOT.RooFormulaVar(
            f"sigma_L_{i}_{outname}", f"sqrt({i+1})*@0", ROOT.RooArgList(sigmaL))
        sigma_R_i = ROOT.RooFormulaVar(
            f"sigma_R_{i}_{outname}", f"sqrt({i+1})*@0", ROOT.RooArgList(sigmaR))
        alpha_L_cb_i = ROOT.RooRealVar(
            f"alpha_L_cb_{i}_{outname}", "CrystalBall alpha", 1.5, 0.5, 5.0)
        alpha_R_cb_i = ROOT.RooRealVar(
            f"alpha_R_cb_{i}_{outname}", "CrystalBall alpha right", 1.5, 0.5, 5.0)
        n_L_cb_i = ROOT.RooRealVar(
            f"n_L_cb_{i}_{outname}", "CrystalBall n", 5.0, 0.01, 50)
        n_R_cb_i = ROOT.RooRealVar(
            f"n_R_cb_{i}_{outname}", "CrystalBall n right", 5.0, 0.01, 50)

        cb_i = ROOT.RooCrystalBall(
            f"cb_{i}_{outname}", f"CB {i} p.e.", x, mu_i, sigma_L_i, sigma_R_i,
            alphaL, nL, alphaR, nR)

        mu_list.append(mu_i)
        sigma_list.append(sigma_L_i)
        sigma_list.append(sigma_R_i)
        cb_list.append(cb_i)
        alpha_list.append(alpha_L_cb_i)
        alpha_list.append(alpha_R_cb_i)
        n_cb_list.append(n_R_cb_i)
        n_cb_list.append(n_L_cb_i)

        pdf_list.add(cb_i)

        if i < npe_max - 1:  # Coefficients only for N-1 peaks
            coef_i = ROOT.RooRealVar(
                f"coef_{i}_{outname}", f"Fraction {i}", poisson_weights[i], 0.0, 1.0
            )
            # coef_i.setConstant(True)  # Fix coefficients for now
            coef_list.add(coef_i)
            coef_refs.append(coef_i)

    total_pdf = ROOT.RooAddPdf(
        "total_pdf_" + outname, "SiPM CB model", pdf_list, coef_list, True
    )

    # --------------------------
    # 5. Background
    # --------------------------
    # tau = ROOT.RooRealVar(
    #    "tau_" + outname, "Background slope", -0.01, -2.0, -1e-3)
    # bkg = ROOT.RooExponential("bkg_" + outname, "Background", x, tau)
    # frac_bkg = ROOT.RooRealVar(
    #    "frac_bkg_" + outname, "Background fraction", 0.0, 0.0, 1.0)
    # frac_bkg.setConstant(True)  # Fix background fraction for now

    # final_pdf = ROOT.RooAddPdf(
    #    "final_pdf_" + outname, "Peaks + Background",
    #    ROOT.RooArgList(bkg, total_pdf),
    #    ROOT.RooArgList(frac_bkg)
    # )
    final_pdf = total_pdf  # No background for now

    # --------------------------
    # 6. Fit
    # --------------------------
    fit_result = final_pdf.fitTo(
        data_hist, ROOT.RooFit.Save(), ROOT.RooFit.Range("fitrange"), ROOT.RooFit.Strategy(0))
    fit_result = final_pdf.fitTo(
        data_hist, ROOT.RooFit.Save(), ROOT.RooFit.Range("fitrange"), ROOT.RooFit.Minimizer(
            "Minuit2", "scan"), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1)
    )
    fit_result = final_pdf.fitTo(
        data_hist, ROOT.RooFit.Save(), ROOT.RooFit.Range(
            "fitrange"),
        ROOT.RooFit.Minimizer("Minuit2", "improve"), ROOT.RooFit.Strategy(
            1), ROOT.RooFit.PrintLevel(-1)
    )
    fit_result = final_pdf.fitTo(
        data_hist, ROOT.RooFit.Save(), ROOT.RooFit.Range(
            "fitrange"),
        ROOT.RooFit.Minimizer("Minuit2", "minimize"), ROOT.RooFit.Strategy(
            1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save()
    )
    fit_result.Print()

    # --------------------------
    # 7. Plot
    # --------------------------
    setTDRStyle()
    frame = x.frame()
    frame.SetTitle("")
    frame.GetXaxis().SetTitle("ADC")
    frame.GetYaxis().SetTitle("Events")
    data_hist.plotOn(frame, ROOT.RooFit.Name(
        "data_hist_" + outname), ROOT.RooFit.Range("plotrange"))
    final_pdf.plotOn(frame, ROOT.RooFit.Name("final_pdf_" + outname),
                     ROOT.RooFit.Range("plotrange"), ROOT.RooFit.NormRange("fitrange"))

    # Plot background
    # final_pdf.plotOn(frame,
    #                 ROOT.RooFit.Components(f"bkg_{outname}"),
    #                 ROOT.RooFit.LineStyle(ROOT.kDashed),
    #                 ROOT.RooFit.LineColor(ROOT.kRed),
    #                 ROOT.RooFit.Name(f"bkg_{outname}"),
    #                 ROOT.RooFit.Range("plotrange"),
    #                 ROOT.RooFit.NormRange("fitrange"))

    # Plot pedestal and peaks
    for i in range(npe_max):
        style = ROOT.kDotted if i == 0 else ROOT.kDashed
        color = ROOT.kGreen if i == 0 else 6
        final_pdf.plotOn(frame,
                         ROOT.RooFit.Components(f"cb_{i}_{outname}"),
                         ROOT.RooFit.LineStyle(style),
                         ROOT.RooFit.LineColor(color),
                         ROOT.RooFit.Name(f"cb_{i}_{outname}"),
                         ROOT.RooFit.Range("plotrange"),
                         ROOT.RooFit.NormRange("fitrange"))

    # --------------------------
    # 8. Legend
    # --------------------------
    legend = ROOT.TLegend(0.6, 0.70, 0.9, 0.9)
    legend.AddEntry(frame.findObject("data_hist_" + outname), "Data", "ep")
    legend.AddEntry(frame.findObject("final_pdf_" + outname), "Fit", "l")
    # legend.AddEntry(frame.findObject(f"bkg_{outname}"), "Background", "l")
    legend.AddEntry(frame.findObject(f"cb_0_{outname}"), "Pedestal CB", "l")
    if npe_max > 1:
        legend.AddEntry(frame.findObject(
            f"cb_1_{outname}"), f"p.e. peak", "l")
    frame.addObject(legend)

    # --------------------------
    # 9. Save Canvas
    # --------------------------
    c = ROOT.TCanvas("c_" + outname, "SiPM Fit", 800, 600)
    c.SetLogy()
    frame.Draw()

    # Draw TLatex annotations
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    latex.DrawLatexNDC(0.60, 0.60, "Pedestal = %.2f #pm %.2f" %
                       (mu_ped.getVal(), mu_ped.getError()))
    latex.DrawLatexNDC(0.60, 0.55, "Gain = %.2f #pm %.2f" %
                       (dpe.getVal(), dpe.getError()))

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outputname = f"{outdir}/{outname}.png"
    c.SaveAs(outputname)

    # Return references to avoid GC
    return mu_ped.getVal(), dpe.getVal()
