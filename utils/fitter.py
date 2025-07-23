import ROOT
from CMSPLOTS.tdrstyle import setTDRStyle
import CMSPLOTS.CMS_lumi as CMS_lumi
import os
import math


def channelFit(hist, outdir, outname, npe_max=3, is3mm=False, runNumber=None):
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

    CMS_lumi.lumi_sqrtS = ""
    CMS_lumi.relPosX = 0.25
    # CMS_lumi.extraText = "Internal"
    CMS_lumi.extraText = ""

    if runNumber is not None:
        CMS_lumi.lumi_13TeV = f"Cosmic Run {runNumber}"

    CMS_lumi.CMS_lumi(c, 4, 0)

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


def eventFit(h, suffix, outdir="plots/fits", addMIP=False, addHE=False, xlabel="Energy [ADCCount]",
             xmin=0, xmax=1000,
             xfitmin=0, xfitmax=450,
             xgausmean=6800, xgausmin=4000, xgausmax=7100,
             wgausmean=2500, wgausmin=500, wgausmax=2500,
             xmipmean=4000, xmipmin=2000, xmipmax=6000,
             wmipmean=1000, wmipmin=100, wmipmax=2000,
             runNumber=None):

    xfitmin = max(xfitmin, xmin)
    xfitmax = min(xfitmax, xmax)

    var = ROOT.RooRealVar("energy_" + suffix, "energy", xmin, xmax, "ADCCount")
    var.setRange("fitrange", xfitmin, xfitmax)
    var.setRange("plotrange", xmin, xmax)

    n_all = h.Integral(0, h.GetNbinsX() + 1)
    # get the number of events in the plot range
    n_events = h.Integral(h.FindBin(xmin), h.FindBin(xmax))
    n_gaus = 0
    n_exp = 0
    n_mip = 0

    datahist = ROOT.RooDataHist("datahist_" + suffix, "datahist",
                                ROOT.RooArgList(var, "argdatahist"), h)

    vmean = ROOT.RooRealVar("vmean_" + suffix, "vmean",
                            xgausmean, xgausmin, xgausmax, "ADCCount")
    vsigma = ROOT.RooRealVar("vsigma_" + suffix, "vsigma",
                             wgausmean, wgausmin, wgausmax, "ADCCount")
    pdf_Gaussian = ROOT.RooGaussian(
        "pdf_gaus_" + suffix, "pdf", var, vmean, vsigma)
    if addMIP:
        frac_gaus = ROOT.RooRealVar(
            "frac_gaus_" + suffix, "frac_gaus", 0.45, 0.2, 0.7)
    elif addHE:
        frac_gaus = ROOT.RooRealVar(
            "frac_gaus_" + suffix, "frac_gaus", 0.99, 0.95, 0.1)

    if addHE:
        vexp = ROOT.RooRealVar(
            "vexp_" + suffix, "vexp", -0.00001, -0.01, 0.0, "1/ADCCount")
        frac_exp = ROOT.RooRealVar(
            "frac_exp_" + suffix, "frac_exp", 0.01, 0.001, 0.1)
        pdf_Exp = ROOT.RooExponential(
            "pdf_exp_" + suffix, "pdf_exp", var, vexp)

    if addMIP:
        vmean_mip = ROOT.RooRealVar("vmean_mip_" + suffix, "vmean_mip",
                                    xmipmean, xmipmin, xmipmax, "ADCCount")
        vwidth_mip = ROOT.RooRealVar("vwidth_mip_" + suffix, "vwidth_mip",
                                     wmipmean, wmipmin, wmipmax, "ADCCount")
        pdf_mip = ROOT.RooGaussian(
            "pdf_mip_" + suffix, "pdf_mip", var, vmean_mip, vwidth_mip)
        frac_mip = ROOT.RooRealVar(
            "frac_mip_" + suffix, "frac_mip", 0.45, 0.2, 0.7)

    if not addMIP and not addHE:
        final_pdf = pdf_Gaussian
    else:
        pdfs = []
        fracs = []
        pdfs.append(pdf_Gaussian)
        fracs.append(frac_gaus)

        if addHE:
            pdfs.append(pdf_Exp)
            fracs.append(frac_exp)

        if addMIP:
            pdfs.append(pdf_mip)
            fracs.append(frac_mip)

        # remove the last fraction
        fracs.pop()

        pdfsArg = ROOT.RooArgList()
        for pdf in pdfs:
            pdfsArg.add(pdf)

        fracsArg = ROOT.RooArgList()
        for frac in fracs:
            fracsArg.add(frac)

        final_pdf = ROOT.RooAddPdf(
            "final_pdf_" + suffix, "pdf", pdfsArg, fracsArg)

    # run the fit
    fit_result = final_pdf.fitTo(
        datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range("fitrange"), ROOT.RooFit.Strategy(0))
    fit_result = final_pdf.fitTo(
        datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range("fitrange"), ROOT.RooFit.Minimizer(
            "Minuit2", "scan"), ROOT.RooFit.Strategy(2), ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1)
    )
    fit_result = final_pdf.fitTo(
        datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range(
            "fitrange"),
        ROOT.RooFit.Minimizer("Minuit2", "improve"), ROOT.RooFit.Strategy(
            1), ROOT.RooFit.PrintLevel(-1)
    )
    fit_result = final_pdf.fitTo(
        datahist, ROOT.RooFit.Save(), ROOT.RooFit.Range(
            "fitrange"),
        ROOT.RooFit.Minimizer("Minuit2", "minimize"), ROOT.RooFit.Strategy(
            1), ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.Save()
    )
    fit_result.Print()

    # make plots
    tsize = 0.045

    ROOT.gStyle.SetPadLeftMargin(0.14)
    ROOT.gStyle.SetPadRightMargin(0.04)
    c = ROOT.TCanvas("c_" + suffix, "c", 800, 1000)
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.05)
    padsize1 = 0.67
    padsize2 = 0.33
    paddown = ROOT.TPad('bottompad_' + suffix,
                        'bottompad', 0.0, 0.0, 1, padsize2)
    padtop = ROOT.TPad('toppad_' + suffix, 'toppad', 0.0, padsize2, 1, 1)
    padtop.SetTopMargin(0.05*0.667/(padsize1*padsize1))
    padtop.SetBottomMargin(0.012/padsize1)
    paddown.SetTopMargin(0.010/padsize2)
    paddown.SetBottomMargin(0.13/padsize2)
    paddown.SetGridy(1)
    c.cd()
    paddown.Draw()
    padtop.Draw()

    padtop.cd()
    padtop.SetLogy()
    frame = var.frame()
    frame.GetXaxis().SetTitle("energy [ADCCount]")
    datahist.plotOn(frame, ROOT.RooFit.Name("datahist_" + suffix))
    # Plot individual components by exact internal names

    # Exponential component (if included)
    if addHE:
        final_pdf.plotOn(frame,
                         ROOT.RooFit.Components(f"pdf_exp_{suffix}"),
                         ROOT.RooFit.LineColor(ROOT.kGreen+2),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFIt.Name(f"pdf_exp_{suffix}"),
                         ROOT.RooFit.Range("plotrange"),
                         ROOT.RooFit.NormRange("fitrange"))

    # Landau component (if included)
    if addMIP:
        final_pdf.plotOn(frame,
                         ROOT.RooFit.Components(f"pdf_mip_{suffix}"),
                         ROOT.RooFit.LineColor(6),
                         ROOT.RooFit.LineStyle(ROOT.kDashed),
                         ROOT.RooFit.Name(f"pdf_mip_{suffix}"),
                         ROOT.RooFit.Range("plotrange"),
                         ROOT.RooFit.NormRange("fitrange"))

    # Gaussian component (always present)
    final_pdf.plotOn(frame,
                     ROOT.RooFit.Components(
                         f"pdf_gaus_{suffix}"),
                     ROOT.RooFit.LineColor(2),
                     ROOT.RooFit.LineStyle(ROOT.kDashed),
                     ROOT.RooFit.Name(f"pdf_gaus_{suffix}"),
                     ROOT.RooFit.Range("plotrange"),
                     ROOT.RooFit.NormRange("fitrange"))

    # Total PDF
    final_pdf.plotOn(frame, ROOT.RooFit.Name(f"final_pdf_{suffix}"))

    legend = ROOT.TLegend(0.6, 0.70, 0.9, 0.9)
    legend.AddEntry(frame.findObject("datahist_" + suffix), "Data", "ep")
    legend.AddEntry(frame.findObject("final_pdf_" + suffix), "Fit", "l")
    legend.AddEntry(frame.findObject("pdf_gaus_" + suffix), "Dark Count", "l")
    if addHE:
        legend.AddEntry(frame.findObject(
            f"pdf_exp_{suffix}"), "HE", "l")
    if addMIP:
        legend.AddEntry(frame.findObject(
            f"pdf_mip_{suffix}"), "MIP", "l")
    legend.SetFillColor(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(tsize)
    legend.SetTextFont(42)
    frame.addObject(legend)

    frame.SetTitle("")
    frame.GetYaxis().SetTitle("Events")
    frame.GetXaxis().SetTitleSize(0.)
    frame.GetXaxis().SetLabelSize(0.)
    frame.GetYaxis().SetRangeUser(1.0, datahist.sumEntries())
    frame.GetYaxis().SetLabelSize(tsize)
    frame.GetYaxis().SetTitleSize(tsize*1.25)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.Draw()

    chi2 = frame.chiSquare()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    # latex.DrawLatexNDC(0.65, 0.80, "#chi^{2}/ndf = %.2f" % (chi2))
    ylabel = 0.65
    latex.DrawLatexNDC(0.60, ylabel, "#mu_{DC} = %.2f #pm %.2f" %
                       (vmean.getVal(), vmean.getError()))
    ylabel -= 0.05
    latex.DrawLatexNDC(0.60, ylabel, "#sigma_{DC} = %.2f #pm %.2f" %
                       (vsigma.getVal(), vsigma.getError()))
    ylabel -= 0.05
    if addHE or addMIP:
        latex.DrawLatexNDC(0.60, ylabel, "f_{DC} = %.2f #pm %.2f" %
                           (frac_gaus.getVal(), frac_gaus.getError()))
        ylabel -= 0.05
        n_gaus = frac_gaus.getVal() * n_events
    if addHE:
        latex.DrawLatexNDC(0.60, ylabel, "exp slope = %.6f #pm %.6f" %
                           (vexp.getVal(), vexp.getError()))
        ylabel -= 0.05
        n_exp = frac_exp.getVal() * n_events
    if addMIP:
        latex.DrawLatexNDC(0.60, ylabel, "#mu_{MIP} = %.2f #pm %.2f" %
                           (vmean_mip.getVal(), vmean_mip.getError()))
        ylabel -= 0.05
        latex.DrawLatexNDC(0.60, ylabel, "#sigma_{MIP} = %.2f #pm %.2f" %
                           (vwidth_mip.getVal(), vwidth_mip.getError()))
        n_mip = n_events - n_gaus - n_exp

    ylabel = 0.85
    latex.DrawLatexNDC(0.20, ylabel, f"# all = {n_all:.0f}")
    ylabel -= 0.05
    latex.DrawLatexNDC(0.20, ylabel, f"# fit range = {n_events:.0f}")
    if n_gaus > 0:
        ylabel -= 0.05
        latex.DrawLatexNDC(0.20, ylabel, f"# Dark Count = {n_gaus:.0f}")
    if n_exp > 0:
        ylabel -= 0.05
        latex.DrawLatexNDC(0.20, ylabel, f"# Shower  = {n_exp:.0f}")
    if n_mip > 0:
        ylabel -= 0.05
        latex.DrawLatexNDC(0.20, ylabel, f"# MIP = {n_mip:.0f}")

    CMS_lumi.lumi_sqrtS = ""
    CMS_lumi.relPosX = 0.25
    # CMS_lumi.extraText = "Internal"
    CMS_lumi.extraText = "Cer" if "Cer" in suffix else "Sci"

    if runNumber is not None:
        CMS_lumi.lumi_13TeV = f"Cosmic Run {runNumber}"

    CMS_lumi.CMS_lumi(padtop, 4, 0)

    # make pull histograms
    hpull = frame.pullHist()
    frame2 = var.frame()
    frame2.SetTitle("")
    frame2.addPlotable(hpull, "P")
    frame2.GetXaxis().SetTitle(xlabel)
    frame2.GetYaxis().SetTitle("Pull")
    frame2.GetYaxis().CenterTitle()
    frame2.GetYaxis().SetNdivisions(505)
    frame2.GetYaxis().SetLabelSize(tsize*2)
    frame2.GetYaxis().SetTitleSize(tsize*2.5)
    frame2.GetYaxis().SetTitleOffset(0.5)
    frame2.GetXaxis().SetTitleOffset(1.2)
    frame2.GetXaxis().SetLabelSize(tsize*2)
    frame2.GetXaxis().SetTitleSize(tsize*2.5)
    frame2.GetYaxis().SetRangeUser(-5, 5)

    paddown.cd()
    frame2.Draw()

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs(f"{outdir}/fit_{suffix}.pdf")
    c.SaveAs(f"{outdir}/fit_{suffix}.png")

    print("nevents: ", h.Integral(0, h.GetNbinsX()+1))

    return f"fit_{suffix}.png"
