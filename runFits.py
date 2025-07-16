import os
import ROOT


def runFit(h, suffix, xmin=0.0, xmax=7500.0, xfitmin=2950.0, xfitmax=3400.0, outdir="plots/fits", addLandau=False, addExp=False, xlabel="Energy [ADCCount]", xgausmax=7100, wgausmax=2500):
    h.GetXaxis().SetRangeUser(xmin, xmax)
    xmean = h.GetMean()
    xrms = h.GetRMS()
    xbins = h.GetNbinsX()

    print("xmean = ", xmean)
    print("xrms = ", xrms)
    print("xbins = ", xbins)

    var = ROOT.RooRealVar("energy_" + suffix, "energy", xmin, xmax, "ADCCount")
    var.setRange("r1", xfitmin, xfitmax)
    datahist = ROOT.RooDataHist("datahist_" + suffix, "datahist",
                                ROOT.RooArgList(var, "argdatahist"), h)

    vmean = ROOT.RooRealVar("vmean_" + suffix, "vmean",
                            6800, 4000, xgausmax, "ADCCount")
    vsigma = ROOT.RooRealVar("vsigma_" + suffix, "vsigma",
                             2500, 500, wgausmax, "ADCCount")
    pdf_Gaussian = ROOT.RooGaussian("pdf_" + suffix, "pdf", var, vmean, vsigma)

    # if addExp and not addLandau:
    #    # If only Exponential is added, set a very small value for vexp
    #    vexp = ROOT.RooRealVar(
    #        "vexp_" + suffix, "vexp", -0.00005, -0.01, 0.0, "1/ADCCount")
    #    pdf_Exp = ROOT.RooExponential(
    #        "pdf_exp_" + suffix, "pdf_exp", var, vexp)
    #    # Combine the Gaussian and Exponential PDFs
    #    frac_gaus = ROOT.RooRealVar(
    #        "frac_gaus_" + suffix, "frac_gaus", 0.95, 0.90, 1.0)
    #    # If only Exponential is added, we don't need to add Landau
    #    pdf = ROOT.RooAddPdf(
    #        "pdf_" + suffix, "pdf", ROOT.RooArgList(pdf_Gaussian, pdf_Exp), ROOT.RooArgList(frac_gaus))

    # elif addLandau:
    #    frac_gaus = ROOT.RooRealVar(
    #        "frac_gaus_" + suffix, "frac_gaus", 0.5, 0.0, 1.0)
    #    vexp = ROOT.RooRealVar(
    #        "vexp_" + suffix, "vexp", -0.0000001, -0.01, 0.0, "1/ADCCount")
    #    pdf_Exp = ROOT.RooExponential(
    #        "pdf_exp_" + suffix, "pdf_exp", var, vexp)
    #    vmean_landau = ROOT.RooRealVar("vmean_landau_" + suffix, "vmean_landau",
    #                                   16000, 13000, 18000, "ADCCount")
    #    vwidth_landau = ROOT.RooRealVar("vwidth_landau_" + suffix, "vwidth_landau",
    #                                    4000, 500, 8000, "ADCCount")
    #    pdf_Landau = ROOT.RooGaussian(
    #        "pdf_landau_" + suffix, "pdf_landau", var, vmean_landau, vwidth_landau)

    #    frac_landau = ROOT.RooRealVar(
    #        "frac_landau_" + suffix, "frac_landau", 0.5, 0.2, 0.8)
    #    pdf = ROOT.RooAddPdf(
    #        "pdf_" + suffix, "pdf", ROOT.RooArgList(pdf_Landau, pdf_Gaussian, pdf_Exp), ROOT.RooArgList(frac_landau, frac_gaus))

    # else:
    #    pdf = pdf_Gaussian

    if addExp or addLandau:
        if not addLandau:
            vexp = ROOT.RooRealVar(
                "vexp_" + suffix, "vexp", -0.00001, -0.01, 0.0, "1/ADCCount")
        else:
            vexp = ROOT.RooRealVar(
                "vexp_" + suffix, "vexp", -0.00003, -0.01, 0.0, "1/ADCCount")
        pdf_Exp = ROOT.RooExponential(
            "pdf_exp_" + suffix, "pdf_exp", var, vexp)
        if not addLandau:
            frac_gaus = ROOT.RooRealVar(
                "frac_gaus_" + suffix, "frac_gaus", 0.99, 0.95, 0.1)
        elif addLandau:
            frac_gaus = ROOT.RooRealVar(
                "frac_gaus_" + suffix, "frac_gaus", 0.45, 0.2, 0.7)
        # combine the Gaussian and Exponential PDFs
        pdf_temp = ROOT.RooAddPdf(
            "pdf_temp_" + suffix, "pdf", ROOT.RooArgList(pdf_Gaussian, pdf_Exp), ROOT.RooArgList(frac_gaus))
    else:
        pdf_temp = pdf_Gaussian

    if addLandau:
        vmean_landau = ROOT.RooRealVar("vmean_landau_" + suffix, "vmean_landau",
                                       16000, 13000, 18000, "ADCCount")
        vwidth_landau = ROOT.RooRealVar("vwidth_landau_" + suffix, "vwidth_landau",
                                        4000, 500, 8000, "ADCCount")
        pdf_Landau = ROOT.RooGaussian(
            "pdf_landau_" + suffix, "pdf_landau", var, vmean_landau, vwidth_landau)

        frac_landau = ROOT.RooRealVar(
            "frac_landau_" + suffix, "frac_landau", 0.45, 0.2, 0.6)
        pdf = ROOT.RooAddPdf(
            "pdf_" + suffix, "pdf", ROOT.RooArgList(pdf_Landau, pdf_Gaussian, pdf_Exp), ROOT.RooArgList(frac_landau, frac_gaus))
    else:
        pdf = pdf_temp

    # run the fit
    pdf.fitTo(datahist, ROOT.RooFit.Range("r1"))

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
    datahist.plotOn(frame)
    # Plot individual components by exact internal names

    # Exponential component (if included)
    if addExp:
        pdf.plotOn(frame,
                   ROOT.RooFit.Components(f"pdf_exp_{suffix}"),
                   ROOT.RooFit.LineColor(ROOT.kGreen+2),
                   ROOT.RooFit.LineStyle(ROOT.kDashed),
                   ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected))

    # Landau component (if included)
    if addLandau:
        pdf.plotOn(frame,
                   ROOT.RooFit.Components(f"pdf_landau_{suffix}"),
                   ROOT.RooFit.LineColor(6),
                   ROOT.RooFit.LineStyle(ROOT.kDashed),
                   ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected))

    # Gaussian component (always present)
    pdf.plotOn(frame,
               ROOT.RooFit.Components(
                   f"pdf_{suffix}" if not addExp else f"pdf_{suffix},pdf_exp_{suffix}"),
               ROOT.RooFit.LineColor(ROOT.kRed),
               ROOT.RooFit.LineStyle(ROOT.kDashed),
               ROOT.RooFit.Normalization(1.0, ROOT.RooAbsReal.RelativeExpected))

    # Total PDF
    pdf.plotOn(frame)

    frame.SetTitle("CaloX Cosmic Run")
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
    latex.SetTextSize(tsize)
    latex.SetTextFont(42)
    # latex.DrawLatexNDC(0.65, 0.80, "#chi^{2}/ndf = %.2f" % (chi2))
    latex.DrawLatexNDC(0.65, 0.75, "#mu = %.2f #pm %.2f" %
                       (vmean.getVal(), vmean.getError()))
    latex.DrawLatexNDC(0.65, 0.70, "#sigma = %.2f #pm %.2f" %
                       (vsigma.getVal(), vsigma.getError()))
    if addExp:
        latex.DrawLatexNDC(0.65, 0.65, "f_gaus = %.2f #pm %.2f" %
                           (frac_gaus.getVal(), frac_gaus.getError()))
        latex.DrawLatexNDC(0.65, 0.60, "exp slope = %.6f #pm %.6f" %
                           (vexp.getVal(), vexp.getError()))
    if addLandau:
        latex.DrawLatexNDC(0.65, 0.55, "f_landau = %.2f #pm %.2f" %
                           (frac_landau.getVal(), frac_landau.getError()))
        latex.DrawLatexNDC(0.65, 0.50, "landau mean = %.2f #pm %.2f" %
                           (vmean_landau.getVal(), vmean_landau.getError()))
        latex.DrawLatexNDC(0.65, 0.45, "landau width = %.2f #pm %.2f" %
                           (vwidth_landau.getVal(), vwidth_landau.getError()))

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

    # pdf.Print("t")


if __name__ == "__main__":
    runNumber = 1040
    filename = f"results/root/Run{runNumber}/fers_energy_sum_subtracted.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run prepareDQMPlots.py first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    hCer = ifile.Get("hist_FERS_CerEnergyHG_subtracted")
    hSci = ifile.Get("hist_FERS_SciEnergyHG_subtracted")

    xmax = 1e5
    xfitmax = 1e5
    xmax_sci = 1e5
    xfit_sci_max = 1e5
    runFit(hCer, f"Run{runNumber}_CerHG", xmin=0.0,
           xmax=xmax, xfitmin=0, xfitmax=xfitmax, outdir=f"results/Run{runNumber}/plots/fits", addExp=True, xlabel="Cer Energy")
    runFit(hSci, f"Run{runNumber}_SciHG", xmin=0.0,
           xmax=xmax_sci, xfitmin=0, xfitmax=xfit_sci_max, outdir=f"results/Run{runNumber}/plots/fits", addLandau=True, addExp=True, xlabel="Sci Energy")

    runNumber = 996
    filename = f"results/root/Run{runNumber}/fers_energy_sum_subtracted.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run prepareDQMPlots.py first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    hCer = ifile.Get("hist_FERS_CerEnergyHG_subtracted")
    hSci = ifile.Get("hist_FERS_SciEnergyHG_subtracted")

    xmax = 1e5
    xfitmax = 1e5
    xmax_sci = 1e5
    xfit_sci_max = 1e5
    runFit(hCer, f"Run{runNumber}_CerHG", xmin=0.0,
           xmax=xmax, xfitmin=0, xfitmax=xfitmax, outdir=f"results/Run{runNumber}/plots/fits", addExp=False, xlabel="Cer Energy", xgausmax=9000, wgausmax=3000)
    runFit(hSci, f"Run{runNumber}_SciHG", xmin=0.0,
           xmax=xmax_sci, xfitmin=0, xfitmax=xfit_sci_max, outdir=f"results/Run{runNumber}/plots/fits", addLandau=False, addExp=False, xlabel="Sci Energy", xgausmax=9000, wgausmax=3000)
