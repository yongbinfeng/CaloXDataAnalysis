import ROOT
import os
import sys
# from CMSPLOTS import CMS_lumi
# from CMSPLOTS import tdrstyle
import CMS_lumi
import tdrstyle


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    RED = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    RED = '\033[1;31m'
    UNDERLINE = '\033[4m'


def DumpHist(hist):
    # dump the number of events in histogram
    print("histogram %s with %d bins" % (hist.GetName, hist.GetNbinsX()))
    for ibin in range(0, hist.GetNbinsX()+2):
        print(" %d bin with %.2f events" % (ibin, hist.GetBinContent(ibin)))


def PrepareHisto(myfile, listhistos, printlist=True):
    if printlist:
        print(listhistos)
    myhisto = myfile.Get("%s" % listhistos[0])
    for idx in range(1, len(listhistos)):
        myhisto.Add(myfile.Get("%s" % listhistos[idx]))
    return myhisto


def RebinHisto(histo, *args):
    if len(args) == 0:
        return histo
    elif len(args) == 1:
        return histo.Rebin(int(args[0]))
    else:
        return histo.Rebin(len(args[0])-1, "Rebinned_%s" % args[1], args[0])


def GetHisto(myfile, lhistos, *args):
    htemp = PrepareHisto(myfile, lhistos[0].split("+"))
    htemp_rebin = RebinHisto(htemp, *args)

    if len(lhistos) > 1:
        htemp_den = PrepareHisto(myfile, lhistos[1].split("+"))
        htemp_den_rebin = RebinHisto(htemp_den, *args)
        # deal with the case with numerator is 0 but denominator is not
        # temporary fix, use the upper error as the error
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
        for ibin in range(1, htemp_rebin.GetNbinsX()+1):
            if htemp_rebin.GetBinContent(ibin) == 0 and htemp_den_rebin.GetBinContent(ibin) != 0:
                htemp_rebin.SetBinError(
                    ibin, ROOT.Math.gamma_quantile_c((1 - 0.6827)/2.0, 1, 1))
                htemp_den_rebin.SetBinError(ibin, ROOT.Math.gamma_quantile_c(
                    (1 - 0.6827)/2.0, htemp_den_rebin.GetBinContent(ibin)+1, 1))
        # divide
        # set the uncertainty to binomial, see line 2800 of https://root.cern.ch/doc/master/TH1_8cxx_source.html
        htemp_rebin.Divide(htemp_rebin, htemp_den_rebin, 1.0, 1.0, "B")
        ##
    return htemp_rebin


def Normalize(th1, option=0):
    if option == 0:
        th1.Scale(1.0/(th1.Integral(0, th1.GetNbinsX()+1)+1e-6))
    elif option == 1:
        th1.Scale(1.0/(th1.Integral()+1e-6))


def ScaleWithWidth(th1):
    th1.Scale(1.0, "width")


def myDivide(num, den):
    # set the ratio to 1 if denominator is 0
    if den == 0:
        # print bcolors.RED, "denominator is zero. Set the ratio to 1", bcolors.ENDC
        return 0.
    else:
        return num/den


def CalculateChi(hobsclone, hexp, doNewman=False, doPearson=False, ignoreHistError=False, printchi2=True):

    if doNewman and doPearson:
        print(bcolors.RED, "can NOT do Newman and Pearson chi2 at the same time!!! Use Pearson instead", bcolors.ENDC)

    for ibin in range(1, hexp.GetNbinsX()+1):
        num = hobsclone.GetBinContent(ibin) - hexp.GetBinContent(ibin)
        if doPearson and ignoreHistError:
            den = ROOT.TMath.Sqrt(hexp.GetBinContent(ibin))
        elif doNewman and ignoreHistError:
            den = ROOT.TMath.Sqrt(hobsclone.GetBinContent(ibin))
        elif doPearson and not ignoreHistError:
            den = hexp.GetBinError(ibin)
        elif doNewman and not ignoreHistError:
            den = hobsclone.GetBinError(ibin)

        # print ibin, hobsclone.GetBinContent(ibin), hexp.GetBinContent(ibin), num, den

        hobsclone.SetBinContent(ibin, myDivide(num, den))
        hobsclone.SetBinError(ibin, 0.)

    if printchi2:
        chi2 = 0
        for ibin in range(1, hobsclone.GetNbinsX()+1):
            chi2 += ROOT.TMath.Power(hobsclone.GetBinContent(ibin), 2.0)
        print("Chi2 for %s is %.2f" % (hobsclone.GetName(), chi2))


def myRead(inputfile):
    """
    Read file, color, label and list of histograms information from the input file
    ignore lines starting with #
    """
    list_input = open("%s" % inputfile)

    files = []
    colors = []
    labels = []
    lhistos = []

    for line in list_input:
        tline = []
        tline = line.rstrip().split()

        if tline[0][0] == '#':  # ignore the line starting with #
            continue

        if len(tline) < 4 or len(tline) > 5:
            print("Input Error in line\n   %s of file %s" % (line, inputfile))
            sys.exit(0)

        files.append(tline[0])
        colors.append(int(tline[1]))
        labels.append(tline[2].replace("_", " "))

        if len(tline) == 4:
            lhistos.append([tline[3]])
        else:
            lhistos.append([tline[3], tline[4]])

    return files, colors, labels, lhistos


def getResolution(th2d, useRMS=False):
    if useRMS:
        return getRMSResolution(th2d)
    else:
        return getQuantileResolution(th2d)


def getQuantileResolution(th2d):
    # return the resolution distribution of th2d, as a format of th1d
    hqup = th2d.QuantilesX(0.84, th2d.GetName()+"_q84")
    hqdn = th2d.QuantilesX(0.16, th2d.GetName()+"_q16")

    hresol = hqup.Clone(th2d.GetName() + "_resolQuantile")
    hresol.Add(hqdn, -1.0)
    hresol.Scale(0.5)

    return hresol


def getRMSResolution(th2d):
    hprof = th2d.ProfileX(th2d.GetName()+"_profX", 1, -1, "s")
    hresol = th2d.ProjectionX(th2d.GetName()+"_resolRMS", 0, -1, "e")
    for ibin in range(1, hresol.GetNbinsX()+1):
        count = hresol.GetBinContent(ibin)
        rms = hprof.GetBinError(ibin)
        rms_err = 0
        if count > 0:
            rms_err = rms / ROOT.TMath.Sqrt(2*count)
        print("counts, ", hresol.GetBinContent(
            ibin), hprof.GetBinContent(ibin))
        print("ibin ", ibin, hresol.GetBinContent(ibin),
              hprof.GetBinContent(ibin), hprof.GetBinError(ibin))
        hresol.SetBinContent(ibin, rms)
        hresol.SetBinError(ibin, rms_err)
    return hresol


def getMedian(th2d):
    # return the median of the th2d
    hmedian = th2d.QuantilesX(0.5, th2d.GetName()+"_q50")
    return hmedian


def getMean(th2d):
    hpmean = th2d.ProfileX("_mean")
    hmean = hpmean.ProjectionX()
    return hmean


def getErrors(hprof, verbose=False):
    herr = hprof.ProjectionX("%s_Error" % hprof.GetName())
    # to clone the binning of input hprof. If there is a better way,
    # this can be changed
    if verbose:
        print("input histogram %s has %d bins" % (
            hprof.GetName(), hprof.GetNbinsX()))
    for ibin in range(1, hprof.GetNbinsX()+1):
        if verbose:
            print("%d bin, with RMS %f" % (ibin, hprof.GetBinError(ibin)))
        herr.SetBinContent(ibin, hprof.GetBinError(ibin))
        herr.SetBinError(ibin, 0)
    for ibin in range(1, hprof.GetNbinsX()+1):
        if verbose:
            print("%d bin, with Content %f" % (ibin, herr.GetBinContent(ibin)))
    return herr


def THStack2TH1(hstack, postfix=""):
    assert isinstance(hstack, ROOT.THStack), "input must be a ROOT.THStack"
    hlist = list(hstack.GetHists())
    h1 = hlist[0].Clone("{}_histogram_{}".format(hstack.GetName(), postfix))
    list(map(h1.Add, hlist[1:]))
    return h1


def AddOverflowsTH1(h1, dolastbin=True):
    assert isinstance(h1, ROOT.TH1), "input must be a ROOT.TH1"
    nbins = h1.GetNbinsX()
    if dolastbin:
        # merge overflow bin to the last bin
        h1.SetBinContent(nbins, h1.GetBinContent(
            nbins)+h1.GetBinContent(nbins+1))
        # treat the uncertainties as uncorrelated
        h1.SetBinError(nbins, ROOT.TMath.Sqrt(
            (h1.GetBinError(nbins))**2 + (h1.GetBinError(nbins+1))**2))
        # clean the old bins
        h1.SetBinContent(nbins + 1, 0)
        h1.SetBinError(nbins + 1, 0)
    else:
        h1.SetBinContent(1, h1.GetBinContent(0)+h1.GetBinContent(1))
        # treat the uncertainties as uncorrelated
        h1.SetBinError(1, ROOT.TMath.Sqrt(
            (h1.GetBinError(1))**2 + (h1.GetBinError(0))**2))
        # clean the old bins
        h1.SetBinContent(0, 0)
        h1.SetBinError(0, 0)


def AddOverflows(hinput, dolastbin=True):
    '''
    move the over/under flow bin to the last/first bin
    '''
    if isinstance(hinput, ROOT.TH1):
        AddOverflowsTH1(hinput, dolastbin)

    if isinstance(hinput, ROOT.THStack):
        # do the AddOverflowsTH1 for all the histograms in THStack
        hlist = list(hinput.GetHists())
        list(map(AddOverflowsTH1, hlist, [dolastbin]*len(hlist)))

    else:
        print("input must be a ROOT.TH1 or ROOT.THStack for Over/Underflows")


def IncludeOverflow2D(h2, doUnderflow=False):
    """
    this might not work for one th2 with only one bin in one direction
    """
    nbinsX = h2.GetNbinsX()
    nbinsY = h2.GetNbinsY()

    for ix in range(1, nbinsX+1):
        h2 = CombineOneBin2D(h2, ix, nbinsY, ix, nbinsY+1)
    for iy in range(1, nbinsY+1):
        h2 = CombineOneBin2D(h2, nbinsX, iy, nbinsX+1, iy)
    # correct the corner
    h2 = CombineOneBin2D(h2, nbinsX, nbinsY, nbinsX+1, nbinsY+1)

    if doUnderflow:
        for ix in range(1, nbinsX+1):
            h2 = CombineOneBin2D(h2, ix, 1, ix, 0)
        for iy in range(1, nbinsY+1):
            h2 = CombineOneBin2D(h2, 1, iy, 0, iy)
        # correct the corner
        h2 = CombineOneBin2D(h2, 1, 1, 0, 0)
        h2 = CombineOneBin2D(h2, 1, nbinsY, 0, nbinsY+1)
        h2 = CombineOneBin2D(h2, nbinsX, 1, nbinsX+1, 0)
    return h2


def Ratio2Diff(hratio, inpercent=True):
    '''
    convert the ratio plot to diff, by substracting 1.0,
    and may show it in percentile
    '''
    for ibin in range(0, hratio.GetNbinsX()+1):
        hratio.SetBinContent(ibin, hratio.GetBinContent(ibin)-1.0)
    if inpercent:
        hratio.Scale(100.0)


def MultiplyH2(h1, h2):
    assert h1.GetNbinsX() == h2.GetNbinsX(
    ), "h1 and h2 have different number of bins in x"
    assert h1.GetNbinsY() == h2.GetNbinsY(
    ), "h1 and h2 have different number of bins in y"

    for ix in range(1, h1.GetNbinsX()+1):
        for iy in range(1, h1.GetNbinsY()+1):
            val1 = h1.GetBinContent(ix, iy)
            val2 = h2.GetBinContent(ix, iy)
            err1 = h1.GetBinError(ix, iy)
            err2 = h2.GetBinError(ix, iy)
            h1.SetBinContent(ix, iy, val1*val2)
            h1.SetBinError(ix, iy, ROOT.TMath.Sqrt(
                (val1*err2)**2 + (val2*err1)**2))
    return h1


def PositiveProtection(h):
    if isinstance(h, ROOT.TH2):
        print("input is a ROOT.TH2", h.GetName())
        PositiveProtection2D(h)
    elif isinstance(h, ROOT.TH1):
        PositiveProtection1D(h)
    else:
        print("input must be a ROOT.TH1 or ROOT.TH2 for PositiveProtection")
        sys.exit(1)


def PositiveProtection1D(h1):
    for ix in range(1, h1.GetNbinsX()+1):
        if h1.GetBinContent(ix) < 0:
            print("WARNING: negative bin content found, set to 0 for bin",
                  ix, h1.GetBinContent(ix), " in histogram ", h1.GetName())
            h1.SetBinContent(ix, 0)
            h1.SetBinError(ix, 0)


def PositiveProtection2D(h2):
    for ix in range(1, h2.GetNbinsX()+1):
        for iy in range(1, h2.GetNbinsY()+1):
            if h2.GetBinContent(ix, iy) < 0:
                print("WARNING: negative bin content found, set to 0 for bin", ix,
                      iy, h2.GetBinContent(ix, iy), " in histogram ", h2.GetName())
                h2.SetBinContent(ix, iy, 0)
                h2.SetBinError(ix, iy, 0)


def IntegralAndError2D(h2s):
    if isinstance(h2s, ROOT.TH2):
        h2s = [h2s]
    val = 0
    err2 = 0
    for h2 in h2s:
        # print("hname ", h2.GetName())
        for ibinx in range(1, h2.GetNbinsX()+1):
            for ibiny in range(1, h2.GetNbinsY()+1):
                val += h2.GetBinContent(ibinx, ibiny)
                err2 += h2.GetBinError(ibinx, ibiny)**2
    return val, ROOT.TMath.Sqrt(err2)


def CombineOneBin2D(h, ix, iy, jx, jy):
    """
    combine the j-th bin to i-th bin, and clean j-th bin
    """
    h.SetBinContent(ix, iy, h.GetBinContent(ix, iy) + h.GetBinContent(jx, jy))
    h.SetBinError(ix, iy, ROOT.TMath.Sqrt(
        h.GetBinError(ix, iy)**2 + h.GetBinError(jx, jy)**2))
    # clean the original bin so that there would not be double counting
    h.SetBinContent(jx, jy, 0)
    h.SetBinError(jx, jy, 0)
    return h


def GetRatioPanel(hs):
    """
    build the ratio error bars used for THStack
    Assuming different histograms are uncorrelated
    """
    hlist = hs.GetHists()
    hratio = None
    for h in hlist:
        if not hratio:
            hratio = h.Clone(h.GetName() + "_ratio")
        else:
            hratio.Add(h)
    hratio.Divide(hratio)
    hratio.SetFillColor(15)
    return hratio


def LHistos2Hist(hs, hname):
    """
    combine a list of hists to one hist
    """
    h_added = None
    for h in hs:
        if not h_added:
            h_added = h.Clone(hname)
        else:
            h_added.Add(h)
    return h_added


def SymmetrizeHisto(h, hSysUp, hname="hSysDown", useRatio=False):
    """
    given one histogram and its up variation,
    prepare the down variation.
    If useRatio is True, then the down variation is c*c/c_up
    else it is 2*c - c_up
    """
    hSysDown = h.Clone(hname)
    for ibin in range(1, h.GetNbinsX()+1):
        val_c = h.GetBinContent(ibin)
        val_up = hSysUp.GetBinContent(ibin)
        if useRatio:
            val_down = val_c*val_c / val_up
        else:
            val_down = 2*val_c - val_up
        hSysDown.SetBinContent(ibin, val_down)
    return hSysDown


def TH2ToTH1s(h2, projY=False, label="X"):
    """
    return a list of h1 from a 2D histogram,
    by default project on X axis
    """
    hs = []
    labels = []
    if not projY:
        for i in range(1, h2.GetNbinsY()+1):
            h1 = h2.ProjectionX(h2.GetName() + f"_projX_{i}", i, i)
            h1.SetLineColor(i)
            h1.SetMarkerColor(i)
            hs.append(h1)
            ymin, ymax = h2.GetYaxis().GetBinLowEdge(i), h2.GetYaxis().GetBinUpEdge(i)
            labels.append(f"{ymin: .2f}<{label}<{ymax: .2f}")
    else:
        for i in range(1, h2.GetNbinsX()+1):
            h1 = h2.ProjectionY(h2.GetName() + f"_projY_{i}", i, i)
            h1.SetLineColor(i)
            h1.SetMarkerColor(i)
            hs.append(h1)
            xmin, xmax = h2.GetXaxis().GetBinLowEdge(i), h2.GetXaxis().GetBinUpEdge(i)
            labels.append(f"{xmin: .2f}<{label}<{xmax: .2f}")
    return hs, labels


def DrawHistos(myhistos, mylabels, xmin, xmax, xlabel, ymin, ymax, ylabel, outputname, dology=True, showratio=False, dologx=False, lheader=None, donormalize=False, binomialratio=False, yrmax=2.0, yrmin=0.0, yrlabel=None, MCOnly=False, leftlegend=False, mycolors=None, legendPos=None, legendNCols=1, linestyles=None, markerstyles=None, showpull=False, doNewman=False, doPearson=False, ignoreHistError=False, ypullmin=-3.99, ypullmax=3.99, drawashist=False, padsize=(2, 0.9, 1.1), setGridx=False, setGridy=False, drawoptions=None, legendoptions=None, ratiooptions=None, dologz=False, doth2=False, ratiobase=0, redrawihist=-1, extraText=None, noCMS=False, noLumi=False, nMaxDigits=None, addOverflow=False, addUnderflow=False, plotdiff=False, hratiopanel=None, doratios=None, hpulls=None, W_ref=600, is5TeV=False, outdir="plots", savepdf=True, zmin=0, zmax=2, extralabels=None, extralheader=None, extraToDraw=None, exlegoffset=0.08):
    """
    draw histograms with the CMS tdr style
    """
    # python feature: for immutable objects, default values are evaluated only once
    # need to be cautious when using mutable objects as default values
    # https://docs.python.org/3/reference/compound_stmts.html#function-definitions
    if mycolors is None:
        mycolors = []
    if legendPos is None:
        legendPos = []
    if linestyles is None:
        linestyles = []
    if markerstyles is None:
        markerstyles = []
    if drawoptions is None:
        drawoptions = []
    if legendoptions is None:
        legendoptions = []
    if ratiooptions is None:
        ratiooptions = []

    # set the tdr style
    tdrstyle.setTDRStyle()
    # not sure why need this...
    ROOT.gStyle.SetErrorX(0.5)

    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetPaintTextFormat(".0f")

    # change the CMS_lumi variables (see CMS_lumi.py)
    # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
    CMS_lumi.lumi_sqrtS = ""
    CMS_lumi.relPosX = 0.25
    # CMS_lumi.extraText = "Internal"
    CMS_lumi.extraText = ""
    if MCOnly:
        CMS_lumi.extraText = "Simulation"
    if isinstance(extraText, str):
        CMS_lumi.extraText = extraText

    if nMaxDigits:
        # print(f"set the maximum number of digits {nMaxDigits}")
        ROOT.TGaxis.SetMaxDigits(nMaxDigits)
    else:
        # default val
        ROOT.TGaxis.SetMaxDigits(5)

    if ymax == None:
        ymax = max([h.GetMaximum() for h in myhistos]) * 1.25
    if ymin == None:
        ymin = min([h.GetMinimum() for h in myhistos]) * 0.75

    H_ref = 500
    W = W_ref
    H = H_ref

    iPeriod = 0

    npads = 1 + showratio + showpull

    if npads == 2:
        canvas = ROOT.TCanvas("c2"+outputname, "c2", 50, 50, W, 600)
        padsize1 = float(padsize[0])/(padsize[0]+padsize[1])
        padsize2 = float(padsize[1])/(padsize[0]+padsize[1])
        padsize3 = 0.
    elif npads == 1:
        canvas = ROOT.TCanvas("c2"+outputname, "c2", 50, 50, W, 600)
        canvas.SetGrid(setGridx, setGridy)
        canvas.SetTicks(1, 1)
        padsize1 = 1.0
        if doth2:
            if W_ref < 1000:
                canvas.SetRightMargin(0.18)
            else:
                canvas.SetRightMargin(0.10)
        padsize2 = 0.
        padsize3 = 0.
        canvas.cd()
    else:
        canvas = ROOT.TCanvas("c2"+outputname, "c2", 50, 50, W, 800)
        padsize1 = float(padsize[0])/(padsize[0]+padsize[1]+padsize[2])
        padsize2 = float(padsize[1])/(padsize[0]+padsize[1]+padsize[2])
        padsize3 = float(padsize[2])/(padsize[0]+padsize[1]+padsize[2])

    if dology:
        # print "set y axis to log scale"
        canvas.SetLogy()
    if dologx:
        # print "set x axis to log scale"
        canvas.SetLogx()
    if dologz:
        # print "set z axis to log scale"
        canvas.SetLogz()

    if npads == 1:
        canvas.SetLeftMargin(0.15)
        canvas.SetBottomMargin(0.13)
        canvas.SetTopMargin(0.06)

    if npads == 2:
        pad1 = ROOT.TPad("pad1" + outputname, "pad1", 0, padsize2, 1, 1)
        pad2 = ROOT.TPad("pad2" + outputname, "pad1", 0, 0, 1, padsize2)
        pad1.SetTopMargin(0.06/padsize1)
        pad1.SetBottomMargin(0.012/padsize1)
        pad1.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetTopMargin(0.010/padsize2)
        pad2.SetBottomMargin(0.13/padsize2)
        pad2.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetGridy(1)
        pad2.SetTicks(1, 1)
        pad1.Draw()
        pad2.Draw()

    if npads == 3:
        pad1 = ROOT.TPad("pad1" + outputname, "pad1", 0, 1-padsize1, 1, 1)
        pad2 = ROOT.TPad("pad2" + outputname, "pad2",
                         0, padsize3, 1, 1-padsize1)
        pad3 = ROOT.TPad("pad3" + outputname, "pad3", 0, 0.,   1, padsize3)
        pad1.SetTopMargin(0.06/(padsize1+padsize3))
        pad1.SetBottomMargin(0.012/padsize1)
        pad1.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetTopMargin(0.010/padsize2)
        pad2.SetBottomMargin(0.02/padsize2)
        pad2.SetLeftMargin(0.15 * (600.0)/W)
        pad2.SetGridy(1)
        pad2.SetTicks(1, 1)
        pad3.SetTopMargin(0.01/padsize3)
        pad3.SetBottomMargin(0.13/padsize3)
        pad3.SetLeftMargin(0.15 * (600.0)/W)
        pad3.SetGridy(1)
        pad3.SetTicks(1, 1)
        pad1.Draw()
        pad2.Draw()
        pad3.Draw()

    if npads > 1:
        pad1.cd()
        pad1.SetGrid(setGridx, setGridy)
        pad1.SetTicks(1, 1)
        if dology:
            pad1.SetLogy()
        if dologx:
            pad1.SetLogx()

    if not doth2:
        h1 = ROOT.TH1F("h1" + outputname, "h1", 80, xmin, xmax)
        h1.SetMinimum(ymin)
        h1.SetMaximum(ymax)
    else:
        h1 = ROOT.TH2F("h2" + outputname, "h2", 80, xmin, xmax, 80, ymin, ymax)
        h1.GetZaxis().SetLabelSize(0.03)
        if zmin != None and zmax != None:
            # print(f"configuring z range to {zmin}, {zmax}")
            h1.GetZaxis().SetRangeUser(zmin, zmax)

    # print "xmin : %f xmax : %f"%(xmin, xmax)

    h1.GetXaxis().SetNdivisions(6, 5, 0)
    h1.GetYaxis().SetNdivisions(6, 5, 0)
    h1.GetYaxis().SetTitle("%s" % ylabel)
    h1.GetYaxis().SetTitleSize(0.050/(padsize1+padsize3))
    h1.GetYaxis().SetLabelSize(0.045/(padsize1+padsize3))
    h1.GetXaxis().SetTitleSize(0.050/(padsize1+padsize3))
    h1.GetXaxis().SetLabelSize(0.045/(padsize1+padsize3))
    h1.GetYaxis().SetTitleOffset(1.35*(padsize1+padsize3)*(600.0/W))
    h1.GetXaxis().SetTitleOffset(1.1*(padsize1+padsize3))

    if showratio or showpull:
        h1.GetXaxis().SetLabelSize(0)
    else:
        h1.GetXaxis().SetTitle("%s" % xlabel)

    h1.Draw()

    x1_l = 0.92
    y1_l = 0.88

    # dx_l = 0.20
    dy_l = 0.22*0.66/padsize1
    dx_l = 0.35*0.66/padsize1
    x0_l = x1_l-dx_l
    y0_l = y1_l-dy_l

    if len(legendPos) == 4:
        x0_l = legendPos[0]
        y0_l = legendPos[1]
        x1_l = legendPos[2]
        y1_l = legendPos[3]

    if npads > 2:
        y1_l -= 0.03
        y0_l -= 0.03

    if leftlegend:
        x0_l -= 0.40
        x1_l -= 0.40

    legend = ROOT.TLegend(x0_l, y0_l, x1_l, y1_l)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.SetFillColor(0)
    legend.SetNColumns(legendNCols)
    if lheader and lheader != "":
        # legend.SetHeader(lheader, "C")
        legend.SetHeader(lheader)

    if extralabels:
        exlegend = ROOT.TLegend(x0_l-exlegoffset, y0_l, x0_l, y1_l)
        exlegend.SetBorderSize(0)
        exlegend.SetTextSize(0.04)
        exlegend.SetTextFont(42)
        exlegend.SetFillColor(0)
        exlegend.SetMargin(0.99)
        if extralheader and extralheader != "":
            exlegend.SetHeader(extralheader)

    myhistos_clone = []
    for ihisto in myhistos:
        # clone the histogram for plotting,
        # so the original histogram is unchanged
        ihcolone = ihisto.Clone("%s_Clone" % ihisto.GetName())
        # ihcolone.SetDirectory(0)
        myhistos_clone.append(ihcolone)

    if drawashist:
        drawoptions = ["HIST" for i in range(len(myhistos_clone))]
        legendoptions = ['L' for i in range(len(myhistos_clone))]

    if drawoptions and not isinstance(drawoptions, list):
        # copy the option n times
        tmp = [drawoptions for i in range(len(myhistos_clone))]
        drawoptions = tmp
    if not isinstance(legendoptions, list):
        tmp = [legendoptions for i in range(len(myhistos_clone))]
        legendoptions = tmp

    # print(("legend options ", legendoptions))

    ileg = 0
    for idx in range(0, len(myhistos_clone)):
        if addOverflow:
            # include overflow bin
            AddOverflows(myhistos_clone[idx])
        if addUnderflow:
            # include underflow bin
            AddOverflows(myhistos_clone[idx], dolastbin=False)
        if idx < len(mycolors):
            myhistos_clone[idx].SetLineColor(mycolors[idx])
            myhistos_clone[idx].SetMarkerColor(mycolors[idx])
        if idx < len(linestyles):
            myhistos_clone[idx].SetLineStyle(linestyles[idx])
        if idx < len(markerstyles):
            myhistos_clone[idx].SetMarkerStyle(markerstyles[idx])
        if isinstance(myhistos_clone[idx], ROOT.TH1):
            myhistos_clone[idx].SetLineWidth(3)
        # myhistos_clone[idx].GetXaxis().SetRangeUser(xmin, xmax)
        if donormalize:
            Normalize(myhistos_clone[idx])
        if idx >= len(drawoptions):
            drawoptions.append("EL")
        if idx >= len(legendoptions):
            if isinstance(myhistos_clone[idx], ROOT.TH1) and myhistos_clone[idx].GetMarkerStyle() != 1:
                legendoptions.append("LP")
            else:
                legendoptions.append("LE")
        if isinstance(myhistos_clone[idx], ROOT.THStack):
            drawoptions[idx] = 'HIST'
            legendoptions[idx] = "F"
        myhistos_clone[idx].Draw(
            " ".join(filter(None, [drawoptions[idx], "same"])))

        exlen = 0
        if extralabels:
            exlen = len(extralabels)

        if ileg < exlen:
            exlegend.AddEntry(myhistos_clone[idx], " ", legendoptions[idx])
            ileg += 1

        elif ileg < len(mylabels) + exlen:
            # number of labels might be different from number of histograms.
            # do the match logically
            if isinstance(myhistos_clone[idx], ROOT.THStack):
                hlist = list(myhistos_clone[idx].GetHists())
                for hist in reversed(hlist):
                    legend.AddEntry(
                        hist, str(mylabels[ileg-exlen]), legendoptions[idx])
                    ileg += 1
            else:
                legend.AddEntry(myhistos_clone[idx], str(
                    mylabels[ileg-exlen]), legendoptions[idx])
                ileg += 1

    # print("draw options ", drawoptions)
    # print("legend options ", legendoptions)

    if showratio and hratiopanel:
        legend.AddEntry(hratiopanel, "Uncertainty", "F")

    if redrawihist >= 0:
        myhistos_clone[redrawihist].Draw(
            " ".join(filter(None, [drawoptions[redrawihist], "same"])))

    iPosX = 0
    plotCMS = True
    if noCMS:
        # do not draw CMS, only lumi info and extraText
        plotCMS = False
        CMS_lumi.cmsText = ""
        # CMS_lumi.extraText = ""
        # iPosX = 11
    if not is5TeV:
        iPeriod = 4
    else:
        iPeriod = 5
    if MCOnly:
        iPeriod = 0
    if showratio or showpull:
        if not noLumi:
            CMS_lumi.CMS_lumi(pad1, iPeriod, iPosX, plotCMS=plotCMS)
        pad1.Update()
        pad1.RedrawAxis()
        # frame = pad1.GetFrame()
        # frame.Draw()
        if len(mylabels):
            legend.Draw()
        if extralabels:
            exlegend.Draw()
        if extraToDraw:
            if not isinstance(extraToDraw, list):
                extraToDraw = [extraToDraw]
            for todraw in extraToDraw:
                todraw.Draw()
        pad1.Update()

    else:
        if not noLumi:
            CMS_lumi.CMS_lumi(canvas, iPeriod, iPosX, plotCMS=plotCMS)
        canvas.cd()
        canvas.Update()
        canvas.RedrawAxis()
        canvas.RedrawAxis("G")
        # frame = canvas.GetFrame()
        # frame.Draw()
        if len(mylabels):
            # print("draw legend", len(mylabels))
            legend.Draw()
        if extralabels:
            exlegend.Draw()
        if extraToDraw:
            if not isinstance(extraToDraw, list):
                extraToDraw = [extraToDraw]
            for todraw in extraToDraw:
                todraw.Draw()
        canvas.Update()

    ##
    # ration plot
    ##
    if showratio or showpull:
        if isinstance(myhistos_clone[ratiobase], ROOT.THStack):
            hden = THStack2TH1(myhistos_clone[ratiobase])
        else:
            hden = myhistos_clone[ratiobase].Clone("hden_%s" % idx)

    if showratio:
        hratios = {}
        for idx in range(0, len(myhistos_clone)):
            if doratios != None and doratios[idx] == False:
                continue
            if idx == ratiobase:
                continue
            if isinstance(myhistos_clone[idx], ROOT.THStack):
                hratios[idx] = THStack2TH1(myhistos_clone[idx])
            else:
                hratios[idx] = myhistos_clone[idx].Clone("hratios_%s" % idx)
            if binomialratio:
                hratios[idx].Divide(hratios[idx], hden, 1.0, 1.0, "b")
            else:
                hratios[idx].Divide(hden)
            if plotdiff:
                Ratio2Diff(hratios[idx])

    if showpull:
        if hpulls == None:
            hpulls = []
            for idx in range(0, len(myhistos_clone)):
                if idx == ratiobase:
                    continue
                hpull = myhistos_clone[idx].Clone("hpullsnew_%s" % idx)
                hpulls.append(hpull)
                CalculateChi(hpull, hden, doNewman=doNewman,
                             doPearson=doPearson, ignoreHistError=ignoreHistError)

    if npads > 1:
        pad2.cd()
        if dologx:
            pad2.SetLogx()

        h2 = ROOT.TH1F("h2", "h2", 80, xmin, xmax)
        h2.GetXaxis().SetTitleSize(0.050/(padsize2+0.3*padsize3))
        h2.GetXaxis().SetLabelSize(0.045/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetTitleSize(0.050/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetLabelSize(0.045/(padsize2+0.3*padsize3))
        h2.GetYaxis().SetTitleOffset(1.35*(padsize2+0.35*padsize3)*(600.0/W))

        h2.GetYaxis().SetNdivisions(8)
        h2.GetYaxis().CenterTitle()
        if yrlabel:
            ytitle = yrlabel
        elif showratio:
            ytitle = "Ratio"
        else:
            ytitle = "Pull"
        h2.GetYaxis().SetTitle(ytitle)
        h2.SetMaximum(yrmax if showratio else ypullmax)
        h2.SetMinimum(yrmin if showratio else ypullmin)

        if npads == 2:
            h2.GetXaxis().SetTitle("%s" % xlabel)
            h2.GetXaxis().SetTitleOffset(1.1)
        else:
            h2.GetXaxis().SetTitleSize(0.0)
            h2.GetXaxis().SetLabelSize(0)

        h2.Draw()

        if showratio:
            if hratiopanel:
                # print("draw hratiopanel")
                # hratiopanel.SetFillColor(15)
                hratiopanel.SetFillColorAlpha(15, 0.5)
                hratiopanel.SetLineColor(0)
                # hratiopanel.SetLineColor(2)
                hratiopanel.SetMarkerSize(0)
                hratiopanel.Draw("E2 same")
                # hratiopanel.Draw("HIST same")
                # for ibin in xrange(1, hratiopanel.GetNbinsX()+1):
                #    print("bin content ", hratiopanel.GetBinContent(ibin), " error ", hratiopanel.GetBinError(ibin))
                h2.Draw("same")
            idx = 0
            for hratio in list(hratios.values()):
                if idx >= len(ratiooptions):
                    if drawashist:
                        ratiooptions.append("HIST")
                    else:
                        ratiooptions.append("")
                hratio.SetFillColor(0)
                hratio.Draw(ratiooptions[idx] + " same")
                idx += 1
        elif showpull:
            for hpull in hpulls:
                hpull.Draw("HIST same")

        pad2.RedrawAxis("G")
        pad2.Update()

    if npads > 2:
        pad3.cd()
        if dologx:
            pad3.SetLogx()

        h3 = ROOT.TH1F("h3", "h3", 80, xmin, xmax)
        h3.GetXaxis().SetTitle("%s" % xlabel)
        h3.GetXaxis().SetTitleOffset(1.1)
        h3.GetXaxis().SetTitleSize(0.050/(padsize3+0.3*padsize2))
        h3.GetXaxis().SetLabelSize(0.045/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetTitleSize(0.050/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetLabelSize(0.045/(padsize3+0.3*padsize2))
        h3.GetYaxis().SetTitleOffset(1.35*(padsize3+0.35*padsize2)*(600.0/W))
        h3.GetYaxis().SetNdivisions(5)
        h3.GetYaxis().CenterTitle()
        # h3.GetYaxis().SetTitle("#frac{Obs. - Exp.}{Err}")
        h3.GetYaxis().SetTitle("Pull")

        h3.SetMaximum(ypullmax)
        h3.SetMinimum(ypullmin)
        h3.Draw()

        for hpull in hpulls:
            hpull.SetFillColorAlpha(2, 0.6)
            hpull.SetLineWidth(1)
            hpull.SetLineColor(1)
            hpull.Draw("HIST same")

        pad3.Update()

    canvas.Update()

    if "/" not in outputname:
        # path not included; by default put to plots/outputname
        outputname = outdir + "/" + outputname

    dirpath = outputname.rpartition('/')[0]
    if not os.path.exists(dirpath):
        print(f"Make the directory {dirpath}")
        os.makedirs(dirpath)

    if savepdf:
        # print("save plot to %s.pdf" % outputname)
        # canvas.Print("%s.C"%outputname)
        canvas.Print("%s.pdf" % outputname)
        # canvas.Print("%s.png" % outputname)
        # canvas.Print("%s.root" % outputname)
    return hratios if showratio else None
