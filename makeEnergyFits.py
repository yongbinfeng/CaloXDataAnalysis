import os
import sys
import ROOT
import array
from utils.fitter import eventFit
from utils.html_generator import generate_html
from configs.plotranges import getRangesForFERSEnergySums, getEventEnergyFitParameters
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos

ROOT.gROOT.SetBatch(True)


def runFit(runNumber, subtractPedestal=False, calibrate=False, clip=False, outdir="results/plots/summary"):
    suffix, _, _, _, _, _, _, xtitle = getRangesForFERSEnergySums(
        subtractPedestal=subtractPedestal, calibrate=calibrate, clip=clip)
    rootdir = f"results/root/Run{runNumber}/"
    filename = f"{rootdir}/fers_energy_sum{suffix}.root"
    if not os.path.exists(filename):
        print(
            f"File {filename} does not exist. Please run makeFERSEnergyPlots.py first.")
        exit(1)

    ifile = ROOT.TFile(filename, "READ")
    hCer = ifile.Get(f"hist_FERS_CerEnergyHG{suffix}")
    hSci = ifile.Get(f"hist_FERS_SciEnergyHG{suffix}")

    args_cer = getEventEnergyFitParameters(
        runNumber, isCer=True, clip=clip)
    args_sci = getEventEnergyFitParameters(
        runNumber, isCer=False, clip=clip)

    plots = []
    output_name, results_cer = eventFit(hCer, f"Run{runNumber}_CerHG{suffix}",
                                        outdir=outdir, xlabel="Cer " + xtitle,
                                        **args_cer)
    plots.append(output_name)
    output_name, results_sci = eventFit(hSci, f"Run{runNumber}_SciHG{suffix}",
                                        outdir=outdir, xlabel="Sci " + xtitle,
                                        **args_sci)
    plots.append(output_name)

    return plots, results_cer, results_sci


def makeOneGraph(np_xs, np_ys, np_xerrs, np_yerrs, np_runs=None, color=ROOT.kRed):
    graph = ROOT.TGraphErrors(
        len(np_xs), np_xs, np_ys, np_xerrs, np_yerrs)
    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)

    extraToDraws = []
    for i in range(len(np_xs)):
        if np_runs is not None and i < len(np_runs):
            x = np_xs[i]
            y = np_ys[i]
            l = ROOT.TLatex(x+5, y, f"{np_runs[i]}")
            l.SetTextSize(0.02)
            l.SetTextColor(color)
            l.SetTextAlign(12)  # Centered text
            extraToDraws.append(l)

    return graph, extraToDraws


def processRuns():
    outdir = "results/plots/summary"
    htmldir = "results/html/summary"
    plots = []

    runinfo = {
        1190: 80.0,
        1221: 20.0,
        1256: 80.0,
        1260: 60.0,
        1261: 30.0,
    }

    gMeans_cer = []
    extras_mean_cer = []
    gMeans_sci = []
    extras_mean_sci = []
    gSigmas_cer = []
    extras_sigma_cer = []
    gSigmas_sci = []
    extras_sigma_sci = []

    for run, energy in runinfo.items():
        plots_run, results_cer, results_sci = runFit(run, subtractPedestal=True,
                                                     calibrate=True, clip=False, outdir=outdir)
        plots += plots_run

        np_energy = array.array("f", [energy])
        np_energy_err = array.array(
            "f", [energy * 0.03])  # assumed 3% uncertainty
        np_cer_mean = array.array("f", [results_cer["mean"][0] / energy])
        np_cer_mean_err = array.array("f", [results_cer["mean"][1] / energy])
        np_sci_mean = array.array("f", [results_sci["mean"][0] / energy])
        np_sci_mean_err = array.array("f", [results_sci["mean"][1] / energy])

        np_cer_sigma = array.array(
            "f", [results_cer["sigma"][0] / results_cer["mean"][0]])
        np_cer_sigma_err = array.array(
            "f", [results_cer["sigma"][1] / results_cer["mean"][0]])
        np_sci_sigma = array.array(
            "f", [results_sci["sigma"][0] / results_sci["mean"][0]])
        np_sci_sigma_err = array.array(
            "f", [results_sci["sigma"][1] / results_sci["mean"][0]])

        gMean_cer, extra_cer = makeOneGraph(
            np_energy, np_cer_mean, np_energy_err, np_cer_mean_err, np_runs=[run], color=ROOT.kRed)
        gMeans_cer.append(gMean_cer)
        extras_mean_cer += extra_cer
        gMean_sci, extra_sci = makeOneGraph(
            np_energy, np_sci_mean, np_energy_err, np_sci_mean_err, np_runs=[run], color=ROOT.kBlue)
        gMeans_sci.append(gMean_sci)
        extras_mean_sci += extra_sci

        gSigma_cer, extra_sigma_cer = makeOneGraph(
            np_energy, np_cer_sigma, np_energy_err, np_cer_sigma_err, np_runs=[run], color=ROOT.kRed)
        gSigmas_cer.append(gSigma_cer)
        extras_sigma_cer += extra_sigma_cer
        gSigma_sci, extra_sigma_sci = makeOneGraph(
            np_energy, np_sci_sigma, np_energy_err, np_sci_sigma_err, np_runs=[run], color=ROOT.kBlue)
        gSigmas_sci.append(gSigma_sci)
        extras_sigma_sci += extra_sigma_sci

    output_name = "energy_cer"
    DrawHistos(gMeans_cer, "", 0, 100, "Energy [GeV]", 0, 35, "# p.e. / GeV", output_name, dology=False, drawoptions=[
               "P"]*len(gMeans_cer), outdir=outdir, runNumber="e^{+}", extraToDraw=extras_mean_cer)
    plots.insert(0, output_name + ".png")
    output_name = "energy_sci"
    DrawHistos(gMeans_sci, "", 0, 100, "Energy [GeV]", 0, 220, "# p.e. / GeV", output_name, dology=False, drawoptions=[
               "P"]*len(gMeans_sci), outdir=outdir, runNumber="e^{+}", extraToDraw=extras_mean_sci)
    plots.insert(1, output_name + ".png")

    output_name = "resol_cer"
    DrawHistos(gSigmas_cer, "", 0, 100, "Energy [GeV]", 0, 0.4, "<#sigma>/<#mu>", output_name, dology=False, drawoptions=[
               "P"]*len(gSigmas_cer), outdir=outdir, runNumber="e^{+}", extraToDraw=extras_sigma_cer)
    plots.insert(2, output_name + ".png")
    output_name = "resol_sci"
    DrawHistos(gSigmas_sci, "", 0, 100, "Energy [GeV]", 0, 0.4, "<#sigma>/<#mu>", output_name, dology=False, drawoptions=[
               "P"]*len(gSigmas_sci), outdir=outdir, runNumber="e^{+}", extraToDraw=extras_sigma_sci)
    plots.insert(3, output_name + ".png")

    output_html = f"{htmldir}/summary/energyfits.html"
    generate_html(plots, outdir, plots_per_row=2,
                  output_html=output_html)
    print(f"Generated HTML file: {output_html}")
    return output_html


if __name__ == "__main__":
    processRuns()
