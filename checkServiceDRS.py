import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import getUpstreamVetoChannel, getPreShowerChannel, getDownStreamMuonChannel
from utils.html_generator import generate_html
from runconfig import runNumber, firstEvent, lastEvent
from utils.utils import loadRDF, loadRDF, preProcessDRSBoards
from configs.plotranges import getServiceDRSProcessedInfoRanges
from selections.selections import checkUpstreamVeto

ROOT.gROOT.SetBatch(True)  # Run in batch mode
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")


def analyzePulse(channels, names):
    rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
    rdf = preProcessDRSBoards(rdf)
    hists = {}

    channel_preshower = getPreShowerChannel(runNumber)
    channel_muon = getDownStreamMuonChannel(runNumber)

    rdf = rdf.Define(f"{channel_preshower}_peak_position",
                     f"ArgMinRange({channel_preshower}_subtractMedian, 100, 400)")
    rdf = rdf.Define(f"{channel_preshower}_peak_value",
                     f"MinRange({channel_preshower}_subtractMedian, 100, 400)")
    rdf = rdf.Define(f"{channel_preshower}_sum",
                     f"SumRange({channel_preshower}_subtractMedian, 100, 400)")

    rdf = rdf.Define(f"{channel_muon}_peak_position",
                     f"ROOT::VecOps::ArgMin({channel_muon}_subtractMedian)")
    rdf = rdf.Define(f"{channel_muon}_peak_value",
                     f"ROOT::VecOps::Min({channel_muon}_subtractMedian)")
    rdf = rdf.Define(f"{channel_muon}_sum",
                     f"ROOT::VecOps::Sum({channel_muon}_subtractMedian)")

    for name, channel in zip(names, channels):
        hists[channel] = {}
        hists[channel]["peak_position"] = rdf.Histo1D(
            (f"{name}_peak_position",
             f"Peak Position {channel};Time Slice;Counts", 128, 0, 1024),
            f"{channel}_peak_position"
        )
        xmin, xmax = getServiceDRSProcessedInfoRanges(name, "peak_value")
        hists[channel]["peak_value"] = rdf.Histo1D(
            (f"{name}_peak_value",
             f"Peak Value {channel};ADC Counts;Counts", 50, xmin, xmax),
            f"{channel}_peak_value"
        )
        xmin, xmax = getServiceDRSProcessedInfoRanges(name, "sum")
        hists[channel]["sum"] = rdf.Histo1D(
            (f"{name}_sum",
             f"Sum {channel};ADC Counts;Counts", 100, xmin, xmax),
            f"{channel}_sum"
        )

    print("Writing histograms to output file...")
    ofile = ROOT.TFile(
        f"results/root/Run{runNumber}/drs_service.root", "RECREATE")
    for name, channel in zip(names, channels):
        for _, hist in hists[channel].items():
            hist.SetDirectory(ofile)
            hist.Write()
    ofile.Close()


def plotPulse(channels, names):
    infile = ROOT.TFile(
        f"results/root/Run{runNumber}/drs_service.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    plots = []
    outdir = f"results/plots/Run{runNumber}/drs_service/"

    for name, channel in zip(names, channels):
        hist = infile.Get(f"{name}_peak_position")
        if not hist:
            print(
                f"Histogram {name}_peak_position not found in {infile.GetName()}")
            return
        DrawHistos([hist], [name], 0, 1024, "Peak Position", 0, None, "Counts",
                   outputname=f"{name}_peak_position", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, runNumber=runNumber,
                   addOverflow=True, addUnderflow=True)
        plots.append(f"{name}_peak_position.png")

        hist = infile.Get(f"{name}_peak_value")
        if not hist:
            print(
                f"Histogram {name}_peak_value not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges(name, "peak_value")
        DrawHistos([hist], [name], xmin, xmax, "Peak Value", 0, None, "Counts",
                   outputname=f"{name}_peak_value", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, runNumber=runNumber,
                   addOverflow=True, addUnderflow=True)
        plots.append(f"{name}_peak_value.png")

        hist = infile.Get(f"{name}_sum")
        if not hist:
            print(f"Histogram {name}_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges(name, "sum")
        DrawHistos([hist], [name], xmin, xmax, "Sum", 0, None, "Counts",
                   outputname=f"{name}_sum", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, runNumber=runNumber,
                   addOverflow=True, addUnderflow=True)
        plots.append(f"{name}_sum.png")

    generate_html(plots, f"results/plots/Run{runNumber}/drs_service/", plots_per_row=3,
                  output_html=f"results/html/Run{runNumber}/drs_service/viewer.html")


if __name__ == "__main__":
    chan_preshower = getPreShowerChannel(runNumber)
    chan_muon = getDownStreamMuonChannel()
    channels = [chan_preshower, chan_muon]
    names = ["preshower", "muon"]

    analyzePulse(channels, names)
    plotPulse(channels, names)
