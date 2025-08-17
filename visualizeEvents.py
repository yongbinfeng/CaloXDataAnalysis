import random
import os
import sys
import ROOT
from utils.channel_map import buildFERSBoards
from utils.utils import loadRDF, calculateEnergySumFERS, vectorizeFERS, calibrateFERSChannels, preProcessDRSBoards
from utils.html_generator import generate_html
from utils.colors import colors
from utils.visualization import visualizeFERSBoards, makeEventDisplay
from selections.selections import vetoMuonCounter, applyUpstreamVeto
from runconfig import runNumber, firstEvent, lastEvent
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos


print(f"Start running {os.path.basename(__file__)}")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

# load the gains and pedestals from SiPM fits
file_gains = f"results/root/Run{runNumber}/valuemaps_gain.json"
file_pedestals = f"results/root/Run{runNumber}/valuemaps_pedestal.json"

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)
rdf, rdf_filterveto = applyUpstreamVeto(rdf, runNumber)

FERSBoards = buildFERSBoards(run=runNumber)

rdf = vectorizeFERS(rdf, FERSBoards)
# define energy sums with different configurations
rdf = calibrateFERSChannels(
    rdf, FERSBoards, file_gains=file_gains, file_pedestals=file_pedestals)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=True, calibrate=True, clip=False)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"


# the values right now only works for 80GeV positrons
def filterBandEvents(rdf):
    suffix = "_subtracted_calibrated"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 2000 && FERS_SciEnergyHG{suffix} < 4000")


def filterDarkCountEvents(rdf):
    suffix = "_subtracted_calibrated"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} < 200")


def filterElePeakEvents(rdf):
    suffix = "_subtracted_calibrated"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 5000 && FERS_SciEnergyHG{suffix} < 7000")


def filterHEEvents(rdf):
    suffix = "_subtracted_calibrated"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 12000")


def fillInValueMap(rdf, event_numbers):
    """
    Fill in the value map for the given event number.
    """
    # this should be lazy action
    # do everything in a loop
    values_FERS_events_lazy = {}
    rdfs_temp = []
    for event_number in event_numbers:
        # filter everything to one event, then use the mean value for readout
        # this achieves "lazy" evaluation
        rdf_temp = rdf.Filter(f"event_n == {event_number}")
        values_FERS = {}
        for _, FERSBoard in FERSBoards.items():
            for chan in FERSBoard:
                channelName_HG = chan.GetHGChannelName()
                variableName = channelName_HG + "_subtracted_calibrated"
                value = rdf_temp.Mean(variableName)

                values_FERS[channelName_HG] = value
        # total energy for the event
        for var in ["Cer", "Sci"]:
            variableName = f"FERS_{var}EnergyHG_subtracted_calibrated"
            value = rdf_temp.Mean(variableName)
            values_FERS[f"FERS_{var}EnergyHG"] = value
        values_FERS_events_lazy[event_number] = values_FERS
        rdfs_temp.append(rdf_temp)

    # read all values
    values_FERS_events = {}
    for event_number, values_FERS in values_FERS_events_lazy.items():
        values_FERS_events[event_number] = {}
        for channelName_HG, value in values_FERS.items():
            value = value.GetValue()
            values_FERS_events[event_number][channelName_HG] = value

    return values_FERS_events


def DrawEventDisplay(values_FERS_events, suffix="MIP", applyScaling=False, zmax=100, zmin=2):
    outdir_plots = f"{plotdir}/EventDisplay_{suffix}/"
    plots = []
    for event_number, values_FERS in values_FERS_events.items():
        suffix_plot = f"_Run{runNumber}_Event{event_number}_Display"
        [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
            FERSBoards, values_FERS, suffix=suffix_plot, useHG=True)

        if applyScaling:
            # scale 6mm channel to 3mm channel size
            h2_Cer_HG_mean.Scale(1/4.0)
            h2_Sci_HG_mean.Scale(1/4.0)

        output_name = f"EventDisplay_Run{runNumber}_Event{event_number}"
        extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.40, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(
            f"E_{{Cer}} = {values_FERS['FERS_CerEnergyHG']:.0f} p.e.")
        makeEventDisplay(h2_Cer_HG_mean, h2_Cer_3mm_HG_mean, output_name +
                         "_Cer", outdir_plots, runNumber, zmin, zmax, isCer=True, extraToDraw=extraToDraw)
        plots.append(output_name + "_Cer.png")

        extraToDraw.Clear()
        extraToDraw.AddText(
            f"E_{{Sci}} = {values_FERS['FERS_SciEnergyHG']:.0f} p.e.")
        # make the event display for Sci
        makeEventDisplay(h2_Sci_HG_mean, h2_Sci_3mm_HG_mean, output_name +
                         "_Sci", outdir_plots, runNumber, zmin, zmax, isCer=False, extraToDraw=extraToDraw)
        plots.append(output_name + "_Sci.png")

    output_html = f"{htmldir}/EventDisplay_{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


rdfs_filtered = {}
rdfs_filtered["ElePeak"] = filterElePeakEvents(rdf)
rdfs_filtered["HE"] = filterHEEvents(rdf)
rdfs_filtered["Band"] = filterBandEvents(rdf)
rdfs_filtered["DarkCount"] = filterDarkCountEvents(rdf)

zranges = {
    "ElePeak": (2, 100),
    "HE": (2, 100),
    "Band": (1, 100),
    "DarkCount": (0, 100)
}

event_numbers = {}
for cat in rdfs_filtered:
    # get the event numbers for each category
    event_numbers[cat] = list(
        rdfs_filtered[cat].Take["unsigned int"]("event_n"))
    event_numbers[cat].sort()  # sort for reproducibility

output_htmls = {}

# randomly select 10 events from MIP and Dark Count events
random.seed(42)  # for reproducibility
for cat in event_numbers:
    selected_events = random.sample(
        event_numbers[cat], min(20, len(event_numbers[cat])))
    print(f"selected {cat} events: {selected_events}")

    values_FERS_events = fillInValueMap(rdfs_filtered[cat], selected_events)
    output_html = DrawEventDisplay(
        values_FERS_events, suffix=cat, applyScaling=True, zmin=zranges[cat][0], zmax=zranges[cat][1])
    output_htmls[cat] = output_html


print("Event display HTML files generated:")
for cat, html_file in output_htmls.items():
    print(f"{cat}: {html_file}")
print("All done!")
