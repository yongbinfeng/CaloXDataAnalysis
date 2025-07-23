import random
import os
import sys
import ROOT
from utils.channel_map import buildFERSBoards
from utils.utils import filterPrefireEvents, loadRDF, calculateEnergySumFERS, vectorizeFERS, calibrateFERSChannels
from utils.html_generator import generate_html
from utils.colors import colors
from utils.visualization import visualizeFERSBoards, makeEventDisplay
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
rdf, rdf_prefilter = filterPrefireEvents(rdf, runNumber)

FERSBoards = buildFERSBoards(run=runNumber)

rdf = vectorizeFERS(rdf, FERSBoards)
# define energy sums with different configurations
rdf = calibrateFERSChannels(
    rdf, FERSBoards, file_gains=file_gains, file_pedestals=file_pedestals)
rdf = calculateEnergySumFERS(
    rdf, FERSBoards, subtractPedestal=True, calibrate=True, clip=True)

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"


# select MIP and Dark Count events
def filterMIPEvents(rdf):
    suffix = "_subtracted_calibrated_clipped"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 200 && FERS_CerEnergyHG{suffix} > 140 && FERS_SciEnergyHG{suffix} < 400 && FERS_CerEnergyHG{suffix} < 300")


def filterDarkCountEvents(rdf):
    suffix = "_subtracted_calibrated_clipped"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} < 200 && FERS_CerEnergyHG{suffix} < 190")


def filterHEEvents(rdf):
    suffix = "_subtracted_calibrated_clipped"
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 150 && FERS_CerEnergyHG{suffix} > 600")


def fillInValueMap(rdf, event_numbers):
    """
    Fill in the value map for the given event number.
    """
    # this should be lazy action
    # do everything in a loop
    values_FERS_events_lazy = {}
    rdfs_temp = []
    for event_number in event_numbers:
        rdf_temp = rdf.Filter(f"event_n == {event_number}")
        values_FERS = {}
        for _, FERSBoard in FERSBoards.items():
            for chan in FERSBoard:
                channelName_HG = chan.GetHGChannelName()
                variableName = channelName_HG + "_subtracted_calibrated_clipped"
                value = rdf_temp.Mean(variableName)

                values_FERS[channelName_HG] = value
        values_FERS_events_lazy[event_number] = values_FERS
        rdfs_temp.append(rdf_temp)

    values_FERS_events = {}
    for event_number, values_FERS in values_FERS_events_lazy.items():
        values_FERS_events[event_number] = {}
        for channelName_HG, value in values_FERS.items():
            value = value.GetValue()
            values_FERS_events[event_number][channelName_HG] = value

    return values_FERS_events


def DrawEventDisplay(values_FERS_events, suffix="MIP"):
    outdir_plots = f"{plotdir}/EventDisplay_{suffix}/"
    plots = []
    for event_number, values_FERS in values_FERS_events.items():
        suffix_plot = f"_Run{runNumber}_Event{event_number}_Display"
        [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
            FERSBoards, values_FERS, suffix=suffix_plot, useHG=True)

        output_name = f"EventDisplay_Run{runNumber}_Event{event_number}"
        makeEventDisplay(h2_Cer_HG_mean, h2_Cer_3mm_HG_mean, output_name +
                         "_Cer", outdir_plots, runNumber, 2, 10, isCer=True)
        plots.append(output_name + "_Cer.png")

        makeEventDisplay(h2_Sci_HG_mean, h2_Sci_3mm_HG_mean, output_name +
                         "_Sci", outdir_plots, runNumber, 2, 10, isCer=False)
        plots.append(output_name + "_Sci.png")

    output_html = f"{htmldir}/EventDisplay_{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


rdfs_filtered = {}
rdfs_filtered["HE"] = filterHEEvents(rdf)
rdfs_filtered["MIP"] = filterMIPEvents(rdf)
rdfs_filtered["DarkCount"] = filterDarkCountEvents(rdf)

event_numbers = {}
event_numbers["HE"] = list(rdfs_filtered["HE"].Take["unsigned int"]("event_n"))
event_numbers["MIP"] = list(
    rdfs_filtered["MIP"].Take["unsigned int"]("event_n"))
event_numbers["DarkCount"] = list(
    rdfs_filtered["DarkCount"].Take["unsigned int"]("event_n"))

output_htmls = {}

# randomly select 10 events from MIP and Dark Count events
random.seed(42)  # for reproducibility
for cat in event_numbers:
    selected_events = random.sample(event_numbers[cat], 20)
    print(f"selected {cat} events: {selected_events}")

    values_FERS_events = fillInValueMap(rdfs_filtered[cat], selected_events)
    output_html = DrawEventDisplay(values_FERS_events, suffix=cat)
    output_htmls[cat] = output_html


print("Event display HTML files generated:")
for cat, html_file in output_htmls.items():
    print(f"{cat}: {html_file}")
print("All done!")
