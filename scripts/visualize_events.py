import random
import os
import sys
from collections import OrderedDict
import ROOT
from channels.channel_map import build_fers_boards, build_drs_boards
from utils.utils import number_to_string
from utils.data_loader import loadRDF, getRunInfo
from variables.fers import vectorizeFERS, calibrateFERSChannels, get_fers_energy_sum
from variables.drs import preProcessDRSBoards
from utils.html_generator import generate_html
from utils.colors import colors
from utils.visualization import visualizeFERSBoards, makeEventDisplay
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection
from configs.plot_ranges import get_drs_plot_ranges
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos


print(f"Start running {os.path.basename(__file__)}")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

args = get_args()
run_number = args.run
firstEvent = args.first_event
lastEvent = args.last_event
jsonFile = args.json_file
btype, benergy = getRunInfo(run_number)

rdf, rdf_org = loadRDF(run_number, firstEvent, lastEvent, jsonFile)
rdf = preProcessDRSBoards(rdf, run_number=run_number)
rdf, rdf_filterveto = applyUpstreamVeto(rdf, run_number)

FERSBoards = build_fers_boards(run=run_number)
DRSBoards = build_drs_boards(run=run_number)

rdf = vectorizeFERS(rdf, FERSBoards)
rdf = get_fers_energy_sum(rdf, FERSBoards, gain="HG")
# define energy sums with different configurations
# rdf = calibrateFERSChannels(
#    rdf, FERSBoards, file_gains=file_gains, file_pedestals=file_pedestals)

rootdir = f"results/root/Run{run_number}/"
plotdir = f"results/plots/Run{run_number}/"
htmldir = f"results/html/Run{run_number}/"

# suffix = "_subtracted_calibrated"
suffix = ""


def filterBandEvents(rdf):
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 2000 && FERS_SciEnergyHG{suffix} < 4000")


def filterDarkCountEvents(rdf):
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} < 200")


def filterElePeakEvents(rdf):
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 4e5 && FERS_SciEnergyHG{suffix} < 5e5")


def filterHEEvents(rdf):
    return rdf.Filter(f"FERS_SciEnergyHG{suffix} > 7e5")


def fillInValueMap(rdf, event_numbers):
    """
    Fill in the value map for the given event number.
    """
    # this should be lazy action
    # do everything in a loop
    rdfs_temp = []
    for event_number in event_numbers:
        # filter everything to one event, then use the mean value for readout
        # this achieves "lazy" evaluation
        rdf_temp = rdf.Filter(f"event_n == {event_number}")
        rdfs_temp.append(rdf_temp)

    values_FERS_events_lazy = {}
    for event_number, rdf_temp in zip(event_numbers, rdfs_temp):
        values_FERS = {}
        for _, FERSBoard in FERSBoards.items():
            for chan in FERSBoard:
                channelName_HG = chan.get_channel_name(gain="HG")
                variableName = f"{channelName_HG}{suffix}"
                value = rdf_temp.Mean(variableName)

                values_FERS[channelName_HG] = value
        # total energy for the event
        for var in ["Cer", "Sci"]:
            variableName = f"FERS_energyHG{suffix}_{var}_sum"
            value = rdf_temp.Mean(variableName)
            values_FERS[f"FERS_energyHG_{var}_sum"] = value
        values_FERS_events_lazy[event_number] = values_FERS

    hists_DRS_events_lazy = {}
    for event_number, rdf_temp in zip(event_numbers, rdfs_temp):
        hists_DRS_board_cer = OrderedDict()
        hists_DRS_board_sci = OrderedDict()
        suffix_plot = f"_Run{run_number}_Event{event_number}_DRS_vs_TS"
        for _, DRSBoard in DRSBoards.items():
            board_no = DRSBoard.board_no
            hists_DRS_board_cer[board_no] = OrderedDict()
            hists_DRS_board_sci[board_no] = OrderedDict()
            for group in [0, 1, 2, 3]:
                hists_DRS_board_cer[board_no][group] = OrderedDict()
                hists_DRS_board_sci[board_no][group] = OrderedDict()
            for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan is None:
                        continue
                    channelName = chan.get_channel_name()
                    group_no = chan.group_no
                    channel_no = chan.channel_no
                    # hist_subtractMedian = rdf_temp.Histo2D((
                    #    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}_subtractMedian_{suffix_plot}",
                    #    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);TS;{var} Variable",
                    #    1024, 0, 1024, 50, ymin, ymax),
                    #    "TS", channelName + "_subtractMedian"
                    # )
                    hist_subtractMedian = rdf_temp.Graph(
                        "TS", channelName + "_blsub")
                    if var == "Cer":
                        hists_DRS_board_cer[board_no][group_no][channel_no] = hist_subtractMedian
                    else:
                        hists_DRS_board_sci[board_no][group_no][channel_no] = hist_subtractMedian
            # add channel 8 (time reference)
            for group in [0, 1, 2, 3]:
                group_no = group
                # channel 8 for time reference
                hist_channel8 = rdf_temp.Graph(
                    "TS", f"DRS_Board{board_no}_Group{group_no}_Channel8_blsub")
                hists_DRS_board_cer[board_no][group_no][8] = hist_channel8
                hists_DRS_board_sci[board_no][group_no][8] = hist_channel8
        hists_DRS_events_lazy[event_number] = [
            hists_DRS_board_cer, hists_DRS_board_sci]

    # read all values
    values_FERS_events = {}
    for event_number, values_FERS in values_FERS_events_lazy.items():
        values_FERS_events[event_number] = {}
        for channelName_HG, value in values_FERS.items():
            value = value.GetValue()
            values_FERS_events[event_number][channelName_HG] = value

    return values_FERS_events, hists_DRS_events_lazy


def DrawEventDisplay(values_FERS_events, suffix="MIP", applyScaling=False, zmax=100, zmin=2):
    outdir_plots = f"{plotdir}/EventDisplay{suffix}/"
    plots = []
    for event_number, values_FERS in values_FERS_events.items():
        suffix_plot = f"_Run{run_number}_Event{event_number}_Display"
        [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
            FERSBoards, values_FERS, suffix=suffix_plot, gain="HG")

        if applyScaling:
            # scale 6mm channel to 3mm channel size
            h2_Cer_HG_mean.Scale(1/4.0)
            h2_Sci_HG_mean.Scale(1/4.0)

        output_name = f"EventDisplay_Run{run_number}_Event{event_number}"
        extraToDraw = ROOT.TPaveText(0.20, 0.85, 0.40, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        unit = "p.e." if "calibrated" in suffix else "ADC"
        extraToDraw.AddText(
            f"E_{{Cer}} = {values_FERS['FERS_energyHG_Cer_sum']:.0f} {unit}")
        makeEventDisplay(h2_Cer_HG_mean, h2_Cer_3mm_HG_mean, output_name +
                         "_Cer", outdir_plots, run_number, zmin, zmax, isCer=True, extraToDraw=extraToDraw)
        plots.append(output_name + "_Cer.png")

        extraToDraw.Clear()
        extraToDraw.AddText(
            f"E_{{Sci}} = {values_FERS['FERS_energyHG_Sci_sum']:.0f} {unit}")
        # make the event display for Sci
        makeEventDisplay(h2_Sci_HG_mean, h2_Sci_3mm_HG_mean, output_name +
                         "_Sci", outdir_plots, run_number, zmin, zmax, isCer=False, extraToDraw=extraToDraw)
        plots.append(output_name + "_Sci.png")

    output_html = f"{htmldir}/EventDisplay/FERS_{suffix}.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def makeDRS2DPlots(hists_DRS_events):
    outdir_plots = f"{plotdir}/DRS_vs_TS{suffix}/"
    colors = {
        0: ROOT.kRed,
        1: ROOT.kBlue,
        2: ROOT.kGreen+2,
        3: ROOT.kMagenta,
        4: ROOT.kCyan+2,
        5: ROOT.kOrange+7,
        6: ROOT.kGray+2,
        7: ROOT.kPink+9,
        8: ROOT.kBlack
    }
    output_htmls_drs = {}
    for event_number, (hists_DRS_board_cer, hists_DRS_board_sci) in hists_DRS_events.items():
        plots = []
        for hists in [hists_DRS_board_cer, hists_DRS_board_sci]:
            output_name = None
            var = "Cer" if hists == hists_DRS_board_cer else "Sci"
            for board, hists_board in hists.items():
                for group, hists_group in hists_board.items():
                    hists_to_draw = []
                    labels = []
                    for channel_no, hist in hists_group.items():
                        hist = hist.GetValue()
                        hist.SetLineColor(colors.get(channel_no, ROOT.kBlack))
                        hist.SetLineWidth(2)
                        hist.SetMarkerStyle(20)
                        hist.SetMarkerSize(0.5)
                        hist.SetMarkerColor(colors.get(channel_no, ROOT.kBlack))
                        hists_to_draw.append(hist)
                        labels.append(f"Ch{number_to_string(channel_no)}")
                    output_name = f"DRS_vs_TS_Run{run_number}_Event{event_number}_Board{board}_{var}_Group{group}"

                    extra_text = ROOT.TPaveText(0.20, 0.77, 0.90, 0.83, "NDC")
                    extra_text.SetTextAlign(11)
                    extra_text.SetFillColorAlpha(0, 0)
                    extra_text.SetBorderSize(0)
                    extra_text.SetTextFont(42)
                    extra_text.SetTextSize(0.04)
                    extra_text.AddText(
                        f"Evt {event_number} B {board} G {group}")

                    DrawHistos(hists_to_draw, labels, 0, 1024, "Time slice", -500, 2000, f"A.U.",
                               output_name,
                               dology=False, drawoptions="LP",
                               outdir=outdir_plots, extra_text="", run_number=run_number, legendNCols=5, legendPos=[0.20, 0.85, 0.90, 0.90], addOverflow=False, addUnderflow=False,
                               extraToDraw=extra_text, usePDF=True, usePNG=False)
                    plots.append(output_name + ".pdf")

        output_html = f"{htmldir}/EventDisplay/DRS_vs_TS{suffix}_event{event_number}.html"
        generate_html(plots, outdir_plots, plots_per_row=4,
                      output_html=output_html)
        output_htmls_drs[event_number] = output_html
    return output_htmls_drs


rdfs_filtered = {}
rdfs_filtered["Inc"] = rdf
# rdfs_filtered["ElePeak"] = filterElePeakEvents(rdf)
# rdfs_filtered["HE"] = filterHEEvents(rdf)
# rdfs_filtered["Band"] = filterBandEvents(rdf)
# rdfs_filtered["DarkCount"] = filterDarkCountEvents(rdf)

zranges = {
    "Inc": (2, 100),
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

random.seed(42)  # for reproducibility
for cat in event_numbers:
    selected_events = random.sample(
        event_numbers[cat], min(2, len(event_numbers[cat])))
    print(f"selected {cat} events: {selected_events}")

    values_FERS_events, values_DRS_events = fillInValueMap(
        rdfs_filtered[cat], selected_events)
    output_html = DrawEventDisplay(
        values_FERS_events, suffix=cat, applyScaling=True, zmin=zranges[cat][0], zmax=zranges[cat][1])
    output_htmls[cat] = output_html
    output_htmls_drs = makeDRS2DPlots(values_DRS_events)
    for event_number, html in output_htmls_drs.items():
        output_htmls[f"{cat}_DRS_Event{event_number}"] = html


print("Event display HTML files generated:")
for cat, html_file in output_htmls.items():
    print(f"{cat}: {html_file}")
print("All done!")
