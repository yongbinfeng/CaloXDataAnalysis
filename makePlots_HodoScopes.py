import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
import json
from results.events import events_interested
from utils.channel_map import get_hodoscope_channels

ROOT.gROOT.SetBatch(True)  # Run in batch mode
# multithread
ROOT.ROOT.EnableImplicitMT(10)

hodoscope_channels = get_hodoscope_channels()


with open("results/hodoscope_noises.json", "r") as json_file:
    hodo_noises = json.load(json_file)


def analyzePeak(infilename):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    histos1D_diff = {}
    histos1D_sum = {}
    histos1D_left = {}
    histos1D_right = {}
    histos2D_left_vs_right = {}
    histos1D_left_peak = {}
    histos1D_right_peak = {}
    for group, channels in hodoscope_channels.items():
        for channel in channels:
            rdf = rdf.Define(f"{channel}_subtracted",
                             f"{channel} - {hodo_noises['hodoscope_' + channel]}")
            # find the minimum index of the pulse shape
            rdf = rdf.Define(f"{channel}_peak_position",
                             f"ROOT::VecOps::ArgMin({channel}_subtracted)")
            rdf = rdf.Define(f"{channel}_peak_value",
                             f"ROOT::VecOps::Min({channel}_subtracted)")

    rdfs_filtered = []
    maps_mean = {}
    for group, channels in hodoscope_channels.items():
        if group != "trigger":
            rdf_filtered = rdf.Filter(
                f"({channels[0]}_peak_value < -100.0 ) && ({channels[1]}_peak_value < -100.0 )"
            )

            for channel in channels:
                # normalize the pulse shape
                rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_sum",
                                                   f"ROOT::VecOps::Sum({channel}_subtracted) + 1e-6")
                rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_norm",
                                                   f"{channel}_subtracted / {channel}_subtracted_sum")
            # calculate the difference between left and right peaks
            rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak",
                                               f"(int){channels[1]}_peak_position - (int){channels[0]}_peak_position")
            rdf_filtered = rdf_filtered.Define(f"{group}_sum_peak",
                                               f"(int){channels[0]}_peak_position + (int){channels[1]}_peak_position")
            rdfs_filtered.append((group, rdf_filtered))

            histos1D_diff[group] = rdf_filtered.Histo1D(
                (f"{group}_delta_peak",
                 f"Delta Peak {group};Peak Position Difference;Counts", 2048, -1024, 1024),
                f"{group}_delta_peak"
            )
            histos1D_sum[group] = rdf_filtered.Histo1D(
                (f"{group}_sum_peak",
                 f"Sum Peak {group};Peak Position Sum;Counts", 2048, 0, 2048),
                f"{group}_sum_peak"
            )
            histos1D_left[group] = rdf_filtered.Histo1D(
                (f"{group}_left_peak",
                 f"Left Peak {group};Peak Position;Counts", 1024, 0, 1024),
                f"{channels[0]}_peak_position"
            )
            histos1D_right[group] = rdf_filtered.Histo1D(
                (f"{group}_right_peak",
                 f"Right Peak {group};Peak Position;Counts", 1024, 0, 1024),
                f"{channels[1]}_peak_position"
            )
            histos2D_left_vs_right[group] = rdf_filtered.Histo2D(
                (f"{group}_left_peak_vs_right_peak",
                 f"Left vs Right Peak {group};Left Peak Position;Right Peak Position;Counts",
                 1024, 0, 1024, 1024, 0, 1024),
                f"{channels[0]}_peak_position",
                f"{channels[1]}_peak_position"
            )
            histos1D_left_peak[group] = rdf_filtered.Histo1D(
                (f"{group}_left_peak_value",
                 f"Left Peak Value {group};Peak Value;Counts", 2048, -2500, 1999),
                f"{channels[0]}_peak_value"
            )
            histos1D_right_peak[group] = rdf_filtered.Histo1D(
                (f"{group}_right_peak_value",
                 f"Right Peak Value {group};Peak Value;Counts", 2048, -2500, 1999),
                f"{channels[1]}_peak_value"
            )

            means = []
            for channel in channels:
                for i in range(0, 1024):
                    mean = rdf_filtered.Define(f"{channel}_subtracted_norm_{i}", f"{channel}_subtracted_norm[{i}]").Mean(
                        f"{channel}_subtracted_norm_{i}")
                    means.append(mean)
                maps_mean[channel] = means

    print("Average channel normalized pulse shapes calculated.")

    # save the means to TH1F histograms
    histos1D_means = {}
    for channel, means in maps_mean.items():
        histo = ROOT.TH1F(
            f"{channel}_means", f"Means of {channel};Time Slice;Mean Value", 1024, 0, 1024)
        for i, mean in enumerate(means):
            print(
                f"Channel {channel}, Time Slice {i}, Mean Value: {mean.GetValue()}")
            histo.SetBinContent(i + 1, mean.GetValue())
        histos1D_means[channel] = histo

    print("Writing histograms to output file...")
    ofile = ROOT.TFile("root/hodoscope_peaks.root", "RECREATE")
    if not ofile or ofile.IsZombie():
        raise RuntimeError(f"Failed to open output file: {ofile}")
    for group, hist in histos1D_diff.items():
        hist.SetDirectory(ofile)
        hist.Write()
        histos1D_sum[group].SetDirectory(ofile)
        histos1D_sum[group].Write()
        histos1D_left[group].SetDirectory(ofile)
        histos1D_left[group].Write()
        histos1D_right[group].SetDirectory(ofile)
        histos1D_right[group].Write()
        histos2D_left_vs_right[group].SetDirectory(ofile)
        histos2D_left_vs_right[group].Write()
        histos1D_left_peak[group].SetDirectory(ofile)
        histos1D_left_peak[group].Write()
        histos1D_right_peak[group].SetDirectory(ofile)
        histos1D_right_peak[group].Write()

    for group, hist in histos1D_means.items():
        histos1D_means[group].SetDirectory(ofile)
        histos1D_means[group].Write()

    ofile.Close()


def analyzeHodoPulse(infilename):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        print(f"Processing event {ievt + 1} of {rdf.Count().GetValue()}")
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        # if ievt > 10:
        #    break
        if evtNumber not in events_interested:
            print(
                f"Skipping event {evtNumber} as it is not in the interested events list.")
            continue
        # dump the pulse shapes for hodoscope channels
        h1_hodos = {}
        for group, channels in hodoscope_channels.items():
            hs_todraw = []
            for channel in channels:
                pulse_shape = rdf.Take["ROOT::VecOps::RVec<float>"](
                    channel).GetValue()[ievt]
                h1_hodos[channel] = ROOT.TH1F(
                    f"pulse_shape_{channel}_Evt{evtNumber}",
                    f"Pulse Shape {channel};Time;Amplitude",
                    1024, 0, 1024
                )
                for i in range(len(pulse_shape)):
                    h1_hodos[channel].Fill(
                        i, pulse_shape[i] - hodo_noises["hodoscope_" + channel])
                hs_todraw.append(h1_hodos[channel])

            if ievt < 100000:
                labels = []
                extraToDraw = None
                if group != "trigger":
                    labels = ["Left", "Right"]
                    peak_left = hs_todraw[0].GetMinimumBin()
                    peak_right = hs_todraw[1].GetMinimumBin()
                    deltaTS = peak_right - peak_left
                    extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
                    extraToDraw.SetTextAlign(11)
                    extraToDraw.SetFillColorAlpha(0, 0)
                    extraToDraw.SetBorderSize(0)
                    extraToDraw.SetTextFont(42)
                    extraToDraw.SetTextSize(0.04)
                    extraToDraw.AddText(
                        f"Event: {evtNumber}, {group}")
                    extraToDraw.AddText(f"Left Peak: {peak_left}")
                    extraToDraw.AddText(f"Right Peak: {peak_right}")
                    extraToDraw.AddText(
                        f"delta Peak: {deltaTS} ts")
                DrawHistos(hs_todraw, labels, 0, 1024, "TS",
                           -2500, 1999, "Amplitude", f"pulse_shape_Hodoscopes_Evt{evtNumber}_{group}", dology=False, mycolors=[2, 4], drawashist=True, extraToDraw=extraToDraw)
    print(f"Events left after filtering: {rdf.Count().GetValue()}")


if __name__ == "__main__":
    # input_file = "root/filtered_events_board1.root"
    # print(f"Processing file: {input_file}")
    # analyzeHodoPulse(input_file)

    input_file = "/Users/yfeng/Desktop/TTU/CaloX/Data/run316_250517140056_converted.root"
    print(f"Processing file: {input_file}")
    analyzePeak(input_file)
