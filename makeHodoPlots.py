import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import buildHodoPosChannels, getUpstreamVetoChannel
from utils.html_generator import generate_html
from runconfig import runNumber, firstEvent, lastEvent
from utils.utils import loadRDF, loadRDF, preProcessDRSBoards
from selections.selections import checkUpstreamVeto

ROOT.gROOT.SetBatch(True)  # Run in batch mode
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")

hodo_pos_channels = buildHodoPosChannels(run=runNumber)


def analyzePeak():
    rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)

    rdf = preProcessDRSBoards(rdf)

    histos1D_diff = {}
    histos1D_diff_realtive = {}
    histos1D_sum = {}
    histos1D_left = {}
    histos1D_right = {}
    histos2D_left_vs_right = {}
    histos1D_left_peak = {}
    histos1D_right_peak = {}

    for group, channels in hodo_pos_channels.items():
        for channel in channels:
            # find the minimum index of the pulse shape
            rdf = rdf.Define(f"{channel}_peak_position",
                             f"ROOT::VecOps::ArgMin({channel}_subtractMedian)")
            rdf = rdf.Define(f"{channel}_peak_value",
                             f"ROOT::VecOps::Min({channel}_subtractMedian)")

    rdf = checkUpstreamVeto(rdf, runNumber)

    rdfs_filtered = []
    maps_mean = {}
    map_means_normalized = {}
    for group, channels in hodo_pos_channels.items():
        rdf_filtered = rdf.Filter(
            f"({channels[0]}_peak_value < -100.0 ) && ({channels[1]}_peak_value < -100.0 )"
        )

        for channel in channels:
            # normalize the pulse shape
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_sum",
                                               f"ROOT::VecOps::Sum({channel}_subtractMedian) + 1e-6")
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_norm",
                                               f"{channel}_subtractMedian / {channel}_subtracted_sum")
        # calculate the difference between left and right peaks
        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak",
                                           f"(int){channels[1]}_peak_position - (int){channels[0]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_sum_peak",
                                           f"(int){channels[0]}_peak_position + (int){channels[1]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak_relative",
                                           f"{group}_delta_peak / ({group}_sum_peak + 1e-6)")
        rdfs_filtered.append((group, rdf_filtered))

        for sel in ["pass_upstream_veto", "pass_NoSel"]:
            cat = f"{group}_{sel}"

            histos1D_diff[cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak",
                 f"Delta Peak {cat};Peak Position Difference;Counts", 256, -1024, 1024),
                f"{group}_delta_peak", sel
            )
            histos1D_diff_realtive[cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak_relative",
                 f"Delta Peak Relative {cat};Peak Position Difference Relative;Counts", 256, -1, 1),
                f"{group}_delta_peak_relative", sel
            )
            histos1D_sum[cat] = rdf_filtered.Histo1D(
                (f"{cat}_sum_peak",
                 f"Sum Peak {cat};Peak Position Sum;Counts", 256, 0, 2048),
                f"{group}_sum_peak", sel
            )
            histos1D_left[cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak",
                 f"Left Peak {cat};Peak Position;Counts", 256, 0, 1024),
                f"{channels[0]}_peak_position", sel
            )
            histos1D_right[cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak",
                 f"Right Peak {cat};Peak Position;Counts", 256, 0, 1024),
                f"{channels[1]}_peak_position", sel
            )
            histos2D_left_vs_right[cat] = rdf_filtered.Histo2D(
                (f"{cat}_left_peak_vs_right_peak",
                 f"Left vs Right Peak {cat};Left Peak Position;Right Peak Position;Counts",
                 256, 0, 1024, 256, 0, 1024),
                f"{channels[0]}_peak_position",
                f"{channels[1]}_peak_position", sel
            )
            histos1D_left_peak[cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak_value",
                 f"Left Peak Value {cat};Peak Value;Counts", 200, -2500, 1999),
                f"{channels[0]}_peak_value", sel
            )
            histos1D_right_peak[cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak_value",
                 f"Right Peak Value {cat};Peak Value;Counts", 200, -2500, 1999),
                f"{channels[1]}_peak_value", sel
            )

        means = []
        means_normalized = []
        for channel in channels:
            for i in range(0, 1024):
                mean_normalized = rdf_filtered.Define(f"{channel}_subtracted_norm_{i}", f"{channel}_subtracted_norm[{i}]").Mean(
                    f"{channel}_subtracted_norm_{i}")
                mean = rdf_filtered.Define(f"{channel}_subtractMedian_{i}", f"{channel}_subtractMedian[{i}]").Mean(
                    f"{channel}_subtractMedian_{i}")
                means.append(mean)
                means_normalized.append(mean_normalized)

            maps_mean[channel] = means
            map_means_normalized[channel] = means_normalized

    print("Average channel normalized pulse shapes calculated.")

    # save the means to TH1F histograms
    histos1D_means = {}
    histos1D_means_normalized = {}
    for channel, means in maps_mean.items():
        histo = ROOT.TH1F(
            f"{channel}_means", f"Means of {channel};Time Slice;Mean Value", 1024, 0, 1024)
        histo_normalized = ROOT.TH1F(
            f"{channel}_means_normalized", f"Means of {channel} normalized;Time Slice;Mean Value", 1024, 0, 1024)
        for i, mean in enumerate(means):
            # print(
            #    f"Channel {channel}, Time Slice {i}, Mean Value: {mean.GetValue()}")
            histo.SetBinContent(i + 1, mean.GetValue())
        for i, mean_normalized in enumerate(map_means_normalized[channel]):
            histo_normalized.SetBinContent(i + 1, mean_normalized.GetValue())
        histos1D_means[channel] = histo
        histos1D_means_normalized[channel] = histo_normalized

    print("Writing histograms to output file...")
    ofile = ROOT.TFile(
        f"results/root/Run{runNumber}/hodoscope_peaks.root", "RECREATE")
    if not ofile or ofile.IsZombie():
        raise RuntimeError(f"Failed to open output file: {ofile}")
    for cat, hist in histos1D_diff.items():
        hist.SetDirectory(ofile)
        hist.Write()
        histos1D_diff_realtive[cat].SetDirectory(ofile)
        histos1D_diff_realtive[cat].Write()
        histos1D_sum[cat].SetDirectory(ofile)
        histos1D_sum[cat].Write()
        histos1D_left[cat].SetDirectory(ofile)
        histos1D_left[cat].Write()
        histos1D_right[cat].SetDirectory(ofile)
        histos1D_right[cat].Write()
        histos2D_left_vs_right[cat].SetDirectory(ofile)
        histos2D_left_vs_right[cat].Write()
        histos1D_left_peak[cat].SetDirectory(ofile)
        histos1D_left_peak[cat].Write()
        histos1D_right_peak[cat].SetDirectory(ofile)
        histos1D_right_peak[cat].Write()

    for group, hist in histos1D_means.items():
        histos1D_means[group].SetDirectory(ofile)
        histos1D_means[group].Write()

    for group, hist in histos1D_means_normalized.items():
        histos1D_means_normalized[group].SetDirectory(ofile)
        histos1D_means_normalized[group].Write()

    ofile.Close()


def plotPeak():
    infile = ROOT.TFile(
        f"results/root/Run{runNumber}/hodoscope_peaks.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    plots = []
    outdir = f"results/plots/Run{runNumber}/HodoPos/"
    for group, channels in hodo_pos_channels.items():
        histos1D_diff = []
        histos1D_diff_realtive = []
        histos1D_sum = []
        histos1D_left = []
        histos1D_right = []
        histos2D_left_vs_right = []
        histos1D_left_peak = []
        histos1D_right_peak = []
        for sel in ["pass_NoSel", "pass_upstream_veto"]:
            cat = f"{group}_{sel}"
            hdiff = infile.Get(f"{cat}_delta_peak")
            hdiff_relat = infile.Get(f"{cat}_delta_peak_relative")
            hsum = infile.Get(f"{cat}_sum_peak")
            hleft = infile.Get(f"{cat}_left_peak")
            hright = infile.Get(f"{cat}_right_peak")
            hleft_vs_right = infile.Get(f"{cat}_left_peak_vs_right_peak")
            hleft_peak_value = infile.Get(f"{cat}_left_peak_value")
            hright_peak_value = infile.Get(f"{cat}_right_peak_value")

            histos1D_diff.append(hdiff)
            histos1D_diff_realtive.append(hdiff_relat)
            histos1D_sum.append(hsum)
            histos1D_left.append(hleft)
            histos1D_right.append(hright)
            histos2D_left_vs_right.append(hleft_vs_right)
            histos1D_left_peak.append(hleft_peak_value)
            histos1D_right_peak.append(hright_peak_value)

        labels = ["No Selection", "Veto Passed"]
        linestyles = [1, 2]

        outputname = f"{group}_diff_peak"
        DrawHistos(histos1D_diff, labels, -250, 250, "Peak Position Difference", 0, None, "Counts",
                   outputname=outputname, outdir=outdir,
                   dology=False, mycolors=[1, 1], drawashist=True, runNumber=runNumber, addOverflow=True, addUnderflow=True, linestyles=linestyles
                   )
        plots.append(outputname + ".png")

        outputname = f"{group}_diff_peak_relative"
        DrawHistos(histos1D_diff_realtive, labels, -0.5, 0.5, "Peak Position Difference Relative", 0, None, "Counts",
                   outputname=outputname, outdir=outdir,
                   dology=False, mycolors=[1, 1], drawashist=True, runNumber=runNumber, addOverflow=True, addUnderflow=True, linestyles=linestyles
                   )
        plots.append(outputname + ".png")

        outputname = f"{group}_sum_peak"
        DrawHistos(
            histos1D_sum, labels, 600, 1200, "Peak Position Sum", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1], drawashist=True, runNumber=runNumber, addOverflow=True, addUnderflow=True, linestyles=linestyles
        )
        plots.append(outputname + ".png")

        outputname = f"{group}_peaks"
        DrawHistos(
            histos1D_left +
            histos1D_right, [
                "Left No Sel", "Left Pass Veto", "Right No Sel", "Right Pass Veto"], 300, 800, "Peak Position", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1, 2, 2], drawashist=True, runNumber=runNumber, addOverflow=True, addUnderflow=True, linestyles=[1, 2, 1, 2]
        )
        plots.append(outputname + ".png")

        for idx, cat in enumerate(["pass_NoSel", "pass_upstream_veto"]):
            outputname = f"{group}_left_vs_right_peak_{cat}"
            DrawHistos(
                [histos2D_left_vs_right[idx]], [
                ], 0, 1024, "Left Peak Position", 0, 1024, "Right Peak Position",
                outputname=outputname, outdir=outdir,
                drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                dology=False, runNumber=runNumber, addOverflow=True, doth2=True,
            )
            plots.append(outputname + ".png")

        outputname = f"{group}_peak_values"
        DrawHistos(
            histos1D_left_peak + histos1D_right_peak, ["Left No Sel", "Left Pass Veto",
                                                       "Right No Sel", "Right Pass Veto"], -1500, 100, "Peak Value", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1, 2, 2], drawashist=True, runNumber=runNumber, addOverflow=True, addUnderflow=True, linestyles=[1, 2, 1, 2]
        )
        plots.append(outputname + ".png")
    generate_html(plots, f"results/plots/Run{runNumber}/HodoPos/", plots_per_row=7,
                  output_html=f"results/html/Run{runNumber}/HodoPos/viewer.html")


def analyzeHodoPulse(infilename):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    plots_toDraw = []

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        print(f"Processing event {ievt + 1} of {rdf.Count().GetValue()}")
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        if ievt > 30:
            break
        # if evtNumber not in events_interested:
        #    print(
        #        f"Skipping event {evtNumber} as it is not in the interested events list.")
        #    continue
        # dump the pulse shapes for hodoscope channels
        h1_hodos = {}
        for group, channels in hodo_pos_channels.items():
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
                output_name = f"pulse_shape_HodoPos_Evt{evtNumber}_{group}"
                DrawHistos(hs_todraw, labels, 0, 1024, "TS",
                           -2500, 1999, "Amplitude", output_name,
                           dology=False, mycolors=[2, 4], drawashist=True, extraToDraw=extraToDraw,
                           outdir=f"results/plots/Run{runNumber}/HodoPos/")
                plots_toDraw.append(
                    f"{output_name}.png")
    print(f"Events left after filtering: {rdf.Count().GetValue()}")
    generate_html(plots_toDraw, f"results/plots/Run{runNumber}/HodoPos/", plots_per_row=5,
                  output_html=f"results/html/Run{runNumber}/HodoPos/viewer.html")


if __name__ == "__main__":
    analyzePeak()
    plotPeak()
    sys.exit(0)
    # input_file = "root/filtered_events_board1.root"
    # print(f"Processing file: {input_file}")
    # analyzeHodoPulse(input_file)

    # input_file = "/Users/yfeng/Desktop/TTU/CaloX/Data/run316_250517140056_converted.root"
    inputfile_name = f"results/root/Run{runNumber}/filtered.root"
    print(f"Processing file: {inputfile_name}")
    analyzeHodoPulse(inputfile_name)
    # analyzePeak(inputfile_name)
