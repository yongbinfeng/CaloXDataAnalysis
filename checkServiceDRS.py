import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from channels.channel_map import getPreShowerChannel, getDownStreamMuonChannel, buildHodoPosChannels, getCerenkovCounters, buildFERSBoards, buildDRSBoards, getDownStreamTTUMuonChannel
from utils.html_generator import generate_html
from utils.dataloader import loadRDF
from variables.drs import preProcessDRSBoards
from variables.fers import vectorizeFERS
from configs.plotranges import getServiceDRSProcessedInfoRanges
from selections.selections import applyUpstreamVeto, applyUpstreamVeto, getServiceDRSSumCutValue
from utils.parser import get_args
from utils.timing import auto_timer  # noqa
auto_timer("Total Execution Time")

ROOT.gROOT.SetBatch(True)  # Run in batch mode
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")

parser = get_args()
runNumber = parser.run
firstEvent = parser.first_event
lastEvent = parser.last_event


def analyzePulse(channels, names):
    rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
    rdf = preProcessDRSBoards(rdf)
    rdf, _ = applyUpstreamVeto(rdf, runNumber)

    FERSBoards = buildFERSBoards(run=runNumber)
    rdf = vectorizeFERS(rdf, FERSBoards)

    hists = {}
    hists2d = {}

    channel_preshower = getPreShowerChannel(runNumber)
    # channel_muon = getDownStreamMuonChannel(runNumber)
    channel_muon = getDownStreamTTUMuonChannel(runNumber)
    channels_cerenkov = getCerenkovCounters(runNumber)

    rdf = rdf.Define(f"{channel_preshower}_peak_position",
                     f"ArgMinRange({channel_preshower}_blsub, 100, 400)")
    rdf = rdf.Define(f"{channel_preshower}_peak_value",
                     f"MinRange({channel_preshower}_blsub, 100, 400)")
    rdf = rdf.Define(f"{channel_preshower}_sum",
                     f"SumRange({channel_preshower}_blsub, 100, 400)")

    rdf = rdf.Define(f"{channel_muon}_peak_position",
                     f"ROOT::VecOps::ArgMin({channel_muon}_blsub)")
    rdf = rdf.Define(f"{channel_muon}_peak_value",
                     f"ROOT::VecOps::Min({channel_muon}_blsub)")
    rdf = rdf.Define(f"{channel_muon}_sum",
                     f"SumRange({channel_muon}_blsub, {channel_muon}_peak_position - 50, {channel_muon}_peak_position + 50)")

    for channel in channels_cerenkov:
        rdf = rdf.Define(f"{channel}_peak_position",
                         f"ROOT::VecOps::ArgMin({channel}_blsub)")
        rdf = rdf.Define(f"{channel}_peak_value",
                         f"ROOT::VecOps::Min({channel}_blsub)")
        rdf = rdf.Define(f"{channel}_sum",
                         f"SumRange({channel}_blsub, 0, 1000)")

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
             f"Sum {channel};ADC Counts;Counts", 500, xmin, xmax),
            f"{channel}_sum"
        )

    for idx1, name1 in enumerate(names):
        for idx2, name2 in enumerate(names):
            if idx2 <= idx1:
                continue
            channel1 = channels[idx1]
            channel2 = channels[idx2]
            hists2d[f"{name1}_vs_{name2}"] = {}
            xmin, xmax = getServiceDRSProcessedInfoRanges(name1, "sum")
            ymin, ymax = getServiceDRSProcessedInfoRanges(name2, "sum")
            hists2d[f"{name1}_vs_{name2}"]["sum2D"] = rdf.Histo2D(
                (f"{name1}_sum_vs_{name2}_sum",
                 f"{name1} vs {name2} Sum;{name1} Sum;{name2} Sum;Counts",
                 500, xmin, xmax, 500, ymin, ymax),
                f"{channel1}_sum",
                f"{channel2}_sum"
            )

    print("Writing histograms to output file...")
    ofile = ROOT.TFile(
        f"results/root/Run{runNumber}/drs_service.root", "RECREATE")
    for name, channel in zip(names, channels):
        for _, hist in hists[channel].items():
            hist.SetDirectory(ofile)
            hist.Write()
    for _, map2d in hists2d.items():
        for _, hist2d in map2d.items():
            hist2d.SetDirectory(ofile)
            hist2d.Write()
    ofile.Close()


def analyzeHodoPeak():
    hodo_pos_channels = buildHodoPosChannels(run=runNumber)
    rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)

    rdf = preProcessDRSBoards(rdf)

    FERSBoards = buildFERSBoards(run=runNumber)
    rdf = vectorizeFERS(rdf, FERSBoards)

    histos1D_diff = {}
    histos1D_diff_realtive = {}
    histos1D_sum = {}
    histos1D_left = {}
    histos1D_right = {}
    histos2D_left_vs_right = {}
    histos1D_left_peak = {}
    histos1D_right_peak = {}
    histos2D_LR_vs_UD = {}

    for group, channels in hodo_pos_channels.items():
        for channel in channels:
            # find the minimum index of the pulse shape
            rdf = rdf.Define(f"{channel}_peak_position",
                             f"ROOT::VecOps::ArgMin({channel}_blsub)")
            rdf = rdf.Define(f"{channel}_peak_value",
                             f"ROOT::VecOps::Min({channel}_blsub)")

    rdf = applyUpstreamVeto(rdf, runNumber, applyCut=False)

    conditions = ["pass_upstream_veto", "pass_NoSel"]

    # rdfs_filtered = []
    rdf_filtered = rdf
    maps_mean = {}
    map_means_normalized = {}
    for group, channels in hodo_pos_channels.items():
        # rdf_filtered = rdf.Filter(
        #    f"({channels[0]}_peak_value < -100.0 ) && ({channels[1]}_peak_value < -100.0 )"
        # )

        for channel in channels:
            # normalize the pulse shape
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_sum",
                                               f"ROOT::VecOps::Sum({channel}_blsub) + 1e-6")
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_norm",
                                               f"{channel}_blsub / {channel}_subtracted_sum")
        # calculate the difference between left and right peaks
        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak",
                                           f"(int){channels[1]}_peak_position - (int){channels[0]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_sum_peak",
                                           f"(int){channels[0]}_peak_position + (int){channels[1]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak_relative",
                                           f"{group}_delta_peak / ({group}_sum_peak + 1e-6)")
        # rdfs_filtered.append((group, rdf_filtered))

        for sel in conditions:
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
                mean = rdf_filtered.Define(f"{channel}_blsub_{i}", f"{channel}_blsub[{i}]").Mean(
                    f"{channel}_blsub_{i}")
                means.append(mean)
                means_normalized.append(mean_normalized)

            maps_mean[channel] = means
            map_means_normalized[channel] = means_normalized

    # LR vs UD
    for sel in conditions:
        for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
            if dwc[:3] not in hodo_pos_channels:
                continue
            cat = f"{dwc}_{sel}"
            histos2D_LR_vs_UD[cat] = rdf_filtered.Histo2D(
                (f"{cat}", f"LR vs UD", 256, -1, 1, 256, -1, 1),
                f"{dwc.split('_vs_')[0]}_delta_peak_relative",
                f"{dwc.split('_vs_')[1]}_delta_peak_relative", sel
            )

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

    for cat, hist in histos2D_LR_vs_UD.items():
        hist.SetDirectory(ofile)
        hist.Write()

    for group, hist in histos1D_means.items():
        histos1D_means[group].SetDirectory(ofile)
        histos1D_means[group].Write()

    for group, hist in histos1D_means_normalized.items():
        histos1D_means_normalized[group].SetDirectory(ofile)
        histos1D_means_normalized[group].Write()

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
                   addOverflow=True, addUnderflow=True, leftlegend=True)
        plots.append(f"{name}_peak_value.png")

        hist = infile.Get(f"{name}_sum")
        if not hist:
            print(f"Histogram {name}_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges(name, "sum")
        ntot = hist.Integral(0, 10000)
        valCut = getServiceDRSSumCutValue(name)
        nhad = hist.Integral(hist.FindBin(valCut), 10000)
        nele = hist.Integral(0, hist.FindBin(valCut))
        extraToDraw = ROOT.TPaveText(0.23, 0.75, 0.55, 0.85, "NDC")
        extraToDraw.SetFillColor(0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"N (sum > {valCut:.2g}): {nhad:.0f}")
        extraToDraw.AddText(f"N (sum < {valCut:.2g}): {nele:.0f}")
        DrawHistos([hist], [name], xmin, xmax, "Sum", 1, None, "Counts",
                   outputname=f"{name}_sum", outdir=outdir,
                   dology=True, mycolors=[1], drawashist=True, runNumber=runNumber,
                   addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw)
        plots.append(f"{name}_sum.png")

    # 2D plots
    for idx1, name1 in enumerate(names):
        for idx2, name2 in enumerate(names):
            if idx2 <= idx1:
                continue
            hist2d = infile.Get(f"{name1}_sum_vs_{name2}_sum")
            if not hist2d:
                print(
                    f"Histogram {name1}_sum_vs_{name2}_sum not found in {infile.GetName()}")
                return
            xmin, xmax = getServiceDRSProcessedInfoRanges(name1, "sum")
            ymin, ymax = getServiceDRSProcessedInfoRanges(name2, "sum")
            valCut1 = getServiceDRSSumCutValue(name1)
            valCut2 = getServiceDRSSumCutValue(name2)
            xPass = hist2d.GetXaxis().FindBin(valCut1)
            yPass = hist2d.GetYaxis().FindBin(valCut2)
            nPP = hist2d.Integral(xPass, 1000, yPass, 1000)
            nPF = hist2d.Integral(xPass, 1000, 0, yPass - 1)
            nFP = hist2d.Integral(0, xPass - 1, yPass, 1000)
            nFF = hist2d.Integral(0, xPass - 1, 0, yPass - 1)
            extraToDraw = ROOT.TPaveText(0.23, 0.20, 0.5, 0.40, "NDC")
            extraToDraw.SetFillColorAlpha(0, 0.0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.03)
            name1_txt = name1.replace("Cerenkov", "Cer")
            name2_txt = name2.replace("Cerenkov", "Cer")
            extraToDraw.AddText(
                f"N ({name1_txt} > {valCut1:.2g}, {name2_txt} > {valCut2:.2g}): {nPP:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} > {valCut1:.2g}, {name2_txt} < {valCut2:.2g}): {nPF:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} < {valCut1:.2g}, {name2_txt} > {valCut2:.2g}): {nFP:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} < {valCut1:.2g}, {name2_txt} < {valCut2:.2g}): {nFF:.0f}")
            DrawHistos([hist2d], "", xmin, xmax, f"{name1} Sum", ymin, ymax, f"{name2} Sum",
                       outputname=f"{name1}_vs_{name2}_sum2D", outdir=outdir,
                       drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                       dology=False, runNumber=runNumber, addOverflow=True, doth2=True,
                       extraToDraw=extraToDraw)
            plots.append(f"{name1}_vs_{name2}_sum2D.png")

    output_html = f"results/html/Run{runNumber}/ServiceDRS/PID.html"
    generate_html(plots, f"results/plots/Run{runNumber}/drs_service/", plots_per_row=3,
                  output_html=output_html)

    return output_html


def plotHodoPeak():
    infile = ROOT.TFile(
        f"results/root/Run{runNumber}/hodoscope_peaks.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    hodo_pos_channels = buildHodoPosChannels(run=runNumber)

    plots = []
    outdir = f"results/plots/Run{runNumber}/DWC/"
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

    # LR vs UD
    plots_summary = []
    for sel in ["pass_NoSel", "pass_upstream_veto"]:
        for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
            cat = f"{dwc}_{sel}"
            hist2d = infile.Get(f"{cat}")
            if hist2d:
                outputname = f"DWC_{dwc}_{sel}"
                DrawHistos(
                    [hist2d], "", -0.4, 0.4, "LR Delta Peak Relative",
                    -0.4, 0.4, "UD Delta Peak Relative",
                    outputname=outputname, outdir=outdir,
                    drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                    dology=False, runNumber=runNumber, addOverflow=True, doth2=True,
                )
            plots_summary.append(outputname + ".png")
    plots_summary.append("NEWLINE")

    plots = plots_summary + plots
    output_html = f"results/html/Run{runNumber}/ServiceDRS/DWC.html"
    generate_html(plots, outdir, plots_per_row=7,
                  output_html=output_html)

    return output_html


if __name__ == "__main__":
    chan_preshower = getPreShowerChannel(runNumber)
    # chan_muon = getDownStreamMuonChannel()
    chan_muon = getDownStreamTTUMuonChannel(runNumber)
    chans_cerenkov = getCerenkovCounters(runNumber)
    channels = [chan_preshower, chan_muon] + chans_cerenkov
    names = ["preshower", "muon"] + ["Cerenkov" +
                                     str(i) for i in range(1, len(chans_cerenkov) + 1)]

    analyzePulse(channels, names)
    analyzeHodoPeak()

    outputs = {}
    outputs["PSD_Muon"] = plotPulse(channels, names)
    outputs["DWC"] = plotHodoPeak()

    for key, value in outputs.items():
        print(f"Output for {key}: {value}")
