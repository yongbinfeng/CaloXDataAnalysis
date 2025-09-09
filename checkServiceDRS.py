import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import getPreShowerChannel, getDownStreamMuonChannel, buildHodoPosChannels, getCerenkovCounters
from utils.html_generator import generate_html
from utils.utils import loadRDF, preProcessDRSBoards
from utils.utils import calculateEnergySumFERS, vectorizeFERS
from utils.channel_map import buildFERSBoards, buildDRSBoards
from configs.plotranges import getServiceDRSProcessedInfoRanges
from selections.selections import checkUpstreamVeto, applyUpstreamVeto, getCC1SumCutValue, getPSDSumCutValue
from utils.parser import get_args
import time

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
    rdf_old = calculateEnergySumFERS(
        rdf, FERSBoards, calibrate=False, subtractPedestal=False, clip=False)
    # rdf = rdf_old.Filter(f"FERS_SciEnergyHG > 7e5")

    hists = {}
    hists2d = {}

    channel_preshower = getPreShowerChannel(runNumber)
    channel_muon = getDownStreamMuonChannel(runNumber)
    channels_cerenkov = getCerenkovCounters(runNumber)

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

    for channel in channels_cerenkov:
        rdf = rdf.Define(f"{channel}_peak_position",
                         f"ROOT::VecOps::ArgMin({channel}_subtractMedian)")
        rdf = rdf.Define(f"{channel}_peak_value",
                         f"ROOT::VecOps::Min({channel}_subtractMedian)")
        rdf = rdf.Define(f"{channel}_sum",
                         f"SumRange({channel}_subtractMedian, 600, 800)")

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

    if "Cerenkov1" in names and "preshower" in names:
        idx_cer1 = names.index("Cerenkov1")
        idx_pre = names.index("preshower")
        channel_cer1 = channels[idx_cer1]
        channel_pre = channels[idx_pre]
        hists2d["Cerenkov1_vs_preshower"] = {}
        xmin, xmax = getServiceDRSProcessedInfoRanges("preshower", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("Cerenkov1", "sum")
        hists2d["Cerenkov1_vs_preshower"]["sum2D"] = rdf.Histo2D(
            ("Cerenkov1_sum_vs_preshower_sum",
             "Cerenkov1 vs Preshower Sum;Preshower Sum;Cerenkov1 Sum;Counts",
             500, xmin, xmax, 500, ymin, ymax),
            f"{channel_pre}_sum",
            f"{channel_cer1}_sum"
        )

    if "muon" in names and "preshower" in names:
        idx_muon = names.index("muon")
        idx_pre = names.index("preshower")
        channel_muon = channels[idx_muon]
        channel_pre = channels[idx_pre]
        hists2d["muon_vs_preshower"] = {}
        xmin, xmax = getServiceDRSProcessedInfoRanges("preshower", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("muon", "sum")
        hists2d["muon_vs_preshower"]["sum2D"] = rdf.Histo2D(
            ("muon_sum_vs_preshower_sum",
             "muon vs preshower Sum;preshower Sum;muon Sum;Counts",
             500, xmin, xmax, 500, ymin, ymax),
            f"{channel_pre}_sum",
            f"{channel_muon}_sum"
        )

    if "muon" in names and "Cerenkov1" in names:
        idx_muon = names.index("muon")
        idx_cer1 = names.index("Cerenkov1")
        channel_muon = channels[idx_muon]
        channel_cer1 = channels[idx_cer1]
        hists2d["muon_vs_Cerenkov1"] = {}
        xmin, xmax = getServiceDRSProcessedInfoRanges("Cerenkov1", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("muon", "sum")
        hists2d["muon_vs_Cerenkov1"]["sum2D"] = rdf.Histo2D(
            ("muon_sum_vs_Cerenkov1_sum",
             "muon vs Cerenkov1 Sum;Cerenkov1 Sum;muon Sum;Counts",
             500, xmin, xmax, 500, ymin, ymax),
            f"{channel_cer1}_sum",
            f"{channel_muon}_sum"
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
    rdf_old = calculateEnergySumFERS(
        rdf, FERSBoards, calibrate=False, subtractPedestal=False, clip=False)
    # rdf = rdf_old.Filter(f"FERS_SciEnergyHG > 7e5")

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
        extraToDraw = None
        if name == "preshower" or name == "Cerenkov1":
            ntot = hist.Integral(0, 10000)
            valCut = getPSDSumCutValue() if name == "preshower" else getCC1SumCutValue()
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
    if "Cerenkov1" in names and "preshower" in names:
        hist2d = infile.Get("Cerenkov1_sum_vs_preshower_sum")
        if not hist2d:
            print(
                f"Histogram Cerenkov1_sum_vs_preshower_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges("preshower", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("Cerenkov1", "sum")
        valCut_psd = getPSDSumCutValue()
        valCut_cc1 = getCC1SumCutValue()
        xPass = hist2d.GetXaxis().FindBin(valCut_psd)
        yPass = hist2d.GetYaxis().FindBin(valCut_cc1)
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
        extraToDraw.AddText(
            f"N (PSD > {valCut_psd:.2g}, CC > {valCut_cc1:.2g}): {nPP:.0f}")
        extraToDraw.AddText(
            f"N (PSD > {valCut_psd:.2g}, CC < {valCut_cc1:.2g}): {nPF:.0f}")
        extraToDraw.AddText(
            f"N (PSD < {valCut_psd:.2g}, CC > {valCut_cc1:.2g}): {nFP:.0f}")
        extraToDraw.AddText(
            f"N (PSD < {valCut_psd:.2g}, CC < {valCut_cc1:.2g}): {nFF:.0f}")
        DrawHistos([hist2d], "", xmin, xmax, "PSD Sum", ymin, ymax, "CC1 Sum",
                   outputname="Cerenkov1_vs_preshower_sum2D", outdir=outdir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, runNumber=runNumber, addOverflow=True, doth2=True, extraToDraw=extraToDraw
                   )
        plots.append("Cerenkov1_vs_preshower_sum2D.png")
    if "muon" in names and "preshower" in names:
        hist2d = infile.Get("muon_sum_vs_preshower_sum")
        if not hist2d:
            print(
                f"Histogram muon_sum_vs_preshower_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges("preshower", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("muon", "sum")
        DrawHistos([hist2d], "", xmin, xmax, "PSD Sum", ymin, ymax, "MuonVeto Sum",
                   outputname="Muon_vs_preshower_sum2D", outdir=outdir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, runNumber=runNumber, addOverflow=True, doth2=True,
                   )
        plots.append("Muon_vs_preshower_sum2D.png")

    if "muon" in names and "Cerenkov1" in names:
        hist2d = infile.Get("muon_sum_vs_Cerenkov1_sum")
        if not hist2d:
            print(
                f"Histogram muon_sum_vs_Cerenkov1_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges("Cerenkov1", "sum")
        ymin, ymax = getServiceDRSProcessedInfoRanges("muon", "sum")
        DrawHistos([hist2d], "", xmin, xmax, "CC1 Sum", ymin, ymax, "MuonVeto Sum",
                   outputname="Muon_vs_Cerenkov1_sum2D", outdir=outdir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, runNumber=runNumber, addOverflow=True, doth2=True,
                   )
        plots.append("Muon_vs_Cerenkov1_sum2D.png")

    output_html = f"results/html/Run{runNumber}/drs_service/viewer.html"
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

    output_html = f"results/html/Run{runNumber}/HodoPos/viewer.html"
    generate_html(plots, f"results/plots/Run{runNumber}/HodoPos/", plots_per_row=7,
                  output_html=output_html)

    return output_html


if __name__ == "__main__":
    start_time = time.time()
    chan_preshower = getPreShowerChannel(runNumber)
    chan_muon = getDownStreamMuonChannel()
    chans_cerenkov = getCerenkovCounters(runNumber)
    channels = [chan_preshower, chan_muon] + chans_cerenkov
    names = ["preshower", "muon"] + ["Cerenkov" +
                                     str(i) for i in range(1, len(chans_cerenkov) + 1)]

    analyzePulse(channels, names)
    # analyzeHodoPeak()

    outputs = {}
    outputs["PSD_Muon"] = plotPulse(channels, names)
    # outputs["hodo"] = plotHodoPeak()

    for key, value in outputs.items():
        print(f"Output for {key}: {value}")

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
