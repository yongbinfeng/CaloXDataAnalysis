import ROOT
from collections import OrderedDict
from plotting.my_function import DrawHistos
from channels.channel_map import build_hodo_pos_channels, get_service_drs_channels
from utils.html_generator import generate_html
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_ranges import getServiceDRSProcessedInfoRanges
from configs.selection_values import get_service_drs_cut
from utils.parser import get_args
from utils.plot_helper import save_hists_to_file
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()
            .apply_selections(flag_only=True))
rdf_org = analysis.get_rdf()
paths = analysis.paths


def analyzePulse(channels):
    # only look at the ones passing the hole veto
    # rdf = rdf_org.Filter("isHoleVetoFired == 0", "Pass Hole Veto")
    rdf = rdf_org

    hists = {}

    for det, channel in channels.items():
        ts_min, ts_max, value_cut, method = get_service_drs_cut(det)
        rdf = rdf.Define(f"{channel}_peak_position",
                         f"ArgMinRange({channel}_blsub, {ts_min}, {ts_max})")
        rdf = rdf.Define(f"{channel}_peak_value",
                         f"MinRange({channel}_blsub, {ts_min}, {ts_max})")
        rdf = rdf.Define(f"{channel}_sum",
                         f"SumRange({channel}_blsub, {ts_min}, {ts_max})")

    # rdf = rdf.Filter(f"{channels['PSD']}_sum > -1e3")
    # rdf = rdf.Filter(
    #    f"{channels['PSD']}_sum < -1e3").Filter(f"{channels['PSD']}_sum >= -5e3")
    # rdf = rdf.Filter(f"{channels['TTUMuonVeto']}_sum > -5e3")

    # create histograms
    for det, channel in channels.items():
        hists[channel] = {}
        hists[channel]["ADC_VS_TS"] = rdf.Histo2D(
            (f"{det}_ADC_vs_TS",
             f"ADC vs Time Slice {channel};Time Slice;ADC Counts;Counts",
             1024, 0, 1024, 500, -3000, 1000),
            "TS",
            f"{channel}_blsub"
        )
        hists[channel]["peak_position"] = rdf.Histo1D(
            (f"{det}_peak_position",
             f"Peak Position {channel};Time Slice;Counts", 128, 0, 1024),
            f"{channel}_peak_position"
        )
        xmin, xmax = getServiceDRSProcessedInfoRanges(
            det, "peak_value")
        hists[channel]["peak_value"] = rdf.Histo1D(
            (f"{det}_peak_value",
             f"Peak Value {channel};ADC Counts;Counts", 50, xmin, xmax),
            f"{channel}_peak_value"
        )
        xmin, xmax = getServiceDRSProcessedInfoRanges(det, "sum")
        hists[channel]["sum"] = rdf.Histo1D(
            (f"{det}_sum",
             f"Sum {channel};ADC Counts;Counts", 500, xmin, xmax),
            f"{channel}_sum"
        )

    for idx1, det1 in enumerate(channels.keys()):
        for idx2, det2 in enumerate(channels.keys()):
            if idx2 <= idx1:
                continue
            channel1 = channels[det1]
            channel2 = channels[det2]
            hists[f"{det1}_vs_{det2}"] = {}
            xmin, xmax = getServiceDRSProcessedInfoRanges(det1, "sum")
            ymin, ymax = getServiceDRSProcessedInfoRanges(det2, "sum")
            hists[f"{det1}_vs_{det2}"]["sum2D"] = rdf.Histo2D(
                (f"{det1}_sum_vs_{det2}_sum",
                 f"{det1} vs {det2} Sum;{det1} Sum;{det2} Sum;Counts",
                 500, xmin, xmax, 500, ymin, ymax),
                f"{channel1}_sum",
                f"{channel2}_sum"
            )

    output_hists = []
    for _, hists_map in hists.items():
        for _, hist in hists_map.items():
            output_hists.append(hist)

    return output_hists


def analyzeHodoPeak():
    rdf = rdf_org
    hodo_pos_channels = build_hodo_pos_channels(run=run_number)

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

    conditions = ["is_HoleVeto_vetoed", "passNone"]

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

    hists_output = []
    for hist_map in [histos1D_diff, histos1D_diff_realtive, histos1D_sum, histos1D_left,
                     histos1D_right, histos2D_left_vs_right,
                     histos1D_left_peak, histos1D_right_peak,
                     histos2D_LR_vs_UD,
                     histos1D_means, histos1D_means_normalized]:
        for _, hist in hist_map.items():
            hists_output.append(hist)

    return hists_output


def plotPulse(channels):
    infile = ROOT.TFile(
        f"{paths['root']}/drs_services_pid.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    output_htmls = []
    plots = []
    outdir = f"{paths['plots']}/drs_service/"

    for det in channels.keys():
        hist = infile.Get(f"{det}_ADC_vs_TS")
        if not hist:
            print(
                f"Histogram {det}_ADC_vs_TS not found in {infile.GetName()}")
            return
        DrawHistos([hist], [], 0, 1024, "Time Slice", -3000, 1000, "ADC Counts",
                   outputname=f"{det}_ADC_vs_TS", outdir=outdir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, run_number=run_number, addOverflow=True, doth2=True)
        plots.append(f"{det}_ADC_vs_TS.png")

        hist = infile.Get(f"{det}_peak_position")
        if not hist:
            print(
                f"Histogram {det}_peak_position not found in {infile.GetName()}")
            return
        DrawHistos([hist], [det], 0, 1024, "Peak Position", 0, None, "Counts",
                   outputname=f"{det}_peak_position", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, run_number=run_number,
                   addOverflow=True, addUnderflow=True)
        plots.append(f"{det}_peak_position.png")

        hist = infile.Get(f"{det}_peak_value")
        if not hist:
            print(
                f"Histogram {det}_peak_value not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges(det, "peak_value")
        DrawHistos([hist], [det], xmin, xmax, "Peak Value", 0, None, "Counts",
                   outputname=f"{det}_peak_value", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, run_number=run_number,
                   addOverflow=True, addUnderflow=True, leftlegend=True)
        plots.append(f"{det}_peak_value.png")

        hist = infile.Get(f"{det}_sum")
        if not hist:
            print(f"Histogram {det}_sum not found in {infile.GetName()}")
            return
        xmin, xmax = getServiceDRSProcessedInfoRanges(det, "sum")
        _, _, value_cut, _ = get_service_drs_cut(det)
        nhad = hist.Integral(hist.FindBin(value_cut), 10000)
        nele = hist.Integral(0, hist.FindBin(value_cut))
        extraToDraw = ROOT.TPaveText(0.23, 0.75, 0.55, 0.85, "NDC")
        extraToDraw.SetFillColor(0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"N (sum > {value_cut:.2g}): {nhad:.0f}")
        extraToDraw.AddText(f"N (sum < {value_cut:.2g}): {nele:.0f}")
        DrawHistos([hist], [det], xmin, xmax, "Sum", 1, None, "Counts",
                   outputname=f"{det}_sum", outdir=outdir,
                   dology=True, mycolors=[1], drawashist=True, run_number=run_number,
                   addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw)
        plots.append(f"{det}_sum.png")

        # make a cdf plot based on sum
        # include over/underflow
        hist.GetXaxis().SetRange(0, hist.GetNbinsX()+1)
        hist_cdf = hist.GetCumulative()
        hist_cdf.Scale(1.0 / hist.Integral())
        DrawHistos([hist_cdf], [det], xmin, xmax, "Sum", 0, 1.3, "Cumulative Fraction",
                   outputname=f"{det}_sum_cdf", outdir=outdir,
                   dology=False, mycolors=[1], drawashist=True, run_number=run_number,
                   addOverflow=False, addUnderflow=False, legendPos=[0.3, 0.80, 0.5, 0.85])
        plots.append(f"{det}_sum_cdf.png")

    output_html = f"{paths['html']}/ServiceDRS/PID.html"
    generate_html(plots, outdir, plots_per_row=5,
                  output_html=output_html)
    output_htmls.append(output_html)

    plots = []

    # 2D plots
    for idx1, det1 in enumerate(channels.keys()):
        for idx2, det2 in enumerate(channels.keys()):
            if idx2 <= idx1:
                continue
            hist2d = infile.Get(f"{det1}_sum_vs_{det2}_sum")
            if not hist2d:
                print(
                    f"Histogram {det1}_sum_vs_{det2}_sum not found in {infile.GetName()}")
                return
            xmin, xmax = getServiceDRSProcessedInfoRanges(det1, "sum")
            ymin, ymax = getServiceDRSProcessedInfoRanges(det2, "sum")
            _, _, value_cut1, _ = get_service_drs_cut(det1)
            _, _, value_cut2, _ = get_service_drs_cut(det2)
            xPass = hist2d.GetXaxis().FindBin(value_cut1)
            yPass = hist2d.GetYaxis().FindBin(value_cut2)
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
            name1_txt = det1.replace("Cerenkov", "Cer")
            name2_txt = det2.replace("Cerenkov", "Cer")
            extraToDraw.AddText(
                f"N ({name1_txt} > {value_cut1:.2g}, {name2_txt} > {value_cut2:.2g}): {nPP:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} > {value_cut1:.2g}, {name2_txt} < {value_cut2:.2g}): {nPF:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} < {value_cut1:.2g}, {name2_txt} > {value_cut2:.2g}): {nFP:.0f}")
            extraToDraw.AddText(
                f"N ({name1_txt} < {value_cut1:.2g}, {name2_txt} < {value_cut2:.2g}): {nFF:.0f}")
            DrawHistos([hist2d], "", xmin, xmax, f"{det1} Sum", ymin, ymax, f"{det2} Sum",
                       outputname=f"{det1}_vs_{det2}_sum2D", outdir=outdir,
                       drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                       dology=False, run_number=run_number, addOverflow=True, doth2=True,
                       extraToDraw=extraToDraw)
            plots.append(f"{det1}_vs_{det2}_sum2D.png")
        plots.append("NEWLINE")

    output_html = f"{paths['html']}/ServiceDRS/PID_correlation.html"
    generate_html(plots, outdir, plots_per_row=6,
                  output_html=output_html)

    output_htmls.append(output_html)

    return output_htmls


def plotHodoPeak():
    infile = ROOT.TFile(
        f"{paths['root']}/hodoscope_peaks.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    hodo_pos_channels = build_hodo_pos_channels(run=run_number)

    plots = []
    outdir = f"{paths['plots']}/DWC/"
    for group, channels in hodo_pos_channels.items():
        histos1D_diff = []
        histos1D_diff_realtive = []
        histos1D_sum = []
        histos1D_left = []
        histos1D_right = []
        histos2D_left_vs_right = []
        histos1D_left_peak = []
        histos1D_right_peak = []
        for sel in ["passNone", "is_HoleVeto_vetoed"]:
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

        labels = ["No Selection", "Pass HoleVeto"]
        linestyles = [1, 2]

        outputname = f"{group}_diff_peak"
        DrawHistos(histos1D_diff, labels, -250, 250, "Peak Position Difference", 0, None, "Counts",
                   outputname=outputname, outdir=outdir,
                   dology=False, mycolors=[1, 1], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, linestyles=linestyles
                   )
        plots.append(outputname + ".png")

        outputname = f"{group}_diff_peak_relative"
        DrawHistos(histos1D_diff_realtive, labels, -0.5, 0.5, "Peak Position Difference Relative", 0, None, "Counts",
                   outputname=outputname, outdir=outdir,
                   dology=False, mycolors=[1, 1], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, linestyles=linestyles
                   )
        plots.append(outputname + ".png")

        outputname = f"{group}_sum_peak"
        DrawHistos(
            histos1D_sum, labels, 600, 1200, "Peak Position Sum", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, linestyles=linestyles
        )
        plots.append(outputname + ".png")

        outputname = f"{group}_peaks"
        DrawHistos(
            histos1D_left +
            histos1D_right, [
                "Left No Sel", "Left Pass Veto", "Right No Sel", "Right Pass Veto"], 300, 800, "Peak Position", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1, 2, 2], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, linestyles=[1, 2, 1, 2]
        )
        plots.append(outputname + ".png")

        for idx, cat in enumerate(["pass_NoSel", "pass_upstream_veto"]):
            outputname = f"{group}_left_vs_right_peak_{cat}"
            DrawHistos(
                [histos2D_left_vs_right[idx]], [
                ], 0, 1024, "Left Peak Position", 0, 1024, "Right Peak Position",
                outputname=outputname, outdir=outdir,
                drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                dology=False, run_number=run_number, addOverflow=True, doth2=True,
            )
            plots.append(outputname + ".png")

        outputname = f"{group}_peak_values"
        DrawHistos(
            histos1D_left_peak + histos1D_right_peak, ["Left No Sel", "Left Pass Veto",
                                                       "Right No Sel", "Right Pass Veto"], -1500, 100, "Peak Value", 0, None, "Counts",
            outputname=outputname, outdir=outdir,
            dology=False, mycolors=[1, 1, 2, 2], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, linestyles=[1, 2, 1, 2]
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
                    dology=False, run_number=run_number, addOverflow=True, doth2=True,
                )
            plots_summary.append(outputname + ".png")
    plots_summary.append("NEWLINE")

    plots = plots_summary + plots
    output_html = f"{paths['html']}/ServiceDRS/DWC.html"
    generate_html(plots, outdir, plots_per_row=7,
                  output_html=output_html)

    return output_html


def main():
    channels = OrderedDict()
    for det in ["HoleVeto", "PSD", "TTUMuonVeto", "Cer474", "Cer519", "Cer537"]:
        channel = get_service_drs_channels(run_number).get(det)
        channels[det] = channel
    hists_pid = analyzePulse(channels)
    hists_hodo = analyzeHodoPeak()

    # save histograms to files
    save_hists_to_file(hists_pid,
                       f"{paths['root']}/drs_services_pid.root")

    save_hists_to_file(hists_hodo,
                       f"{paths['root']}/hodoscope_peaks.root")

    # make plots
    outputs = {}
    outputs["PID"] = plotPulse(channels)
    outputs["DWC"] = plotHodoPeak()

    for key, value in outputs.items():
        print(f"Output for {key}: {value}")


if __name__ == "__main__":
    main()
