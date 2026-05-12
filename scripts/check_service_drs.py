"""
Service DRS Check Script.

Analyzes service DRS channels for particle identification and hodoscope peak timing.
"""

import ROOT
from collections import OrderedDict
from channels.channel_map import build_hodo_pos_channels, get_service_drs_channels, get_mcp_channels
from configs.plot_config import get_service_drs_processed_info_ranges
from configs.selection_config import get_service_drs_cut
from core.analysis_manager import CaloXAnalysisManager
from utils.parser import get_args
from utils.plot_helper import save_hists_to_file
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from plotting.calox_plot_helper import create_pave_text
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()
            .apply_hole_veto(flag_only=True)
            )
paths = analysis.paths
rdf_org = analysis.get_rdf()

# Common styles
STYLE_1D = PlotStyle(dology=False, drawoptions="HIST",
                     mycolors=[1])
STYLE_1D_LOG = PlotStyle(dology=True, drawoptions="HIST",
                         mycolors=[1])
STYLE_2D_LOG = PlotStyle(dology=False, dologz=True,
                         drawoptions="COLz", zmin=1, zmax=None)
STYLE_COMPARE = PlotStyle(dology=False, drawoptions="HIST", mycolors=[
                          1, 1])


def analyze_pulse(channels):
    """Analyze pulse shapes for PID channels using precomputed DRS variables."""
    rdf = rdf_org
    hists = {}

    # Use precomputed columns from process_drs_channels / process_mcp_channels:
    #   {channel}_TS_cfd    — CFD crossing time slice (float)
    #   {channel}_peak_value — waveform peak amplitude
    #   {channel}_energy     — CFD energy integral
    for det, channel in channels.items():
        # MCP columns are keyed by det name; regular service DRS by channel name
        # col = det if det.startswith("MCP") else channel
        col = channel
        hists[channel] = {}
        wf_ymin, wf_ymax = get_service_drs_processed_info_ranges(
            det, "waveform")
        hists[channel]["ADC_VS_TS"] = rdf.Histo2D(
            (f"{det}_ADC_vs_TS", f"ADC vs Time Slice {channel};Time Slice;ADC Counts;Counts",
             1024, 0, 1024, 500, wf_ymin, wf_ymax),
            "ts", f"{channel}_blsub")
        hists[channel]["ADC_VS_TS_prof"] = rdf.Profile1D(
            (f"{det}_ADC_vs_TS_prof", f"Mean ADC vs Time Slice {channel};Time Slice;Mean ADC Counts",
             1024, 0, 1024, wf_ymin, wf_ymax),
            "ts", f"{channel}_blsub")
        hists[channel]["cfd_ts"] = rdf.Histo1D(
            (f"{det}_cfd_ts",
             f"CFD TS {channel};CFD TS;Counts", 128, 0, 1024),
            f"{col}_TS_cfd")
        xmin, xmax = get_service_drs_processed_info_ranges(det, "peak_value")
        hists[channel]["peak_value_vs_cfd_ts"] = rdf.Histo2D(
            (f"{det}_peak_value_vs_cfd_ts",
             f"Peak Value vs CFD TS {channel};CFD TS;Peak Value;Counts",
             128, 0, 1024, 100, xmin, xmax),
            f"{col}_TS_cfd", f"{col}_peak_value")
        hists[channel]["peak_value"] = rdf.Histo1D(
            (f"{det}_peak_value",
             f"Peak Value {channel};ADC Counts;Counts", 50, xmin, xmax),
            f"{col}_peak_value")
        xmin, xmax = get_service_drs_processed_info_ranges(det, "sum")
        hists[channel]["energy"] = rdf.Histo1D(
            (f"{det}_energy", f"Energy {channel};Energy (ADC);Counts", 500, xmin, xmax),
            f"{col}_energy")

    # 2D correlation histograms
    det_list = list(channels.keys())
    for idx1, det1 in enumerate(det_list):
        for idx2, det2 in enumerate(det_list):
            if idx2 <= idx1:
                continue
            special_dets = {"KT1", "KT2", "T3", "T4"}
            if (det1 in special_dets) ^ (det2 in special_dets):
                continue  # Skip correlations between trigger and non-trigger detectors
            channel1, channel2 = channels[det1], channels[det2]
            col1 = det1 if det1.startswith("MCP") else channel1
            col2 = det2 if det2.startswith("MCP") else channel2
            _, _, _, method1 = get_service_drs_cut(det1)
            _, _, _, method2 = get_service_drs_cut(det2)
            var1 = "peak_value" if method1 == "PeakValue" else "energy"
            var2 = "peak_value" if method2 == "PeakValue" else "energy"
            xmin, xmax = get_service_drs_processed_info_ranges(
                det1, "peak_value" if method1 == "PeakValue" else "sum")
            ymin, ymax = get_service_drs_processed_info_ranges(
                det2, "peak_value" if method2 == "PeakValue" else "sum")
            hists[f"{det1}_vs_{det2}"] = {}
            hists[f"{det1}_vs_{det2}"]["corr2D"] = rdf.Histo2D(
                (f"{det1}_{var1}_vs_{det2}_{var2}",
                 f"{det1} vs {det2};{det1} {var1};{det2} {var2};Counts",
                 500, xmin, xmax, 500, ymin, ymax),
                f"{col1}_{var1}", f"{col2}_{var2}")

    output_hists = []
    for _, hists_map in hists.items():
        for _, hist in hists_map.items():
            output_hists.append(hist)

    return output_hists


def analyze_hodo_peak():
    """Analyze hodoscope peak timing."""
    rdf = rdf_org
    hodo_pos_channels = build_hodo_pos_channels(run=run_number)

    histos = {key: {} for key in [
        "diff", "diff_relative", "sum", "left", "right",
        "left_vs_right", "left_peak", "right_peak", "LR_vs_UD"
    ]}

    # for group, channels in hodo_pos_channels.items():
    #    for channel in channels:
    #        rdf = rdf.Define(f"{channel}_peak_position",
    #                         f"ArgMaxRange({channel}_blsub, 0, -1)")
    #        rdf = rdf.Define(f"{channel}_peak_value",
    #                         f"MaxRange({channel}_blsub, 0, -1)")

    conditions = ["is_HoleVeto_vetoed", "passNone"]
    rdf_filtered = rdf

    for group, channels in hodo_pos_channels.items():
        for channel in channels:
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_sum",
                                               f"ROOT::VecOps::Sum({channel}_blsub) + 1e-6")
            rdf_filtered = rdf_filtered.Define(f"{channel}_subtracted_norm",
                                               f"{channel}_blsub / {channel}_subtracted_sum")

        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak",
                                           f"(int){channels[1]}_peak_position - (int){channels[0]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_sum_peak",
                                           f"(int){channels[0]}_peak_position + (int){channels[1]}_peak_position")
        rdf_filtered = rdf_filtered.Define(f"{group}_delta_peak_relative",
                                           f"{group}_delta_peak / ({group}_sum_peak + 1e-6)")

        for sel in conditions:
            cat = f"{group}_{sel}"
            histos["diff"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak", f"Delta Peak {cat}", 256, -1024, 1024),
                f"{group}_delta_peak", sel)
            histos["diff_relative"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak_relative",
                 f"Delta Peak Relative {cat}", 256, -1, 1),
                f"{group}_delta_peak_relative", sel)
            histos["sum"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_sum_peak", f"Sum Peak {cat}", 256, 0, 2048),
                f"{group}_sum_peak", sel)
            histos["left"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak", f"Left Peak {cat}", 256, 0, 1024),
                f"{channels[0]}_peak_position", sel)
            histos["right"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak", f"Right Peak {cat}", 256, 0, 1024),
                f"{channels[1]}_peak_position", sel)
            histos["left_vs_right"][cat] = rdf_filtered.Histo2D(
                (f"{cat}_left_peak_vs_right_peak", f"Left vs Right Peak {cat}",
                 256, 0, 1024, 256, 0, 1024),
                f"{channels[0]}_peak_position", f"{channels[1]}_peak_position", sel)
            histos["left_peak"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak_value",
                 f"Left Peak Value {cat}", 200, -2500, 1999),
                f"{channels[0]}_peak_value", sel)
            histos["right_peak"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak_value",
                 f"Right Peak Value {cat}", 200, -2500, 1999),
                f"{channels[1]}_peak_value", sel)

    # LR vs UD
    for sel in conditions:
        for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
            if dwc[:3] not in hodo_pos_channels:
                continue
            cat = f"{dwc}_{sel}"
            histos["LR_vs_UD"][cat] = rdf_filtered.Histo2D(
                (f"{cat}", "LR vs UD", 256, -1, 1, 256, -1, 1),
                f"{dwc.split('_vs_')[0]}_delta_peak_relative",
                f"{dwc.split('_vs_')[1]}_delta_peak_relative", sel)

    hists_output = []
    for hist_map in histos.values():
        for hist in hist_map.values():
            hists_output.append(hist)

    return hists_output


def plot_pulse(channels, suffix="services"):
    """Plot PID pulse analysis results."""
    infile_name = f"{paths['root']}/drs_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile_name}")

    output_htmls = []

    with PlotManager(paths['root'], paths['plots'], paths['html'], run_number) as pm:
        pm.set_output_dir(f"drs_{suffix}")

        for det in channels.keys():
            wf_ymin, wf_ymax = get_service_drs_processed_info_ranges(
                det, "waveform")
            pv_xmin, pv_xmax = get_service_drs_processed_info_ranges(
                det, "peak_value")
            en_xmin, en_xmax = get_service_drs_processed_info_ranges(
                det, "sum")
            _, _, value_cut, cut_method = get_service_drs_cut(det)

            # ADC vs TS (2D + profile)
            hist = infile.Get(f"{det}_ADC_vs_TS")
            if hist:
                pm.plot_2d(hist, f"{det}_ADC_vs_TS", "Time Slice", (0, 1024),
                           "ADC Counts", (wf_ymin, wf_ymax), style=STYLE_2D_LOG)
            prof = infile.Get(f"{det}_ADC_vs_TS_prof")
            if prof:
                prof_px = prof.ProjectionX()
                pm.plot_1d(prof_px, f"{det}_ADC_vs_TS_prof", "Time Slice", (0, 1024),
                           "Mean ADC Counts", (wf_ymin/10.0, wf_ymax/5.0), style=STYLE_1D)

            # CFD TS
            hist = infile.Get(f"{det}_cfd_ts")
            if hist:
                pm.plot_1d(hist, f"{det}_cfd_ts", "CFD TS", (0, 1024),
                           style=STYLE_1D)

            # Peak value vs CFD TS
            hist = infile.Get(f"{det}_peak_value_vs_cfd_ts")
            if hist:
                pm.plot_2d(hist, f"{det}_peak_value_vs_cfd_ts",
                           "CFD TS", (0, 1024), "Peak Value", (pv_xmin, pv_xmax),
                           style=STYLE_2D_LOG)

            # Peak value (cut annotation shown here when method is PeakValue)
            hist = infile.Get(f"{det}_peak_value")
            if hist:
                if cut_method == "PeakValue":
                    nhad = hist.Integral(hist.FindBin(value_cut), 10000)
                    nele = hist.Integral(0, hist.FindBin(value_cut))
                    ntot = nhad + nele + 1e-6
                    pave = create_pave_text(0.23, 0.75, 0.55, 0.85)
                    pave.SetFillColor(0)
                    pave.AddText(
                        f"N (peak > {value_cut:.2g}): {nhad:.0f} ({nhad/ntot:.1%})")
                    pave.AddText(
                        f"N (peak < {value_cut:.2g}): {nele:.0f} ({nele/ntot:.1%})")
                    line = ROOT.TLine(
                        value_cut, 0, value_cut, hist.GetMaximum())
                    line.SetLineColor(ROOT.kRed)
                    line.SetLineWidth(2)
                    line.SetLineStyle(ROOT.kDashed)
                    pm.plot_1d(hist, f"{det}_peak_value", "Peak Value", (pv_xmin, pv_xmax),
                               yrange=(1, None), style=STYLE_1D_LOG,
                               extraToDraw=[pave, line])
                else:
                    pm.plot_1d(hist, f"{det}_peak_value", "Peak Value", (pv_xmin, pv_xmax),
                               style=STYLE_1D, leftlegend=True)

            # Energy (cut annotation shown here when method is Sum/Energy)
            hist = infile.Get(f"{det}_energy")
            if hist:
                if cut_method != "PeakValue":
                    nhad = hist.Integral(hist.FindBin(value_cut), 10000)
                    nele = hist.Integral(0, hist.FindBin(value_cut))
                    ntot = nhad + nele + 1e-6
                    pave = create_pave_text(0.23, 0.75, 0.55, 0.85)
                    pave.SetFillColor(0)
                    pave.AddText(
                        f"N (energy > {value_cut:.2g}): {nhad:.0f} ({nhad/ntot:.1%})")
                    pave.AddText(
                        f"N (energy < {value_cut:.2g}): {nele:.0f} ({nele/ntot:.1%})")
                    line = ROOT.TLine(
                        value_cut, 0, value_cut, hist.GetMaximum())
                    line.SetLineColor(ROOT.kRed)
                    line.SetLineWidth(2)
                    line.SetLineStyle(ROOT.kDashed)
                    pm.plot_1d(hist, f"{det}_energy", "Energy (ADC)", (en_xmin, en_xmax),
                               yrange=(1, None), style=STYLE_1D_LOG,
                               extraToDraw=[pave, line])
                else:
                    pm.plot_1d(hist, f"{det}_energy", "Energy (ADC)", (en_xmin, en_xmax),
                               yrange=(1, None), style=STYLE_1D_LOG)

                # CDF plot
                hist.GetXaxis().SetRange(0, hist.GetNbinsX() + 1)
                hist_cdf = hist.GetCumulative()
                hist_cdf.Scale(1.0 / hist.Integral())
                line_cdf = ROOT.TLine(value_cut, 0, value_cut, 1)
                line_cdf.SetLineColor(ROOT.kRed)
                line_cdf.SetLineWidth(2)
                line_cdf.SetLineStyle(ROOT.kDashed)
                pm.plot_1d(hist_cdf, f"{det}_energy_cdf", "Energy (ADC)", (en_xmin, en_xmax),
                           ylabel="Cumulative Fraction", yrange=(0, 1.3),
                           style=STYLE_1D, addOverflow=False, addUnderflow=False,
                           legendPos=[0.3, 0.80, 0.5, 0.85], extraToDraw=line_cdf)

        intro_text = """This page shows the pulse shape analysis for the service DRS channels for particle identification.
No selection is applied unless specified."""

        output_htmls.append(pm.generate_html(
            f"ServiceDRS/{suffix}.html", plots_per_row=7, intro_text=intro_text, title=f"{suffix.upper()} Analysis"))

        # 2D correlation plots
        # pm.reset_plots()
        det_list = list(channels.keys())

        trigger_dets = ["KT1", "KT2", "T3", "T4"]

        # remove trigger_dets from det_list for correlation plotting
        det_list = [det for det in det_list if det not in trigger_dets]

        corr_categories = [("PID", det_list)] if suffix == "mcp" else [
            ("PID", det_list), ("Trigger", trigger_dets)]
        for cat, tmp_list in corr_categories:
            pm.reset_plots()
            for idx1, det1 in enumerate(tmp_list):
                for idx2, det2 in enumerate(tmp_list):
                    if idx2 <= idx1:
                        continue

                    _, _, value_cut1, method1 = get_service_drs_cut(det1)
                    _, _, value_cut2, method2 = get_service_drs_cut(det2)
                    var1 = "peak_value" if method1 == "PeakValue" else "energy"
                    var2 = "peak_value" if method2 == "PeakValue" else "energy"
                    xmin, xmax = get_service_drs_processed_info_ranges(
                        det1, "peak_value" if method1 == "PeakValue" else "sum")
                    ymin, ymax = get_service_drs_processed_info_ranges(
                        det2, "peak_value" if method2 == "PeakValue" else "sum")

                    hist2d = infile.Get(f"{det1}_{var1}_vs_{det2}_{var2}")
                    if not hist2d:
                        continue

                    xPass = hist2d.GetXaxis().FindBin(value_cut1)
                    yPass = hist2d.GetYaxis().FindBin(value_cut2)
                    nPP = hist2d.Integral(xPass, 1000, yPass, 1000)
                    nPF = hist2d.Integral(xPass, 1000, 0, yPass - 1)
                    nFP = hist2d.Integral(0, xPass - 1, yPass, 1000)
                    nFF = hist2d.Integral(0, xPass - 1, 0, yPass - 1)

                    pave = create_pave_text(0.23, 0.20, 0.5, 0.40)
                    pave.SetTextSize(0.03)
                    name1 = det1.replace("Cerenkov", "Cer")
                    name2 = det2.replace("Cerenkov", "Cer")
                    pave.AddText(
                        f"N ({name1} > {value_cut1:.2g}, {name2} > {value_cut2:.2g}): {nPP:.0f}")
                    pave.AddText(
                        f"N ({name1} > {value_cut1:.2g}, {name2} < {value_cut2:.2g}): {nPF:.0f}")
                    pave.AddText(
                        f"N ({name1} < {value_cut1:.2g}, {name2} > {value_cut2:.2g}): {nFP:.0f}")
                    pave.AddText(
                        f"N ({name1} < {value_cut1:.2g}, {name2} < {value_cut2:.2g}): {nFF:.0f}")

                    line1 = ROOT.TLine(value_cut1, ymin, value_cut1, ymax)
                    line2 = ROOT.TLine(xmin, value_cut2, xmax, value_cut2)
                    for line in [line1, line2]:
                        line.SetLineWidth(2)
                        line.SetLineStyle(ROOT.kDashed)
                        line.SetLineColor(ROOT.kRed)

                    pm.plot_2d(hist2d, f"{det1}_vs_{det2}_corr2D",
                               f"{det1} {var1}", (xmin, xmax),
                               f"{det2} {var2}", (ymin, ymax),
                               style=STYLE_2D_LOG,
                               extraToDraw=[pave, line1, line2])

                pm.add_newline()

            intro_text = f"""This page shows the correlation plots of the service DRS channels for {cat}.
No selection is applied unless specified."""

            output_htmls.append(pm.generate_html(
                f"ServiceDRS/{cat}_correlation_{suffix}.html", plots_per_row=5, intro_text=intro_text))

    infile.Close()
    return output_htmls


def analyze_mcp_timing_diff(channels_mcp):
    """Book event-by-event MCP CFD time difference histograms w.r.t. the first MCP."""
    rdf = rdf_org
    hists = []

    dets = list(channels_mcp.keys())
    ref_det = dets[0]
    ref_col = channels_mcp[ref_det]

    for det, channel in channels_mcp.items():
        if det == ref_det:
            continue
        diff_col = f"{channel}_minus_{ref_col}_cfd_diff"
        rdf = rdf.Define(
            diff_col, f"{channel}_TS_cfd_ref - {ref_col}_TS_cfd_ref")
        hists.append(rdf.Histo1D(
            (f"{det}_cfd_diff_vs_{ref_det}",
             f"{det} - {ref_det};#Delta t_{{CFD,ref}} [TS];Counts",
             400, -10, 10),
            diff_col))

    return hists, ref_det


def plot_mcp_timing_diff(channels_mcp, ref_det):
    """Plot MCP CFD timing difference histograms in a dedicated HTML page."""
    infile_name = f"{paths['root']}/drs_mcp_timing_diff.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile_name}")

    dets = [det for det in channels_mcp if det != ref_det]

    with PlotManager(paths['root'], paths['plots'], paths['html'], run_number) as pm:
        pm.set_output_dir("drs_mcp_timing_diff")

        for det in dets:
            hist = infile.Get(f"{det}_cfd_diff_vs_{ref_det}")
            if not hist:
                continue
            hist.Fit("gaus", "Q")
            fit = hist.GetFunction("gaus")
            fit.SetLineColor(ROOT.kRed)
            fit.SetLineWidth(2)
            pave = create_pave_text(0.60, 0.75, 0.90, 0.88)
            pave.AddText(f"Mean: {fit.GetParameter(1):.2f} TS")
            pave.AddText(f"Sigma: {fit.GetParameter(2):.2f} TS")
            pm.plot_1d(hist, f"{det}_cfd_diff_vs_{ref_det}",
                       f"#Delta t_{{CFD,ref}} ({det} - {ref_det}) [TS]", (-10, 10),
                       yrange=(0.1, 200),
                       ylabel="Counts", style=STYLE_1D,
                       extraToDraw=[fit, pave])

        intro_text = (f"Event-by-event MCP CFD time differences with respect to "
                      f"{ref_det}. No selection applied.\n This should be the timing resolution of the system, (DRS+Method) for the same MCP; DRS+Method+MCP for different MCPs.")
        output_html = pm.generate_html(
            "ServiceDRS/mcp_timing_diff.html",
            plots_per_row=3,
            title="MCP CFD Timing Differences",
            intro_text=intro_text)

    infile.Close()
    return output_html


def plot_hodo_peak():
    """Plot hodoscope peak analysis results."""
    infile = ROOT.TFile(f"{paths['root']}/hodoscope_peaks.root", "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file")

    hodo_pos_channels = build_hodo_pos_channels(run=run_number)

    with PlotManager(paths['root'], paths['plots'], paths['html'], run_number) as pm:
        pm.set_output_dir("DWC")

        for group, channels in hodo_pos_channels.items():
            labels = ["No Selection", "Pass HoleVeto"]
            linestyles = [1, 2]

            # Collect histograms
            histos = {key: [] for key in ["diff", "diff_rel", "sum", "left", "right",
                                          "left_vs_right", "left_peak", "right_peak"]}

            for sel in ["passNone", "is_HoleVeto_vetoed"]:
                cat = f"{group}_{sel}"
                histos["diff"].append(infile.Get(f"{cat}_delta_peak"))
                histos["diff_rel"].append(
                    infile.Get(f"{cat}_delta_peak_relative"))
                histos["sum"].append(infile.Get(f"{cat}_sum_peak"))
                histos["left"].append(infile.Get(f"{cat}_left_peak"))
                histos["right"].append(infile.Get(f"{cat}_right_peak"))
                histos["left_vs_right"].append(
                    infile.Get(f"{cat}_left_peak_vs_right_peak"))
                histos["left_peak"].append(
                    infile.Get(f"{cat}_left_peak_value"))
                histos["right_peak"].append(
                    infile.Get(f"{cat}_right_peak_value"))

            style_compare = PlotStyle(dology=False, drawoptions="HIST", mycolors=[1, 1],
                                      linestyles=linestyles)

            pm.plot_1d(histos["diff"], f"{group}_diff_peak",
                       "Peak Position Difference", (-250, 250),
                       legends=labels, style=style_compare)

            pm.plot_1d(histos["diff_rel"], f"{group}_diff_peak_relative",
                       "Peak Position Difference Relative", (-0.5, 0.5),
                       legends=labels, style=style_compare)

            pm.plot_1d(histos["sum"], f"{group}_sum_peak",
                       "Peak Position Sum", (600, 1200),
                       legends=labels, style=style_compare)

            pm.plot_1d(
                histos["left"] + histos["right"],
                f"{group}_peaks",
                "Peak Position", (300, 800),
                legends=["Left No Sel", "Left Pass Veto",
                         "Right No Sel", "Right Pass Veto"],
                style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[1, 1, 2, 2],
                                linestyles=[1, 2, 1, 2]))

            for idx, cat in enumerate(["pass_NoSel", "pass_upstream_veto"]):
                if histos["left_vs_right"][idx]:
                    pm.plot_2d(histos["left_vs_right"][idx],
                               f"{group}_left_vs_right_peak_{cat}",
                               "Left Peak Position", (0, 1024),
                               "Right Peak Position", (0, 1024),
                               style=STYLE_2D_LOG)

            pm.plot_1d(
                histos["left_peak"] + histos["right_peak"],
                f"{group}_peak_values",
                "Peak Value", (-1500, 100),
                legends=["Left No Sel", "Left Pass Veto",
                         "Right No Sel", "Right Pass Veto"],
                style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[1, 1, 2, 2],
                                linestyles=[1, 2, 1, 2]))

        # LR vs UD summary
        plots_summary = []
        for sel in ["pass_NoSel", "pass_upstream_veto"]:
            for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
                cat = f"{dwc}_{sel}"
                hist2d = infile.Get(cat)
                if hist2d:
                    pm.plot_2d(hist2d, f"DWC_{dwc}_{sel}",
                               "LR Delta Peak Relative", (-0.4, 0.4),
                               "UD Delta Peak Relative", (-0.4, 0.4),
                               style=STYLE_2D_LOG, prepend=True)

        pm.add_newline()

        output_html = pm.generate_html(
            "ServiceDRS/DWC.html", plots_per_row=7, title="DWC Positions")

    infile.Close()
    return output_html


def main():
    channels = OrderedDict()
    for det in ["HoleVeto", "PSD", "TTUMuonVeto", "Cer474", "Cer519", "Cer537", "KT1", "KT2", "T3", "T4"]:
        channel = get_service_drs_channels(run_number).get(det)
        channels[det] = channel

    channels_mcp = get_mcp_channels(run_number)

    hists_pid = analyze_pulse(channels)
    # hists_hodo = analyze_hodo_peak()
    hists_mcp = analyze_pulse(channels_mcp)
    hists_mcp_diff, ref_det = analyze_mcp_timing_diff(channels_mcp)

    print("Triggering all computations...")
    # type: ignore[attr-defined]
    ROOT.RDF.RunGraphs(hists_pid + hists_mcp + hists_mcp_diff)

    # Save histograms
    save_hists_to_file(hists_pid, f"{paths['root']}/drs_services.root")
    # save_hists_to_file(hists_hodo, f"{paths['root']}/hodoscope_peaks.root")
    save_hists_to_file(hists_mcp, f"{paths['root']}/drs_mcp.root")
    save_hists_to_file(
        hists_mcp_diff, f"{paths['root']}/drs_mcp_timing_diff.root")

    # Make plots
    outputs = {}
    outputs["PID"] = plot_pulse(channels, suffix="services")
    # outputs["DWC"] = plot_hodo_peak()
    outputs["MCP"] = plot_pulse(channels_mcp, suffix="mcp")
    outputs["MCP_timing_diff"] = plot_mcp_timing_diff(channels_mcp, ref_det)

    for key, value in outputs.items():
        if isinstance(value, list):
            for item in value:
                print(f"Output for {key}: {item}")
        else:
            print(f"Output for {key}: {value}")


if __name__ == "__main__":
    main()
