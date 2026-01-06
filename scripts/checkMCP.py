import re
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from channels.channel_map import getMCPChannels
from utils.html_generator import generate_html
from utils.dataloader import loadRDF
from variables.drs import preProcessDRSBoards
from utils.parser import get_args
from utils.auto_compile import auto_compile
from utils.timing import auto_timer  # noqa
auto_timer("Total Execution Time")

ROOT.gROOT.SetBatch(True)  # Run in batch mode
ROOT.ROOT.EnableImplicitMT(10)


def main():
    parser = get_args()
    runNumber = parser.run
    firstEvent = parser.first_event
    lastEvent = parser.last_event
    jsonFile = parser.json_file

    TSmin = 500
    TSmax = 700
    RelTSmin = -320
    RelTSmax = -130
    DiffRelTSmin = -5
    DiffRelTSmax = 5
    DiffRelTSmin_US = -8
    DiffRelTSmax_US = 10

    # jitter offset
    value_diffcorrs = [0, 0, 0, -16]

    outputs_html = {}

    rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent, jsonFile)
    rdf = preProcessDRSBoards(rdf, runNumber=runNumber)

    # get the MCP peak TS
    map_mcp_channels = getMCPChannels(runNumber)

    # get time reference channel TS
    channels_TF = [
        "DRS_Board0_Group3_Channel8",
        "DRS_Board1_Group3_Channel8",
        "DRS_Board2_Group3_Channel8",
        "DRS_Board3_Group3_Channel8",
    ]

    for channel in channels_TF:
        rdf = rdf.Define(f"{channel}_PeakTS",
                         f"ArgMinRange({channel}_blsub, 0, 1024, -70.0)")

    for det, channels in map_mcp_channels.items():
        for idx, channel in enumerate(channels):
            rdf = rdf.Define(f"{channel}_PeakTS",
                             f"ArgMinRange({channel}_blsub, {TSmin}, {TSmax}, -200.0)")
            # define the relative peak TS with respect to the reference channel
            channel_TS = re.sub(r"_Channel[0-7]", "_Channel8", channel)
            rdf = rdf.Define(
                f"{channel}_RelPeakTS", f"(int){channel}_PeakTS - (int){channel_TS}_PeakTS")

            # define the difference of relative peak TS with respect to the first channel
            rdf = rdf.Define(
                f"{channel}_DiffRelPeakTS", f"(int){channel}_RelPeakTS - (int){channels[0]}_RelPeakTS - {value_diffcorrs[idx]}")

    for det, channels in map_mcp_channels.items():
        for idx, channel in enumerate(channels):
            # define the difference of relative peak TS with respect to the first channel of US
            rdf = rdf.Define(
                f"{channel}_DiffRelPeakTS_US", f"(int){channel}_RelPeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[idx]}")

    rdf_prefilter = rdf
    condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
    condition += f" && {map_mcp_channels['DS'][0]}_RelPeakTS > -350 && {map_mcp_channels['DS'][0]}_RelPeakTS < -100"
    condition += f" && {map_mcp_channels['US'][0]}_PeakTS > {TSmin} && {map_mcp_channels['US'][0]}_PeakTS < {TSmax}"
    rdf = rdf.Filter(condition,
                     "Pre-filter on MCP US channel 0 Peak TS")

    # plot the time reference channel TS
    hists_tf = []
    for channel in channels_TF:
        h1 = rdf.Histo1D(
            (f"h_{channel}_PeakTS", f"{channel} Peak TS;Peak TS;Counts",
             500, 500, 1000),
            f"{channel}_PeakTS"
        )
        hists_tf.append(h1)

    # plot the MCP peak TS
    hists_1d = {}
    hists_2d = {}
    for det, channels in map_mcp_channels.items():
        hists_1d[det] = []
        hists_2d[det] = []
        for idx, channel in enumerate(channels):
            h1 = rdf.Histo1D(
                (f"h_{channel}_PeakTS", f"{channel} Peak TS;Peak TS;Counts",
                 TSmax - TSmin, TSmin, TSmax),
                f"{channel}_PeakTS"
            )
            hists_1d[det].append(h1)

            if idx > 0:
                h2 = rdf.Histo2D((f"h2_{channel}_PeakTS_vs_{channels[0]}_PeakTS", f"{channel} Peak TS vs {channels[0]} Peak TS;{channels[0]} Peak TS;{channel} Peak TS",
                                 TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
                                 f"{channels[0]}_PeakTS", f"{channel}_PeakTS"
                                 )
                hists_2d[det].append(h2)

    hists_2d_mcps = []
    for idx in range(min(len(map_mcp_channels["US"]), len(map_mcp_channels["DS"]))):
        channel_us = map_mcp_channels["US"][idx]
        channel_ds = map_mcp_channels["DS"][idx]
        h2 = rdf.Histo2D((f"h2_{channel_ds}_PeakTS_vs_{channel_us}_PeakTS", f"{channel_ds} Peak TS vs {channel_us} Peak TS;{channel_us} Peak TS;{channel_ds} Peak TS",
                         TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
                         f"{channel_us}_PeakTS", f"{channel_ds}_PeakTS")
        hists_2d_mcps.append(h2)

    # plot the relative MCP peak TS
    hists_rel_1d = {}
    hists_rel_2d = {}
    for det, channels in map_mcp_channels.items():
        hists_rel_1d[det] = []
        hists_rel_2d[det] = []
        for idx, channel in enumerate(channels):
            h1 = rdf.Histo1D(
                (f"h_{channel}_RelPeakTS", f"{channel} Relative Peak TS;Relative Peak TS;Counts",
                 RelTSmax - RelTSmin, RelTSmin, RelTSmax),
                f"{channel}_RelPeakTS"
            )
            hists_rel_1d[det].append(h1)

            if idx > 0:
                h2 = rdf.Histo2D((f"h2_{channel}_RelPeakTS_vs_{channels[0]}_RelPeakTS", f"{channel} Relative Peak TS vs {channels[0]} Relative Peak TS;{channels[0]} Relative Peak TS;{channel} Relative Peak TS",
                                 RelTSmax - RelTSmin, RelTSmin, RelTSmax, RelTSmax - RelTSmin, RelTSmin, RelTSmax),
                                 f"{channels[0]}_RelPeakTS", f"{channel}_RelPeakTS"
                                 )
                hists_rel_2d[det].append(h2)

    hists_2d_rel_mcps = []
    for idx in range(min(len(map_mcp_channels["US"]), len(map_mcp_channels["DS"]))):
        channel_us = map_mcp_channels["US"][idx]
        channel_ds = map_mcp_channels["DS"][idx]
        h2 = rdf.Histo2D((f"h2_{channel_ds}_RelPeakTS_vs_{channel_us}_RelPeakTS", f"{channel_ds} Relative Peak TS vs {channel_us} Relative Peak TS;{channel_us} Relative Peak TS;{channel_ds} Relative Peak TS",
                         RelTSmax - RelTSmin, RelTSmin, RelTSmax, RelTSmax - RelTSmin, RelTSmin, RelTSmax),
                         f"{channel_us}_RelPeakTS", f"{channel_ds}_RelPeakTS")
        hists_2d_rel_mcps.append(h2)

    # plot the diffrence of relative MCP peak TS
    hists_diff_rel_1d = {}
    for det, channels in map_mcp_channels.items():
        hists_diff_rel_1d[det] = []
        for idx, channel in enumerate(channels):
            if idx == 0:
                continue
            h1 = rdf.Histo1D(
                (f"h_{channel}_DiffRelPeakTS", f"{channel} Difference of Relative Peak TS;Difference of Relative Peak TS;Counts",
                 DiffRelTSmax - DiffRelTSmin, DiffRelTSmin, DiffRelTSmax),
                f"{channel}_DiffRelPeakTS"
            )
            hists_diff_rel_1d[det].append(h1)

    hists_diff_rel_1d_us = {}
    for det, channels in map_mcp_channels.items():
        hists_diff_rel_1d_us[det] = []
        for idx, channel in enumerate(channels):
            h1 = rdf.Histo1D(
                (f"h_{channel}_DiffRelPeakTS_US", f"{channel} Difference of Relative Peak TS wrt US;Difference of Relative Peak TS wrt US;Counts",
                 DiffRelTSmax_US - DiffRelTSmin_US, DiffRelTSmin_US, DiffRelTSmax_US),
                f"{channel}_DiffRelPeakTS_US"
            )
            hists_diff_rel_1d_us[det].append(h1)

    #
    # make plots
    #
    labels = ["Board0", "Board1", "Board2", "Board3"]
    legendPos = [0.20, 0.85, 0.90, 0.9]
    output_dir = f"results/plots/Run{runNumber}/MCP"
    plots = []
    output_name = f"Time_Reference_Peak_TS"
    DrawHistos(hists_tf, labels, 750, 900, "Peak TS", 0, None, "Counts",
               outputname=output_name, outdir=output_dir,
               dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
               addOverflow=True, addUnderflow=True, legendNCols=4, legendPos=legendPos)
    plots.append(f"{output_name}.png")
    # Special marker to indicate a new line in the HTML
    plots.append("NEWLINE")

    for det in map_mcp_channels.keys():
        hists = hists_1d[det]
        outputname = f"{det}_peak_ts"
        DrawHistos(hists, labels, TSmin, TSmax, "Peak TS", 0, None, "Counts",
                   outputname=outputname, outdir=output_dir,
                   dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
                   addOverflow=False, addUnderflow=False, legendNCols=4, legendPos=legendPos)
        plots.append(f"{outputname}.png")

        for idx, hist2d in enumerate(hists_2d[det]):
            outputname = f"{det}_peak_ts_{channels[idx+1]}_vs_{channels[0]}"
            DrawHistos([hist2d], [], TSmin, TSmax, f"Peak TS Board0", TSmin, TSmax, f"Peak TS Board{idx+1}",
                       outputname=outputname,
                       outdir=output_dir,
                       drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                       dology=False, runNumber=runNumber,
                       addOverflow=True, doth2=True)
            plots.append(f"{outputname}.png")

    for idx, hist2d in enumerate(hists_2d_mcps):
        channel_us = map_mcp_channels["US"][idx]
        channel_ds = map_mcp_channels["DS"][idx]
        outputname = f"MCP_peak_ts_{channel_ds}_vs_{channel_us}"
        DrawHistos([hist2d], [], TSmin, TSmax, f"Peak TS US", TSmin, TSmax, f"Peak TS DS",
                   outputname=outputname,
                   outdir=output_dir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, runNumber=runNumber,
                   addOverflow=True, doth2=True)
        plots.append(f"{outputname}.png")

    output_html = f"results/html/Run{runNumber}/ServiceDRS/MCPPeakTS.html"
    outputs_html["abs"] = generate_html(
        plots, output_dir, plots_per_row=4, output_html=output_html)

    plots = []

    for det in map_mcp_channels.keys():
        hists = hists_rel_1d[det]
        outputname = f"{det}_relative_peak_ts"
        extraToDraw = ROOT.TPaveText(0.20, 0.63, 0.50, 0.83, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColor(0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(
            f"B0: {hists[0].GetMean():.1f} #pm {hists[0].GetRMS():.1f}")
        extraToDraw.AddText(
            f"B1: {hists[1].GetMean():.1f} #pm {hists[1].GetRMS():.1f}")
        extraToDraw.AddText(
            f"B2: {hists[2].GetMean():.1f} #pm {hists[2].GetRMS():.1f}")
        extraToDraw.AddText(
            f"B3: {hists[3].GetMean():.1f} #pm {hists[3].GetRMS():.1f}")
        DrawHistos(hists, labels, RelTSmin, RelTSmax, "Relative Peak TS", 0, hists[0].GetMaximum()*1.5, "Counts",
                   outputname=outputname, outdir=output_dir,
                   dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
                   addOverflow=False, addUnderflow=False, legendNCols=4, legendPos=legendPos, extraToDraw=extraToDraw)
        plots.append(f"{outputname}.png")

        for idx, hist2d in enumerate(hists_rel_2d[det]):
            outputname = f"{det}_relative_peak_ts_{channels[idx+1]}_vs_{channels[0]}"
            DrawHistos([hist2d], [], RelTSmin, RelTSmax, f"Relative Peak TS Board0", RelTSmin, RelTSmax, f"Relative Peak TS Board{idx+1}",
                       outputname=outputname,
                       outdir=output_dir,
                       drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                       dology=False, runNumber=runNumber,
                       addOverflow=True, doth2=True)
            plots.append(f"{outputname}.png")

    for idx, hist2d in enumerate(hists_2d_rel_mcps):
        channel_us = map_mcp_channels["US"][idx]
        channel_ds = map_mcp_channels["DS"][idx]
        outputname = f"MCP_relative_peak_ts_{channel_ds}_vs_{channel_us}"
        DrawHistos([hist2d], [], RelTSmin, RelTSmax, f"Relative Peak TS US", RelTSmin, RelTSmax, f"Relative Peak TS DS",
                   outputname=outputname,
                   outdir=output_dir,
                   drawoptions="COLz", zmin=1, zmax=None, dologz=True,
                   dology=False, runNumber=runNumber,
                   addOverflow=True, doth2=True)
        plots.append(f"{outputname}.png")

    for det in map_mcp_channels.keys():
        hists = hists_diff_rel_1d[det]
        outputname = f"{det}_difference_relative_peak_ts_wrt_board0"
        extraToDraw = ROOT.TPaveText(0.20, 0.63, 0.50, 0.83, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColor(0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        for idx, hist in enumerate(hists):
            extraToDraw.AddText(
                f"B{idx+1}: {hist.GetMean():.1f} #pm {hist.GetRMS():.1f}")
        DrawHistos(hists, labels[1:], DiffRelTSmin, DiffRelTSmax, "Difference of Relative Peak TS", 0, hists[0].GetMaximum()*1.5, "Counts",
                   outputname=outputname, outdir=output_dir,
                   dology=False, mycolors=[2, 3, 4], drawashist=True, runNumber=runNumber,
                   addOverflow=False, addUnderflow=False, legendNCols=3, legendPos=legendPos, extraToDraw=extraToDraw)
        plots.append(f"{outputname}.png")

    # Special marker to indicate a new line in the HTML
    plots.append("NEWLINE")

    for det in map_mcp_channels.keys():
        hists = hists_diff_rel_1d_us[det]
        outputname = f"{det}_difference_relative_peak_ts_wrt_US_board0"
        extraToDraw = ROOT.TPaveText(0.20, 0.63, 0.50, 0.83, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColor(0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        for idx, hist in enumerate(hists):
            extraToDraw.AddText(
                f"B{idx}: {hist.GetMean():.1f} #pm {hist.GetRMS():.1f}")
        DrawHistos(hists, labels, DiffRelTSmin_US, DiffRelTSmax_US, "Difference of Relative Peak TS wrt US", 0, hists[1].GetMaximum()*1.5, "Counts",
                   outputname=outputname, outdir=output_dir,
                   dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
                   addOverflow=False, addUnderflow=False, legendNCols=4, legendPos=legendPos, extraToDraw=extraToDraw)
        plots.append(f"{outputname}.png")

    output_html = f"results/html/Run{runNumber}/ServiceDRS/MCPPeakTS_relative.html"
    outputs_html["rel"] = generate_html(
        plots, output_dir, plots_per_row=4, output_html=output_html)

    for key, value in outputs_html.items():
        print(f"MCP {key} plots html: {value}")


if __name__ == "__main__":
    main()
