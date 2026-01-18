"""
MCP (Micro-Channel Plate) Check Script.

Analyzes MCP peak timing for time calibration
"""

import re
import ROOT
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from utils.parser import get_args
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from plotting.calox_plot_helper import create_pave_text
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

# Common styles
STYLE_MULTI_BOARD = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[1, 2, 3, 4],
    legendNCols=4,
    legendPos=[0.20, 0.85, 0.90, 0.9]
)

STYLE_2D_LOG = PlotStyle(
    dology=False,
    dologz=True,
    drawoptions="COLZ",
    zmin=1,
    zmax=None
)


def main():
    args = get_args()
    run_number = args.run

    # Configuration
    TSmin, TSmax = 500, 700
    RelTSmin, RelTSmax = -320, -130
    DiffRelTSmin, DiffRelTSmax = -5, 5
    DiffRelTSmin_US, DiffRelTSmax_US = -8, 10

    # Jitter offset
    value_diffcorrs = [0, 0, 0, -16]

    outputs_html = {}

    analysis = (CaloXAnalysisManager(args)
                .prepare()
                .apply_hole_veto(flag_only=True)
                )
    rdf = analysis.get_rdf()

    map_mcp_channels = get_mcp_channels(run_number)

    # Time reference channels
    channels_TF = [
        "DRS_Board0_Group3_Channel8",
        "DRS_Board1_Group3_Channel8",
        "DRS_Board2_Group3_Channel8",
        "DRS_Board3_Group3_Channel8",
    ]

    # Define time reference peak TS
    for channel in channels_TF:
        rdf = rdf.Define(f"{channel}_PeakTS",
                         f"ArgMinRange({channel}_blsub, 0, 1024, -70.0)")

    # Define MCP peak TS and related variables
    for det, channels in map_mcp_channels.items():
        for idx, channel in enumerate(channels):
            rdf = rdf.Define(f"{channel}_PeakTS",
                             f"ArgMinRange({channel}_blsub, {TSmin}, {TSmax}, -200.0)")
            channel_TS = re.sub(r"_Channel[0-7]", "_Channel8", channel)
            rdf = rdf.Define(f"{channel}_RelPeakTS",
                             f"(int){channel}_PeakTS - (int){channel_TS}_PeakTS")
            rdf = rdf.Define(f"{channel}_DiffRelPeakTS",
                             f"(int){channel}_RelPeakTS - (int){channels[0]}_RelPeakTS - {value_diffcorrs[idx]}")

    for det, channels in map_mcp_channels.items():
        for idx, channel in enumerate(channels):
            rdf = rdf.Define(
                f"{channel}_DiffRelPeakTS_US",
                f"(int){channel}_RelPeakTS - (int){map_mcp_channels['US'][0]}_RelPeakTS - {value_diffcorrs[idx]}")

    # Apply pre-filter
    condition = f"{map_mcp_channels['US'][0]}_RelPeakTS > -350 && {map_mcp_channels['US'][0]}_RelPeakTS < -100"
    condition += f" && {map_mcp_channels['DS'][0]}_RelPeakTS > -350 && {map_mcp_channels['DS'][0]}_RelPeakTS < -100"
    condition += f" && {map_mcp_channels['US'][0]}_PeakTS > {TSmin} && {map_mcp_channels['US'][0]}_PeakTS < {TSmax}"
    rdf = rdf.Filter(condition, "Pre-filter on MCP US channel 0 Peak TS")

    # Book histograms
    hists_tf = [rdf.Histo1D(
        (f"h_{ch}_PeakTS", f"{ch} Peak TS;Peak TS;Counts", 500, 500, 1000),
        f"{ch}_PeakTS") for ch in channels_TF]

    hists_1d, hists_2d = {}, {}
    for det, channels in map_mcp_channels.items():
        hists_1d[det] = [rdf.Histo1D(
            (f"h_{ch}_PeakTS", f"{ch} Peak TS;Peak TS;Counts",
             TSmax - TSmin, TSmin, TSmax),
            f"{ch}_PeakTS") for ch in channels]
        hists_2d[det] = [rdf.Histo2D(
            (f"h2_{channels[i]}_PeakTS_vs_{channels[0]}_PeakTS", "",
             TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
            f"{channels[0]}_PeakTS", f"{channels[i]}_PeakTS") for i in range(1, len(channels))]

    hists_2d_mcps = [rdf.Histo2D(
        (f"h2_{map_mcp_channels['DS'][i]}_PeakTS_vs_{map_mcp_channels['US'][i]}_PeakTS", "",
         TSmax - TSmin, TSmin, TSmax, TSmax - TSmin, TSmin, TSmax),
        f"{map_mcp_channels['US'][i]}_PeakTS", f"{map_mcp_channels['DS'][i]}_PeakTS")
        for i in range(min(len(map_mcp_channels["US"]), len(map_mcp_channels["DS"])))]

    hists_rel_1d, hists_rel_2d = {}, {}
    for det, channels in map_mcp_channels.items():
        hists_rel_1d[det] = [rdf.Histo1D(
            (f"h_{ch}_RelPeakTS", f"{ch} Relative Peak TS",
             RelTSmax - RelTSmin, RelTSmin, RelTSmax),
            f"{ch}_RelPeakTS") for ch in channels]
        hists_rel_2d[det] = [rdf.Histo2D(
            (f"h2_{channels[i]}_RelPeakTS_vs_{channels[0]}_RelPeakTS", "",
             RelTSmax - RelTSmin, RelTSmin, RelTSmax, RelTSmax - RelTSmin, RelTSmin, RelTSmax),
            f"{channels[0]}_RelPeakTS", f"{channels[i]}_RelPeakTS") for i in range(1, len(channels))]

    hists_2d_rel_mcps = [rdf.Histo2D(
        (f"h2_{map_mcp_channels['DS'][i]}_RelPeakTS_vs_{map_mcp_channels['US'][i]}_RelPeakTS", "",
         RelTSmax - RelTSmin, RelTSmin, RelTSmax, RelTSmax - RelTSmin, RelTSmin, RelTSmax),
        f"{map_mcp_channels['US'][i]}_RelPeakTS", f"{map_mcp_channels['DS'][i]}_RelPeakTS")
        for i in range(min(len(map_mcp_channels["US"]), len(map_mcp_channels["DS"])))]

    hists_diff_rel_1d = {det: [rdf.Histo1D(
        (f"h_{ch}_DiffRelPeakTS", f"{ch} Diff Rel Peak TS",
         DiffRelTSmax - DiffRelTSmin, DiffRelTSmin, DiffRelTSmax),
        f"{ch}_DiffRelPeakTS") for ch in channels[1:]]  # Skip first channel
        for det, channels in map_mcp_channels.items()}

    hists_diff_rel_1d_us = {det: [rdf.Histo1D(
        (f"h_{ch}_DiffRelPeakTS_US", f"{ch} Diff Rel Peak TS wrt US",
         DiffRelTSmax_US - DiffRelTSmin_US, DiffRelTSmin_US, DiffRelTSmax_US),
        f"{ch}_DiffRelPeakTS_US") for ch in channels]
        for det, channels in map_mcp_channels.items()}

    # Make plots
    labels = ["Board0", "Board1", "Board2", "Board3"]
    output_dir = f"results/plots/Run{run_number}/MCP"

    with PlotManager("", output_dir, f"results/html/Run{run_number}", run_number) as pm:
        pm.set_output_dir("")

        # Time reference plot
        pm.plot_1d(
            hists_tf,
            "Time_Reference_Peak_TS",
            "Peak TS",
            (750, 900),
            legends=labels,
            style=STYLE_MULTI_BOARD
        )
        pm.add_newline()

        # Absolute peak TS plots
        for det, channels in map_mcp_channels.items():
            pm.plot_1d(
                hists_1d[det],
                f"{det}_peak_ts",
                "Peak TS",
                (TSmin, TSmax),
                legends=labels,
                style=STYLE_MULTI_BOARD,
                addOverflow=False,
                addUnderflow=False
            )

            for idx, hist2d in enumerate(hists_2d[det]):
                pm.plot_2d(
                    hist2d,
                    f"{det}_peak_ts_{channels[idx+1]}_vs_{channels[0]}",
                    "Peak TS Board0",
                    (TSmin, TSmax),
                    f"Peak TS Board{idx+1}",
                    (TSmin, TSmax),
                    style=STYLE_2D_LOG
                )

        for idx, hist2d in enumerate(hists_2d_mcps):
            pm.plot_2d(
                hist2d,
                f"MCP_peak_ts_{map_mcp_channels['DS'][idx]}_vs_{map_mcp_channels['US'][idx]}",
                "Peak TS US",
                (TSmin, TSmax),
                "Peak TS DS",
                (TSmin, TSmax),
                style=STYLE_2D_LOG
            )

        outputs_html["abs"] = pm.generate_html(
            "ServiceDRS/MCPPeakTS.html", plots_per_row=4)

        # Relative peak TS plots
        pm.reset_plots()

        for det, channels in map_mcp_channels.items():
            hists = hists_rel_1d[det]

            pave = create_pave_text(0.20, 0.63, 0.50, 0.83)
            for i, h in enumerate(hists):
                pave.AddText(f"B{i}: {h.GetMean():.1f} #pm {h.GetRMS():.1f}")

            pm.plot_1d(
                hists,
                f"{det}_relative_peak_ts",
                "Relative Peak TS",
                (RelTSmin, RelTSmax),
                yrange=(0, hists[0].GetMaximum() * 1.5),
                legends=labels,
                style=STYLE_MULTI_BOARD,
                addOverflow=False,
                addUnderflow=False,
                extraToDraw=pave
            )

            for idx, hist2d in enumerate(hists_rel_2d[det]):
                pm.plot_2d(
                    hist2d,
                    f"{det}_relative_peak_ts_{channels[idx+1]}_vs_{channels[0]}",
                    "Relative Peak TS Board0",
                    (RelTSmin, RelTSmax),
                    f"Relative Peak TS Board{idx+1}",
                    (RelTSmin, RelTSmax),
                    style=STYLE_2D_LOG
                )

        for idx, hist2d in enumerate(hists_2d_rel_mcps):
            pm.plot_2d(
                hist2d,
                f"MCP_relative_peak_ts_{map_mcp_channels['DS'][idx]}_vs_{map_mcp_channels['US'][idx]}",
                "Relative Peak TS US",
                (RelTSmin, RelTSmax),
                "Relative Peak TS DS",
                (RelTSmin, RelTSmax),
                style=STYLE_2D_LOG
            )

        # Difference of relative peak TS
        for det, channels in map_mcp_channels.items():
            hists = hists_diff_rel_1d[det]
            if not hists:
                continue

            pave = create_pave_text(0.20, 0.63, 0.50, 0.83)
            for idx, hist in enumerate(hists):
                pave.AddText(
                    f"B{idx+1}: {hist.GetMean():.1f} #pm {hist.GetRMS():.1f}")

            pm.plot_1d(
                hists,
                f"{det}_difference_relative_peak_ts_wrt_board0",
                "Difference of Relative Peak TS",
                (DiffRelTSmin, DiffRelTSmax),
                yrange=(0, hists[0].GetMaximum() * 1.5),
                legends=labels[1:],
                style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[2, 3, 4],
                                legendNCols=3, legendPos=[0.20, 0.85, 0.90, 0.9]),
                addOverflow=False,
                addUnderflow=False,
                extraToDraw=pave
            )

        pm.add_newline()

        # Difference wrt US
        for det, channels in map_mcp_channels.items():
            hists = hists_diff_rel_1d_us[det]

            pave = create_pave_text(0.20, 0.63, 0.50, 0.83)
            for idx, hist in enumerate(hists):
                pave.AddText(
                    f"B{idx}: {hist.GetMean():.1f} #pm {hist.GetRMS():.1f}")

            pm.plot_1d(
                hists,
                f"{det}_difference_relative_peak_ts_wrt_US_board0",
                "Difference of Relative Peak TS wrt US",
                (DiffRelTSmin_US, DiffRelTSmax_US),
                yrange=(0, hists[1].GetMaximum() * 1.5),
                legends=labels,
                style=STYLE_MULTI_BOARD,
                addOverflow=False,
                addUnderflow=False,
                extraToDraw=pave
            )

        outputs_html["rel"] = pm.generate_html(
            "ServiceDRS/MCPPeakTS_relative.html", plots_per_row=4)

    for key, value in outputs_html.items():
        print(f"MCP {key} plots html: {value}")


if __name__ == "__main__":
    main()
