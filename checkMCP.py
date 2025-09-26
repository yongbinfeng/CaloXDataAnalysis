import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from channels.channel_map import getPreShowerChannel, getDownStreamMuonChannel, buildHodoPosChannels, getCerenkovCounters, buildFERSBoards, buildDRSBoards, getDownStreamTTUMuonChannel, getMCPChannels
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

TSmin = 450
TSmax = 650


rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf = preProcessDRSBoards(rdf)

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
    for channel in channels:
        rdf = rdf.Define(f"{channel}_PeakTS",
                         f"ArgMinRange({channel}_blsub, {TSmin}, {TSmax}, -70.0)")

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

output_dir = f"results/plots/Run{runNumber}/MCP"
plots = []
DrawHistos(hists_tf, [], 500, 1000, "Peak TS", 0, None, "Counts",
           outputname="Time_Reference_Peak_TS", outdir=output_dir,
           dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
           addOverflow=True, addUnderflow=True)
plots.append("Time_Reference_Peak_TS.png")
plots.append("\n")


for det in map_mcp_channels.keys():
    hists = hists_1d[det]
    outputname = f"{det}_peak_ts"
    DrawHistos(hists, [], TSmin, TSmax, "Peak TS", 0, None, "Counts",
               outputname=outputname, outdir=output_dir,
               dology=False, mycolors=[1, 2, 3, 4], drawashist=True, runNumber=runNumber,
               addOverflow=True, addUnderflow=True)
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

output_html = f"results/html/Run{runNumber}/MCP/index.html"
generate_html(plots, output_dir, plots_per_row=4, output_html=output_html)
