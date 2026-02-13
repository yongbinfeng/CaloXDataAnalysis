import json
import ROOT
from channels.channel_map import get_mcp_channels
from core.analysis_manager import CaloXAnalysisManager
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle, STYLE_2D_COLZ
from plotting.calox_plot_helper import create_pave_text
from plotting.my_function import LHistos2Hist
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.timing import auto_timer
from utils.utils import number_to_string
from variables.drs import calibrateDRSPeakTS
from channels.channel_map import build_fers_boards, build_drs_boards
from utils.data_loader import getRunInfo

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

ifilename = f"results/root/Run{run_number}/drs_vs_ts_calibrated.root"
ifile = ROOT.TFile.Open(ifilename, "READ")

drsboards = build_drs_boards(run=run_number)

beam_type, beam_energy = getRunInfo(run_number)


xmax = 14
xmin = -14
ymax = 10
ymin = -10

h3_outer = ROOT.TH3D("h3_outer", "Outer Tower 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin), ymin, ymax)
h3_center = ROOT.TH3D("h3_center", "Center Tower 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin) * 4, ymin, ymax)
h3_center_quartz = ROOT.TH3D("h3_center_quartz", "Center Tower Quartz 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin) * 4, ymin, ymax)
h3_center_plastic = ROOT.TH3D("h3_center_plastic", "Center Tower Plastic 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin) * 4, ymin, ymax)

# get the list of profiled histograms
hprofs_DRS_VS_TS = []

hlist = ifile.GetListOfKeys()
hnames = [key.GetName() for key in hlist]

for _, drsboard in drsboards.items():
    board_no = drsboard.board_no
    if board_no > 3:
        continue
    for i_tower_x, i_tower_y in drsboard.get_list_of_towers():
        sTowerX = number_to_string(i_tower_x)
        sTowerY = number_to_string(i_tower_y)

        channelNames = {}
        for var in ["Cer"]:
            chan_DRS = drsboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=(var == "Cer"))
            if chan_DRS is None:
                continue

            hname = f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}"
            if hname not in hnames:
                continue

            isQuartz = (chan_DRS.isQuartz)

            if abs(i_tower_x) <= 2 and abs(i_tower_y) <= 2:
                hprof = ifile.Get(hname)
                for ibin in range(1, hprof.GetNbinsX()+1):
                    drs_value = hprof.GetBinContent(ibin)
                    h3_center.Fill(i_tower_x, ibin, i_tower_y, drs_value)
                if isQuartz:
                    for ibin in range(1, hprof.GetNbinsX()+1):
                        drs_value = hprof.GetBinContent(ibin)
                        h3_center_quartz.Fill(
                            i_tower_x, ibin, i_tower_y, drs_value)
                else:
                    for ibin in range(1, hprof.GetNbinsX()+1):
                        drs_value = hprof.GetBinContent(ibin)
                        h3_center_plastic.Fill(
                            i_tower_x, ibin, i_tower_y, drs_value)
            else:
                hprof = ifile.Get(hname)
                for ibin in range(1, hprof.GetNbinsX()+1):
                    drs_value = hprof.GetBinContent(ibin)
                    h3_outer.Fill(i_tower_x, ibin, i_tower_y, drs_value)

plots = []


def Draw3DHist(h3, output_name, ts_min=230, ts_max=250, adc_min=20, var="Cer"):
    from plotting.tdrstyle import setTDRStyle
    setTDRStyle()
    c1 = ROOT.TCanvas("c1", "3D Histogram", 800, 600)
    c1.cd()

    c1.SetTheta(15.64)
    c1.SetPhi(106.68)

    # Set axis title offsets to prevent overlapping in 3D
    h3.GetXaxis().SetTitle("Tower X")
    h3.GetYaxis().SetTitle("Time [x 200ps]")
    h3.GetZaxis().SetTitle("Tower Y")
    h3.GetXaxis().SetTitleOffset(2.0)
    h3.GetYaxis().SetTitleOffset(2.0)
    h3.GetZaxis().SetTitleOffset(1.5)

    h3.GetXaxis().SetRangeUser(-2, 2)
    h3.GetZaxis().SetRangeUser(-2, 2)
    h3.GetYaxis().SetRangeUser(ts_min, ts_max)
    h3.SetMinimum(adc_min)

    h3.Draw("BOX2")

    btypes = {
        "pion": "#pi^{+}",
        "pions": "#pi^{+}",
        "pi+": "#pi^{+}",
        "positron": "e^{+}",
        "positrons": "e^{+}",
        "e+": "e^{+}",
        "mu+": "#mu^{+}",
    }

    pave = ROOT.TPaveText(0.15, 0.73, 0.4, 0.88, "NDC")
    pave.SetFillColor(0)          # White background
    pave.SetFillStyle(0)         # Transparent
    pave.SetLineColor(1)          # Black border
    pave.SetBorderSize(0)         # Thin border
    pave.SetTextAlign(12)         # Left-adjusted, vertically centered
    pave.SetTextFont(42)          # Standard Helvetica (non-bold)
    pave.SetTextSize(0.04)        # Medium size
    pave.AddText(f"Run {run_number}")
    pave.AddText(f"{btypes[beam_type]} at {beam_energy} GeV")
    pave.AddText(f"{var}")
    pave.Draw()

    # 5. Export to PDF
    c1.Modified()
    c1.Update()
    c1.Print(f"{output_name}.pdf")

    plots.append(output_name)


Draw3DHist(h3_center_plastic,
           f"results/plots/Run{run_number}/drs_timing_3D_center_plastic", var="Plastic")
Draw3DHist(h3_center_quartz,
           f"results/plots/Run{run_number}/drs_timing_3D_center_quartz", var="Quartz")

ofilename = f"results/root/Run{run_number}/drs_timing_3D.root"
ofile = ROOT.TFile.Open(ofilename, "RECREATE")
h3_outer.Write()
h3_center.Write()
h3_center_quartz.Write()
h3_center_plastic.Write()
ofile.Close()
