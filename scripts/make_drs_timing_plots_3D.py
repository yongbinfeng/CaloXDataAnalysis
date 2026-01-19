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

auto_timer("Total Execution Time")

# setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

ifilename = f"results/root/Run{run_number}/drs_vs_ts_calibrated.root"

ifile = ROOT.TFile.Open(ifilename, "READ")


xmax = 14
xmin = -14
ymax = 10
ymin = -10

h3_outer = ROOT.TH3D("h3_outer", "Outer Tower 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin), ymin, ymax)
h3_center = ROOT.TH3D("h3_center", "Center Tower 3D Histogram", int(
    xmax-xmin), xmin, xmax, 400, 0, 400, int(ymax - ymin) * 4, ymin, ymax)

# get the list of profiled histograms
hprofs_DRS_VS_TS = []

hlist = ifile.GetListOfKeys()
for key in hlist:
    # name format: prof_DRS_vs_TS_Cer_0p5_4p5
    name = key.GetName()
    if "prof_DRS_vs_TS_Cer" not in name:
        continue
    tower_ix = float(name.split("_")[-2].replace("p", ".").replace("m", "-"))
    tower_iy = float(name.split("_")[-1].replace("p", ".").replace("m", "-"))

    if abs(tower_ix) <= 2 and abs(tower_iy) <= 2:
        # center tower
        hprof = ifile.Get(name)
        for ibin in range(1, hprof.GetNbinsX()+1):
            # ts = hprof.GetXaxis().GetBinCenter(ibin)
            drs_value = hprof.GetBinContent(ibin)
            h3_center.Fill(tower_ix, ibin, tower_iy, drs_value)
    else:
        # outer tower
        hprof = ifile.Get(name)
        for ibin in range(1, hprof.GetNbinsX()+1):
            # ts = hprof.GetXaxis().GetBinCenter(ibin)
            drs_value = hprof.GetBinContent(ibin)
            h3_outer.Fill(tower_ix, ibin, tower_iy, drs_value)

ofilename = f"results/root/Run{run_number}/drs_timing_3D.root"
ofile = ROOT.TFile.Open(ofilename, "RECREATE")
h3_outer.Write()
h3_center.Write()
ofile.Close()
