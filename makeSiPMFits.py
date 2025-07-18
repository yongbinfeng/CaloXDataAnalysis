import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.html_generator import generate_html
from utils.fitter import runFit
from runconfig import runNumber
from utils.visualization import visualizeFERSBoards
from utils.channel_map import buildFERSBoards
from utils.utils import number2string
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)  # Disable interactive mode for batch processing

FERSBoards = buildFERSBoards(runNumber)

infile = f"results/root/Run{runNumber}/fers_all_channels_1D.root"
outdir = f"results/plots/Run{runNumber}/fit_plots"
ifile = ROOT.TFile(infile)
hists = ifile.GetListOfKeys()
hists_3mm = OrderedDict()
hists_6mm = OrderedDict()

for _, FERSBoard in FERSBoards.items():
    if FERSBoard.Is3mm():
        hists = hists_3mm
    else:
        hists = hists_6mm
    for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
        sTowerX = number2string(iTowerX)
        sTowerY = number2string(iTowerY)
        for var in ["Cer", "Sci"]:
            hname = f"hist_FERS_Board{FERSBoard.boardNo}_{var}_{sTowerX}_{sTowerY}"
            chan = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=(var == "Cer"))
            if chan is None:
                print(f"Warning: Channel not found for {hname}")
                continue

            channelName = chan.GetHGChannelName()

            print("hname =", hname, "channelName =", channelName)

            hists[channelName] = ifile.Get(hname).Clone(f"{hname}_clone")


# for h in hists:
#    hname = h.GetName()
#    if "hist_FERS_Board11" in hname or "hist_FERS_Board3" in hname:
#        hists_3mm.append(hname)
#    else:
#        hists_6mm.append(hname)

plots = []
valuemaps_pedestal = {}
valuemaps_gain = {}
for channelName, h in hists_6mm.items():
    print(f"Processing histogram: {channelName}")
    mu_ped, dpe = runFit(h, outdir, channelName, is3mm=False, npe_max=4)
    plots.append(channelName + ".png")
    valuemaps_pedestal[channelName] = mu_ped
    valuemaps_gain[channelName] = dpe

for channelName, h in hists_3mm.items():
    print(f"Processing histogram: {channelName}")
    mu_ped, dpe = runFit(h, outdir, channelName, is3mm=True, npe_max=4)
    plots.append(channelName + ".png")
    valuemaps_pedestal[channelName] = mu_ped
    valuemaps_gain[channelName] = dpe

generate_html(plots, outdir, plots_per_row=4,
              output_html=f"results/html/Run{runNumber}/sipm_fit_results/fits.html")

plots = []
[h2_Cer_pedestal, h2_Cer_3mm_pedestal], [h2_Sci_pedestal, h2_Sci_3mm_pedestal] = visualizeFERSBoards(
    FERSBoards, valuemaps_pedestal, suffix=f"Run{runNumber}_fitted_pedestal", useHG=True)
[h2_Cer_gain, h2_Cer_3mm_gain], [h2_Sci_gain, h2_Sci_3mm_gain] = visualizeFERSBoards(
    FERSBoards, valuemaps_gain, suffix=f"Run{runNumber}_fitted_gain", useHG=True)

output_name = f"FERS_Boards_Run{runNumber}_Stats_pedestal"
xmax = 14
xmin = -14
ymax = 10
ymin = -10
W_ref = 1000
H_ref = 1100
outdir_plots = f"results/plots/Run{runNumber}/fit_plots"
DrawHistos([h2_Cer_pedestal, h2_Cer_3mm_pedestal], "", xmin, xmax, "iX", ymin,
           ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
           outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=120, zmax=160)
plots.append(output_name + "_Cer.png")
DrawHistos([h2_Sci_pedestal, h2_Sci_3mm_pedestal], "", xmin, xmax, "iX", ymin,
           ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
           outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=120, zmax=160)
plots.append(output_name + "_Sci.png")
output_name = f"FERS_Boards_Run{runNumber}_Stats_gain"
DrawHistos([h2_Cer_gain, h2_Cer_3mm_gain], "", xmin, xmax, "iX", ymin,
           ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
           outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=50, zmax=80)
plots.append(output_name + "_Cer.png")
DrawHistos([h2_Sci_gain, h2_Sci_3mm_gain], "", xmin, xmax, "iX", ymin,
           ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
           outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=50, zmax=80)
plots.append(output_name + "_Sci.png")

generate_html(plots, outdir_plots, plots_per_row=2,
              output_html=f"results/html/Run{runNumber}/sipm_fit_results/summary.html")
