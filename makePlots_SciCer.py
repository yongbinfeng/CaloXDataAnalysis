import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import build_map_Cer_Sci

ROOT.gROOT.SetBatch(True)  # Run in batch mode

map_Cer_Sci = build_map_Cer_Sci()

# 2D histograms
infile = ROOT.TFile("root/fers_all_channels.root", "READ")
hists2d = {}
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        hist_name = f"hist_board{board}_Cer{iCer}_vs_Sci{iSci}"
        hist = infile.Get(hist_name)
        hist_zoom = infile.Get(f"{hist_name}_zoom")

        hists2d[hist_name] = hist
        DrawHistos([hist], f"", 0, 1500, "Cer", 0, 9000, "Sci",
                   f"Cer_Sci_Board{board}_Cer{iCer}_vs_Sci{iSci}", dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True)
        DrawHistos([hist_zoom], f"", 0, 500, "Cer", 0, 1000, "Sci",
                   f"Cer_Sci_Board{board}_Cer{iCer}_vs_Sci{iSci}_zoom", dology=False, drawoptions="COLZ", doth2=True, zmin=0, zmax=100)

# 1D histograms
hists1d = {}
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        hist_C_name = f"hist_board{board}_ch{iCer}"
        hist_C = infile.Get(hist_C_name)
        hists1d[hist_C_name] = hist_C
        hist_S_name = f"hist_board{board}_ch{iSci}"
        hist_S = infile.Get(hist_S_name)
        hists1d[hist_S_name] = hist_S
        DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                   f"Energy_Board{board}_Channel{iCer}", dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True)
