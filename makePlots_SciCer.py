import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos
from utils.channel_map import build_map_Cer_Sci, build_map_ixy_DRSVar

print("Start running makePlots_SciCer.py")
# batch mode
ROOT.gROOT.SetBatch(True)

map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()

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
for board in [1, 2, 3, 4, 5]:
    for iCer, iSci in map_Cer_Sci.items():
        hist_C_name = f"hist_board{board}_ch{iCer}"
        hist_C = infile.Get(hist_C_name)
        hists1d[hist_C_name] = hist_C
        hist_S_name = f"hist_board{board}_ch{iSci}"
        hist_S = infile.Get(hist_S_name)
        hists1d[hist_S_name] = hist_S

        extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(
            f"Board: {board}")
        extraToDraw.AddText(f"Cer Channel: {iCer}")
        extraToDraw.AddText(f"Sci Channel: {iSci}")

        DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                   f"Energy_Board{board}_Channel{iCer}", dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, extraToDraw=extraToDraw)

# 1D histograms for DRS variables
hists1d_DRS = {}
infile_DRS = ROOT.TFile("root/drs_all_channels.root", "READ")
for (ix, iy), var_cer in map_ixy_DRSVar_Cer.items():
    var_sci = map_ixy_DRSVar_Sci[(ix, iy)]
    hist_cer_name = f"hist_{var_cer}"
    hist_cer = infile_DRS.Get(hist_cer_name)
    hists1d_DRS[hist_cer_name] = hist_cer
    hist_sci_name = f"hist_{var_sci}"
    hist_sci = infile_DRS.Get(hist_sci_name)
    hists1d_DRS[hist_sci_name] = hist_sci
    DrawHistos([hist_cer, hist_sci], ["Cer", "Sci"], 1400, 2500, "DRS Output", 1, 1e10, "Counts",
               f"DRS_Variable_iX{ix}_iY{iy}", dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True)
