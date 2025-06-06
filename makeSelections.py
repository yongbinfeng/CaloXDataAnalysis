import sys
import ROOT
from utils.channel_map import build_map_Cer_Sci, build_map_ixy_DRSVar

print("Start running makeSelections.py")

map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)

# Open the input ROOT file
ifile = "/Users/yfeng/Desktop/TTU/CaloX/Data/run316_250517140056_converted.root"
infile = ROOT.TFile(ifile, "READ")
rdf = ROOT.RDataFrame("EventTree", infile)

tree = infile.Get("EventTree")
for board in range(0, 6):
    for ch in range(0, 64):
        rdf = rdf.Define(
            f"FERS_Board{board}_energyHG_{ch}",
            f"FERS_Board{board}_energyHG[{ch}]"
        )


hists1d = []
for board in [1]:
    for ch in range(0, 64):
        hist = rdf.Histo1D((
            f"hist_board{board}_ch{ch}",
            f"FERS Board {board + 1} - Channel {ch};Energy HG;Counts",
            3000, 0, 9000),
            f"FERS_Board{board}_energyHG_{ch}",
        )
        hists1d.append(hist)

hists2d = []
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        hist = rdf.Histo2D((
            f"hist_board{board}_Cer{iCer}_vs_Sci{iSci}",
            f"CER {iCer} vs SCI {iSci};CER Energy HG;SCI Energy HG",
            500, 0, 9000, 500, 0, 9000),
            f"FERS_Board1_energyHG_{iCer}",
            f"FERS_Board1_energyHG_{iSci}"
        )
        hist_zoomed = rdf.Histo2D((
            f"hist_board{board}_Cer{iCer}_vs_Sci{iSci}_zoom",
            f"CER {iCer} vs SCI {iSci} (zoomed);CER Energy HG;SCI Energy HG",
            500, 0, 1000, 1000, 0, 2000),
            f"FERS_Board1_energyHG_{iCer}",
            f"FERS_Board1_energyHG_{iSci}"
        )
        hists2d.append(hist)
        hists2d.append(hist_zoomed)

hists1d_DRS = []
for var in list(map_ixy_DRSVar_Cer.values()) + list(map_ixy_DRSVar_Sci.values()):
    hist = rdf.Histo1D((
        f"hist_DRS_{var}",
        f"DRS Variable {var};Counts;Events",
        1500, 1000, 2500),
        var
    )
    hists1d_DRS.append(hist)

    # filter some events for displays and analysis
rdfs_temp = []
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        requirement = (
            f"FERS_Board{board}_energyHG_{iCer} > 1000 && "
            f"FERS_Board{board}_energyHG_{iSci} > 500"
        )
        rdf_temp = rdf.Filter(requirement)
        rdfs_temp.append(rdf_temp)

print("Save histograms and filtered RDataFrames")

# Save histograms to an output ROOT file
outfile = ROOT.TFile("root/fers_all_channels.root", "RECREATE")
for hist in hists1d:
    hist.Write()
for hist in hists2d:
    hist.Write()
outfile.Close()

outfile_DRS = ROOT.TFile("root/fers_all_DRS_variables.root", "RECREATE")
for hist in hists1d_DRS:
    hist.SetDirectory(outfile_DRS)
    hist.Write()
outfile_DRS.Close()

# save the filtered RDataFrames
for i, rdf_temp in enumerate(rdfs_temp):
    if i < 2:
        print(
            f"Events left after filtering for board 1, CER {iCer}, SCI {iSci}: {rdf_temp.Count().GetValue()}")
        rdf_temp.Snapshot(
            "EventTree", f"root/filtered_events_board1_cersci_{i}.root")
