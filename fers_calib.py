import sys
import ROOT
from utils.channel_map import build_map_Cer_Sci

map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)

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
            300, 0, 9000),
            f"FERS_Board{board}_energyHG_{ch}",
        )
        hists1d.append(hist)

hists2d = []
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        hist = rdf.Histo2D((
            f"hist_board{board}_Cer{iCer}_vs_Sci{iSci}",
            f"CER {iCer} vs SCI {iSci};CER Energy HG;SCI Energy HG",
            300, 0, 9000, 300, 0, 9000),
            f"FERS_Board1_energyHG_{iCer}",
            f"FERS_Board1_energyHG_{iSci}"
        )
        hists2d.append(hist)

# Save histograms to an output ROOT file
outfile = ROOT.TFile("root/fers_all_channels.root", "RECREATE")
for hist in hists1d:
    hist.Write()
for hist in hists2d:
    hist.Write()
outfile.Close()
