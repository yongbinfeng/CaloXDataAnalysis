import sys
import ROOT

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

hist_list = []
for board in range(0, 6):
    for ch in range(0, 64):
        hist = rdf.Histo1D((
            f"hist_board{board}_ch{ch}",
            f"FERS Board {board + 1} - Channel {ch};Energy HG;Counts",
            1000, 0, 1000),
            f"FERS_Board{board}_energyHG_{ch}",
        )
        hist_list.append(hist)

# Save histograms to an output ROOT file
outfile = ROOT.TFile("root/fers_all_channels.root", "RECREATE")
for hist in hist_list:
    hist.Write()
outfile.Close()
