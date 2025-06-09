import sys
import ROOT
from utils.channel_map import build_map_Cer_Sci, build_map_ixy_DRSVar, get_hodoscope_channels, build_map_FERSs_ixy

print("Start running makeSelections.py")

map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()
hodoscope_channels = get_hodoscope_channels()
map_ixy_FERSs = build_map_FERSs_ixy()

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
for board in [1, 2, 3, 4, 5]:
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
        f"hist_{var}",
        f"DRS Variable {var};Counts;Events",
        1500, 1000, 2500),
        var
    )
    hists1d_DRS.append(hist)

hists2d_DRS = []
for (ix, iy), var_cer in map_ixy_DRSVar_Cer.items():
    var_sci = map_ixy_DRSVar_Sci[(ix, iy)]
    hist = rdf.Histo2D((
        f"hist_DRS_iX{ix}_iY{iy}",
        f"DRS Variables iX {ix} vs iY {iy};CER {var_cer};SCI {var_sci}",
        500, 1000, 2500, 500, 1000, 2500),
        var_cer,
        var_sci
    )
    hists2d_DRS.append(hist)

hists1d_hodos = []
for group, channels in hodoscope_channels.items():
    for channel in channels:
        hist = rdf.Histo1D((
            f"hist_hodoscope_{channel}",
            f"Hodoscope {group} - Channel {channel};Amplitude;Counts",
            3000, 0, 3000),
            channel
        )
        hists1d_hodos.append(hist)

# filter some events for displays and analysis
requirement = ""
idx = 0
for board in [1]:
    for iCer, iSci in map_Cer_Sci.items():
        if idx != 0:
            requirement += " || "
        requirement += f"(FERS_Board{board}_energyHG_{iCer} > 1000 && FERS_Board{board}_energyHG_{iSci} > 500)"
        idx += 1
print(f"Filtering events with requirement: {requirement}")
rdf_filtered = rdf.Filter(requirement)

print("Save histograms and filtered RDataFrames")

# Save histograms to an output ROOT file
outfile = ROOT.TFile("root/fers_all_channels.root", "RECREATE")
for hist in hists1d:
    hist.Write()
for hist in hists2d:
    hist.Write()
outfile.Close()

outfile_DRS = ROOT.TFile("root/drs_all_channels.root", "RECREATE")
for hist in hists1d_DRS:
    hist.SetDirectory(outfile_DRS)
    hist.Write()
for hist in hists2d_DRS:
    hist.SetDirectory(outfile_DRS)
    hist.Write()
outfile_DRS.Close()

outfile_hodos = ROOT.TFile("root/hodoscope_all_channels.root", "RECREATE")
for hist in hists1d_hodos:
    hist.SetDirectory(outfile_hodos)
    hist.Write()
outfile_hodos.Close()


# save the filtered RDataFrames
print(
    f"Events left after filtering for board 1, CER {iCer}, SCI {iSci}: {rdf_filtered.Count().GetValue()}")
rdf_filtered.Snapshot(
    "EventTree", f"root/filtered_events_board1.root")
