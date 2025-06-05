import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import build_map_Cer_Sci, build_map_ixy_FERS1, build_map_ixy_DRS

ROOT.gROOT.SetBatch(True)  # Run in batch mode


map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)
map_ixy_DRS = build_map_ixy_DRS()
print("DRS mapping ", map_ixy_DRS)
map_ixy_FERS1 = build_map_ixy_FERS1()
print("FERS mapping ", map_ixy_FERS1)


def make_event_displays(infilename, prefix=""):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    hists_eventdisplay = []

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        hist2d_Cer = ROOT.TH2F(
            f"event_display_Evt{ievt}_Cer",
            f"Event Display {ievt};X;Y (Cherenkov)",
            4, 0, 4, 8, 0, 8
        )
        hist2d_Sci = ROOT.TH2F(
            f"event_display_Evt{ievt}_Sci",
            f"Event Display {ievt};X;Y (Scintillator)",
            4, 0, 4, 8, 0, 8
        )
        for iCer, iSci in map_Cer_Sci.items():
            # get the energy for CER and SCI
            e_Cer = rdf.Take["unsigned short"](
                f"FERS_Board1_energyHG_{iCer}").GetValue()[ievt]
            e_Sci = rdf.Take["unsigned short"](
                f"FERS_Board1_energyHG_{iSci}").GetValue()[ievt]

            # get the ixy for CER and SCI
            ix_Cer, iy_Cer = map_ixy_FERS1[iCer]
            ix_Sci, iy_Sci = map_ixy_FERS1[iSci]

            # fill the histograms
            hist2d_Cer.Fill(ix_Cer, iy_Cer, e_Cer)
            hist2d_Sci.Fill(ix_Sci, iy_Sci, e_Sci)
        hists_eventdisplay.append(hist2d_Cer)
        hists_eventdisplay.append(hist2d_Sci)

        DrawHistos([hist2d_Cer], f"", 0, 4, "iX",
                   0, 8, "iY", f"{prefix}_event_display_Evt{ievt}_Cer", dology=False, drawoptions=["COLZ,text"], zmin=200.0, zmax=3000.0, doth2=True)
        DrawHistos([hist2d_Sci], f"", 0, 4, "iX",
                   0, 8, "iY", f"{prefix}_event_display_Evt{ievt}_Sci", dology=False, drawoptions=["COLZ,text"], zmin=200.0, zmax=9000.0, doth2=True)
    print(f"Events left after filtering: {rdf.Count().GetValue()}")

    # Save event display histograms
    ofilename = infilename.replace(".root", "_event_display.root")
    ofile = ROOT.TFile(ofilename, "RECREATE")
    if not ofile or ofile.IsZombie():
        raise RuntimeError(f"Failed to create output file: {ofilename}")
    print(f"Saving event display histograms to {ofilename}")
    for hist in hists_eventdisplay:
        hist.SetDirectory(ofile)
        hist.Write()


if __name__ == "__main__":
    files = [
        "root/filtered_events_board1_cersci_0.root",
        "root/filtered_events_board1_cersci_1.root",
        "root/filtered_events_board1_cersci_2.root",
        "root/filtered_events_board1_cersci_3.root",
        "root/filtered_events_board1_cersci_4.root",
    ]
    for idx, input_file in enumerate(files):
        print(f"Processing file: {input_file}")
        make_event_displays(input_file, prefix=f"{idx}")
