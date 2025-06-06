import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import build_map_Cer_Sci, build_map_FERS1_ixy, build_map_ixy_DRSVar
import json

ROOT.gROOT.SetBatch(True)  # Run in batch mode


map_Cer_Sci = build_map_Cer_Sci()
print("Map of CER to SCI channels in FERS1:", map_Cer_Sci)
map_FERS1_ixy = build_map_FERS1_ixy()
print("Map of FERS1 channels to (ix, iy):", map_FERS1_ixy)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()
print("Map of DRS variable names to (ix, iy) (CER):", map_ixy_DRSVar_Cer)
print("Map of DRS variable names to (ix, iy) (SCI):", map_ixy_DRSVar_Sci)

# read the noise from json file
with open("results/drs_noises.json", "r") as json_file:
    drs_noises = json.load(json_file)


def make_event_displays(infilename, prefix=""):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    hists_eventdisplay = []
    hists_pulse_shapes = []

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        hist2d_Cer = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Cer",
            f"Event Display {evtNumber};X;Y (Cherenkov)",
            4, -0.5, 3.5, 8, -0.5, 7.5
        )
        hist2d_Sci = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            4, -0.5, 3.5, 8, -0.5, 7.5
        )
        for iCer, iSci in map_Cer_Sci.items():
            # get the energy for CER and SCI
            e_Cer = rdf.Take["unsigned short"](
                f"FERS_Board1_energyHG_{iCer}").GetValue()[ievt]
            e_Sci = rdf.Take["unsigned short"](
                f"FERS_Board1_energyHG_{iSci}").GetValue()[ievt]

            # get the ixy for CER and SCI
            ix_Cer, iy_Cer = map_FERS1_ixy[iCer]
            ix_Sci, iy_Sci = map_FERS1_ixy[iSci]
            # they should be exactly the same
            assert ix_Cer == ix_Sci, f"ix mismatch: {ix_Cer} != {ix_Sci}"
            assert iy_Cer == iy_Sci, f"iy mismatch: {iy_Cer} != {iy_Sci}"

            # fill the histograms
            hist2d_Cer.Fill(ix_Cer, iy_Cer, e_Cer)
            hist2d_Sci.Fill(ix_Sci, iy_Sci, e_Sci)

            # print the DRS variable names and their pulse shapes
            varname_Cer = map_ixy_DRSVar_Cer[(ix_Cer, iy_Cer)]
            varname_Sci = map_ixy_DRSVar_Sci[(ix_Sci, iy_Sci)]
            pulse_shape_Cer = rdf.Take["ROOT::VecOps::RVec<float>"](
                f"{varname_Cer}").GetValue()[ievt]
            pulse_shape_Sci = rdf.Take["ROOT::VecOps::RVec<float>"](
                f"{varname_Sci}").GetValue()[ievt]

            h1_Cer = ROOT.TH1F(
                f"pulse_shape_Cer_Evt{evtNumber}_iX{ix_Cer}_iY{iy_Cer}",
                f"Pulse Shape {varname_Cer} (CER);Time;Amplitude",
                1024, 0, 1024
            )
            h1_Sci = ROOT.TH1F(
                f"pulse_shape_Sci_Evt{evtNumber}_iX{ix_Sci}_iY{iy_Sci}",
                f"Pulse Shape {varname_Sci} (SCI);Time;Amplitude",
                1024, 0, 1024
            )
            for i in range(len(pulse_shape_Cer)):
                h1_Cer.Fill(i, pulse_shape_Cer[i] - drs_noises[varname_Cer])
            for i in range(len(pulse_shape_Sci)):
                h1_Sci.Fill(i, pulse_shape_Sci[i] - drs_noises[varname_Sci])

            hists_pulse_shapes.append(h1_Cer)
            hists_pulse_shapes.append(h1_Sci)

            DrawHistos([h1_Cer, h1_Sci], ["Cer", "Sci"], 0, 1024, "TS",
                       0, 30, "Amplitude", f"{prefix}_pulse_shape_Evt{evtNumber}_iX{ix_Cer}_iY{iy_Cer}", dology=False, mycolors=[2, 4], drawashist=True)

        hists_eventdisplay.append(hist2d_Cer)
        hists_eventdisplay.append(hist2d_Sci)

        DrawHistos([hist2d_Cer], f"", -0.5, 3.5, "iX",
                   -0.5, 7.5, "iY", f"{prefix}_event_display_Evt{evtNumber}_Cer", dology=False, drawoptions=["COLZ,text"], zmin=200.0, zmax=3000.0, doth2=True)
        DrawHistos([hist2d_Sci], f"", -0.5, 3.5, "iX",
                   -0.5, 7.5, "iY", f"{prefix}_event_display_Evt{evtNumber}_Sci", dology=False, drawoptions=["COLZ,text"], zmin=200.0, zmax=9000.0, doth2=True)
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
    ofile.Close()

    # Save pulse shape histograms
    ofilename_pulse_shapes = infilename.replace(".root", "_pulse_shapes.root")
    ofile = ROOT.TFile(ofilename_pulse_shapes, "RECREATE")
    if not ofile or ofile.IsZombie():
        raise RuntimeError(
            f"Failed to create output file: {ofilename_pulse_shapes}")
    print(f"Saving pulse shape histograms to {ofilename_pulse_shapes}")
    for hist in hists_pulse_shapes:
        hist.SetDirectory(ofile)
        hist.Write()
    ofile.Close()


if __name__ == "__main__":
    files = [
        "root/filtered_events_board1_cersci_0.root",
        # "root/filtered_events_board1_cersci_1.root",
        # "root/filtered_events_board1_cersci_2.root",
        # "root/filtered_events_board1_cersci_3.root",
        # "root/filtered_events_board1_cersci_4.root",
    ]
    for idx, input_file in enumerate(files):
        print(f"Processing file: {input_file}")
        make_event_displays(input_file, prefix=f"{idx}")
