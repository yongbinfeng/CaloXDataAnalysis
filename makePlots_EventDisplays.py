import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import build_map_Cer_Sci, build_map_FERS1_ixy, build_map_ixy_DRSVar, build_map_FERSs_ixy, get_hodoscope_channels, hodoTS2iX
import json
from results.events import events_interested

ROOT.gROOT.SetBatch(True)  # Run in batch mode


map_Cer_Sci = build_map_Cer_Sci()
print("Map of CER to SCI channels in FERS1:", map_Cer_Sci)
map_FERS1_ixy = build_map_FERS1_ixy()
print("Map of FERS1 channels to (ix, iy):", map_FERS1_ixy)
map_FERSs_ixy = build_map_FERSs_ixy()
print("Map of FERS channels to (ix, iy):", map_FERSs_ixy)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()
print("Map of DRS variable names to (ix, iy) (CER):", map_ixy_DRSVar_Cer)
print("Map of DRS variable names to (ix, iy) (SCI):", map_ixy_DRSVar_Sci)


hodoscope_channels = get_hodoscope_channels()

with open("results/hodoscope_noises.json", "r") as json_file:
    hodo_noises = json.load(json_file)

# read the noise from json file
with open("results/drs_noises.json", "r") as json_file:
    drs_noises = json.load(json_file)

with open("results/fers_noises.json", "r") as json_file:
    fers_noises = json.load(json_file)


def make_event_displays(infilename, onlyFERS1=True):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    hists_eventdisplay = []
    hists_pulse_shapes = []

    boards = [1] if onlyFERS1 else [1, 2, 3, 4, 5]
    ix_min = -0.5 if onlyFERS1 else -16.5
    iy_min = -0.5 if onlyFERS1 else -4.5
    W_ref = 800 if onlyFERS1 else 1400

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        # if ievt > 100:
        #    break
        print(f"Processing event {ievt + 1} of {rdf.Count().GetValue()}")
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        if evtNumber not in events_interested:
            print(
                f"Skipping event {evtNumber} as it is not in the interested events list.")
            continue
        hist2d_Cer = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Cer",
            f"Event Display {evtNumber};X;Y (Cherenkov)",
            int(3.5 - ix_min), ix_min, 3.5, int(7.5 - iy_min), iy_min, 7.5
        )
        hist2d_Sci = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(3.5 - ix_min), ix_min, 3.5, int(7.5 - iy_min), iy_min, 7.5
        )
        for iCer, iSci in map_Cer_Sci.items():
            # get the energy for CER and SCI
            for board in boards:
                e_Cer = rdf.Take["unsigned short"](
                    f"FERS_Board{board}_energyHG_{iCer}").GetValue()[ievt]
                e_Sci = rdf.Take["unsigned short"](
                    f"FERS_Board{board}_energyHG_{iSci}").GetValue()[ievt]

                # get the ixy for CER and SCI
                ix_Cer_FERS, iy_Cer_FERS = map_FERSs_ixy[f"Board{board}"][iCer]
                ix_Sci_FERS, iy_Sci_FERS = map_FERSs_ixy[f"Board{board}"][iSci]

                # fill the histograms
                hist2d_Cer.Fill(ix_Cer_FERS, iy_Cer_FERS, int(e_Cer -
                                fers_noises[f"board{board}_ch{iCer}"]))
                hist2d_Sci.Fill(ix_Sci_FERS, iy_Sci_FERS, int(e_Sci -
                                fers_noises[f"board{board}_ch{iSci}"]))

            # print the DRS variable names and their pulse shapes
            ix_Cer_FERS1, iy_Cer_FERS1 = map_FERS1_ixy[iCer]
            ix_Sci_FERS1, iy_Sci_FERS1 = map_FERS1_ixy[iSci]
            varname_Cer = map_ixy_DRSVar_Cer[(ix_Cer_FERS1, iy_Cer_FERS1)]
            varname_Sci = map_ixy_DRSVar_Sci[(ix_Sci_FERS1, iy_Sci_FERS1)]
            pulse_shape_Cer = rdf.Take["ROOT::VecOps::RVec<float>"](
                f"{varname_Cer}").GetValue()[ievt]
            pulse_shape_Sci = rdf.Take["ROOT::VecOps::RVec<float>"](
                f"{varname_Sci}").GetValue()[ievt]

            h1_Cer = ROOT.TH1F(
                f"pulse_shape_Cer_Evt{evtNumber}_iX{ix_Cer_FERS1}_iY{iy_Cer_FERS1}",
                f"Pulse Shape {varname_Cer} (CER);Time;Amplitude",
                1024, 0, 1024
            )
            h1_Sci = ROOT.TH1F(
                f"pulse_shape_Sci_Evt{evtNumber}_iX{ix_Sci_FERS1}_iY{iy_Sci_FERS1}",
                f"Pulse Shape {varname_Sci} (SCI);Time;Amplitude",
                1024, 0, 1024
            )
            for i in range(len(pulse_shape_Cer)):
                h1_Cer.Fill(i, pulse_shape_Cer[i] - drs_noises[varname_Cer])
            for i in range(len(pulse_shape_Sci)):
                h1_Sci.Fill(i, pulse_shape_Sci[i] - drs_noises[varname_Sci])

            hists_pulse_shapes.append(h1_Cer)
            hists_pulse_shapes.append(h1_Sci)

            if ievt < 1000000:
                peak_Cer = h1_Cer.GetMaximumBin()
                peak_Sci = h1_Sci.GetMaximumBin()
                deltaTS = peak_Sci - peak_Cer

                extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Event: {evtNumber}, iX: {ix_Cer_FERS1}, iY: {iy_Cer_FERS1}")
                extraToDraw.AddText(f"Cer Peak: {peak_Cer}")
                extraToDraw.AddText(f"Sci Peak: {peak_Sci}")
                extraToDraw.AddText(
                    f"delta Peak: {deltaTS} ts")
                extraToDraw.AddText(f"delta T: {deltaTS * 0.2: .2f} ns")

                DrawHistos([h1_Cer, h1_Sci], ["Cer", "Sci"], 0, 1024, "TS",
                           -5, 30, "Amplitude", f"pulse_shape_Evt{evtNumber}_iX{ix_Cer_FERS1}_iY{iy_Cer_FERS1}", dology=False, mycolors=[2, 4], drawashist=True, extraToDraw=extraToDraw)

        hists_eventdisplay.append(hist2d_Cer)
        hists_eventdisplay.append(hist2d_Sci)

        # analyze hodoscope channels
        h1_hodos = {}
        for group, channels in hodoscope_channels.items():
            if group != "TopX" and group != "BottomX":
                continue
            print(f"Analyzing hodoscope group: {group}")
            hs_todraw = []
            for channel in channels:
                pulse_shape = rdf.Take["ROOT::VecOps::RVec<float>"](
                    channel).GetValue()[ievt]
                h1 = ROOT.TH1F(
                    f"pulse_shape_{channel}_Evt{evtNumber}",
                    f"Pulse Shape {channel};Time;Amplitude",
                    1024, 0, 1024
                )
                print("make pulse shape for channel:", channel)
                for i in range(len(pulse_shape)):
                    h1.Fill(
                        i, pulse_shape[i] - hodo_noises["hodoscope_" + channel])
                hs_todraw.append(h1)
            h1_hodos[group] = hs_todraw

        peak_TopX_left = h1_hodos['TopX'][0].GetMinimumBin()
        peak_TopX_right = h1_hodos['TopX'][1].GetMinimumBin()
        deltaTS_TopX = peak_TopX_right - peak_TopX_left
        peak_BottomX_left = h1_hodos['BottomX'][0].GetMinimumBin()
        peak_BottomX_right = h1_hodos['BottomX'][1].GetMinimumBin()
        deltaTS_BottomX = peak_BottomX_right - peak_BottomX_left

        extraToDraw = ROOT.TPaveText(0.01, 0.55, 0.12, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(
            f"Event: {evtNumber}")
        # extraToDraw.AddText(f"TopX Left Peak: {peak_TopX_left}")
        # extraToDraw.AddText(f"TopX Right Peak: {peak_TopX_right}")
        # extraToDraw.AddText(f"BottomX Left Peak: {peak_BottomX_left}")
        # extraToDraw.AddText(f"BottomX Right Peak: {peak_BottomX_right}")
        extraToDraw.AddText(f"TopX Diff: {deltaTS_TopX}")
        extraToDraw.AddText(f"BottomX Diff: {deltaTS_BottomX}")
        extraToDraw.AddText(f"Top iX: {hodoTS2iX(deltaTS_TopX)}")
        extraToDraw.AddText(f"Bottom iX: {hodoTS2iX(deltaTS_BottomX)}")

        if ievt < 100000:
            DrawHistos([hist2d_Cer], f"", ix_min, 3.5, "iX",
                       iy_min, 7.5, "iY", f"event_display_Evt{evtNumber}_Cer", dology=False, drawoptions=["COLZ,text"], zmin=50.0, zmax=3000.0, doth2=True, W_ref=W_ref, extraToDraw=extraToDraw)
            DrawHistos([hist2d_Sci], f"", ix_min, 3.5, "iX",
                       iy_min, 7.5, "iY", f"event_display_Evt{evtNumber}_Sci", dology=False, drawoptions=["COLZ,text"], zmin=50.0, zmax=8000.0, doth2=True, W_ref=W_ref, extraToDraw=extraToDraw)
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
        "root/filtered_events_board1.root",
    ]
    for idx, input_file in enumerate(files):
        print(f"Processing file: {input_file}")
        make_event_displays(input_file, onlyFERS1=False)
