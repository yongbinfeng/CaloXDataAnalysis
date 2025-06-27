#!/usr/bin/env python3

import sys
import os
import time

# ensure CMSPLOTS is on the path
sys.path.append("CMSPLOTS")  # noqa

import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile
from utils.html_generator import generate_html

ROOT.gROOT.SetBatch(True)  # Run in batch mode

runNumber = 624
outdir = f"plots/Run{runNumber}/"


def makeEventDisplays(infilename, random_per_block=None, block_size=100,
                      nth_interval=None):
    """Generate event and pulse-shape displays.

    Parameters
    ----------
    infilename : str
        Input ROOT file containing the event tree.
    random_per_block : int or None, optional
        Number of random events to include in every ``block_size`` events.
    block_size : int, optional
        Block size used for ``random_per_block``.
    nth_interval : int or None, optional
        Additionally include every ``nth_interval``-th event.
    """
    start_time = time.time()
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Grab the TTree instead of RDataFrame
    tree = infile.Get("EventTree")
    if not tree:
        raise RuntimeError("Failed to get EventTree from file")

    DRSBoards = buildDRSBoards(run=runNumber)
    FERSBoards = buildFERSBoards(run=runNumber)

    hists_eventdisplay = []
    hists_pulse_shapes = []
    plots_eventdisplay = []
    plots_pulse_shapes = []

    W_ref = 1000
    iTowerX_max = 5.5
    iTowerX_min = -12.5
    iTowerY_max = 10.5
    iTowerY_min = -9.5

    nevents = tree.GetEntries()

    indices = set()
    if nth_interval and nth_interval > 0:
        indices.update(range(0, nevents, nth_interval))
    if random_per_block and block_size and random_per_block > 0:
        import random
        for start in range(0, nevents, block_size):
            end = min(start + block_size, nevents)
            count = min(random_per_block, end - start)
            indices.update(random.sample(range(start, end), count))
    if not indices:
        nproc = min(nevents, 100)
        indices.update(range(nproc))

    for ievt in sorted(indices):
        print(f"Processing event {ievt + 1} of {nevents}")
        tree.GetEntry(ievt)
        evtNumber = int(tree.event_n)

        # Build your 2D histograms
        hist2d_Cer = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Cer",
            f"Event Display {evtNumber};X;Y (Cherenkov)",
            int(iTowerX_max - iTowerX_min), iTowerX_min, iTowerX_max,
            int(iTowerY_max - iTowerY_min), iTowerY_min, iTowerY_max
        )
        hist2d_Sci = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(iTowerX_max - iTowerX_min), iTowerX_min, iTowerX_max,
            int(iTowerY_max - iTowerY_min), iTowerY_min, iTowerY_max
        )
        hist2d_Cer_3mm = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Cer_3mm",
            f"Event Display {evtNumber};X;Y (Cherenkov)",
            int(iTowerX_max - iTowerX_min), iTowerX_min, iTowerX_max,
            (int(iTowerY_max - iTowerY_min))*4, iTowerY_min, iTowerY_max
        )
        hist2d_Sci_3mm = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci_3mm",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(iTowerX_max - iTowerX_min), iTowerX_min, iTowerX_max,
            (int(iTowerY_max - iTowerY_min))*4, iTowerY_min, iTowerY_max
        )

        # Fill them by reading each branch directly from the tree
                # --- REPLACE your old FERS loop with this ---
        for _, FERSBoard in FERSBoards.items():
            for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
                # get the two channels (Cer & Sci) for this tower
                cer_chan = FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=True)
                sci_chan = FERSBoard.GetChannelByTower(iTowerX, iTowerY, isCer=False)

                # these vectors live in a single branch per board:
                branch_name = f"FERS_Board{cer_chan.boardNo}_energyHG"
                energies = getattr(tree, branch_name)

                # index into the vector by channelNo
                e_Cer = energies[cer_chan.channelNo]
                e_Sci = energies[sci_chan.channelNo]

                # now fill your 2D histos
                if FERSBoard.Is3mm():
                    hist2d_Cer_3mm.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci_3mm.Fill(iTowerX, iTowerY, int(e_Sci))
                else:
                    hist2d_Cer.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci.Fill(iTowerX, iTowerY, int(e_Sci))


        # Prepare the TPaveText overlays
        extraToDraw = ROOT.TPaveText(0.01, 0.73, 0.12, 0.88, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Run: {runNumber}")
        extraToDraw.AddText(f"Event: {evtNumber}")

        extraToDraw_Cer = extraToDraw.Clone("extraToDraw_2")
        extraToDraw_Cer.AddText("Cerenkov")
        extraToDraw_Sci = extraToDraw.Clone("extraToDraw_3")
        extraToDraw_Sci.AddText("Scintillator")

        # Draw and save event displays
        outdir_eventdisplay = f"{outdir}/event_display/"
        DrawHistos(
            [hist2d_Cer, hist2d_Cer_3mm], "", iTowerX_min, iTowerX_max, "iX",
            iTowerY_min, iTowerY_max, "iY", f"event_display_Evt{evtNumber}_Cer",
            dology=False, drawoptions=["COLZ,text", "COLZ,text"],
            zmin=50.0, zmax=3000.0, doth2=True, W_ref=W_ref,
            extraToDraw=extraToDraw_Cer, outdir=outdir_eventdisplay
        )
        DrawHistos(
            [hist2d_Sci, hist2d_Sci_3mm], "", iTowerX_min, iTowerX_max, "iX",
            iTowerY_min, iTowerY_max, "iY", f"event_display_Evt{evtNumber}_Sci",
            dology=False, drawoptions=["COLZ,text", "COLZ,text"],
            zmin=50.0, zmax=8000.0, doth2=True, W_ref=W_ref,
            extraToDraw=extraToDraw_Sci, outdir=outdir_eventdisplay
        )
        plots_eventdisplay += [
            f"event_display_Evt{evtNumber}_Cer.png",
            f"event_display_Evt{evtNumber}_Sci.png"
        ]

        # Now build pulse shapes from DRS boards
        for _, DRSBoard in DRSBoards.items():
            boardNo = DRSBoard.boardNo
            for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
                sTowerX = number2string(iTowerX)
                sTowerY = number2string(iTowerY)

                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer")
                    )
                    if chan is None:
                        print(f"Warning: Channel not found for tower ({iTowerX}, {iTowerY}) "
                              f"on board {boardNo} for {var} channel.")
                        continue

                    varname = chan.GetChannelName()
                    pulse_shape = getattr(tree, varname)

                    h1 = ROOT.TH1F(
                        f"pulse_shape_Evt{evtNumber}_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                        f"Pulse Shape {var} (Board {boardNo}, Tower {sTowerX}, {sTowerY});Time;Amplitude",
                        1024, 0, 1024
                    )
                    for i in range(len(pulse_shape)):
                        h1.Fill(i, pulse_shape[i])
                    hists_pulse_shapes.append(h1)

                # annotate and save
                peak_Cer = hists_pulse_shapes[-2].GetMaximumBin()
                peak_Sci = hists_pulse_shapes[-1].GetMaximumBin()
                deltaTS = peak_Sci - peak_Cer

                extraToDraw2 = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
                extraToDraw2.SetTextAlign(11)
                extraToDraw2.SetFillColorAlpha(0, 0)
                extraToDraw2.SetBorderSize(0)
                extraToDraw2.SetTextFont(42)
                extraToDraw2.SetTextSize(0.04)
                extraToDraw2.AddText(f"Event: {evtNumber}, iTowerX: {iTowerX}, iTowerY: {iTowerY}")
                extraToDraw2.AddText(f"Cer Peak: {peak_Cer}")
                extraToDraw2.AddText(f"Sci Peak: {peak_Sci}")
                extraToDraw2.AddText(f"delta Peak: {deltaTS} ts")
                extraToDraw2.AddText(f"delta T: {deltaTS * 0.2: .2f} ns")

                outdir_pulse_shapes = f"{outdir}/pulse_shapes/"
                DrawHistos(
                    hists_pulse_shapes[-2:], ["Cer", "Sci"], 0, 1024, "TS",
                    1400, 2200, "Amplitude",
                    f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}",
                    dology=False, mycolors=[2, 4], drawashist=True,
                    extraToDraw=extraToDraw2, outdir=outdir_pulse_shapes
                )
                plots_pulse_shapes.append(
                    f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}.png"
                )

    print(f"Events left after filtering: {nevents}")

    # Generate HTML viewers
    generate_html(
        plots_eventdisplay,
        outdir_eventdisplay,
        plots_per_row=2,
        output_html="html/event_display/viewer.html",
        random_per_block=random_per_block,
        block_size=block_size,
        nth_interval=nth_interval,
        group_size=2,
    )
    generate_html(
        plots_pulse_shapes,
        outdir_pulse_shapes,
        output_html="html/pulse_shapes/viewer.html",
        random_per_block=random_per_block,
        block_size=block_size,
        nth_interval=nth_interval,
    )

    # Save histograms out to ROOT files
    ofilename = infilename.replace(".root", "_event_display.root")
    ofile = ROOT.TFile(ofilename, "RECREATE")
    print(f"Saving event display histograms to {ofilename}")
    for hist in hists_eventdisplay:
        hist.SetDirectory(ofile)
        hist.Write()
    ofile.Close()

    ofilename_ps = infilename.replace(".root", "_pulse_shapes.root")
    ofile_ps = ROOT.TFile(ofilename_ps, "RECREATE")
    print(f"Saving pulse shape histograms to {ofilename_ps}")
    for hist in hists_pulse_shapes:
        hist.SetDirectory(ofile_ps)
        hist.Write()
    ofile_ps.Close()

    print(f"Processing completed in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate event displays")
    parser.add_argument("input_file", help="Input ROOT file")
    parser.add_argument("--random-per-block", type=int, default=None,
                        dest="random_per_block",
                        help="Number of random events to display per block")
    parser.add_argument("--block-size", type=int, default=100,
                        dest="block_size",
                        help="Block size used for random sampling")
    parser.add_argument("--nth-interval", type=int, default=None,
                        dest="nth_interval",
                        help="Additionally display every nth event")

    args = parser.parse_args()

    print(f"Processing file: {args.input_file}")
    makeEventDisplays(
        args.input_file,
        random_per_block=args.random_per_block,
        block_size=args.block_size,
        nth_interval=args.nth_interval,
    )