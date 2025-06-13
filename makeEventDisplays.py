import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
import os
from utils.channel_map import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile
from utils.html_generator import generate_html
import time


ROOT.gROOT.SetBatch(True)  # Run in batch mode

runNumber = 583
outdir = f"plots/Run{runNumber}/"


def makeEventDisplays(infilename):
    start_time = time.time()
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    DRSBoards = buildDRSBoards(run=runNumber)
    FERSBoards = buildFERSBoards(run=runNumber)

    hists_eventdisplay = []
    hists_pulse_shapes = []
    plots_eventdisplay = []
    plots_pulse_shapes = []

    W_ref = 1000
    iTowerX_max = 11.5
    iTowerX_min = -12.5
    iTowerY_max = 15.5
    iTowerY_min = -8.5

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        if ievt > 30:
            break
        print(f"Processing event {ievt + 1} of {rdf.Count().GetValue()}")
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        # if evtNumber not in events_interested:
        #    print(
        #        f"Skipping event {evtNumber} as it is not in the interested events list.")
        #    continue
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
            int(iTowerY_max - iTowerY_min)*4, iTowerY_min, iTowerY_max
        )
        hist2d_Sci_3mm = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci_3mm",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(iTowerX_max - iTowerX_min), iTowerX_min, iTowerX_max,
            int(iTowerY_max - iTowerY_min)*4, iTowerY_min, iTowerY_max
        )
        for _, FERSBoard in FERSBoards.items():
            for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
                chan_Cer = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=True)
                chan_Sci = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=False)
                var_Cer = chan_Cer.GetHGChannelName()
                var_Sci = chan_Sci.GetHGChannelName()
                e_Cer = rdf.Take["unsigned short"](var_Cer).GetValue()[ievt]
                e_Sci = rdf.Take["unsigned short"](var_Sci).GetValue()[ievt]
                # print(
                #    f"Event {evtNumber}, Tower ({iTowerX}, {iTowerY}): Cerenkov: {e_Cer}, Scintillator: {e_Sci}")

                if FERSBoard.Is3mm():
                    iX, iY = hist2d_Cer_3mm.GetXaxis().FindBin(iTowerX), \
                        hist2d_Cer_3mm.GetYaxis().FindBin(iTowerY)
                    hist2d_Cer_3mm.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci_3mm.Fill(iTowerX, iTowerY, int(e_Sci))
                else:
                    iX, iY = hist2d_Cer.GetXaxis().FindBin(iTowerX), \
                        hist2d_Cer.GetYaxis().FindBin(iTowerY)
                    hist2d_Cer.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci.Fill(iTowerX, iTowerY, int(e_Sci))

        extraToDraw = ROOT.TPaveText(0.01, 0.73, 0.12, 0.88, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Run: {runNumber}")
        extraToDraw.AddText(
            f"Event: {evtNumber}")

        extraToDraw_Cer = extraToDraw.Clone("extraToDraw_2")
        extraToDraw_Cer.AddText("Cerenkov")
        extraToDraw_Sci = extraToDraw.Clone("extraToDraw_3")
        extraToDraw_Sci.AddText("Scintillator")

        if ievt < 100000:
            outdir_eventdisplay = f"{outdir}/event_display/"
            DrawHistos([hist2d_Cer, hist2d_Cer_3mm], f"", iTowerX_min, iTowerX_max, "iX",
                       iTowerY_min, iTowerY_max, "iY", f"event_display_Evt{evtNumber}_Cer",
                       dology=False, drawoptions=["COLZ,text", "COLZ,text"],
                       zmin=50.0, zmax=3000.0, doth2=True, W_ref=W_ref, extraToDraw=extraToDraw_Cer,
                       outdir=outdir_eventdisplay)
            DrawHistos([hist2d_Sci, hist2d_Sci_3mm], f"", iTowerX_min, iTowerX_max, "iX",
                       iTowerY_min, iTowerY_max, "iY", f"event_display_Evt{evtNumber}_Sci",
                       dology=False, drawoptions=["COLZ,text", "COLZ,text"],
                       zmin=50.0, zmax=8000.0, doth2=True, W_ref=W_ref, extraToDraw=extraToDraw_Sci,
                       outdir=outdir_eventdisplay)
            plots_eventdisplay.append(
                f"event_display_Evt{evtNumber}_Cer.png")
            plots_eventdisplay.append(
                f"event_display_Evt{evtNumber}_Sci.png")

        # make pulses from DRS
        for _, DRSBoard in DRSBoards.items():
            boardNo = DRSBoard.boardNo
            for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
                sTowerX = number2string(iTowerX)
                sTowerY = number2string(iTowerY)

                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer"))
                    varname = chan.GetChannelName()
                    pulse_shape = rdf.Take["ROOT::VecOps::RVec<float>"](
                        varname).GetValue()[ievt]

                    h1 = ROOT.TH1F(
                        f"pulse_shape_Evt{evtNumber}_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                        f"Pulse Shape {var} (Board {boardNo}, Tower {sTowerX}, {sTowerY});Time;Amplitude",
                        1024, 0, 1024
                    )
                    for i in range(len(pulse_shape)):
                        h1.Fill(i, pulse_shape[i])
                    hists_pulse_shapes.append(h1)

                if ievt < 1000000:
                    peak_Cer = hists_pulse_shapes[-2].GetMaximumBin()
                    peak_Sci = hists_pulse_shapes[-1].GetMaximumBin()
                    deltaTS = peak_Sci - peak_Cer

                    extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
                    extraToDraw.SetTextAlign(11)
                    extraToDraw.SetFillColorAlpha(0, 0)
                    extraToDraw.SetBorderSize(0)
                    extraToDraw.SetTextFont(42)
                    extraToDraw.SetTextSize(0.04)
                    extraToDraw.AddText(
                        f"Event: {evtNumber}, iTowerX: {iTowerX}, iTowerY: {iTowerY}")
                    extraToDraw.AddText(f"Cer Peak: {peak_Cer}")
                    extraToDraw.AddText(f"Sci Peak: {peak_Sci}")
                    extraToDraw.AddText(
                        f"delta Peak: {deltaTS} ts")
                    extraToDraw.AddText(f"delta T: {deltaTS * 0.2: .2f} ns")

                    outdir_pulse_shapes = f"{outdir}/pulse_shapes/"
                    DrawHistos(hists_pulse_shapes[-2:], ["Cer", "Sci"], 0, 1024, "TS",
                               1400, 2200, "Amplitude",
                               f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}",
                               dology=False, mycolors=[2, 4], drawashist=True, extraToDraw=extraToDraw,
                               outdir=outdir_pulse_shapes)
                    plots_pulse_shapes.append(
                        f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}.png")

    print(f"Events left after filtering: {rdf.Count().GetValue()}")

    generate_html(plots_eventdisplay, outdir_eventdisplay, plots_per_row=2,
                  output_html=f"html/event_display/viewer.html")
    generate_html(plots_pulse_shapes, outdir_pulse_shapes,
                  output_html=f"html/pulse_shapes/viewer.html")

    # Save event display histograms
    ofilename = infilename.replace(".root", "_event_display.root")
    ofile = ROOT.TFile(ofilename, "RECREATE")
    print(f"Saving event display histograms to {ofilename}")
    for hist in hists_eventdisplay:
        hist.SetDirectory(ofile)
        hist.Write()
    ofile.Close()

    # Save pulse shape histograms
    ofilename_pulse_shapes = infilename.replace(".root", "_pulse_shapes.root")
    ofile = ROOT.TFile(ofilename_pulse_shapes, "RECREATE")
    print(f"Saving pulse shape histograms to {ofilename_pulse_shapes}")
    for hist in hists_pulse_shapes:
        hist.SetDirectory(ofile)
        hist.Write()
    ofile.Close()

    print(f"Processing completed in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    input_file = f"root/Run{runNumber}/filtered.root"
    print(f"Processing file: {input_file}")
    makeEventDisplays(input_file)
