import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
from utils.channel_map import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile, processDRSBoards
from utils.html_generator import generate_html
from runNumber import runNumber
import time


ROOT.gROOT.SetBatch(True)  # Run in batch mode

outdir = f"plots/Run{runNumber}/"

xmax = 14
xmin = -14
ymax = 10
ymin = -10
W_ref = 1000
H_ref = 1100


def makeEventDisplays(infilename):
    start_time = time.time()
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    DRSBoards = buildDRSBoards(run=runNumber)
    FERSBoards = buildFERSBoards(run=runNumber)

    rdf = processDRSBoards(rdf, DRSBoards)

    hists_eventdisplay = []
    hists_pulse_shapes = {}
    plots_eventdisplay = []
    plots_pulse_shapes = []

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        if ievt > 20:
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
            int(xmax - xmin), xmin, xmax,
            int(ymax - ymin), ymin, ymax
        )
        hist2d_Sci = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(xmax - xmin), xmin, xmax,
            int(ymax - ymin), ymin, ymax
        )
        hist2d_Cer_3mm = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Cer_3mm",
            f"Event Display {evtNumber};X;Y (Cherenkov)",
            int(xmax - xmin), xmin, xmax,
            int(ymax - ymin) * 4, ymin, ymax
        )
        hist2d_Sci_3mm = ROOT.TH2F(
            f"event_display_Evt{evtNumber}_Sci_3mm",
            f"Event Display {evtNumber};X;Y (Scintillator)",
            int(xmax - xmin), xmin, xmax,
            int(ymax - ymin) * 4, ymin, ymax
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
                    hist2d_Cer_3mm.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci_3mm.Fill(iTowerX, iTowerY, int(e_Sci))
                else:
                    hist2d_Cer.Fill(iTowerX, iTowerY, int(e_Cer))
                    hist2d_Sci.Fill(iTowerX, iTowerY, int(e_Sci))

        extraToDraw = ROOT.TPaveText(0.01, 0.73, 0.12, 0.88, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(
            f"Event: {evtNumber}")

        hist2d_Cer.SetMarkerSize(0.55)
        hist2d_Sci.SetMarkerSize(0.55)
        hist2d_Cer_3mm.SetMarkerSize(0.55)
        hist2d_Sci_3mm.SetMarkerSize(0.55)

        if ievt < 100000:
            outdir_eventdisplay = f"{outdir}/event_display/"
            DrawHistos([hist2d_Cer, hist2d_Cer_3mm], f"", xmin, xmax, "iX",
                       ymin, ymax, "iY", f"event_display_Evt{evtNumber}_Cer",
                       dology=False, drawoptions=["COL,text", "COL,text"],
                       zmin=0.0, zmax=3000.0, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=extraToDraw, extraText="Cer", runNumber=runNumber,
                       outdir=outdir_eventdisplay)
            DrawHistos([hist2d_Sci, hist2d_Sci_3mm], f"", xmin, xmax, "iX",
                       ymin, ymax, "iY", f"event_display_Evt{evtNumber}_Sci",
                       dology=False, drawoptions=["COL,text", "COL,text"],
                       zmin=0.0, zmax=8000.0, doth2=True, W_ref=W_ref, extraToDraw=extraToDraw, H_ref=H_ref, extraText="Sci", runNumber=runNumber,
                       outdir=outdir_eventdisplay)
            plots_eventdisplay.append(
                f"event_display_Evt{evtNumber}_Cer.png")
            plots_eventdisplay.append(
                f"event_display_Evt{evtNumber}_Sci.png")

        # make pulses from DRS
        for _, DRSBoard in DRSBoards.items():
            for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
                sTowerX = number2string(iTowerX)
                sTowerY = number2string(iTowerY)

                if (sTowerX, sTowerY) not in hists_pulse_shapes:
                    hists_pulse_shapes[(sTowerX, sTowerY)] = {}
                    hists_pulse_shapes[(sTowerX, sTowerY)]['Cer'] = None
                    hists_pulse_shapes[(sTowerX, sTowerY)]['Sci'] = None

                for var in ["Cer", "Sci"]:
                    chan = DRSBoard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer"))
                    if chan is None:
                        continue
                    varname = chan.GetChannelName() + "_subtractMedian"
                    pulse_shape = rdf.Take["ROOT::VecOps::RVec<float>"](
                        varname).GetValue()[ievt]

                    h1 = ROOT.TH1F(
                        f"pulse_shape_Evt{evtNumber}_{var}_{sTowerX}_{sTowerY}",
                        f"Pulse Shape {var} (Tower {sTowerX}, {sTowerY});Time;Amplitude",
                        1024, 0, 1024
                    )
                    for i in range(len(pulse_shape)):
                        h1.Fill(i, pulse_shape[i])
                    hists_pulse_shapes[sTowerX, sTowerY][var] = h1

        for (sTowerX, sTowerY) in hists_pulse_shapes:
            hist_pulse_shape_cer = hists_pulse_shapes[sTowerX, sTowerY]['Cer']
            hist_pulse_shape_sci = hists_pulse_shapes[sTowerX, sTowerY]['Sci']

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"Event: {evtNumber}, iTowerX: {iTowerX}, iTowerY: {iTowerY}")

            to_draws = []
            labels = []
            colors = []
            if hist_pulse_shape_cer:
                peak_Cer = hist_pulse_shape_cer.GetMaximumBin()
                extraToDraw.AddText(f"Cer Peak: {peak_Cer}")
                to_draws.append(hist_pulse_shape_cer)
                labels.append("Cer")
                colors.append(2)  # Red color for Cer
            if hist_pulse_shape_sci:
                peak_Sci = hist_pulse_shape_sci.GetMaximumBin()
                extraToDraw.AddText(f"Sci Peak: {peak_Sci}")
                to_draws.append(hist_pulse_shape_sci)
                labels.append("Sci")
                colors.append(4)  # Blue color for Sci
            if hist_pulse_shape_cer and hist_pulse_shape_sci:
                deltaTS = peak_Sci - peak_Cer
                extraToDraw.AddText(
                    f"delta Peak: {deltaTS} ts")
                extraToDraw.AddText(f"delta T: {deltaTS * 0.2: .2f} ns")

            outdir_pulse_shapes = f"{outdir}/pulse_shapes/"
            DrawHistos(to_draws, labels, 0, 1024, "TS",
                       -5, 30, "Amplitude",
                       f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}",
                       dology=False, mycolors=colors, drawashist=True, extraToDraw=extraToDraw,
                       outdir=outdir_pulse_shapes, runNumber=runNumber)
            plots_pulse_shapes.append(
                f"pulse_shape_Evt{evtNumber}_iTowerX{sTowerX}_iTowerY{sTowerY}.png")

    print(f"Events left after filtering: {rdf.Count().GetValue()}")

    generate_html(plots_eventdisplay, outdir_eventdisplay, plots_per_row=2,
                  output_html=f"html/Run{runNumber}/event_display/viewer.html")
    generate_html(plots_pulse_shapes, outdir_pulse_shapes,
                  output_html=f"html/Run{runNumber}/pulse_shapes/viewer.html")

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
    for hists in hists_pulse_shapes.values():
        for hist in hists.values():
            if hist is not None:
                hist.SetDirectory(ofile)
                hist.Write()
    ofile.Close()

    print(f"Processing completed in {time.time() - start_time:.2f} seconds.")


if __name__ == "__main__":
    input_file = f"root/Run{runNumber}/filtered.root"
    print(f"Processing file: {input_file}")
    makeEventDisplays(input_file)
