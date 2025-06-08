import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from CMSPLOTS.myFunction import DrawHistos
import json
from results.events import events_interested
from utils.channel_map import get_hodoscope_channels

ROOT.gROOT.SetBatch(True)  # Run in batch mode

hodoscope_channels = get_hodoscope_channels()


def analyzeHodoPulse(infilename):
    infile = ROOT.TFile(infilename, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Failed to open input file: {infile}")

    # Create an RDataFrame from the EventTree
    rdf = ROOT.RDataFrame("EventTree", infile)

    # print how many events are left after filtering
    for ievt in range(0, rdf.Count().GetValue()):
        print(f"Processing event {ievt + 1} of {rdf.Count().GetValue()}")
        evtNumber = rdf.Take["unsigned int"]("event_n").GetValue()[ievt]
        if ievt > 10:
            break
        if evtNumber not in events_interested:
            print(
                f"Skipping event {evtNumber} as it is not in the interested events list.")
            continue
        # dump the pulse shapes for hodoscope channels
        h1_hodos = {}
        for group, channels in hodoscope_channels.items():
            hs_todraw = []
            for channel in channels:
                pulse_shape = rdf.Take["ROOT::VecOps::RVec<float>"](
                    channel).GetValue()[ievt]
                h1_hodos[channel] = ROOT.TH1F(
                    f"pulse_shape_{channel}_Evt{evtNumber}",
                    f"Pulse Shape {channel};Time;Amplitude",
                    1024, 0, 1024
                )
                for i in range(len(pulse_shape)):
                    h1_hodos[channel].Fill(i, pulse_shape[i])
                hs_todraw.append(h1_hodos[channel])

            if ievt < 100000:
                if group == "trigger":
                    labels = ["Trigger"]
                else:
                    labels = ["Left", "Right"]
                DrawHistos(hs_todraw, labels, 0, 1024, "TS",
                           -200, 3000, "Amplitude", f"pulse_shape_Hodoscopes_Evt{evtNumber}_{group}", dology=False, mycolors=[2, 4], drawashist=True)
    print(f"Events left after filtering: {rdf.Count().GetValue()}")


if __name__ == "__main__":
    files = [
        "root/filtered_events_board1.root",
    ]
    for idx, input_file in enumerate(files):
        print(f"Processing file: {input_file}")
        analyzeHodoPulse(input_file)
