import re
import sys
import os
import ROOT
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTriggerChannels
from utils.utils import number2string, getDataFile, processDRSBoards
from runNumber import runNumber
import time

print("Start running prepareDQMPlots.py")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)

# Open the input ROOT file
ifile = getDataFile(runNumber)

suffix = f"run{runNumber}"
infile = ROOT.TFile(ifile, "READ")
rdf = ROOT.RDataFrame("EventTree", infile)

# Get total number of entries
n_entries = rdf.Count().GetValue()
nEvents = int(n_entries)
print(f"Total number of events: {nEvents} in run {runNumber}")

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)

# FRES board outputs
# define variables as RDF does not support reading vectors
# with indices directly
for _, FERSBoard in FERSBoards.items():
    boardNo = FERSBoard.boardNo
    for channel in FERSBoard:
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_energyHG_{channel.channelNo}",
            f"FERS_Board{boardNo}_energyHG[{channel.channelNo}]")
        rdf = rdf.Define(
            f"FERS_Board{boardNo}_energyLG_{channel.channelNo}",
            f"FERS_Board{boardNo}_energyLG[{channel.channelNo}]"
        )

rdf = processDRSBoards(rdf, DRSBoards)


def makeFERS1DPlots():
    hists1d_FERS = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist = rdf.Histo1D((
                    f"hist_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS Board {boardNo} - {var} iTowerX {sTowerX} iTowerY {sTowerY};{var} Energy HG;Counts",
                    1000, 0, 9000),
                    chan.GetHGChannelName()
                )
                hists1d_FERS.append(hist)

    return hists1d_FERS


def trackFERSPlots():
    hists2d_FERS_vs_Event = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist = rdf.Histo2D((
                    f"hist_FERS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}",
                    f"FERS Board {boardNo} - Event vs {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};Event;{var} Energy HG",
                    100, 0, nEvents, 1000, 0, 9000),
                    "event_n", chan.GetHGChannelName()
                )
                hists2d_FERS_vs_Event.append(hist)
    return hists2d_FERS_vs_Event


def makeFERS2DPlots():
    hists2d_FERS = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            chan_Cer = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            chan_Sci = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)

            iCer = chan_Cer.channelNo
            iSci = chan_Sci.channelNo
            # high gain
            hist = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_vs_Sci_{sTowerX}_{sTowerY}",
                f"CER {iCer} vs SCI {iSci} in iTowerX {sTowerX} iTowerY {sTowerY};CER Energy HG;SCI Energy HG",
                300, 0, 9000, 300, 0, 9000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hist_zoomed = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_vs_Sci_{sTowerX}_{sTowerY}_zoom",
                f"CER {iCer} vs SCI {iSci} in iTowerX {sTowerX} iTowerY {sTowerY} (zoomed);CER Energy HG;SCI Energy HG",
                300, 0, 1000, 200, 0, 2000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hists2d_FERS.append(hist)
            hists2d_FERS.append(hist_zoomed)

            # high gain vs low gain for Sci
            hist_sci_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"SCI {iSci} HG vs LG;SCI Energy HG;SCI Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_Sci.GetHGChannelName(),
                chan_Sci.GetLGChannelName()
            )
            hists2d_FERS.append(hist_sci_hg_vs_lg)
            hist_cer_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"CER {iCer} HG vs LG;CER Energy HG;CER Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_Cer.GetHGChannelName(),
                chan_Cer.GetLGChannelName()
            )
            hists2d_FERS.append(hist_cer_hg_vs_lg)
    return hists2d_FERS


def makeDRS1DPlots():
    hists1d_DRS = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))

                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                value_mean = stats[channelName]['mean']
                hist = rdf.Histo1D((
                    f"hist_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Board {boardNo} - {var} iTowerX {sTowerX} iTowerY {sTowerY};{var} Variable;Counts",
                    200, value_mean - 100, value_mean + 100),
                    channelName
                )
                hists1d_DRS.append(hist)
    return hists1d_DRS


def makeDRS2DPlots():
    hists2d_DRS_vs_TS = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                # mean_value = stats[channelName]['mean']
                # hist = rdf.Histo2D((
                #    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}",
                #    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};TS;{var} Variable",
                #    1024, 0, 1024, 200, mean_value - 100, mean_value + 100),
                #    "TS", chan.GetChannelName()
                # )
                hist_subtractMedian = rdf.Histo2D((
                    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}_subtractMedian",
                    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);TS;{var} Variable",
                    1024, 0, 1024, 400, -100, 300),
                    "TS", channelName + "_subtractMedian"
                )
                # hists2d_DRS_vs_TS.append(hist)
                hists2d_DRS_vs_TS.append(hist_subtractMedian)
    return hists2d_DRS_vs_TS


def trackDRSPlots():
    hists2d_DRS_vs_Event = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                mean_value = stats[channelName]['mean']
                hist = rdf.Histo2D((
                    f"hist_DRS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}",
                    f"DRS Board {boardNo} Mean - Event vs {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};Event;{var} Variable",
                    int(nEvents/100), 0, nEvents, 200, mean_value - 100, mean_value + 100),
                    "event_n", channelName + "_mean"
                )
                hists2d_DRS_vs_Event.append(hist)
    return hists2d_DRS_vs_Event


def compareDRSChannels(channels_to_compare):
    """
    Compare trigger channels with DRS channels.
    """
    hists_trigger = []
    for chan_name in channels_to_compare:
        hist = rdf.Histo2D((
            f"hist_{chan_name}",
            f"{chan_name};TS;DRS values",
            1024, 0, 1024,
            200, 500, 2500),
            "TS", chan_name
        )
        hist_subtractMedian = rdf.Histo2D((
            f"hist_{chan_name}_subtractMedian",
            f"{chan_name} (subtract median);TS;DRS values",
            1024, 0, 1024,
            200, -1500, 500),
            "TS", chan_name + "_subtractMedian"
        )
        hists_trigger.append(hist)
        hists_trigger.append(hist_subtractMedian)
    return hists_trigger


if __name__ == "__main__":
    start_time = time.time()

    hists1d_FERS = makeFERS1DPlots()
    # hists2d_FERS = makeFERS2DPlots()
    # hists2d_FERS_vs_Event = trackFERSPlots()

    # hists1d_DRS = makeDRS1DPlots()
    hists2d_DRS_vs_TS = makeDRS2DPlots()
    # hists2d_DRS_vs_Event = trackDRSPlots()

    trigger_channels = buildTriggerChannels(run=runNumber)
    hists2d_trigger = compareDRSChannels(trigger_channels)

    print("Save histograms")

    rootdir = f"root/Run{runNumber}"
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)

    # Save histograms to an output ROOT file
    outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_1D.root", "RECREATE")
    for hist in hists1d_FERS:
        hist.Write()
    outfile.Close()
    # outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_2D.root", "RECREATE")
    # for hist in hists2d_FERS:
    #    hist.Write()
    # outfile.Close()
    # outfile = ROOT.TFile(
    #    f"{rootdir}/fers_all_channels_2D_vs_event.root", "RECREATE")
    # for hist in hists2d_FERS_vs_Event:
    #    hist.Write()
    # outfile.Close()
    #
    # outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_1D.root", "RECREATE")
    # for hist in hists1d_DRS:
    #    hist.SetDirectory(outfile_DRS)
    #    hist.Write()
    # outfile_DRS.Close()
    outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_2D.root", "RECREATE")
    for hist in hists2d_DRS_vs_TS:
        hist.SetDirectory(outfile_DRS)
        hist.Write()
    outfile_DRS.Close()
    # outfile_DRS = ROOT.TFile(
    #    f"{rootdir}/drs_all_channels_2D_vs_event.root", "RECREATE")
    # for hist in hists2d_DRS_vs_Event:
    #    hist.SetDirectory(outfile_DRS)
    #    hist.Write()
    # outfile_DRS.Close()

    outfile_trigger = ROOT.TFile(
        f"{rootdir}/trigger_channels.root", "RECREATE")
    for hist in hists2d_trigger:
        hist.SetDirectory(outfile_trigger)
        hist.Write()
    outfile_trigger.Close()

    time_taken = time.time() - start_time
    print(f"Finished running script in {time_taken:.2f} seconds")
