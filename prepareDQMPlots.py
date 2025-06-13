import os
import ROOT
from utils.channel_map_new import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile
import time

print("Start running prepareDQMPlots.py")

runNumber = 583

# multi-threading support
ROOT.ROOT.EnableImplicitMT(5)

# Open the input ROOT file
ifile = getDataFile(runNumber)

suffix = f"run{runNumber}"
infile = ROOT.TFile(ifile, "READ")
rdf = ROOT.RDataFrame("EventTree", infile)

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
                    3000, 0, 9000),
                    chan.GetHGChannelName()
                )
                hists1d_FERS.append(hist)

    return hists1d_FERS


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
                500, 0, 9000, 500, 0, 9000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hist_zoomed = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_vs_Sci_{sTowerX}_{sTowerY}_zoom",
                f"CER {iCer} vs SCI {iSci} in iTowerX {sTowerX} iTowerY {sTowerY} (zoomed);CER Energy HG;SCI Energy HG",
                500, 0, 1000, 1000, 0, 2000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hists2d_FERS.append(hist)
            hists2d_FERS.append(hist_zoomed)

            # high gain vs low gain for Sci
            hist_sci_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"SCI {iSci} HG vs LG;SCI Energy HG;SCI Energy LG",
                500, 0, 9000, 500, 0, 3000),
                chan_Sci.GetHGChannelName(),
                chan_Sci.GetLGChannelName()
            )
            hists2d_FERS.append(hist_sci_hg_vs_lg)
            hist_cer_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"CER {iCer} HG vs LG;CER Energy HG;CER Energy LG",
                500, 0, 9000, 500, 0, 3000),
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

                hist = rdf.Histo1D((
                    f"hist_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Board {boardNo} - {var} iTowerX {sTowerX} iTowerY {sTowerY};{var} Variable;Counts",
                    1500, 1000, 2500),
                    chan.GetChannelName()
                )
                hists1d_DRS.append(hist)
    return hists1d_DRS


def makeDRS2DPlots():
    hists2d_DRS = []
    for _, DRSBoard in DRSBoards.items():
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            chan_Cer = DRSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            chan_Sci = DRSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)

            hist = rdf.Histo2D((
                f"hist_DRS_Board{boardNo}_Cer_vs_Sci_{sTowerX}_{sTowerY}",
                f"CER {chan_Cer.channelNo} vs SCI {chan_Sci.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};CER Variable;SCI Variable",
                500, 1000, 2500, 500, 1000, 2500),
                chan_Cer.GetChannelName(),
                chan_Sci.GetChannelName()
            )
            hists2d_DRS.append(hist)
    return hists2d_DRS


if __name__ == "__main__":
    start_time = time.time()

    hists1d_FERS = makeFERS1DPlots()
    hists2d_FERS = makeFERS2DPlots()
    hists1d_DRS = makeDRS1DPlots()
    hists2d_DRS = makeDRS2DPlots()

    print("Save histograms")

    rootdir = f"root/Run{runNumber}"
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)

    # Save histograms to an output ROOT file
    outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_1D.root", "RECREATE")
    for hist in hists1d_FERS:
        hist.Write()
    outfile.Close()
    outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_2D.root", "RECREATE")
    for hist in hists2d_FERS:
        hist.Write()
    outfile.Close()
    #
    outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_1D.root", "RECREATE")
    for hist in hists1d_DRS:
        hist.SetDirectory(outfile_DRS)
        hist.Write()
    outfile_DRS.Close()
    outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_2D.root", "RECREATE")
    for hist in hists2d_DRS:
        hist.SetDirectory(outfile_DRS)
        hist.Write()
    outfile_DRS.Close()

    time_taken = time.time() - start_time
    print(f"Finished running script in {time_taken:.2f} seconds")
