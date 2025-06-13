import os
import ROOT
from utils.channel_map_new import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile
import time

start_time = time.time()

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


requirements = ""
for _, FERSBoard in FERSBoards.items():
    if FERSBoard.Is6mm():
        continue
    # require on the 3mm boards
    for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
        chan_Sci = FERSBoard.GetChannelByTower(
            iTowerX, iTowerY, isCer=False)
        chan_Cer = FERSBoard.GetChannelByTower(
            iTowerX, iTowerY, isCer=True)
        var_Sci = chan_Sci.GetHGChannelName()
        var_Cer = chan_Cer.GetHGChannelName()

        requirements += f"({var_Sci} > 1000 && {var_Cer} > 1000) || "

# remove the last ' || '
requirements = requirements[:-4]
requirements = f"({requirements})"
print(f"Requirements: {requirements}")

rdf_filtered = rdf.Filter(
    requirements, "Filter FERS boards with energy > 1000")

print(f"Number of events after filtering: {rdf_filtered.Count().GetValue()}")
# snapshot the filtered RDF
outfile_name = f"root/Run{runNumber}/filtered.root"
if not os.path.exists(os.path.dirname(outfile_name)):
    os.makedirs(os.path.dirname(outfile_name))
rdf_filtered.Snapshot("EventTree", outfile_name)
print(f"Filtered data saved to {outfile_name}")

print(f"Time taken: {time.time() - start_time:.2f} seconds")
