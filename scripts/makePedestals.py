import json
from utils.parser import get_args
from configs.plotranges import getDRSPlotRanges, getServiceDRSPlotRanges
from utils.dataloader import getRunInfo, loadRDF
from variables.fers import vectorizeFERS, getFERSEnergyMax
from variables.drs import preProcessDRSBoards, getDRSStats
from utils.utils import number2string
from channels.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels, getUpstreamVetoChannel, getDownStreamMuonChannel, getServiceDRSChannels
import ROOT
import os
from utils.timing import auto_timer  # noqa
auto_timer("Total Execution Time")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

debugDRS = False

runNumber = 1259
firstEvent = 10
lastEvent = -1

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)

fersboards = buildFERSBoards(run=runNumber)
rdf = vectorizeFERS(rdf, fersboards)


def collectFERSPedestals():
    pedestals = {}
    for fersboard in fersboards.values():
        for iTowerX, iTowerY in fersboard.GetListOfTowers():
            for var in ["Cer", "Sci"]:
                for gain in ["HG", "LG"]:
                    chan = fersboard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer"))
                    channelName = chan.GetChannelName(gain=gain)

                    pedestals[channelName] = rdf.Mean(channelName)
    return pedestals


pedestals = collectFERSPedestals()
pedestalValues = {channelName: pedestal.GetValue()
                  for channelName, pedestal in pedestals.items()}

# save to json file
output_dir = "data/fers/"
os.makedirs(output_dir, exist_ok=True)
file_pedestals = os.path.join(
    output_dir, f"FERS_pedestals_run{runNumber}.json")
with open(file_pedestals, 'w') as f:
    json.dump(pedestalValues, f, indent=4)
print(f"Pedestals saved to {file_pedestals}")
