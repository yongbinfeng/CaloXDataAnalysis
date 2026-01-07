import json
from utils.parser import get_args
from configs.plot_ranges import get_drs_plot_ranges, get_service_drs_plot_ranges
from utils.data_loader import getRunInfo, loadRDF
from variables.fers import vectorizeFERS, get_fers_energy_max
from variables.drs import preProcessDRSBoards, get_drs_stats
from utils.utils import number_to_string
from channels.channel_map import build_drs_boards, build_fers_boards, build_time_reference_channels, build_hodo_trigger_channels, build_hodo_pos_channels, get_upstream_veto_channel, get_downstream_muon_channel, get_service_drs_channels
import ROOT
import os
from utils.timing import auto_timer  # noqa
auto_timer("Total Execution Time")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

debug_drs = False

run_number = 1259
firstEvent = 10
lastEvent = -1

rdf, rdf_org = loadRDF(run_number, firstEvent, lastEvent)

fersboards = build_fers_boards(run=run_number)
rdf = vectorizeFERS(rdf, fersboards)


def collectFERSPedestals():
    pedestals = {}
    for fersboard in fersboards.values():
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            for var in ["Cer", "Sci"]:
                for gain in ["HG", "LG"]:
                    chan = fersboard.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    channelName = chan.get_channel_name(gain=gain)

                    pedestals[channelName] = rdf.Mean(channelName)
    return pedestals


pedestals = collectFERSPedestals()
pedestalValues = {channelName: pedestal.GetValue()
                  for channelName, pedestal in pedestals.items()}

# save to json file
output_dir = "data/fers/"
os.makedirs(output_dir, exist_ok=True)
file_pedestals = os.path.join(
    output_dir, f"FERS_pedestals_run{run_number}.json")
with open(file_pedestals, 'w') as f:
    json.dump(pedestalValues, f, indent=4)
print(f"Pedestals saved to {file_pedestals}")
