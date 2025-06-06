import ROOT
from utils.channel_map import build_map_Cer_Sci, build_map_ixy_DRSVar
import json

print("Start running toyCalibration_DRS.py")
# batch mode
ROOT.gROOT.SetBatch(True)

map_Cer_Sci = build_map_Cer_Sci()
print(map_Cer_Sci)
map_ixy_DRSVar_Cer, map_ixy_DRSVar_Sci = build_map_ixy_DRSVar()


def FindPeakPosition(hist):
    """Find the peak position of a histogram."""
    max_bin = hist.GetMaximumBin()
    return hist.GetBinCenter(max_bin)


infile_name = "root/drs_all_channels.root"
infile = ROOT.TFile(infile_name, "READ")

noises_map = {}

for hist_name in infile.GetListOfKeys():
    # filter out 2d histograms
    hist_name = hist_name.GetName()
    if not hist_name.startswith("hist_DRS_Board"):
        continue
    hist = infile.Get(hist_name)

    noise = FindPeakPosition(hist)
    print(f"DRS Noise for {hist_name}: {noise}")
    branch_name = hist_name.replace("hist_", "")
    noises_map[branch_name] = noise

with open("results/drs_noises.json", "w") as json_file:
    json.dump(noises_map, json_file)
