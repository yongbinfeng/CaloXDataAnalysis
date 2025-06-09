import ROOT
import json

print("Start running toyCalibration.py")
# batch mode
ROOT.gROOT.SetBatch(True)


def FindPeakPosition(hist):
    """Find the peak position of a histogram."""
    max_bin = hist.GetMaximumBin()
    return hist.GetBinCenter(max_bin)


def FindPeakPosition_FERS(hist, xmin):
    """Find the peak position of a histogram."""
    val_max = hist.GetMinimum()
    max_bin = -1
    for bin in range(hist.FindBin(xmin), hist.GetNbinsX() + 1):
        if hist.GetBinContent(bin) > val_max:
            val_max = hist.GetBinContent(bin)
            max_bin = bin
    return hist.GetBinCenter(max_bin) if max_bin != -1 else None


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

infile_name = "root/hodoscope_all_channels.root"
infile = ROOT.TFile(infile_name, "READ")
noises_map = {}
for hist_name in infile.GetListOfKeys():
    # filter out 2d histograms
    hist_name = hist_name.GetName()
    if not hist_name.startswith("hist_hodoscope_DRS"):
        continue
    hist = infile.Get(hist_name)

    noise = FindPeakPosition(hist)
    print(f"Hodoscope Noise for {hist_name}: {noise}")
    branch_name = hist_name.replace("hist_", "")
    noises_map[branch_name] = noise
with open("results/hodoscope_noises.json", "w") as json_file:
    json.dump(noises_map, json_file)

# FRES noise
infile_name = "root/fers_all_channels.root"
infile = ROOT.TFile(infile_name, "READ")
noises_map = {}
for hist_name in infile.GetListOfKeys():
    # filter out 2d histograms
    hist_name = hist_name.GetName()
    if not hist_name.startswith("hist_board"):
        continue
    if "_vs_" in hist_name:  # skip 2D histograms
        continue
    hist = infile.Get(hist_name)

    # Sometimes FERS seems to have a spike around 140 - 160, remove that
    noise = FindPeakPosition_FERS(hist, 180)
    print(f"FERS Noise for {hist_name}: {noise}")
    branch_name = hist_name.replace("hist_", "")
    noises_map[branch_name] = noise
with open("results/fers_noises.json", "w") as json_file:
    json.dump(noises_map, json_file)


print("Finished running toyCalibration.py")
