import ROOT
import sys
sys.path.append("CMSPLOTS")
import numpy as np
from utils.html_generator import generate_html
from runconfig import runNumber

root_file = "/Users/elegantuniverse/hodoReco/fers_dummy_data.root"
tree_name = "EventTree"
threshold = 1000

rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"

heatmap_entries = []

def eventProcess():
    global heatmap_entries
    f = ROOT.TFile.Open(root_file, "READ")
    t = f.Get(tree_name)
    
    for i, ev in enumerate(t):
        b1 = getattr(ev, "FERS_Board1_energyHG")
        b2 = getattr(ev, "FERS_Board2_energyHG")
        
        for ix in range(64):
            if b1[ix] < threshold:
                continue
            for iy in range(64):
                if b2[iy] < threshold:
                    continue
                heatmap_entries.append((ix, iy)) 

def buildHeatMapEntries():
    heatmap_matrix = np.zeros((64, 64)) 

    for x, y in heatmap_entries:
        heatmap_matrix[y, x] += 1  

    return heatmap_matrix

def get_mode(hist):
    best_bin, best_count = 1, hist.GetBinContent(1)
    for b in range(2, hist.GetNbinsX() + 1):
        cnt = hist.GetBinContent(b)
        if cnt > best_count:
            best_count, best_bin = cnt, b
    return best_bin - 1

def makeHeatmaps():
    #blank for now

def generateHTML():
    plots = ["hodoHeatmap.png"]
    outdir_plots = f"{plotdir}/Conditions_vs_Event"
    output_html = f"{htmldir}/Conditions_vs_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html

def main():
    eventProcess()
    makeHeatmaps()
    generateHTML()
    
if __name__ == "__main__":
    main()
