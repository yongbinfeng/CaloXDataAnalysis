import os
import sys
import ROOT
import json
from collections import OrderedDict
from channels.channel_map import build_fers_boards, build_drs_boards, get_mcp_channels
from utils.data_loader import loadRDF, getRunInfo
from variables.drs import preProcessDRSBoards, calibrateDRSPeakTS
from utils.html_generator import generate_html
from selections.selections import vetoMuonCounter, applyUpstreamVeto, applyPSDSelection, applyCC1Selection
from utils.parser import get_args
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos, LHistos2Hist
from utils.timing import auto_timer
from utils.utils import number_to_string

# runs = [1416, 1412, 1410]
runs = [1442, 1441, 1439, 1437]

channels = {"m1p5_1p125": "Quartz",
            "m1p5_1p375": "Plastic",
            "m0p5_1p875": "Quartz",
            "m0p5_1p625": "Plastic"
            }

bins = {
    "Quartz": [-70, -55],
    "Plastic": [-70, -50],
}

ROOT.gROOT.SetBatch(True)

hprofs = {}
legends = []
for run_number in runs:
    btype, benergy = getRunInfo(run_number)
    legends.append(f"{benergy} GeV")
    hprofs[run_number] = {}
    rootfile = f"results2/root/Run{run_number}/drs_vs_ts_calibrated.root"
    infile = ROOT.TFile.Open(rootfile, "READ")

    for channel in channels.keys():
        hprof = infile.Get(f"prof_DRS_vs_TS_Cer_{channel}")

        hist = hprof.ProjectionX().Clone(
            f"prof_DRS_vs_TS_Cer_{channel}_Run{run_number}")
        hist.SetDirectory(0)
        # normalize
        binmin = hist.FindBin(-70)
        binmax = hist.FindBin(-55)
        hist.Scale(1.0 / hist.Integral(binmin, binmax))
        hprofs[run_number][channel] = hist

for channel, fibertype in channels.items():
    output_name = f"prof_DRS_vs_TS_Cer_{channel}_Combined"
    DrawHistos([hprofs[run_number][channel] for run_number in runs], legends, -70, -55, "Calibrated TS", 0, None, "DRS ADC",
               output_name,
               dology=False, drawoptions="HIST", mycolors=[ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2], addOverflow=False, addUnderflow=False,
               outdir="results2/plots/CombinedRuns", legendPos=[0.30, 0.85, 0.90, 0.90], legendNCols=3, extra_text=f"{fibertype}", lumitext=btype)
