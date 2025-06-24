import sys
import os
import ROOT
from utils.channel_map import buildDRSBoards, buildFERSBoards
from utils.utils import number2string, getDataFile, getBranchStats
import time
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.html_generator import generate_html
from runNumber import runNumber

# multi-threading support
ROOT.gROOT.SetBatch(True)  # Disable interactive mode
ROOT.ROOT.EnableImplicitMT(5)

suffix = f"run{runNumber}"

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)

FERS_min = 100
FERS_max = 9e3
FERS_LG_max = 2e3
DRS_min = -100
DRS_max = 1e3
DRS_LG_max = 2e3


def prepareFERSDRSPlots():
    ifile = getDataFile(runNumber)
    infile = ROOT.TFile(ifile, "READ")
    rdf = ROOT.RDataFrame("EventTree", infile)

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

    # get the mean of DRS outputs per channel
    ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<float> clipToZero(const ROOT::VecOps::RVec<float>& vec) {
        ROOT::VecOps::RVec<float> out;
        for (float v : vec) {
            if (v < 5.0f) v = 0.0f;  // clip to zero if below threshold
            out.push_back(v);
        }
        return out;
    }
    """)
    ROOT.gInterpreter.Declare("""
    #include "ROOT/RVec.hxx"
    #include <algorithm>
    
    float compute_median(ROOT::RVec<float> vec) {
        if (vec.empty()) return -9999;
        std::sort(vec.begin(), vec.end());
        size_t n = vec.size();
        if (n % 2 == 0)
            return 0.5 * (vec[n / 2 - 1] + vec[n / 2]);
        else
            return vec[n / 2];
    }
    """)
    ROOT.gInterpreter.Declare("""
float SumRange(const ROOT::VecOps::RVec<float>& v, size_t i, size_t j) {
    if (i >= v.size() || j > v.size() || i >= j) return 0.0;
    return std::accumulate(v.begin() + i, v.begin() + j, 0.0f);
}
""")
    ROOT.gInterpreter.Declare("""
                              float MaxRange(const ROOT::VecOps::RVec<float>& v, size_t i, size_t j) {
    if (i >= v.size() || j > v.size() || i >= j) return 0.0;
    float maxVal = v[i];
    for (size_t k = i + 1; k < j; ++k)
        if (v[k] > maxVal) maxVal = v[k];
    return maxVal;
}
""")
    # get the mean of DRS outputs per channel
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for channel in DRSBoard:
            varname = channel.GetChannelName()
            rdf = rdf.Define(
                f"{varname}_median",
                f"compute_median({varname})"
            )
            rdf = rdf.Define(
                f"{varname}_subtractMedian",
                f"{varname} - {varname}_median"
            )
            rdf = rdf.Define(
                f"{varname}_subtractMedian_positive",
                f"clipToZero({varname}_subtractMedian)"
            )
            rdf = rdf.Define(
                f"{varname}_sum",
                f"SumRange({varname}_subtractMedian_positive, 0, 200)"
            )

    # correlate  FERS and DRS outputs
    h2s_FERS_VS_DRS = []
    h2s_FERSLG_VS_DRS = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue
                chan_FERS = None
                for fersNo, FERSBoard in FERSBoards.items():
                    chan_FERS = FERSBoard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer"))
                    if chan_FERS is not None:
                        break
                if chan_FERS is None:
                    print(
                        f"Warning: FERS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                h2_FERS_VS_DRS = rdf.Histo2D((
                    f"hist_FERS_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS vs DRS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, DRS_min, DRS_max, 100, FERS_min, FERS_max
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetHGChannelName()
                )
                h2s_FERS_VS_DRS.append(h2_FERS_VS_DRS)

                h2_FERSLG_VS_DRS = rdf.Histo2D((
                    f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS LG vs DRS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, DRS_min, DRS_LG_max, 100, FERS_min, FERS_LG_max
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetLGChannelName()
                )
                h2s_FERSLG_VS_DRS.append(h2_FERSLG_VS_DRS)

    # sum of FERS and DRS outputs
    h2s_FERS_VS_DRS_sum = []
    h2s_FERSLG_VS_DRS_sum = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for var in ["Cer", "Sci"]:
            h2sum = ROOT.TH2F(
                f"hist_FERS_VS_DRS_Board{boardNo}_{var}_sum",
                f"FERS vs DRS energy correlation for Board{boardNo}, {var}",
                100, DRS_min, DRS_max, 100, FERS_min, FERS_max
            )
            for h2 in h2s_FERS_VS_DRS:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum.Add(h2.GetValue())
            h2s_FERS_VS_DRS_sum.append(h2sum)

            h2sum_LG = ROOT.TH2F(
                f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_sum",
                f"FERS LG vs DRS energy correlation for Board{boardNo}, {var}",
                100, DRS_min, DRS_LG_max, 100, FERS_min, FERS_LG_max
            )
            for h2 in h2s_FERSLG_VS_DRS:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum_LG.Add(h2.GetValue())
            h2s_FERSLG_VS_DRS_sum.append(h2sum_LG)

            # Save the histograms to a ROOT file
    rootdir = f"root/Run{runNumber}"
    output_file = ROOT.TFile(os.path.join(
        rootdir, f"checkFERSDRS_{suffix}.root"), "RECREATE")
    for h2 in h2s_FERS_VS_DRS:
        h2.Write()
    for h2 in h2s_FERS_VS_DRS_sum:
        h2.Write()
    for h2 in h2s_FERSLG_VS_DRS:
        h2.Write()
    for h2 in h2s_FERSLG_VS_DRS_sum:
        h2.Write()
    output_file.Close()
    print(f"Histograms saved to {output_file.GetName()}")

    # snapshot DRS and FERS board 10
    variables = ["event_n"]
    variables += [f"{varname}_subtractMedian_positive" for _, DRSBoard in DRSBoards.items()
                  for channel in DRSBoard for varname in [channel.GetChannelName()]]
    variables += [f"{varname}" for _, DRSBoard in DRSBoards.items()
                  for channel in DRSBoard for varname in [channel.GetChannelName()]]
    variables += [f"FERS_Board{FERSBoard.boardNo}_energyHG_{channel.channelNo}"
                  for channel in FERSBoard]
    variables += [f"FERS_Board{FERSBoard.boardNo}_energyLG_{channel.channelNo}"
                  for channel in FERSBoard]
    rdf.Snapshot("DRSBoards", os.path.join(
        rootdir, f"DRSBoards_{suffix}.root"), variables)


def makeFERSDRSPlots():
    inputfile_name = os.path.join(
        f"root/Run{runNumber}", f"checkFERSDRS_{suffix}.root")
    input_file = ROOT.TFile(inputfile_name, "READ")
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open file {inputfile_name}")
        return
    print(f"Opened file {inputfile_name} successfully")

    plots = []
    outdir_plots = f"plots/Run{runNumber}/checkFERSDRS"

    # correlate  FERS and DRS outputs
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue
                chan_FERS = None
                for fersNo, FERSBoard in FERSBoards.items():
                    chan_FERS = FERSBoard.GetChannelByTower(
                        iTowerX, iTowerY, isCer=(var == "Cer"))
                    if chan_FERS is not None:
                        break
                if chan_FERS is None:
                    print(
                        f"Warning: FERS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                h2_name = f"hist_FERS_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist = input_file.Get(h2_name)
                if not hist:
                    print(f"Warning: Histogram {h2_name} not found in file")
                    continue

                output_name = f"FERS_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", DRS_min, DRS_max, "DRS Integral", FERS_min, FERS_max, f"FERS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber)
                plots.append(output_name + ".png")

                h2_LG_name = f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist_LG = input_file.Get(h2_LG_name)
                if not hist_LG:
                    print(f"Warning: Histogram {h2_LG_name} not found in file")
                    continue

                output_name_LG = f"FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist_LG], "", DRS_min, DRS_LG_max, "DRS Integral", FERS_min, FERS_LG_max, f"FERS Output",
                           output_name_LG,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber)
                plots.append(output_name_LG + ".png")

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for var in ["Cer", "Sci"]:
            h2sum_name = f"hist_FERS_VS_DRS_Board{boardNo}_{var}_sum"
            hist_sum = input_file.Get(h2sum_name)
            if not hist_sum:
                print(f"Warning: Histogram {h2sum_name} not found in file")
                continue

            output_name = f"FERS_VS_DRS_Board{boardNo}_{var}_sum_vs_Event"
            DrawHistos([hist_sum], "", DRS_min, DRS_max, "DRS Integral", FERS_min, FERS_max, f"FERS Output",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                       outdir=outdir_plots, extraText=var, runNumber=runNumber)
            plots.append(output_name + ".png")

            h2sum_LG_name = f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_sum"
            hist_sum_LG = input_file.Get(h2sum_LG_name)
            if not hist_sum_LG:
                print(f"Warning: Histogram {h2sum_LG_name} not found in file")
                continue

            output_name_LG = f"FERSLG_VS_DRS_Board{boardNo}_{var}_sum_vs_Event"
            DrawHistos([hist_sum_LG], "", DRS_min, DRS_LG_max, "DRS Integral", FERS_min, FERS_LG_max, f"FERS Output",
                       output_name_LG,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                       outdir=outdir_plots, extraText=var, runNumber=runNumber)
            plots.append(output_name_LG + ".png")

    generate_html(plots, outdir_plots,
                  output_html=f"html/Run{runNumber}/checkFERSDRS/view.html")


# snapshot DRSBoards and FERSBoards
# variables = [f"{varname}_subtracted" for _, DRSBoard in DRSBoards.items()
#             for channel in DRSBoard for varname in [channel.GetChannelName()]]
# rdf.Snapshot("DRSBoards", os.path.join(
#    rootdir, f"DRSBoards_{suffix}.root"), variables)

if __name__ == "__main__":
    start_time = time.time()
    prepareFERSDRSPlots()
    makeFERSDRSPlots()
    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")
