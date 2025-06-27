import sys
import os
import ROOT
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTriggerChannels, buildHodoTriggerChannels
from utils.utils import number2string, getDataFile, processDRSBoards
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
trigger_channels = buildTriggerChannels(run=runNumber)
hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)

FERS_min = 100
FERS_max = 9e3
FERS_LG_max = 2e3
DRS_min = -100
DRS_max = 1e3
DRS_LG_max = 2e3


def findTrigFireTime(rdf, channels):
    ROOT.gInterpreter.Declare("""
    size_t findTrigFireTime(const ROOT::VecOps::RVec<float>& vec, float val_min) {
        for (size_t i = 0; i < vec.size(); ++i) {
            if (vec[i] < val_min / 2.0) {
                return i;  // return the index of the first value below the threshold
            }
        }
        return -1;  // return -1 if no value is below the threshold
    }
    """)
    for channel in channels:
        rdf = rdf.Define("TrigMin_" + channel,
                         f"ROOT::VecOps::Min({channel}_subtractMedian)")
        # find the half minimum TS
        rdf = rdf.Define("TrigHalfMin_" + channel,
                         f"findTrigFireTime({channel}_subtractMedian, TrigMin_{channel})")
    return rdf


def prepareFERSDRSPlots():
    ifile = getDataFile(runNumber)
    infile = ROOT.TFile(ifile, "READ")
    rdf_temp = ROOT.RDataFrame("EventTree", infile)

    rdf = rdf_temp.Filter("event_n > 1")  # filter out the two events

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

    rdf = processDRSBoards(rdf)

    # get the mean of DRS outputs per channel
    ROOT.gInterpreter.Declare("""
    ROOT::VecOps::RVec<float> clipToZero(const ROOT::VecOps::RVec<float>& vec) {
        ROOT::VecOps::RVec<float> out;
        for (float v : vec) {
            if (fabs(v) < 3.0f) v = 0.0f;  // clip to zero if below threshold
            out.push_back(v);
        }
        return out;
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
                f"{varname}_subtractMedian_positive",
                f"clipToZero({varname}_subtractMedian)"
            )
            rdf = rdf.Define(
                f"{varname}_sum",
                f"SumRange({varname}_subtractMedian_positive, 0, 400)"
            )

    rdf = findTrigFireTime(rdf, trigger_channels + hodo_trigger_channels)
    rdf = rdf.Define("TriggerHodo_deltaT",
                     f"TrigHalfMin_{hodo_trigger_channels[1]} - TrigHalfMin_{hodo_trigger_channels[0]}")

    rdf = rdf.Define(
        "passTime", f"TrigHalfMin_{hodo_trigger_channels[0]} > 200 && TrigMin_{hodo_trigger_channels[1]} < -200")
    # rdf = rdf.Define("passTime", "1.0")

    # correlate  FERS and DRS outputs
    h2s_FERS_VS_DRS = []
    h2s_FERSLG_VS_DRS = []
    h2s_DRSOverFERS_VS_HodoUp = []
    h2s_DRSOverFERSLG_VS_HodoUp = []
    vars_ratio = []
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
                for _, FERSBoard in FERSBoards.items():
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
                    chan_FERS.GetHGChannelName(),
                    "passTime"
                )
                h2s_FERS_VS_DRS.append(h2_FERS_VS_DRS)

                h2_FERSLG_VS_DRS = rdf.Histo2D((
                    f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS LG vs DRS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, DRS_min, DRS_LG_max, 100, FERS_min, FERS_LG_max
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetLGChannelName(),
                    "passTime"
                )
                h2s_FERSLG_VS_DRS.append(h2_FERSLG_VS_DRS)

                rdf = rdf.Define(f"DRSOverFERS_{var}_{sTowerX}_{sTowerY}",
                                 f"std::min(0.19f, {chan_DRS.GetChannelName()}_sum / {chan_FERS.GetHGChannelName()})")
                rdf = rdf.Define(f"DRSOverFERSLG_{var}_{sTowerX}_{sTowerY}",
                                 f"std::min(0.19f, {chan_DRS.GetChannelName()}_sum / {chan_FERS.GetLGChannelName()})")

                vars_ratio.append(f"DRSOverFERS_{var}_{sTowerX}_{sTowerY}")
                vars_ratio.append(f"DRSOverFERSLG_{var}_{sTowerX}_{sTowerY}")

                h2_DRSOverFERS_VS_HodoUp = rdf.Histo2D((
                    f"hist_DRSOverFERS_VS_HodoUp_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS over FERS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var} vs Hodo Up",
                    100, 0, 0.2, 1025, -1, 1024
                ),
                    f"DRSOverFERS_{var}_{sTowerX}_{sTowerY}",
                    f"TrigHalfMin_{hodo_trigger_channels[0]}",
                    "passTime"  # only include events where Hodo Up trigger is passed
                )
                h2_DRSOverFERSLG_VS_HodoUp = rdf.Histo2D((
                    f"hist_DRSOverFERSLG_VS_HodoUp_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS LG over FERS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var} vs Hodo Up",
                    100, 0, 0.2, 1025, -1, 1024
                ),
                    f"DRSOverFERSLG_{var}_{sTowerX}_{sTowerY}",
                    f"TrigHalfMin_{hodo_trigger_channels[0]}",
                    "passTime"  # only include events where Hodo Up trigger is passed
                )
                h2s_DRSOverFERS_VS_HodoUp.append(h2_DRSOverFERS_VS_HodoUp)
                h2s_DRSOverFERSLG_VS_HodoUp.append(h2_DRSOverFERSLG_VS_HodoUp)

    # trigger fire time
    hists_trig_fire_time = []
    for channel in trigger_channels + hodo_trigger_channels:
        h_trig_fire_time = rdf.Histo1D((
            f"hist_TrigFireTime_{channel}",
            f"Trigger fire time for {channel}",
            1025, -1, 1024
        ), f"TrigHalfMin_{channel}")
        hists_trig_fire_time.append(h_trig_fire_time)
    h_trig_fire_deltaT = rdf.Histo1D((
        "hist_TrigFireTime_deltaT",
        "Hodo trigger delta T",
        200, -100, 100
    ), "TriggerHodo_deltaT")
    hists_trig_fire_time.append(h_trig_fire_deltaT)

    # sum of FERS and DRS outputs
    h2s_FERS_VS_DRS_sum = []
    h2s_FERSLG_VS_DRS_sum = []
    h2s_DRSOverFERS_VS_HodoUp_sum = []
    h2s_DRSOverFERSLG_VS_HodoUp_sum = []
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

            h2sum_DRSOverFERS = ROOT.TH2F(
                f"hist_DRSOverFERS_VS_HodoUp_{var}_sum",
                f"DRS over FERS energy correlation for Board{boardNo}, {var}",
                100, 0, 0.2, 1025, -1, 1024
            )
            for h2 in h2s_DRSOverFERS_VS_HodoUp:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum_DRSOverFERS.Add(h2.GetValue())
            h2s_DRSOverFERS_VS_HodoUp_sum.append(h2sum_DRSOverFERS)

            h2sum_DRSOverFERSLG = ROOT.TH2F(
                f"hist_DRSOverFERSLG_VS_HodoUp_{var}_sum",
                f"DRS LG over FERS energy correlation for Board{boardNo}, {var}",
                100, 0, 0.2, 1025, -1, 1024
            )
            for h2 in h2s_DRSOverFERSLG_VS_HodoUp:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum_DRSOverFERSLG.Add(h2.GetValue())
            h2s_DRSOverFERSLG_VS_HodoUp_sum.append(h2sum_DRSOverFERSLG)

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
    for h2 in h2s_DRSOverFERS_VS_HodoUp:
        h2.Write()
    for h2 in h2s_DRSOverFERSLG_VS_HodoUp:
        h2.Write()
    for h2 in h2s_DRSOverFERS_VS_HodoUp_sum:
        h2.Write()
    for h2 in h2s_DRSOverFERSLG_VS_HodoUp_sum:
        h2.Write()
    output_file.Close()
    print(f"Histograms saved to {output_file.GetName()}")

    output_file = ROOT.TFile(os.path.join(
        rootdir, f"checkFERSDRS_trigFireTime_{suffix}.root"), "RECREATE")
    for h in hists_trig_fire_time:
        h.Write()
    output_file.Close()
    print(f"Trigger fire time histograms saved to {output_file.GetName()}")

    # snapshot DRS and FERS board 10
    if 0:
        variables = ["event_n"]
        variables += [f"{varname}_subtractMedian_positive" for _, DRSBoard in DRSBoards.items()
                      for channel in DRSBoard for varname in [channel.GetChannelName()]]
        variables += [f"{varname}_sum" for _, DRSBoard in DRSBoards.items()
                      for channel in DRSBoard for varname in [channel.GetChannelName()]]
        variables += [f"{varname}" for _, DRSBoard in DRSBoards.items()
                      for channel in DRSBoard for varname in [channel.GetChannelName()]]
        variables += [f"FERS_Board{FERSBoard.boardNo}_energyHG_{channel.channelNo}"
                      for channel in FERSBoard]
        variables += [f"FERS_Board{FERSBoard.boardNo}_energyLG_{channel.channelNo}"
                      for channel in FERSBoard]
        variables += [f"TrigMin_{channel}" for channel in trigger_channels +
                      hodo_trigger_channels]
        variables += [
            f"TrigHalfMin_{channel}" for channel in trigger_channels + hodo_trigger_channels]
        variables += vars_ratio
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

                output_name = f"FERS_VS_DRS_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist], "", DRS_min, DRS_max, "DRS Integral", FERS_min, FERS_max, f"FERS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
                plots.append(output_name + ".png")

                h2_LG_name = f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist_LG = input_file.Get(h2_LG_name)
                if not hist_LG:
                    print(f"Warning: Histogram {h2_LG_name} not found in file")
                    continue

                output_name_LG = f"FERSLG_VS_DRS_{var}_{sTowerX}_{sTowerY}_vs_Event"
                DrawHistos([hist_LG], "", DRS_min, DRS_LG_max, "DRS Integral", FERS_min, FERS_LG_max, f"FERS Output",
                           output_name_LG,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
                plots.append(output_name_LG + ".png")

                h2_name = f"hist_DRSOverFERS_VS_HodoUp_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist = input_file.Get(h2_name)
                if not hist:
                    print(f"Warning: Histogram {h2_name} not found in file")
                    continue
                output_name = f"DRSOverFERS_{var}_{sTowerX}_{sTowerY}_vs_HodoUp"
                DrawHistos([hist], "", 0, 0.2, "DRS / FERS",
                           -1, 1025, "Hodo Up Trigger Fire Time (TS)",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=0, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
                plots.append(output_name + ".png")
                h2_LG_name = f"hist_DRSOverFERSLG_VS_HodoUp_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist_LG = input_file.Get(h2_LG_name)
                if not hist_LG:
                    print(f"Warning: Histogram {h2_LG_name} not found in file")
                    continue
                output_name_LG = f"DRSOverFERSLG_{var}_{sTowerX}_{sTowerY}_vs_HodoUp"
                DrawHistos([hist_LG], "", 0, 0.2, "DRS / FERS LG ",
                           -1, 1025, "Hodo Up Trigger Fire Time (TS)",
                           output_name_LG,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=0, zmax=2e3, dologz=True,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
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
                       outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
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
                       outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
            plots.append(output_name_LG + ".png")

            h2sum_DRSOverFERS_name = f"hist_DRSOverFERS_VS_HodoUp_{var}_sum"
            hist_DRSOverFERS_sum = input_file.Get(h2sum_DRSOverFERS_name)
            if not hist_DRSOverFERS_sum:
                print(
                    f"Warning: Histogram {h2sum_DRSOverFERS_name} not found in file")
                continue
            output_name = f"DRSOverFERS_Board{boardNo}_{var}_sum_vs_HodoUp"
            DrawHistos([hist_DRSOverFERS_sum], "", 0, 0.2, "DRS / FERS",
                       -1, 1025, "Hodo Up Trigger Fire Time (TS)",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=0, zmax=2e3, dologz=True,
                       outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
            plots.append(output_name + ".png")

            h2sum_DRSOverFERSLG_name = f"hist_DRSOverFERSLG_VS_HodoUp_{var}_sum"
            hist_DRSOverFERSLG_sum = input_file.Get(h2sum_DRSOverFERSLG_name)
            if not hist_DRSOverFERSLG_sum:
                print(
                    f"Warning: Histogram {h2sum_DRSOverFERSLG_name} not found in file")
                continue
            output_name_LG = f"DRSOverFERSLG_Board{boardNo}_{var}_sum_vs_HodoUp"
            DrawHistos([hist_DRSOverFERSLG_sum], "", 0, 0.2, "DRS / FERS LG ",
                       -1, 1025, "Hodo Up Trigger Fire Time (TS)",
                       output_name_LG,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=0, zmax=2e3, dologz=True,
                       outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
            plots.append(output_name_LG + ".png")
    generate_html(plots, outdir_plots,
                  output_html=f"html/Run{runNumber}/checkFERSDRS/view.html")

    # trigger fire time plots
    inputfile_name = os.path.join(
        f"root/Run{runNumber}", f"checkFERSDRS_trigFireTime_{suffix}.root")
    input_file = ROOT.TFile(inputfile_name, "READ")
    if not input_file or input_file.IsZombie():
        print(f"Error: Could not open file {inputfile_name}")
        return
    print(f"Opened file {inputfile_name} successfully")
    plots_trig_fire_time = []
    outdir_trig_fire_time = f"plots/Run{runNumber}/checkFERSDRS/trigFireTime"
    for channel in trigger_channels:
        hist_name = f"hist_TrigFireTime_{channel}"
        hist = input_file.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in file")
            continue

        output_name = f"TrigFireTime_{channel}_vs_Event"
        DrawHistos([hist], "", -1, 1024, "Trigger Fire Time (TS)", 0, 0.04, "Count",
                   output_name,
                   dology=False,
                   outdir=outdir_trig_fire_time, runNumber=runNumber, mycolors=[1], donormalize=True)
        plots_trig_fire_time.append(output_name + ".png")

    # trigger fire time for hodo up and down
    hist_up_name = f"hist_TrigFireTime_{hodo_trigger_channels[0]}"
    hist_up = input_file.Get(hist_up_name)
    hist_down_name = f"hist_TrigFireTime_{hodo_trigger_channels[1]}"
    hist_down = input_file.Get(hist_down_name)
    output_name = f"HodoTrigger_{hodo_trigger_channels[0]}_{hodo_trigger_channels[1]}_vs_Event"
    DrawHistos([hist_up, hist_down], "", -1, 1024, "Trigger Fire Time (TS)", 0, 0.04, "Count",
               output_name,
               dology=False,
               outdir=outdir_trig_fire_time, runNumber=runNumber, mycolors=[1, 2], donormalize=True)
    plots_trig_fire_time.append(output_name + ".png")

    # trigger fire time for hodo delta T
    hist_deltaT = input_file.Get("hist_TrigFireTime_deltaT")
    if not hist_deltaT:
        print(f"Warning: Histogram hist_TrigFireTime_deltaT not found in file")
    else:
        output_name_deltaT = "HodoTrigger_deltaT_vs_Event"
        DrawHistos([hist_deltaT], "", -50, 100, "Hodo Trigger Delta T (TS)", 0, 0.04, "Count",
                   output_name_deltaT,
                   dology=False,
                   outdir=outdir_trig_fire_time, runNumber=runNumber, mycolors=[1], donormalize=True)
        plots_trig_fire_time.append(output_name_deltaT + ".png")
    generate_html(plots_trig_fire_time, outdir_trig_fire_time,
                  output_html=f"html/Run{runNumber}/checkFERSDRS/trigFireTime.html")


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
