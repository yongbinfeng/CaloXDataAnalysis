import os
import ROOT
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels, getUpstreamVetoChannel, getDownStreamMuonChannel, getServiceDRSChannels
from utils.utils import number2string, preProcessDRSBoards, loadRDF, vectorizeFERS, prepareDRSStats, getFERSBoardMax
from configs.plotranges import getDRSPlotRanges, getServiceDRSPlotRanges
from runconfig import runNumber, firstEvent, lastEvent
import time
import sys

print("Start running prepareDQMPlots.py")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

debugDRS = False

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)

# Get total number of entries
n_entries = rdf.Count().GetValue()
nEvents = int(n_entries)
nbins_Event = min(max(int(nEvents / 100), 1), 500)
print(f"Total number of events to process: {nEvents} in run {runNumber}")

rdf = vectorizeFERS(rdf, FERSBoards)
rdf = preProcessDRSBoards(rdf, debug=debugDRS)
rdf = prepareDRSStats(rdf, DRSBoards, 0, 1000, 9)

rdf = getFERSBoardMax(rdf, FERSBoards)


def monitorConditions():
    # monitor the V, I, T conditions
    hists2d_Condition_vs_Event = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_SipmHV = rdf.Histo2D((
            f"hist_FERS_Board{boardNo}_SipmHV_vs_Event",
            f"FERS Board {boardNo} - SipmHV vs Event;Event;SipmHV (V)",
            nbins_Event, 0, nEvents, 40, 26, 29),
            "event_n", f"FERS_Board{boardNo}_SipmHV"
        )
        hist_SipmI = rdf.Histo2D((
            f"hist_FERS_Board{boardNo}_SipmI_vs_Event",
            f"FERS Board {boardNo} - SipmI vs Event;Event;SipmI (mA)",
            nbins_Event, 0, nEvents, 50, 0.02, 0.2),
            "event_n", f"FERS_Board{boardNo}_SipmI"
        )
        hist_TempDET = rdf.Histo2D((
            f"hist_FERS_Board{boardNo}_TempDET_vs_Event",
            f"FERS Board {boardNo} - TempDET vs Event;Event;TempDET (C)",
            nbins_Event, 0, nEvents, 100, 10, 30),
            "event_n", f"FERS_Board{boardNo}_TempDET"
        )
        hist_TempFPGA = rdf.Histo2D((
            f"hist_FERS_Board{boardNo}_TempFPGA_vs_Event",
            f"FERS Board {boardNo} - TempFPGA vs Event;Event;TempFPGA (C)",
            nbins_Event, 0, nEvents, 100, 30, 50),
            "event_n", f"FERS_Board{boardNo}_TempFPGA"
        )
        hists2d_Condition_vs_Event.append(hist_SipmHV)
        hists2d_Condition_vs_Event.append(hist_SipmI)
        hists2d_Condition_vs_Event.append(hist_TempDET)
        hists2d_Condition_vs_Event.append(hist_TempFPGA)

    return hists2d_Condition_vs_Event


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


def collectFERSPedestals(hists1d_FERS):
    pedestals = {}
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                channelName_HG = chan.GetHGChannelName()

                hname = f"hist_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}"
                hist = None
                for h in hists1d_FERS:
                    if h.GetName() == hname:
                        hist = h
                        break
                if hist is None:
                    print(
                        f"Warning: Histogram {hname} not found in hists1d_FERS")
                    pedestals[var][channelName_HG] = None
                    continue
                pedestal = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())
                pedestals[channelName_HG] = pedestal
    return pedestals


def collectFERSStats():
    from configs.plotranges import getFERSSaturationValue
    saturation_value = getFERSSaturationValue()
    stats = {}
    # mean, max,
    # and how frequent the saturation value is reached
    for _, FERSBoard in FERSBoards.items():
        for chan in FERSBoard:
            channelName_HG = chan.GetHGChannelName()
            stats[channelName_HG] = (
                rdf.Mean(channelName_HG),
                rdf.Max(channelName_HG),
                rdf.Filter(f"{channelName_HG} >= {saturation_value}").Count()
            )
            channelName_LG = chan.GetLGChannelName()
            stats[channelName_LG] = (
                rdf.Mean(channelName_LG),
                rdf.Max(channelName_LG),
                rdf.Filter(f"{channelName_LG} >= {saturation_value}").Count()
            )

    return stats


def makeFERSMaxValueHists():
    """
    get the max FERS readout per board and per event 
    """
    hists_board_cer_HG_max = []
    hists_board_cer_LG_max = []
    hists_board_sci_HG_max = []
    hists_board_sci_LG_max = []
    nbins = 100
    xmin = 7500
    xmax = 8500
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_board_cer_HG_max = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_energy_cer_HG_max",
            f"FERS Board {boardNo} - CER Energy HG max;CER Energy HG max per board;Counts",
            nbins, xmin, xmax),
            f"FERS_Board{boardNo}_energy_cer_HG_max"
        )
        hist_board_cer_LG_max = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_energy_cer_LG_max",
            f"FERS Board {boardNo} - CER Energy LG max;CER Energy LG max per board;Counts",
            nbins, xmin, xmax),
            f"FERS_Board{boardNo}_energy_cer_LG_max"
        )
        hist_board_sci_HG_max = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_energy_sci_HG_max",
            f"FERS Board {boardNo} - SCI Energy HG max;SCI Energy HG max per board;Counts",
            nbins, xmin, xmax),
            f"FERS_Board{boardNo}_energy_sci_HG_max"
        )
        hist_board_sci_LG_max = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_energy_sci_LG_max",
            f"FERS Board {boardNo} - SCI Energy LG max;SCI Energy LG max per board;Counts",
            nbins, xmin, xmax),
            f"FERS_Board{boardNo}_energy_sci_LG_max"
        )
        hists_board_cer_HG_max.append(hist_board_cer_HG_max)
        hists_board_cer_LG_max.append(hist_board_cer_LG_max)
        hists_board_sci_HG_max.append(hist_board_sci_HG_max)
        hists_board_sci_LG_max.append(hist_board_sci_LG_max)

    hist_cer_HG_max = rdf.Histo1D((
        "hist_FERS_energy_cer_HG_max",
        "FERS - CER Energy HG max;CER Energy HG max;Counts",
        nbins, xmin, xmax),
        "FERS_energy_cer_HG_max"
    )
    hist_cer_LG_max = rdf.Histo1D((
        "hist_FERS_energy_cer_LG_max",
        "FERS - CER Energy LG max;CER Energy LG max;Counts",
        nbins, xmin, xmax),
        "FERS_energy_cer_LG_max"
    )
    hist_sci_HG_max = rdf.Histo1D((
        "hist_FERS_energy_sci_HG_max",
        "FERS - SCI Energy HG max;SCI Energy HG max;Counts",
        nbins, xmin, xmax),
        "FERS_energy_sci_HG_max"
    )
    hist_sci_LG_max = rdf.Histo1D((
        "hist_FERS_energy_sci_LG_max",
        "FERS - SCI Energy LG max;SCI Energy LG max;Counts",
        nbins, xmin, xmax),
        "FERS_energy_sci_LG_max"
    )
    return hists_board_cer_HG_max + hists_board_cer_LG_max + hists_board_sci_HG_max + hists_board_sci_LG_max + [hist_cer_HG_max, hist_cer_LG_max, hist_sci_HG_max, hist_sci_LG_max]


def trackFERSPlots():
    hists2d_FERS_vs_Event = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = FERSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist = rdf.Histo2D((
                    f"hist_FERS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}",
                    f"FERS Board {boardNo} - Event vs {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};Event;{var} Energy HG",
                    nbins_Event, 0, nEvents, 1000, 0, 9000),
                    "event_n", chan.GetHGChannelName()
                )
                hists2d_FERS_vs_Event.append(hist)
    return hists2d_FERS_vs_Event


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
                300, 0, 9000, 300, 0, 9000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hist_zoomed = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_vs_Sci_{sTowerX}_{sTowerY}_zoom",
                f"CER {iCer} vs SCI {iSci} in iTowerX {sTowerX} iTowerY {sTowerY} (zoomed);CER Energy HG;SCI Energy HG",
                300, 0, 1000, 200, 0, 2000),
                chan_Cer.GetHGChannelName(),
                chan_Sci.GetHGChannelName()
            )
            hists2d_FERS.append(hist)
            hists2d_FERS.append(hist_zoomed)

            # high gain vs low gain for Sci
            hist_sci_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"SCI {iSci} HG vs LG;SCI Energy HG;SCI Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_Sci.GetHGChannelName(),
                chan_Sci.GetLGChannelName()
            )
            hists2d_FERS.append(hist_sci_hg_vs_lg)
            hist_cer_hg_vs_lg = rdf.Histo2D((
                f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}_hg_vs_lg",
                f"CER {iCer} HG vs LG;CER Energy HG;CER Energy LG",
                300, 0, 9000, 300, 0, 3000),
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

                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                value_mean = stats[channelName]['mean']
                hist = rdf.Histo1D((
                    f"hist_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Board {boardNo} - {var} iTowerX {sTowerX} iTowerY {sTowerY};{var} Variable;Counts",
                    200, value_mean - 100, value_mean + 100),
                    channelName
                )
                hists1d_DRS.append(hist)
    return hists1d_DRS


def makeDRS2DPlots(debug=False):
    hists2d_DRS_vs_TS = []
    if debug:
        hists2d_DRS_vs_RTSpos = []
        hists2d_DRS_vs_RTSneg = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                ymin, ymax = getDRSPlotRanges(
                    subtractMedian=True, isAmplified=chan.isAmplified)
                hist_subtractMedian = rdf.Histo2D((
                    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}_subtractMedian",
                    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);TS;{var} Variable",
                    1024, 0, 1024, 50, ymin, ymax),
                    "TS", channelName + "_subtractMedian"
                )
                hists2d_DRS_vs_TS.append(hist_subtractMedian)

                if debug:
                    hist_subtractMedian_RTSpos = rdf.Histo2D((
                        f"hist_DRS_Board{boardNo}_{var}_vs_RTSpos_{sTowerX}_{sTowerY}_subtractMedian",
                        f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);RTS pos;{var} Variable",
                        1024, 0, 1024, 50, ymin, ymax),
                        f"RTS_pos_{channelName}", channelName +
                        "_subtractMedian"
                    )
                    hist_subtractMedian_RTSneg = rdf.Histo2D((
                        f"hist_DRS_Board{boardNo}_{var}_vs_RTSneg_{sTowerX}_{sTowerY}_subtractMedian",
                        f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);RTS neg;{var} Variable",
                        1024, 0, 1024, 50, ymin, ymax),
                        f"RTS_neg_{channelName}", channelName +
                        "_subtractMedian"
                    )
                    hists2d_DRS_vs_RTSpos.append(hist_subtractMedian_RTSpos)
                    hists2d_DRS_vs_RTSneg.append(hist_subtractMedian_RTSneg)
    if debug:
        return hists2d_DRS_vs_TS, hists2d_DRS_vs_RTSpos, hists2d_DRS_vs_RTSneg
    return hists2d_DRS_vs_TS


def trackDRSPlots():
    hists2d_DRS_vs_Event = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                # Get the channel for CER or SCI
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan is None:
                    continue
                channelName = chan.GetChannelName()
                mean_value = stats[channelName]['mean']
                hist = rdf.Histo2D((
                    f"hist_DRS_Board{boardNo}_{var}_vs_Event_{sTowerX}_{sTowerY}",
                    f"DRS Board {boardNo} Mean - Event vs {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};Event;{var} Variable",
                    nbins_Event, 0, nEvents, 200, mean_value - 100, mean_value + 100),
                    "event_n", channelName + "_mean"
                )
                hists2d_DRS_vs_Event.append(hist)
    return hists2d_DRS_vs_Event


def compareDRSChannels(channels_to_compare, isServiceDRS=False):
    hists_trigger = []
    for chan_name in channels_to_compare:
        ymin = -2500
        ymax = 500
        if isServiceDRS:
            ymin, ymax = getServiceDRSPlotRanges(
                chan_name, subtractMedian=True)
        hist_subtractMedian = rdf.Histo2D((
            f"hist_{chan_name}_subtractMedian",
            f"{chan_name} (subtract median);TS;DRS values",
            1024, 0, 1024,
            300, ymin, ymax),
            "TS", chan_name + "_subtractMedian"
        )
        hists_trigger.append(hist_subtractMedian)
    return hists_trigger


def checkFERSvsDRSSum():
    xymax = {
        "Cer": (1000, 4000),
        "Sci": (9000, 9000)
    }
    xymax_LG = {
        "Cer": (1000, 1000),
        "Sci": (9000, 9000)
    }

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
                    100, 0, xymax[var][0], 100, 0, xymax[var][1]
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetHGChannelName(),
                )
                h2s_FERS_VS_DRS.append(h2_FERS_VS_DRS)

                h2_FERSLG_VS_DRS = rdf.Histo2D((
                    f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS LG vs DRS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, 0, xymax_LG[var][0], 100, 0, xymax_LG[var][1]
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetLGChannelName(),
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
                100, 0, xymax[var][0], 100, 0, xymax[var][1]
            )
            for h2 in h2s_FERS_VS_DRS:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum.Add(h2.GetValue())
            h2s_FERS_VS_DRS_sum.append(h2sum)

            h2sum_LG = ROOT.TH2F(
                f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_sum",
                f"FERS LG vs DRS energy correlation for Board{boardNo}, {var}",
                100, 0, xymax_LG[var][0], 100, 0, xymax_LG[var][1]
            )
            for h2 in h2s_FERSLG_VS_DRS:
                if f"Board{boardNo}_{var}" in h2.GetName():
                    h2sum_LG.Add(h2.GetValue())
            h2s_FERSLG_VS_DRS_sum.append(h2sum_LG)
    return h2s_FERS_VS_DRS_sum + h2s_FERSLG_VS_DRS_sum + h2s_FERS_VS_DRS + h2s_FERSLG_VS_DRS


def checkDRSPeakvsFERS():
    h2s_DRSPeak_vs_FERS = []
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

                _, ymax = getDRSPlotRanges(
                    subtractMedian=True, isAmplified=chan_DRS.isAmplified)
                h2_DRSPeak_vs_FERS = rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"DRS Peak vs FERS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, 0, 9000, 100, 0, ymax
                ),
                    chan_FERS.GetHGChannelName(),
                    f"{chan_DRS.GetChannelName()}_peak",
                )
                h2s_DRSPeak_vs_FERS.append(h2_DRSPeak_vs_FERS)

    return h2s_DRSPeak_vs_FERS


def checkDRSPeakTS():
    h1s_DRSPeakTS = {}
    h1s_DRSPeakTS["Cer"] = []
    h1s_DRSPeakTS["Sci"] = []
    h2s_DRSPeakTS_Cer_vs_Sci = []

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            channelNames = {}
            for var in ["Cer", "Sci"]:
                chan_DRS = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan_DRS is None:
                    print(
                        f"Warning: DRS Channel not found for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}")
                    continue

                channelName = chan_DRS.GetChannelName()
                channelNames[var] = channelName

                h1_DRS_PeakTS = rdf.Histo1D((
                    f"hist_DRS_PeakTS_Board{boardNo}_peakTS_{sTowerX}_{sTowerY}_{var}",
                    f"DRS Peak TS for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var};Peak TS;Counts",
                    1000, 0, 1000),
                    channelName + "_peakTS"
                )
                h1s_DRSPeakTS[var].append(h1_DRS_PeakTS)

            if len(channelNames) < 2:
                print(
                    f"Warning: Not enough channels found for Board{boardNo}, Tower({sTowerX}, {sTowerY})")
                continue

            h2_DRSPeak_Cer_vs_Sci = rdf.Histo2D((
                f"hist_DRSPeak_Cer_vs_Sci_Board{boardNo}_{sTowerX}_{sTowerY}",
                f"DRS Peak TS - CER vs SCI for Board{boardNo}, Tower({sTowerX}, {sTowerY});CER Peak TS;SCI Peak TS",
                1000, 0, 1000, 1000, 0, 1000),
                channelNames["Cer"] + "_peakTS",
                channelNames["Sci"] + "_peakTS"
            )
            h2s_DRSPeakTS_Cer_vs_Sci.append(h2_DRSPeak_Cer_vs_Sci)
    return h1s_DRSPeakTS["Cer"], h1s_DRSPeakTS["Sci"], h2s_DRSPeakTS_Cer_vs_Sci


if __name__ == "__main__":
    start_time = time.time()

    hists_conditions = monitorConditions()

    hists1d_FERS = makeFERS1DPlots()

    # hists2d_FERS = makeFERS2DPlots()
    # hists2d_FERS_vs_Event = trackFERSPlots()

    # hists1d_DRS = makeDRS1DPlots()
    hists2d_DRS_vs_RTSpos = None
    hists2d_DRS_vs_RTSneg = None
    if debugDRS:
        hists2d_DRS_vs_TS, hists2d_DRS_vs_RTSpos, hists2d_DRS_vs_RTSneg = makeDRS2DPlots(
            debug=True)
    else:
        hists2d_DRS_vs_TS = makeDRS2DPlots(debug=False)
    # hists2d_DRS_vs_Event = trackDRSPlots()

    hists2d_DRSPeak_vs_FERS = checkDRSPeakvsFERS()

    hists1d_DRSPeakTS_Cer, hists1d_DRSPeakTS_Sci, hists2d_DRSPeakTS_Cer_vs_Sci = checkDRSPeakTS()

    time_reference_channels = buildTimeReferenceChannels(run=runNumber)
    hists2d_time_reference = compareDRSChannels(time_reference_channels)

    hists2d_hodo_trigger = None
    # hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
    # hists2d_hodo_trigger = compareDRSChannels(hodo_trigger_channels)

    service_drs_channels = getServiceDRSChannels(run=runNumber)
    hists2d_service_drs = compareDRSChannels(
        service_drs_channels, isServiceDRS=True)

    upstream_veto_channel = getUpstreamVetoChannel(run=runNumber)
    hists2d_veto = compareDRSChannels([upstream_veto_channel])

    downstream_muon_channel = getDownStreamMuonChannel(run=runNumber)
    hists2d_muon = compareDRSChannels([downstream_muon_channel])

    hodo_pos_channels = buildHodoPosChannels(run=runNumber)
    channels = [channel for channels in hodo_pos_channels.values()
                for channel in channels]
    hists2d_hodo_pos = compareDRSChannels(channels)

    hists2d_FERS_vs_DRSs = checkFERSvsDRSSum()

    hists_FERS_max = makeFERSMaxValueHists()

    stats = collectFERSStats()

    print("\033[94mSave results\033[0m")

    rootdir = f"results/root/Run{runNumber}"
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)

    # dump stats into a json file
    if 'stats' in locals() and stats:
        import json
        stats_results = {}
        for channelName, (mean, max_value, sat_frequency) in stats.items():
            stats_results[channelName] = (mean.GetValue(), max_value.GetValue(), float(
                sat_frequency.GetValue()) / nEvents)
        with open(f"{rootdir}/fers_stats.json", "w") as f:
            json.dump(stats_results, f, indent=4)

    # Save histograms to an output ROOT file
    if 'hists_conditions' in locals() and hists_conditions:
        outfile = ROOT.TFile(f"{rootdir}/conditions_vs_event.root", "RECREATE")
        for hist in hists_conditions:
            hist.SetDirectory(outfile)
            hist.Write()
        outfile.Close()

    if 'hists1d_FERS' in locals() and hists1d_FERS:
        outfile = ROOT.TFile(
            f"{rootdir}/fers_all_channels_1D.root", "RECREATE")
        for hist in hists1d_FERS:
            hist.Write()
        outfile.Close()
    # outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_2D.root", "RECREATE")
    # for hist in hists2d_FERS:
    #    hist.Write()
    # outfile.Close()
    # outfile = ROOT.TFile(
    #    f"{rootdir}/fers_all_channels_2D_vs_event.root", "RECREATE")
    # for hist in hists2d_FERS_vs_Event:
    #    hist.Write()
    # outfile.Close()
    #
    # outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_1D.root", "RECREATE")
    # for hist in hists1d_DRS:
    #    hist.SetDirectory(outfile_DRS)
    #    hist.Write()
    # outfile_DRS.Close()
    outfile_DRS = ROOT.TFile(f"{rootdir}/drs_vs_TS.root", "RECREATE")
    for hist in hists2d_DRS_vs_TS:
        hist.SetDirectory(outfile_DRS)
        hist.Write()
    outfile_DRS.Close()
    if debugDRS:
        outfile_DRS_RTSpos = ROOT.TFile(
            f"{rootdir}/drs_vs_RTSpos.root", "RECREATE")
        for hist in hists2d_DRS_vs_RTSpos:
            hist.SetDirectory(outfile_DRS_RTSpos)
            hist.Write()
        outfile_DRS_RTSpos.Close()

        outfile_DRS_RTSneg = ROOT.TFile(
            f"{rootdir}/drs_vs_RTSneg.root", "RECREATE")
        for hist in hists2d_DRS_vs_RTSneg:
            hist.SetDirectory(outfile_DRS_RTSneg)
            hist.Write()
        outfile_DRS_RTSneg.Close()
    # outfile_DRS = ROOT.TFile(
    #    f"{rootdir}/drs_all_channels_2D_vs_event.root", "RECREATE")
    # for hist in hists2d_DRS_vs_Event:
    #    hist.SetDirectory(outfile_DRS)
    #    hist.Write()
    # outfile_DRS.Close()

    outfile_DRSPeak_vs_FERS = ROOT.TFile(
        f"{rootdir}/drs_peak_vs_fers.root", "RECREATE")
    for hist in hists2d_DRSPeak_vs_FERS:
        hist.SetDirectory(outfile_DRSPeak_vs_FERS)
        hist.Write()
    outfile_DRSPeak_vs_FERS.Close()

    outfile_DRSPeakTS = ROOT.TFile(f"{rootdir}/drs_peak_ts.root", "RECREATE")
    for hist in hists1d_DRSPeakTS_Cer:
        hist.SetDirectory(outfile_DRSPeakTS)
        hist.Write()
    for hist in hists1d_DRSPeakTS_Sci:
        hist.SetDirectory(outfile_DRSPeakTS)
        hist.Write()
    for hist in hists2d_DRSPeakTS_Cer_vs_Sci:
        hist.SetDirectory(outfile_DRSPeakTS)
        hist.Write()
    outfile_DRSPeakTS.Close()

    outfile_time_reference = ROOT.TFile(
        f"{rootdir}/time_reference_channels.root", "RECREATE")
    for hist in hists2d_time_reference:
        hist.SetDirectory(outfile_time_reference)
        hist.Write()
    outfile_time_reference.Close()

    if 'hists2d_hodo_trigger' in locals() and hists2d_hodo_trigger:
        outfile_hodo_trigger = ROOT.TFile(
            f"{rootdir}/hodo_trigger_channels.root", "RECREATE")
        for hist in hists2d_hodo_trigger:
            hist.SetDirectory(outfile_hodo_trigger)
            hist.Write()
        outfile_hodo_trigger.Close()

    if 'hists2d_service_drs' in locals() and hists2d_service_drs:
        outfile_service_drs = ROOT.TFile(
            f"{rootdir}/service_drs_channels.root", "RECREATE")
        for hist in hists2d_service_drs:
            hist.SetDirectory(outfile_service_drs)
            hist.Write()
        outfile_service_drs.Close()

    if 'hists2d_veto' in locals() and hists2d_veto:
        outfile_veto = ROOT.TFile(
            f"{rootdir}/upstream_veto_channel.root", "RECREATE")
        for hist in hists2d_veto:
            hist.SetDirectory(outfile_veto)
            hist.Write()
        outfile_veto.Close()

    if 'hists2d_muon' in locals() and hists2d_muon:
        outfile_muon = ROOT.TFile(
            f"{rootdir}/downstream_muon_channel.root", "RECREATE")
        for hist in hists2d_muon:
            hist.SetDirectory(outfile_muon)
            hist.Write()
        outfile_muon.Close()

    outfile_hodo_pos = ROOT.TFile(
        f"{rootdir}/hodo_pos_channels.root", "RECREATE")
    for hist in hists2d_hodo_pos:
        hist.SetDirectory(outfile_hodo_pos)
        hist.Write()
    outfile_hodo_pos.Close()

    outfile_FERS_DRS = ROOT.TFile(
        f"{rootdir}/fers_vs_drs.root", "RECREATE")
    for hist in hists2d_FERS_vs_DRSs:
        hist.SetDirectory(outfile_FERS_DRS)
        hist.Write()
    outfile_FERS_DRS.Close()

    outfile_FERS_max = ROOT.TFile(
        f"{rootdir}/fers_max_values.root", "RECREATE")
    for hist in hists_FERS_max:
        hist.SetDirectory(outfile_FERS_max)
        hist.Write()
    outfile_FERS_max.Close()

    time_taken = time.time() - start_time
    print(f"Finished running script in {time_taken:.2f} seconds")
