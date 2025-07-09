import os
import ROOT
from utils.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels
from utils.utils import number2string, processDRSBoards, filterPrefireEvents, loadRDF, calculateEnergySumFERS, getDRSSum, vectorizeFERS
from runconfig import runNumber, firstEvent, lastEvent
import time

print("Start running prepareDQMPlots.py")

# multi-threading support
ROOT.ROOT.EnableImplicitMT(10)
ROOT.gSystem.Load("utils/functions_cc.so")  # Load the compiled C++ functions

rdf, rdf_org = loadRDF(runNumber, firstEvent, lastEvent)
rdf, rdf_prefilter = filterPrefireEvents(rdf)

DRSBoards = buildDRSBoards(run=runNumber)
FERSBoards = buildFERSBoards(run=runNumber)

# Get total number of entries
n_entries = rdf.Count().GetValue()
nEvents = int(n_entries)
nbins_Event = min(max(int(nEvents / 100), 1), 500)
print(f"Total number of events to process: {nEvents} in run {runNumber}")

rdf = vectorizeFERS(rdf, FERSBoards)
rdf = calculateEnergySumFERS(rdf, FERSBoards)
rdf = processDRSBoards(rdf)
rdf = getDRSSum(rdf, DRSBoards) 


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


def makeFERSEnergySumPlots():
    hists_FERS_EnergySum = []
    for _, FERSBoard in FERSBoards.items():
        boardNo = FERSBoard.boardNo
        hist_CerEnergyHG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_CerEnergyHG",
            f"FERS Board {boardNo} - CER Energy HG;CER Energy HG;Counts",
            500, 0, 64000),
            f"FERS_Board{boardNo}_CerEnergyHG"
        )
        hist_CerEnergyLG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_CerEnergyLG",
            f"FERS Board {boardNo} - CER Energy LG;CER Energy LG;Counts",
            500, 0, 64000),
            f"FERS_Board{boardNo}_CerEnergyLG"
        )
        hist_SciEnergyHG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_SciEnergyHG",
            f"FERS Board {boardNo} - SCI Energy HG;SCI Energy HG;Counts",
            500, 0, 64000),
            f"FERS_Board{boardNo}_SciEnergyHG"
        )
        hist_SciEnergyLG_Board = rdf.Histo1D((
            f"hist_FERS_Board{boardNo}_SciEnergyLG",
            f"FERS Board {boardNo} - SCI Energy LG;SCI Energy LG;Counts",
            500, 0, 64000),
            f"FERS_Board{boardNo}_SciEnergyLG"
        )
        hists_FERS_EnergySum.append(hist_CerEnergyHG_Board)
        hists_FERS_EnergySum.append(hist_CerEnergyLG_Board)
        hists_FERS_EnergySum.append(hist_SciEnergyHG_Board)
        hists_FERS_EnergySum.append(hist_SciEnergyLG_Board)

    # per-event energy sum
    hist_CerEnergyHG = rdf.Histo1D((
        "hist_FERS_CerEnergyHG",
        "FERS - CER Energy HG;CER Energy HG;Counts",
        500, 0, 8e5),
        "FERS_CerEnergyHG"
    )
    hist_CerEnergyLG = rdf.Histo1D((
        "hist_FERS_CerEnergyLG",
        "FERS - CER Energy LG;CER Energy LG;Counts",
        500, 0, 8e5),
        "FERS_CerEnergyLG"
    )
    hist_SciEnergyHG = rdf.Histo1D((
        "hist_FERS_SciEnergyHG",
        "FERS - SCI Energy HG;SCI Energy HG;Counts",
        500, 0, 8e5),
        "FERS_SciEnergyHG"
    )
    hist_SciEnergyLG = rdf.Histo1D((
        "hist_FERS_SciEnergyLG",
        "FERS - SCI Energy LG;SCI Energy LG;Counts",
        500, 0, 8e5),
        "FERS_SciEnergyLG"
    )
    hists_FERS_EnergySum.append(hist_CerEnergyHG)
    hists_FERS_EnergySum.append(hist_CerEnergyLG)
    hists_FERS_EnergySum.append(hist_SciEnergyHG)
    hists_FERS_EnergySum.append(hist_SciEnergyLG)

    return hists_FERS_EnergySum


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
                    1000, 0, 9000),
                    chan.GetHGChannelName()
                )
                hists1d_FERS.append(hist)

    return hists1d_FERS


def collectFERSStats():
    stats = {}
    for _, FERSBoard in FERSBoards.items():
        for chan in FERSBoard:
            channelName_HG = chan.GetHGChannelName()
            stats[channelName_HG] = (
                rdf.Mean(channelName_HG),
                rdf.Max(channelName_HG),
            )
            channelName_LG = chan.GetLGChannelName()
            stats[channelName_LG] = (
                rdf.Mean(channelName_LG),
                rdf.Max(channelName_LG),
            )

    return stats


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


def makeDRS2DPlots():
    hists2d_DRS_vs_TS = []
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
                # mean_value = stats[channelName]['mean']
                # hist = rdf.Histo2D((
                #    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}",
                #    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY};TS;{var} Variable",
                #    1024, 0, 1024, 200, mean_value - 100, mean_value + 100),
                #    "TS", chan.GetChannelName()
                # )
                hist_subtractMedian = rdf.Histo2D((
                    f"hist_DRS_Board{boardNo}_{var}_vs_TS_{sTowerX}_{sTowerY}_subtractMedian",
                    f"DRS Board {boardNo} - {var} {chan.channelNo} in iTowerX {sTowerX} iTowerY {sTowerY} (subtract median);TS;{var} Variable",
                    1024, 0, 1024, 400, -200, 600),
                    "TS", channelName + "_subtractMedian"
                )
                # hists2d_DRS_vs_TS.append(hist)
                hists2d_DRS_vs_TS.append(hist_subtractMedian)
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


def compareDRSChannels(channels_to_compare):
    hists_trigger = []
    for chan_name in channels_to_compare:
        # hist = rdf.Histo2D((
        #    f"hist_{chan_name}",
        #    f"{chan_name};TS;DRS values",
        #    1024, 0, 1024,
        #    200, 500, 2500),
        #    "TS", chan_name
        # )
        hist_subtractMedian = rdf.Histo2D((
            f"hist_{chan_name}_subtractMedian",
            f"{chan_name} (subtract median);TS;DRS values",
            1024, 0, 1024,
            300, -2500, 500),
            "TS", chan_name + "_subtractMedian"
        )
        # hists_trigger.append(hist)
        hists_trigger.append(hist_subtractMedian)
    return hists_trigger

def checkFERSvsDRSSum():
    FERS_min = 100
    FERS_max = 9e3
    FERS_LG_max = 2e3
    DRS_min = -100
    DRS_max = 1e3
    DRS_LG_max = 2e3
    
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
                    100, DRS_min, DRS_max, 100, FERS_min, FERS_max
                ),
                    f"{chan_DRS.GetChannelName()}_sum",
                    chan_FERS.GetHGChannelName(),
                )
                h2s_FERS_VS_DRS.append(h2_FERS_VS_DRS)

                h2_FERSLG_VS_DRS = rdf.Histo2D((
                    f"hist_FERSLG_VS_DRS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}",
                    f"FERS LG vs DRS energy correlation for Board{boardNo}, Tower({sTowerX}, {sTowerY}), {var}",
                    100, DRS_min, DRS_LG_max, 100, FERS_min, FERS_LG_max
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
    return  h2s_FERS_VS_DRS_sum + h2s_FERSLG_VS_DRS_sum + h2s_FERS_VS_DRS + h2s_FERSLG_VS_DRS


if __name__ == "__main__":
    start_time = time.time()

    hists_conditions = monitorConditions()

    hists_FERS_EnergySum = makeFERSEnergySumPlots()

    hists1d_FERS = makeFERS1DPlots()
    # hists2d_FERS = makeFERS2DPlots()
    # hists2d_FERS_vs_Event = trackFERSPlots()

    # hists1d_DRS = makeDRS1DPlots()
    hists2d_DRS_vs_TS = makeDRS2DPlots()
    # hists2d_DRS_vs_Event = trackDRSPlots()

    time_reference_channels = buildTimeReferenceChannels(run=runNumber)
    hists2d_time_reference = compareDRSChannels(time_reference_channels)

    hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
    hists2d_hodo_trigger = compareDRSChannels(hodo_trigger_channels)

    hodo_pos_channels = buildHodoPosChannels(run=runNumber)
    channels = [channel for channels in hodo_pos_channels.values()
                for channel in channels]
    hists2d_hodo_pos = compareDRSChannels(channels)
    
    hists2d_FERS_vs_DRSs = checkFERSvsDRSSum()

    stats = collectFERSStats()

    print("\033[94mSave results\033[0m")

    rootdir = f"root/Run{runNumber}"
    if not os.path.exists(rootdir):
        os.makedirs(rootdir)

    # dump stats into a json file
    import json
    stats_results = {}
    for channelName, (mean, max_value) in stats.items():
        stats_results[channelName] = (mean.GetValue(), max_value.GetValue())
    with open(f"{rootdir}/fers_stats.json", "w") as f:
        json.dump(stats_results, f, indent=4)

    # Save histograms to an output ROOT file
    outfile = ROOT.TFile(f"{rootdir}/conditions_vs_event.root", "RECREATE")
    for hist in hists_conditions:
        hist.SetDirectory(outfile)
        hist.Write()
    outfile.Close()

    outfile = ROOT.TFile(f"{rootdir}/fers_energy_sum.root", "RECREATE")
    for hist in hists_FERS_EnergySum:
        hist.SetDirectory(outfile)
        hist.Write()
    outfile.Close()

    outfile = ROOT.TFile(f"{rootdir}/fers_all_channels_1D.root", "RECREATE")
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
    outfile_DRS = ROOT.TFile(f"{rootdir}/drs_all_channels_2D.root", "RECREATE")
    for hist in hists2d_DRS_vs_TS:
        hist.SetDirectory(outfile_DRS)
        hist.Write()
    outfile_DRS.Close()
    # outfile_DRS = ROOT.TFile(
    #    f"{rootdir}/drs_all_channels_2D_vs_event.root", "RECREATE")
    # for hist in hists2d_DRS_vs_Event:
    #    hist.SetDirectory(outfile_DRS)
    #    hist.Write()
    # outfile_DRS.Close()

    outfile_time_reference = ROOT.TFile(
        f"{rootdir}/time_reference_channels.root", "RECREATE")
    for hist in hists2d_time_reference:
        hist.SetDirectory(outfile_time_reference)
        hist.Write()
    outfile_time_reference.Close()

    outfile_hodo_trigger = ROOT.TFile(
        f"{rootdir}/hodo_trigger_channels.root", "RECREATE")
    for hist in hists2d_hodo_trigger:
        hist.SetDirectory(outfile_hodo_trigger)
        hist.Write()
    outfile_hodo_trigger.Close()

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

    time_taken = time.time() - start_time
    print(f"Finished running script in {time_taken:.2f} seconds")
