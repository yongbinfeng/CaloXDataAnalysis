import sys
sys.path.append("CMSPLOTS")  # noqa
import ROOT
from myFunction import DrawHistos, LHistos2Hist
from channels.channel_map import buildDRSBoards, buildFERSBoards, buildTimeReferenceChannels, buildHodoTriggerChannels, buildHodoPosChannels, getUpstreamVetoChannel, getDownStreamMuonChannel, getServiceDRSChannels
from utils.utils import number2string, round_up_to_1eN
from utils.html_generator import generate_html
from utils.visualization import visualizeFERSBoards
from channels.validateMap import DrawFERSBoards, DrawDRSBoards
from utils.colors import colors
from configs.plotranges import getDRSPlotRanges, getServiceDRSPlotRanges
from utils.parser import get_args
from utils.timing import auto_timer
auto_timer("Total Execution Time")

print("Start running script")
ROOT.gROOT.SetBatch(True)

args = get_args()
runNumber = args.run

DRSBoards = buildDRSBoards(run=runNumber)
fersboards = buildFERSBoards(run=runNumber)
time_reference_channels = buildTimeReferenceChannels(run=runNumber)
hodo_trigger_channels = buildHodoTriggerChannels(run=runNumber)
hodo_pos_channels = buildHodoPosChannels(run=runNumber)
upstream_veto_channel = getUpstreamVetoChannel(run=runNumber)
downstream_muon_channel = getDownStreamMuonChannel(run=runNumber)
service_drs_channels = getServiceDRSChannels(run=runNumber)


rootdir = f"results/root/Run{runNumber}/"
plotdir = f"results/plots/Run{runNumber}/"
htmldir = f"results/html/Run{runNumber}/"


def makeConditionsPlots():
    plots = []
    outdir_plots = f"{plotdir}/Conditions_VS_Event"
    infile_name = f"{rootdir}/conditions_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofiles_SipmHV = []
    hprofiles_SipmI = []
    hprofiles_TempDET = []
    hprofiles_TempFPGA = []

    legends = []

    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        hprof_SipmHV_name = f"hprof_{fersboard.GetSipmHVName()}_VS_Event"
        hprof_SipmI_name = f"hprof_{fersboard.GetSipmIName()}_VS_Event"
        hprof_TempDET_name = f"hprof_{fersboard.GetTempDETName()}_VS_Event"
        hprof_TempFPGA_name = f"hprof_{fersboard.GetTempFPGAName()}_VS_Event"

        hprof_SipmHV = infile.Get(hprof_SipmHV_name)
        hprof_SipmI = infile.Get(hprof_SipmI_name)
        hprof_TempDET = infile.Get(hprof_TempDET_name)
        hprof_TempFPGA = infile.Get(hprof_TempFPGA_name)

        if not (hprof_SipmHV and hprof_SipmI and hprof_TempDET and hprof_TempFPGA):
            print(
                f"Warning: Profiles {hprof_SipmHV_name}, {hprof_SipmI_name}, {hprof_TempDET_name}, or {hprof_TempFPGA_name} not found in {infile_name}")
            continue

        hprofiles_SipmHV.append(hprof_SipmHV)
        hprofiles_SipmI.append(hprof_SipmI)
        hprofiles_TempDET.append(hprof_TempDET)
        hprofiles_TempFPGA.append(hprof_TempFPGA)

        legends.append(str(boardNo))

    # Draw the profiles
    nEvents = hprofiles_SipmHV[0].GetXaxis().GetXmax()
    legendPos = [0.3, 0.7, 0.9, 0.9]
    output_name = "Conditions_SipmHV_VS_Event"
    DrawHistos(hprofiles_SipmHV, legends, 0, nEvents, "Event", 26, 30, "Voltage (V)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(0, output_name + ".png")

    output_name = "Conditions_SipmI_VS_Event"
    DrawHistos(hprofiles_SipmI, legends, 0, nEvents, "Event", 0.0, 0.3, "Current (mA)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(1, output_name + ".png")

    output_name = "Conditions_TempDET_VS_Event"
    DrawHistos(hprofiles_TempDET, legends, 0, nEvents, "Event", 14, 35, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(2, output_name + ".png")

    output_name = "Conditions_TempFPGA_VS_Event"
    DrawHistos(hprofiles_TempFPGA, legends, 0, nEvents, "Event", 32, 50, "Temperature (C)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(3, output_name + ".png")

    output_html = f"{htmldir}/Conditions_VS_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeFERSSumPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_EnergySum_VS_Event"
    infile_name = f"{rootdir}/fers_energysum_vs_event.root"
    infile = ROOT.TFile(infile_name, "READ")

    hprofiles_Cer_HG_sum = []
    hprofiles_Sci_HG_sum = []
    hprofiles_Cer_LG_sum = []
    hprofiles_Sci_LG_sum = []

    legends = []

    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        hprof_Cer_HG_sum_name = f"hprof_{fersboard.GetEnergySumName(gain='HG', isCer=True)}_VS_Event"
        hprof_Sci_HG_sum_name = f"hprof_{fersboard.GetEnergySumName(gain='HG', isCer=False)}_VS_Event"
        hprof_Cer_LG_sum_name = f"hprof_{fersboard.GetEnergySumName(gain='LG', isCer=True)}_VS_Event"
        hprof_Sci_LG_sum_name = f"hprof_{fersboard.GetEnergySumName(gain='LG', isCer=False)}_VS_Event"

        hprof_Cer_HG_sum = infile.Get(hprof_Cer_HG_sum_name)
        hprof_Sci_HG_sum = infile.Get(hprof_Sci_HG_sum_name)
        hprof_Cer_LG_sum = infile.Get(hprof_Cer_LG_sum_name)
        hprof_Sci_LG_sum = infile.Get(hprof_Sci_LG_sum_name)

        if not (hprof_Cer_HG_sum and hprof_Sci_HG_sum and hprof_Cer_LG_sum and hprof_Sci_LG_sum):
            print(
                f"Warning: Profiles {hprof_Cer_HG_sum_name}, {hprof_Sci_HG_sum_name}, {hprof_Cer_LG_sum_name}, or {hprof_Sci_LG_sum_name} not found in {infile_name}")
            continue

        hprofiles_Cer_HG_sum.append(hprof_Cer_HG_sum)
        hprofiles_Sci_HG_sum.append(hprof_Sci_HG_sum)
        hprofiles_Cer_LG_sum.append(hprof_Cer_LG_sum)
        hprofiles_Sci_LG_sum.append(hprof_Sci_LG_sum)

        legends.append(str(boardNo))

    nEvents = hprofiles_Cer_HG_sum[0].GetXaxis().GetXmax()
    legendPos = [0.3, 0.7, 0.9, 0.9]
    output_name = "FERS_Cer_HG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Cer_HG_sum, legends, 0, nEvents, "Event", 0, 5e4, "Cer FERS Sum HG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(0, output_name + ".png")

    output_name = "FERS_Sci_HG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Sci_HG_sum, legends, 0, nEvents, "Event", 0, 1.6e5, "Sci FERS Sum HG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(1, output_name + ".png")

    output_name = "FERS_Cer_LG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Cer_LG_sum, legends, 0, nEvents, "Event", 0, 1.4e4, "Cer FERS Sum LG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(2, output_name + ".png")

    output_name = "FERS_Sci_LG_EnergySum_VS_Event"
    DrawHistos(hprofiles_Sci_LG_sum, legends, 0, nEvents, "Event", 0, 5e4, "Sci FERS Sum LG (ADC)",
               output_name,
               dology=False, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=3, legendPos=legendPos)
    plots.insert(3, output_name + ".png")

    # total sum
    hprof_Cer_HG_sum = infile.Get(
        f"hprof_{fersboards.GetEnergySumName(gain='HG', isCer=True)}_VS_Event")
    hprof_Sci_HG_sum = infile.Get(
        f"hprof_{fersboards.GetEnergySumName(gain='HG', isCer=False)}_VS_Event")
    hprof_Cer_LG_sum = infile.Get(
        f"hprof_{fersboards.GetEnergySumName(gain='LG', isCer=True)}_VS_Event")
    hprof_Sci_LG_sum = infile.Get(
        f"hprof_{fersboards.GetEnergySumName(gain='LG', isCer=False)}_VS_Event")

    DrawHistos([hprof_Cer_HG_sum, hprof_Sci_HG_sum], ["Cer", "Sci"], 0, nEvents, "Event", 0, 7e5, "FERS Total Sum HG (ADC)",
               "FERS_Total_HG_EnergySum_VS_Event",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=2, legendPos=[0.6, 0.85, 0.9, 0.9])
    plots.insert(0, "FERS_Total_HG_EnergySum_VS_Event.png")
    DrawHistos([hprof_Cer_LG_sum, hprof_Sci_LG_sum], ["Cer", "Sci"], 0, nEvents, "Event", 0, 2e5, "FERS Total Sum LG (ADC)",
               "FERS_Total_LG_EnergySum_VS_Event",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=2, legendPos=[0.6, 0.85, 0.9, 0.9])
    plots.insert(1, "FERS_Total_LG_EnergySum_VS_Event.png")

    output_html = f"{htmldir}/FERS_EnergySum_VS_Event/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html, plots_per_row=6)

    return output_html


def makeFERS1DPlots():
    plots = []

    infile_name = f"{rootdir}fers_all_channels_1d.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/FERS_1D"
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        for iTowerX, iTowerY in fersboard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hist_C_name = f"hist_FERS_Board{boardNo}_Cer_{sTowerX}_{sTowerY}"
            hist_S_name = f"hist_FERS_Board{boardNo}_Sci_{sTowerX}_{sTowerY}"
            hist_C = infile.Get(hist_C_name)
            hist_S = infile.Get(hist_S_name)
            if not hist_C or not hist_S:
                print(
                    f"Warning: Histograms {hist_C_name} or {hist_S_name} not found in {infile_name}")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"Board: {fersboard.boardNo}")
            extraToDraw.AddText(f"Tower X: {iTowerX}")
            extraToDraw.AddText(f"Tower Y: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: {fersboard.GetChannelByTower(iTowerX, iTowerY, isCer=True).channelNo}")
            extraToDraw.AddText(
                f"Sci Channel: {fersboard.GetChannelByTower(iTowerX, iTowerY, isCer=False).channelNo}")

            output_name = f"Energy_Board{boardNo}_iTowerX{sTowerX}_iTowerY{sTowerY}"
            DrawHistos([hist_C, hist_S], ["Cer", "Sci"], 0, 1000, "Energy HG", 1, 1e5, "Counts",
                       output_name,
                       dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)

            plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_1D/index.html"
    generate_html(plots, outdir_plots,
                  output_html=output_html)
    return output_html


def makeFERSStatsPlots(includePedestals=False):
    plots = []
    outdir_plots = f"{plotdir}/FERS_Stats"
    # load the json file
    import json
    infile_name = f"{rootdir}/fers_stats.json"
    with open(infile_name, "r") as f:
        stats = json.load(f)

    if includePedestals:
        infile_name_HG = f"{rootdir}/fers_pedestals_hg.json"
        infile_name_LG = f"{rootdir}/fers_pedestals_lg.json"
        with open(infile_name_HG, "r") as f:
            pedestals_HG = json.load(f)
        with open(infile_name_LG, "r") as f:
            pedestals_LG = json.load(f)

    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100
    valuemaps_HG_mean = {}
    valuemaps_HG_max = {}
    valuemaps_HG_satfreq = {}
    valuemaps_HG_pedestal = {}
    valuemaps_LG_mean = {}
    valuemaps_LG_max = {}
    valuemaps_LG_satfreq = {}
    valuemaps_LG_pedestal = {}

    for channelName, (vmean, vmax, VSatfreq) in stats.items():
        if "energyHG" in channelName:
            valuemaps_HG_mean[channelName] = vmean
            valuemaps_HG_max[channelName] = vmax
            valuemaps_HG_satfreq[channelName] = VSatfreq
            valuemaps_HG_pedestal[channelName] = pedestals_HG.get(
                channelName, None) if includePedestals else 0.
        elif "energyLG" in channelName:
            valuemaps_LG_mean[channelName] = vmean
            valuemaps_LG_max[channelName] = vmax
            valuemaps_LG_satfreq[channelName] = VSatfreq
            valuemaps_LG_pedestal[channelName] = pedestals_LG.get(
                channelName, None) if includePedestals else 0.

    [h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], [h2_Sci_HG_mean, h2_Sci_3mm_HG_mean] = visualizeFERSBoards(
        fersboards, valuemaps_HG_mean, suffix=f"Run{runNumber}_HG_mean", gain="HG")
    [h2_Cer_HG_max, h2_Cer_3mm_HG_max], [h2_Sci_HG_max, h2_Sci_3mm_HG_max] = visualizeFERSBoards(
        fersboards, valuemaps_HG_max, suffix=f"Run{runNumber}_HG_max", gain="HG")
    [h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], [h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq] = visualizeFERSBoards(
        fersboards, valuemaps_HG_satfreq, suffix=f"Run{runNumber}_HG_satfreq", gain="HG")
    [h2_Cer_HG_pedestal, h2_Cer_3mm_HG_pedestal], [h2_Sci_HG_pedestal, h2_Sci_3mm_HG_pedestal] = visualizeFERSBoards(
        fersboards, valuemaps_HG_pedestal, suffix=f"Run{runNumber}_HG_pedestal", gain="HG")
    [h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], [h2_Sci_LG_mean, h2_Sci_3mm_LG_mean] = visualizeFERSBoards(
        fersboards, valuemaps_LG_mean, suffix=f"Run{runNumber}_LG_mean", gain="LG")
    [h2_Cer_LG_max, h2_Cer_3mm_LG_max], [h2_Sci_LG_max, h2_Sci_3mm_LG_max] = visualizeFERSBoards(
        fersboards, valuemaps_LG_max, suffix=f"Run{runNumber}_LG_max", gain="LG")
    [h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], [h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq] = visualizeFERSBoards(
        fersboards, valuemaps_LG_satfreq, suffix=f"Run{runNumber}_LG_satfreq", gain="LG")
    [h2_Cer_LG_pedestal, h2_Cer_3mm_LG_pedestal], [h2_Sci_LG_pedestal, h2_Sci_3mm_LG_pedestal] = visualizeFERSBoards(
        fersboards, valuemaps_LG_pedestal, suffix=f"Run{runNumber}_LG_pedestal", gain="LG")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_mean"
    DrawHistos([h2_Cer_HG_mean, h2_Cer_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_mean, h2_Sci_3mm_HG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_max"
    DrawHistos([h2_Cer_HG_max, h2_Cer_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_max, h2_Sci_3mm_HG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_satfreq"
    DrawHistos([h2_Cer_HG_satfreq, h2_Cer_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_satfreq, h2_Sci_3mm_HG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_HG_pedestal"
    DrawHistos([h2_Cer_HG_pedestal, h2_Cer_3mm_HG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_HG_pedestal, h2_Sci_3mm_HG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_mean"
    DrawHistos([h2_Cer_LG_mean, h2_Cer_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_mean, h2_Sci_3mm_LG_mean], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_max"
    DrawHistos([h2_Cer_LG_max, h2_Cer_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_max, h2_Sci_3mm_LG_max], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=8000)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_satfreq"
    DrawHistos([h2_Cer_LG_satfreq, h2_Cer_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_satfreq, h2_Sci_3mm_LG_satfreq], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=0, zmax=1, nTextDigits=2)
    plots.append(output_name + "_Sci.png")

    output_name = f"FERS_Boards_Run{runNumber}_Stats_LG_pedestal"
    DrawHistos([h2_Cer_LG_pedestal, h2_Cer_3mm_LG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=runNumber, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Cer.png")
    DrawHistos([h2_Sci_LG_pedestal, h2_Sci_3mm_LG_pedestal], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir_plots, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=runNumber, zmin=100, zmax=300, nTextDigits=0)
    plots.append(output_name + "_Sci.png")

    output_html = f"{htmldir}/FERS_Stats/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def makeFERSMaxValuePlots():
    plots = []
    xmin = 7500
    xmax = 8500
    outdir_plots = f"{plotdir}/FERS_MaxValues"
    infile_name = f"{rootdir}/fers_max_values.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_board_cer_HG_max = []
    hists_board_cer_LG_max = []
    hists_board_sci_HG_max = []
    hists_board_sci_LG_max = []
    legends = []
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        hist_board_cer_HG_max = infile.Get(
            f'hist_{fersboard.GetEnergyMaxName(gain="HG", isCer=True)}')
        hist_board_sci_HG_max = infile.Get(
            f'hist_{fersboard.GetEnergyMaxName(gain="HG", isCer=False)}')
        hist_board_cer_LG_max = infile.Get(
            f'hist_{fersboard.GetEnergyMaxName(gain="LG", isCer=True)}')
        hist_board_sci_LG_max = infile.Get(
            f'hist_{fersboard.GetEnergyMaxName(gain="LG", isCer=False)}')
        hists_board_cer_HG_max.append(hist_board_cer_HG_max)
        hists_board_sci_HG_max.append(hist_board_sci_HG_max)
        hists_board_cer_LG_max.append(hist_board_cer_LG_max)
        hists_board_sci_LG_max.append(hist_board_sci_LG_max)
        legends.append(str(boardNo))

    hist_cer_HG_max = infile.Get(
        f'hist_{fersboards.GetEnergyMaxName(gain="HG", isCer=True)}')
    hist_sci_HG_max = infile.Get(
        f'hist_{fersboards.GetEnergyMaxName(gain="HG", isCer=False)}')
    hist_cer_LG_max = infile.Get(
        f'hist_{fersboards.GetEnergyMaxName(gain="LG", isCer=True)}')
    hist_sci_LG_max = infile.Get(
        f'hist_{fersboards.GetEnergyMaxName(gain="LG", isCer=False)}')

    output_name = "FERS_Boards_CerEnergyHG_max"
    DrawHistos(hists_board_cer_HG_max, legends, xmin, xmax, f"HG Cer Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_SciEnergyHG_max"
    DrawHistos(hists_board_sci_HG_max, legends, xmin, xmax, f"HG Sci Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_CerEnergyLG_max"
    DrawHistos(hists_board_cer_LG_max, legends, xmin, xmax, f"LG Cer Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")
    output_name = "FERS_Boards_SciEnergyLG_max"
    DrawHistos(hists_board_sci_LG_max, legends, xmin, xmax, f"LG Sci Max (Board)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=colors, addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendNCols=5, legendPos=[0.20, 0.75, 0.90, 0.9])
    plots.append(output_name + ".png")

    output_name = "FERS_EnergyHG_max"
    VSat = 8000
    frac_sci = hist_sci_HG_max.Integral(hist_sci_HG_max.FindBin(VSat), 100000) / \
        hist_sci_HG_max.Integral(0, 100000)
    frac_cer = hist_cer_HG_max.Integral(hist_cer_HG_max.FindBin(VSat), 100000) / \
        hist_cer_HG_max.Integral(0, 100000)
    extraToDraw = ROOT.TPaveText(0.20, 0.63, 0.90, 0.72, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
    extraToDraw.AddText(f"Sat Frac Cer : {frac_cer:.3f}")
    DrawHistos([hist_cer_HG_max, hist_sci_HG_max],
               ["Cer",
                   "Sci"], xmin, xmax, f"HG Max (All Boards)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.30, 0.75, 0.50, 0.9], extraToDraw=extraToDraw)
    plots.append(output_name + ".png")
    output_name = "FERS_EnergyLG_max"
    frac_sci = hist_sci_LG_max.Integral(hist_sci_LG_max.FindBin(VSat), 100000) / \
        hist_sci_LG_max.Integral(0, 100000)
    frac_cer = hist_cer_LG_max.Integral(hist_cer_LG_max.FindBin(VSat), 100000) / \
        hist_cer_LG_max.Integral(0, 100000)
    extraToDraw = ROOT.TPaveText(0.20, 0.55, 0.90, 0.63, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Sat Frac Sci : {frac_sci:.3f}")
    extraToDraw.AddText(f"Sat Frac Cer : {frac_cer:.3f}")
    DrawHistos([hist_cer_LG_max, hist_sci_LG_max],
               ["Cer",
                   "Sci"], xmin, xmax, f"LG Max (All Boards)", 1, None, "Events",
               output_name,
               dology=True, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=True,
               outdir=outdir_plots, runNumber=runNumber, legendPos=[0.30, 0.75, 0.50, 0.9], extraToDraw=extraToDraw)
    plots.append(output_name + ".png")

    output_html = f"{htmldir}/FERS_MaxValues/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html

# 2D FERS histograms, hg VS lg


def makeFERS2DPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_2D"
    infile_name = f"{rootdir}/fers_all_channels_2d.root"
    infile = ROOT.TFile(infile_name, "READ")
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        for iTowerX, iTowerY in fersboard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = fersboard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_VS_lg"
                hist = infile.Get(hist_name)

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {fersboard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_hg_VS_lg"
                DrawHistos([hist], f"", 0, 9000, "HG", 0, 1500, "LG",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True, extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERS_2D/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# FERS output VS event
def trackFERSPlots():
    plots = []
    outdir_plots = f"{plotdir}/FERS_VS_Event"
    infile_name = f"{rootdir}/fers_all_channels_2D_VS_event.root"
    infile = ROOT.TFile(infile_name, "READ")
    for fersboard in fersboards.values():
        boardNo = fersboard.boardNo
        for iTowerX, iTowerY in fersboard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = fersboard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_FERS_Board{boardNo}_{var}_VS_Event_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"Board: {fersboard.boardNo}")
                extraToDraw.AddText(f"Tower X: {iTowerX}")
                extraToDraw.AddText(f"Tower Y: {iTowerY}")
                extraToDraw.AddText(f"{var} Channel: {chan.channelNo}")

                nEvents = hist.GetXaxis().GetXmax()

                output_name = f"FERS_Board{boardNo}_{var}_{sTowerX}_{sTowerY}_VS_Event"
                DrawHistos([hist], "", 0, nEvents, "Event", 1, 1e5, f"{var} Energy HG",
                           output_name,
                           dology=True, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
                plots.append(output_name + ".png")
    output_html = f"{htmldir}/FERS_VS_Event/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# DRS VS TS
def makeDRSVSTSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRS_VS_TS"
    infile_name = f"{rootdir}/drs_vs_ts.root"
    infile = ROOT.TFile(infile_name, "READ")
    for _, DRSBoard in DRSBoards.items():
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if chan is None:
                    print(
                        f"Warning: No channel found for Board {DRSBoard.boardNo}, Tower ({iTowerX}, {iTowerY}), Var {var}")
                    continue
                channelName = chan.GetChannelName(blsub=True)
                hist_name = f"hist_{channelName}_VS_TS"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    continue

                output_name = f"DRS_{var}_VS_TS_{sTowerX}_{sTowerY}"
                plots.append(output_name + ".png")

                ymin_tmp, ymax_tmp = getDRSPlotRanges(
                    subtractMedian=True, isAmplified=chan.isAmplified)

                extraToDraw = ROOT.TPaveText(0.20, 0.75, 0.60, 0.90, "NDC")
                extraToDraw.SetTextAlign(11)
                extraToDraw.SetFillColorAlpha(0, 0)
                extraToDraw.SetBorderSize(0)
                extraToDraw.SetTextFont(42)
                extraToDraw.SetTextSize(0.04)
                extraToDraw.AddText(
                    f"B: {DRSBoard.boardNo}, G: {chan.groupNo}, C: {chan.channelNo}")
                extraToDraw.AddText(f"iTowerX: {iTowerX}")
                extraToDraw.AddText(f"iTowerY: {iTowerY}")

                DrawHistos([hist], "", 0, 1024, "Time Slice", ymin_tmp, ymax_tmp, f"DRS Output",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                           extraToDraw=extraToDraw,
                           outdir=outdir_plots, extraText=var, runNumber=runNumber, addOverflow=True)
    output_html = f"{htmldir}/DRS_VS_TS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeakTS"
    infile_name = f"{rootdir}/drspeakts.root"
    infile = ROOT.TFile(infile_name, "READ")
    hists_Cer = []
    hists_Sci = []
    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hists = {}
            channelNos = {}
            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                hist_name = f"hist_DRSPeakTS_{var}_{sTowerX}_{sTowerY}"
                hist = infile.Get(hist_name)

                if not hist:
                    print(
                        f"Warning: Histogram {hist_name} not found in {infile_name}")
                    hists[var] = None
                    channelNos[var] = -1
                else:
                    hists[var] = hist
                    channelNos[var] = chan.channelNo

                output_name = f"hist_DRSPeakTS_{sTowerX}_{sTowerY}"

            if not hists["Cer"] or not hists["Sci"]:
                print(
                    f"Warning: Histograms for Cer or Sci not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY})")
                continue

            extraToDraw = ROOT.TPaveText(0.20, 0.65, 0.60, 0.90, "NDC")
            extraToDraw.SetTextAlign(11)
            extraToDraw.SetFillColorAlpha(0, 0)
            extraToDraw.SetBorderSize(0)
            extraToDraw.SetTextFont(42)
            extraToDraw.SetTextSize(0.04)
            extraToDraw.AddText(
                f"B: {DRSBoard.boardNo}, G: {chan.groupNo}")
            extraToDraw.AddText(f"Tower X: {iTowerX}")
            extraToDraw.AddText(f"Tower Y: {iTowerY}")
            extraToDraw.AddText(
                f"Cer Channel: {channelNos['Cer']}")
            extraToDraw.AddText(
                f"Sci Channel: {channelNos['Sci']}")

            hists_Cer.append(hists["Cer"])
            hists_Sci.append(hists["Sci"])

            DrawHistos([hists["Cer"], hists["Sci"]], ["Cer", "Sci"], 400, 600, "Peak TS", 1, None, "Counts",
                       output_name,
                       dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False, extraToDraw=extraToDraw,
                       outdir=outdir_plots, runNumber=runNumber)
            plots.append(output_name + ".png")

    # summary plots
    hist_Cer_Combined = LHistos2Hist(hists_Cer, "hist_DRSPeakTS_Cer_Combined")
    hist_Sci_Combined = LHistos2Hist(hists_Sci, "hist_DRSPeakTS_Sci_Combined")
    DrawHistos([hist_Cer_Combined, hist_Sci_Combined], ["Cer", "Sci"], 400, 600, "Peak TS", 1, None, "Counts",
               "DRS_PeakTS_Combined",
               dology=False, drawoptions="HIST", mycolors=[2, 4], addOverflow=True, addUnderflow=False,
               outdir=outdir_plots, runNumber=runNumber)
    plots.insert(0, "DRS_PeakTS_Combined.png")

    output_html = f"{htmldir}/DRSPeakTS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakTSCerVSSciPlots():
    plots = []
    hists = []
    outdir_plots = f"{plotdir}/DRSPeakTSCerVSSci"
    infile_name = f"{rootdir}/drspeakts.root"
    infile = ROOT.TFile(infile_name, "READ")

    # Create a dashed diagonal line from (0,0) to (1000,1000)
    diagonal_line = ROOT.TLine(0, 0, 1000, 1000)
    diagonal_line.SetLineStyle(2)  # 2 = dashed
    diagonal_line.SetLineWidth(1)
    diagonal_line.SetLineColor(ROOT.kRed)
    extraToDraw = diagonal_line

    for _, DRSBoard in DRSBoards.items():
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            hist_name = f"hist_DRSPeakTS_Cer_VS_Sci_{sTowerX}_{sTowerY}"
            hist = infile.Get(hist_name)
            output_name = hist_name[5:]

            if not hist:
                print(
                    f"Warning: Histogram {hist_name} not found in {infile_name}")
                continue

            hists.append(hist)

            DrawHistos([hist], "", 400, 600, "Sci Peak TS", 400, 600, f"Cer Peak TS",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
                       outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)

            plots.append(output_name + ".png")

    # summary plots
    hcombined = LHistos2Hist(hists, "hist_DRSPeakTS_Cer_VS_Sci_Combined")
    output_name = "DRS_PeakTS_Cer_VS_Sci_Combined"
    DrawHistos([hcombined], "", 400, 600, "Sci Peak TS", 400, 600, f"Cer Peak TS",
               output_name,
               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e2, dologz=True,
               outdir=outdir_plots, addOverflow=False, runNumber=runNumber, extraToDraw=extraToDraw)
    plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/DRSPeakTSCerVSSci/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


# time reference
def compareTimeReferencePlots(doSubtractMedian=False):
    suffix = ""
    ymin = 500
    ymax = 2500
    if doSubtractMedian:
        suffix = "_subtractMedian"
        ymin = -2500
        ymax = 500
    plots = []
    infile_name = f"{rootdir}/time_reference_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/TimeReference"
    for chan_name in time_reference_channels:
        hist_name = f"hist_{chan_name}{suffix}"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"TimeReference_{chan_name}{suffix}"
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/TimeReference{suffix}/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)
    return output_html


def compareServiceDRSPlots():
    ymin = 500
    ymax = 2500
    plots = []
    infile_name = f"{rootdir}/service_drs_channels.root"
    infile = ROOT.TFile(infile_name, "READ")
    outdir_plots = f"{plotdir}/ServiceDRS"

    for chan_name in service_drs_channels:
        hist_name = f"hist_{chan_name}_blsub"
        hist = infile.Get(hist_name)
        if not hist:
            print(f"Warning: Histogram {hist_name} not found in {infile_name}")
            continue
        extraToDraw = ROOT.TPaveText(0.20, 0.70, 0.60, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"{chan_name}")
        output_name = f"ServiceDRS_{chan_name}"

        ymin, ymax = getServiceDRSPlotRanges(
            chan_name, subtractMedian=True)
        DrawHistos([hist], "", 0, 1024, "Time Slice", ymin, ymax, "Counts",
                   output_name,
                   dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=1e4, dologz=True,
                   extraToDraw=extraToDraw,
                   outdir=outdir_plots, addOverflow=True, runNumber=runNumber)
        plots.append(output_name + ".png")

    output_html = f"{htmldir}/ServiceDRS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=2,
                  output_html=output_html)

    return output_html


def makeDRSSumVSFERSPlots():
    """
    Check if the sum of FERS and DRS energies are consistent.
    """
    plots = []
    outdir_plots = f"{plotdir}/DRSSum_VS_FERS"
    infile_name = f"{rootdir}/drssum_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")

    xymax = {
        "Cer": (20000, 8500),
        "Sci": (30000, 8500)
    }
    xymax_LG = {
        "Cer": (20000, 2000),
        "Sci": (30000, 4000)
    }

    hists = []
    for _, DRSBoard in DRSBoards.items():
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                for gain in ["FERS", "FERSLG"]:
                    histname = f"hist_DRSSum_VS_{gain}_{var}_{sTowerX}_{sTowerY}"
                    hist = infile.Get(histname)
                    if not hist:
                        print(
                            f"Warning: Histogram {histname} not found in {infile_name}")
                        continue

                    zmax = hist.Integral(0, 10000, 0, 10000)
                    zmax = round_up_to_1eN(zmax)

                    tmp = xymax[var] if gain == "FERS" else xymax_LG[var]

                    hists.append(hist)

                    output_name = histname.replace("hist_", "")
                    DrawHistos([hist], "", 0, tmp[1], gain, 0, tmp[0], "DRSSum",
                               output_name,
                               dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
                               outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")
                    plots.append(output_name + ".png")

    # summary plots
    for var in ["Cer", "Sci"]:
        for gain in ["FERS", "FERSLG"]:
            histname = f"hist_DRSSum_VS_{gain}_{var}_Combined"
            hists_to_combine = [
                h for h in hists if f"_{gain}_{var}_" in h.GetName()]
            hist_combined = LHistos2Hist(hists_to_combine, histname)
            print("Combining histograms for", histname, " is ",
                  hist_combined, " len ", len(hists_to_combine))
            output_name = histname.replace("hist_", "")

            zmax = hist_combined.Integral(0, 10000, 0, 10000)
            zmax = round_up_to_1eN(zmax)

            output_name = histname.replace("hist_", "")

            tmp = xymax[var] if gain == "FERS" else xymax_LG[var]
            DrawHistos([hist_combined], "", 0, tmp[1], gain, 0, tmp[0], "DRSSum",
                       output_name,
                       dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=zmax, dologz=True,
                       outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")
            plots.insert(0, output_name + ".png")

    output_html = f"{htmldir}/DRSSum_VS_FERS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


def makeDRSPeakVSFERSPlots():
    plots = []
    outdir_plots = f"{plotdir}/DRSPeak_VS_FERS"
    infile_name = f"{rootdir}/drspeak_vs_fers.root"
    infile = ROOT.TFile(infile_name, "READ")

    for _, DRSBoard in DRSBoards.items():
        boardNo = DRSBoard.boardNo
        for iTowerX, iTowerY in DRSBoard.GetListOfTowers():
            sTowerX = number2string(iTowerX)
            sTowerY = number2string(iTowerY)

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.GetChannelByTower(
                    iTowerX, iTowerY, isCer=(var == "Cer"))
                if not chan:
                    print(
                        f"Warning: Channel not found for Board {boardNo}, Tower ({iTowerX}, {iTowerY}), Var {var}")
                    continue
                _, ymax = getDRSPlotRanges(
                    subtractMedian=True, isAmplified=chan.isAmplified)
                histname = f"hist_DRSPeak_VS_FERS_{var}_{sTowerX}_{sTowerY}"
                output_name = histname.replace("hist_", "")
                plots.append(output_name + ".png")

                hist = infile.Get(histname)
                if not hist:
                    print(
                        f"Warning: Histogram {histname} not found in {infile_name}")
                    continue

                DrawHistos([hist], "", 0, 9000, "FERS ADC", 0, ymax, "DRS Peak",
                           output_name,
                           dology=False, drawoptions="COLZ", doth2=True, zmin=1, zmax=None, dologz=True,
                           outdir=outdir_plots, addOverflow=True, runNumber=runNumber, extraText=f"{var}")

    output_html = f"{htmldir}/DRSPeak_VS_FERS/index.html"
    generate_html(plots, outdir_plots, plots_per_row=4,
                  output_html=output_html)
    return output_html


if __name__ == "__main__":
    output_htmls = {}

    output_htmls["conditions plots"] = makeConditionsPlots()
    output_htmls["fers sum"] = makeFERSSumPlots()

    # validate DRS and FERS boards
    output_htmls["fers mapping"] = DrawFERSBoards(run=runNumber)
    output_htmls["drs mapping"] = DrawDRSBoards(run=runNumber)

    output_htmls["fers 1D"] = makeFERS1DPlots()
    output_htmls["fers stats"] = makeFERSStatsPlots(includePedestals=True)

    output_htmls["drs vs ts"] = makeDRSVSTSPlots()
    output_htmls["drs peak ts"] = makeDRSPeakTSPlots()

    output_htmls["drs peak ts cer vs ts"] = makeDRSPeakTSCerVSSciPlots()

    output_htmls["drs services"] = compareServiceDRSPlots()

    # output_htmls["time reference"] = compareTimeReferencePlots(True)

    output_htmls["fers max values"] = makeFERSMaxValuePlots()

    output_htmls["drs sum vs fers"] = makeDRSSumVSFERSPlots()

    output_htmls["drs peak VS fers"] = makeDRSPeakVSFERSPlots()

    print("\n\n\n")
    print("*" * 30)
    for key, value in output_htmls.items():
        print(f"✅ {key} plots can be viewed at: {value}")

    print("All plots generated successfully.")
