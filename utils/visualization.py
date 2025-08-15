import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos


def visualizeFERSBoards(fers_boards, valuemaps=None, suffix="", useHG=True):
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    h2_Cer = ROOT.TH2D(f"h_FERSBoards_Cer_{suffix}", f"FERS Board Cer Channels for {suffix}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Sci = ROOT.TH2D(f"h_FERSBoards_Sci_{suffix}", f"FERS Board Sci Channels for {suffix}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Cer_3mm = ROOT.TH2D(f"h_FERSBoards_Cer_3mm_{suffix}", f"FERS Board Cer Channels for {suffix} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Sci_3mm = ROOT.TH2D(f"h_FERSBoards_Sci_3mm_{suffix}", f"FERS Board Sci Channels for {suffix} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

    for _, FERSBoard in fers_boards.items():
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            channel_Cer = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            channel_Sci = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)

            channel_Cer_name = channel_Cer.GetHGChannelName(
            ) if useHG else channel_Cer.GetLGChannelName()
            channel_Sci_name = channel_Sci.GetHGChannelName(
            ) if useHG else channel_Sci.GetLGChannelName()

            # default encoding based on board and channel numbers
            cer_encoded = channel_Cer.boardNo * 100 + channel_Cer.channelNo
            if cer_encoded == 0:
                cer_encoded = 0.001
            sci_encoded = channel_Sci.boardNo * 100 + channel_Sci.channelNo
            if sci_encoded == 0:
                sci_encoded = 0.001

            if valuemaps and channel_Cer_name in valuemaps:
                cer_encoded = valuemaps[channel_Cer_name]
            if valuemaps and channel_Sci_name in valuemaps:
                sci_encoded = valuemaps[channel_Sci_name]
            if FERSBoard.Is3mm():
                h2_Cer_3mm.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                                cer_encoded)
                h2_Sci_3mm.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                                sci_encoded)
            else:
                h2_Cer.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                            cer_encoded)
                h2_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                            sci_encoded)

    h2_Cer.SetMarkerSize(0.60)
    h2_Sci.SetMarkerSize(0.60)
    h2_Cer_3mm.SetMarkerSize(0.60)
    h2_Sci_3mm.SetMarkerSize(0.60)

    return [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm]


def visualizeDRSBoards(drs_boards, suffix=""):
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10

    h2_DRS_Cer = ROOT.TH2D(f"h_DRSBoards_Cer_{suffix}", f"DRS Board Channels for {suffix}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_DRS_Cer_3mm = ROOT.TH2D(f"h_DRSBoards_Cer_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

    h2_DRS_Sci = ROOT.TH2D(f"h_DRSBoards_Sci_{suffix}", f"DRS Board Channels for {suffix}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_DRS_Sci_3mm = ROOT.TH2D(f"h_DRSBoards_Sci_3mm_{suffix}", f"DRS Board Channels for {suffix} (3mm)",
                               int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)

    for _, DRSBoard in drs_boards.items():
        for channel_Sci in DRSBoard.GetSciChannels():
            sci_encoded = channel_Sci.boardNo * 100 + \
                channel_Sci.groupNo * 10 + channel_Sci.channelNo
            if sci_encoded == 0:
                sci_encoded = 0.001
            if not channel_Sci.is6mm:
                h2_DRS_Sci_3mm.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                                    sci_encoded)
            else:
                h2_DRS_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                                sci_encoded)
        for channel_Cer in DRSBoard.GetCerChannels():
            cer_encoded = channel_Cer.boardNo * 100 + \
                channel_Cer.groupNo * 10 + channel_Cer.channelNo
            if cer_encoded == 0:
                cer_encoded = 0.001
            if not channel_Cer.is6mm:
                h2_DRS_Cer_3mm.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                                    cer_encoded)
            else:
                h2_DRS_Cer.Fill(channel_Cer.iTowerX,
                                channel_Cer.iTowerY, cer_encoded)

    h2_DRS_Cer.SetMarkerSize(0.60)
    h2_DRS_Sci.SetMarkerSize(0.60)
    h2_DRS_Cer_3mm.SetMarkerSize(0.60)
    h2_DRS_Sci_3mm.SetMarkerSize(0.60)

    return [h2_DRS_Cer, h2_DRS_Cer_3mm], [h2_DRS_Sci, h2_DRS_Sci_3mm]


def makeEventDisplay(h2_6mm, h2_3mm, output_name, outdir, runNumber, zmin, zmax, isCer=False):
    xmax = 14
    xmin = -14
    ymax = 10
    ymin = -10
    W_ref = 1000
    H_ref = 1100
    if isCer:
        extraText = "Cer"
    else:
        extraText = "Sci"
    DrawHistos([h2_6mm, h2_3mm], "", xmin, xmax, "iX", ymin, ymax, "iY", output_name, dology=False, drawoptions=["col,text", "col,text"],
               outdir=outdir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText=extraText, runNumber=runNumber, zmin=zmin, zmax=zmax)
