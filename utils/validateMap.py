import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.channel_map import buildFERSBoards, buildDRSBoards
from utils.html_generator import generate_html

ROOT.gROOT.SetBatch(True)  # Disable interactive mode

xmax = 14
xmin = -14
ymax = 10
ymin = -10
W_ref = 1000
H_ref = 1100


def DrawFERSBoards(run=316):
    """
    Draws the FERS boards for a given run in a TH2D
    """
    fers_boards = buildFERSBoards(run)
    h2_Cer = ROOT.TH2D(f"h_FERSBoards_Run{run}_Cer", f"FERS Board Cer Channels for Run {run}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Sci = ROOT.TH2D(f"h_FERSBoards_Run{run}_Sci", f"FERS Board Sci Channels for Run {run}",
                       int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)
    h2_Cer_3mm = ROOT.TH2D(f"h_FERSBoards_Run{run}_Cer_3mm", f"FERS Board Cer Channels for Run {run} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Sci_3mm = ROOT.TH2D(f"h_FERSBoards_Run{run}_Sci_3mm", f"FERS Board Sci Channels for Run {run} (3mm)",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin) * 4, ymin, ymax)
    h2_Cer_3mm.SetMarkerSize(0.70)
    h2_Sci_3mm.SetMarkerSize(0.70)
    bboxes = []
    colors = [
        ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1,
        ROOT.kCyan+2, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kViolet+1, ROOT.kAzure+1,
        ROOT.kSpring+3, ROOT.kTeal+2, ROOT.kPink+2, ROOT.kGray+1, ROOT.kOrange+7
    ]

    for _, FERSBoard in fers_boards.items():
        towerx_min = 99
        towerx_max = -99
        towery_min = 99
        towery_max = -99
        for iTowerX, iTowerY in FERSBoard.GetListOfTowers():
            channel_Cer = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=True)
            channel_Sci = FERSBoard.GetChannelByTower(
                iTowerX, iTowerY, isCer=False)
            if FERSBoard.Is3mm():
                h2_Cer_3mm.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                                channel_Cer.channelNo)
                h2_Sci_3mm.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                                channel_Sci.channelNo)
                extra_x_side = 0.5
                extra_y_side = 0.5 / 4
            else:
                h2_Cer.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                            channel_Cer.channelNo)
                h2_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                            channel_Sci.channelNo)
                extra_x_side = 0.5
                extra_y_side = 0.5
            towerx_min = min(towerx_min, channel_Cer.iTowerX)
            towerx_max = max(towerx_max, channel_Cer.iTowerX)
            towery_min = min(towery_min, channel_Cer.iTowerY)
            towery_max = max(towery_max, channel_Cer.iTowerY)

        bbox = ROOT.TBox(towerx_min - extra_x_side, towery_min - extra_y_side,
                         towerx_max + extra_x_side, towery_max + extra_y_side)
        bbox.SetFillColorAlpha(colors[len(bboxes)], 0.1)
        bboxes.append(bbox)

    output_dir = f"plots/Run{run}/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes, extraText="Cer", runNumber=run)
    DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes, extraText="Sci", runNumber=run)

    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"],
        output_dir,
        2,
        output_html=f"html/Run{run}/ChannelMaps/FERSBoards/view.html"
    )


def DrawDRSBoards(run=316):
    """
    Draws the DRS boards for a given run in a TH2D
    """
    drs_boards = buildDRSBoards(run)
    h2_DRS_Cer = ROOT.TH2D(f"h_DRSBoards_Run{run}_Cer", f"DRS Board Channels for Run {run}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)

    h2_DRS_Sci = ROOT.TH2D(f"h_DRSBoards_Run{run}_Sci", f"DRS Board Channels for Run {run}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)

    bboxes = []
    colors = [
        ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1,
        ROOT.kCyan+2, ROOT.kOrange+1, ROOT.kYellow+2, ROOT.kViolet+1, ROOT.kAzure+1,
        ROOT.kSpring+3, ROOT.kTeal+2, ROOT.kPink+2, ROOT.kGray+1, ROOT.kOrange+7
    ]

    for _, DRSBoard in drs_boards.items():
        towerx_min = 99
        towerx_max = -99
        towery_min = 99
        towery_max = -99
        for channel_Sci in DRSBoard.GetSciChannels():
            sci_encoded = channel_Sci.boardNo * 100 + \
                channel_Sci.groupNo * 10 + channel_Sci.channelNo
            if sci_encoded == 0:
                sci_encoded = 0.001
            h2_DRS_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                            sci_encoded)
        for channel_Cer in DRSBoard.GetCerChannels():
            cer_encoded = channel_Cer.boardNo * 100 + \
                channel_Cer.groupNo * 10 + channel_Cer.channelNo
            if cer_encoded == 0:
                cer_encoded = 0.001
            h2_DRS_Cer.Fill(channel_Cer.iTowerX,
                            channel_Cer.iTowerY, cer_encoded)
            towerx_min = min(towerx_min, channel_Cer.iTowerX)
            towerx_max = max(towerx_max, channel_Cer.iTowerX)
            towery_min = min(towery_min, channel_Cer.iTowerY)
            towery_max = max(towery_max, channel_Cer.iTowerY)

        bbox = ROOT.TBox(towerx_min - 0.5, towery_min - 0.5,
                         towerx_max + 0.5, towery_max + 0.5)
        bbox.SetFillColorAlpha(colors[len(bboxes)], 0.1)
        bboxes.append(bbox)

    output_dir = f"plots/Run{run}/ChannelMaps/"
    output_name = f"DRS_Boards_Run{run}"
    h2_DRS_Cer.SetMarkerSize(0.60)
    h2_DRS_Sci.SetMarkerSize(0.60)
    DrawHistos([h2_DRS_Cer], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text1,colz"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=run, zmax=300, zmin=0, ncolors=2)
    DrawHistos([h2_DRS_Sci], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text1,colz"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=run, zmax=300, zmin=0, ncolors=2)

    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"],
        output_dir,
        2,
        output_html=f"html/Run{run}/ChannelMaps/DRSBoards/view.html"
    )


if __name__ == "__main__":
    # Example usage
    for run in [316, 571, 624, 662, 685]:
        DrawFERSBoards(run=run)
        DrawDRSBoards(run=run)

    print("Mapping plots generated for FERS and DRS boards.")
