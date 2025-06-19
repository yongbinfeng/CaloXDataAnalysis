import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.channel_map import buildFERSBoards, buildDRSBoards
from utils.html_generator import generate_html

ROOT.gROOT.SetBatch(True)  # Disable interactive mode

xmax = 15.5
xmin = -12.5
ymax = 12.5
ymin = -7.5
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

    extraToDraw = ROOT.TPaveText(0.01, 0.78, 0.12, 0.88, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Run: {run}")
    extraToDraw_Cer = extraToDraw.Clone("extraToDraw_2")
    extraToDraw_Cer.AddText("Cerenkov")
    extraToDraw_Sci = extraToDraw.Clone("extraToDraw_3")
    extraToDraw_Sci.AddText("Scintillator")

    output_dir = f"plots/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes + [extraToDraw_Cer])
    DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text0", "text0"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes + [extraToDraw_Sci])

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
            h2_DRS_Sci.Fill(channel_Sci.iTowerX, channel_Sci.iTowerY,
                            float(f"{channel_Sci.groupNo}.{channel_Sci.channelNo}"))
        for channel_Cer in DRSBoard.GetCerChannels():
            h2_DRS_Cer.Fill(channel_Cer.iTowerX, channel_Cer.iTowerY,
                            float(f"{channel_Cer.groupNo}.{channel_Cer.channelNo}"))
            towerx_min = min(towerx_min, channel_Cer.iTowerX)
            towerx_max = max(towerx_max, channel_Cer.iTowerX)
            towery_min = min(towery_min, channel_Cer.iTowerY)
            towery_max = max(towery_max, channel_Cer.iTowerY)

        bbox = ROOT.TBox(towerx_min - 0.5, towery_min - 0.5,
                         towerx_max + 0.5, towery_max + 0.5)
        bbox.SetFillColorAlpha(colors[len(bboxes)], 0.1)
        bboxes.append(bbox)

    extraToDraw = ROOT.TPaveText(0.01, 0.78, 0.12, 0.88, "NDC")
    extraToDraw.SetTextAlign(11)
    extraToDraw.SetFillColorAlpha(0, 0)
    extraToDraw.SetBorderSize(0)
    extraToDraw.SetTextFont(42)
    extraToDraw.SetTextSize(0.04)
    extraToDraw.AddText(f"Run: {run}")
    extraToDraw_Cer = extraToDraw.Clone("extraToDraw_2")
    extraToDraw_Cer.AddText("Cerenkov")
    extraToDraw_Sci = extraToDraw.Clone("extraToDraw_3")
    extraToDraw_Sci.AddText("Scintillator")

    output_dir = f"plots/ChannelMaps/"
    output_name = f"DRS_Boards_Run{run}"
    DrawHistos([h2_DRS_Cer], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text1"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes + [extraToDraw_Cer], nTextDigits=1)
    DrawHistos([h2_DRS_Sci], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text1"] + ["same"] * len(bboxes),
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraToDraw=bboxes + [extraToDraw_Sci], nTextDigits=1)

    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"],
        output_dir,
        2,
        output_html=f"html/Run{run}/ChannelMaps/DRSBoards/view.html"
    )


if __name__ == "__main__":
    # Example usage
    run_number = 662
    DrawFERSBoards(run=run_number)
    DrawDRSBoards(run=run_number)

    print("Mapping plots generated for FERS and DRS boards.")
