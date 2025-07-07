import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.channel_map import buildFERSBoards, buildDRSBoards
from utils.html_generator import generate_html
from utils.visualization import visualizeFERSBoards

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

    [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm] = visualizeFERSBoards(
        fers_boards, suffix=f"Run{run}")

    output_dir = f"plots/Run{run}/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=run, ncolors=16, zmin=0, zmax=1600)
    DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=run, ncolors=16, zmin=0, zmax=1600)

    output_html = f"html/Run{run}/ChannelMaps/FERSBoards/view.html"
    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"],
        output_dir,
        2,
        output_html=output_html
    )
    return output_html


def DrawDRSBoards(run=316):
    """
    Draws the DRS boards for a given run in a TH2D
    """
    drs_boards = buildDRSBoards(run)
    h2_DRS_Cer = ROOT.TH2D(f"h_DRSBoards_Run{run}_Cer", f"DRS Board Channels for Run {run}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)

    h2_DRS_Sci = ROOT.TH2D(f"h_DRSBoards_Run{run}_Sci", f"DRS Board Channels for Run {run}",
                           int(xmax - xmin), xmin, xmax, int(ymax - ymin), ymin, ymax)

    for _, DRSBoard in drs_boards.items():
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

    output_dir = f"plots/Run{run}/ChannelMaps/"
    output_name = f"DRS_Boards_Run{run}"
    h2_DRS_Cer.SetMarkerSize(0.60)
    h2_DRS_Sci.SetMarkerSize(0.60)
    DrawHistos([h2_DRS_Cer], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text1,col"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=run, zmax=1600, zmin=0, ncolors=16)
    DrawHistos([h2_DRS_Sci], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text1,col"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=run, zmax=1600, zmin=0, ncolors=16)

    output_html = f"html/Run{run}/ChannelMaps/DRSBoards/view.html"
    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"],
        output_dir,
        2,
        output_html=output_html
    )
    return output_html


if __name__ == "__main__":
    # Example usage
    output_htmls = {}
    for run in [316, 571, 624, 662, 685]:
        output_htmls[f"fers mapping {run}"] = DrawFERSBoards(run=run)
        output_htmls[f"drs mapping {run}"] = DrawDRSBoards(run=run)

    for key, value in output_htmls.items():
        print(f"{key}: {value}")

    print("Mapping plots generated for FERS and DRS boards.")
