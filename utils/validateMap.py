import ROOT
import sys
sys.path.append("CMSPLOTS")  # noqa
from myFunction import DrawHistos
from utils.channel_map import buildFERSBoards, buildDRSBoards
from utils.html_generator import generate_html
from utils.visualization import visualizeFERSBoards, visualizeDRSBoards

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

    output_dir = f"results/plots/Run{run}/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=run, ncolors=16, zmin=0, zmax=1600)
    DrawHistos([h2_Sci, h2_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=run, ncolors=16, zmin=0, zmax=1600)

    output_html = f"results/html/Run{run}/ChannelMaps/FERSBoards/index.html"
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

    output_dir = f"results/plots/Run{run}/ChannelMaps/"
    output_name = f"DRS_Boards_Run{run}"
    [h2_DRS_Cer, h2_DRS_Cer_3mm], [h2_DRS_Sci, h2_DRS_Sci_3mm] = visualizeDRSBoards(
        drs_boards, suffix=f"Run{run}")

    DrawHistos([h2_DRS_Cer, h2_DRS_Cer_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Cer", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Cer", runNumber=run, zmax=1600, zmin=0, ncolors=16)
    DrawHistos([h2_DRS_Sci, h2_DRS_Sci_3mm], "", xmin, xmax, "iX", ymin,
               ymax, "iY", output_name + "_Sci", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=W_ref, H_ref=H_ref, extraText="Sci", runNumber=run, zmax=1600, zmin=0, ncolors=16)

    output_html = f"results/html/Run{run}/ChannelMaps/DRSBoards/index.html"
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
