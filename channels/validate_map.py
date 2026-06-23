import ROOT
from plotting.my_function import DrawHistos
from channels.channel_map import build_fers_boards, build_drs_boards
from utils.html_generator import generate_html
from utils.visualization import (
    visualizeFERSBoards, visualizeDRSBoards,
    FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, FERS_YMIN_DISPLAY, FERS_YMAX_DISPLAY,
    FERS_W_REF, FERS_H_REF,
)

ROOT.gROOT.SetBatch(True)  # Disable interactive mode


def DrawFERSBoards(run_number=316):
    """
    Draws the FERS boards for a given run_number in a TH2D
    """
    fersboards = build_fers_boards(run_number)

    suffix = f"Run{run_number}"

    [h2_Cer, h2_Cer_3mm], [h2_Sci, h2_Sci_3mm] = visualizeFERSBoards(
        fersboards, suffix=suffix, quartzOnly=0)
    [h2_Cer_quartz, h2_Cer_3mm_quartz], _ = visualizeFERSBoards(
        fersboards, suffix=suffix + "_quartz", quartzOnly=1)
    [h2_Cer_plastic, h2_Cer_3mm_plastic], _ = visualizeFERSBoards(
        fersboards, suffix=suffix + "_plastic", quartzOnly=2)

    output_dir = f"results/plots/Run{run_number}/ChannelMaps/"
    output_name = f"FERS_Boards_Run{run_number}"
    DrawHistos([h2_Cer, h2_Cer_3mm], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Cer", run_number=run_number, ncolors=16, zmin=0, zmax=1600)
    DrawHistos([h2_Sci, h2_Sci_3mm], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Sci", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Sci", run_number=run_number, ncolors=16, zmin=0, zmax=1600)
    DrawHistos([h2_Cer_quartz, h2_Cer_3mm_quartz], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer_quartz", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Quartz", run_number=run_number, ncolors=16, zmin=0, zmax=1600)
    DrawHistos([h2_Cer_plastic, h2_Cer_3mm_plastic], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer_plastic", dology=False, drawoptions=["col,text", "col,text"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Plastic", run_number=run_number, ncolors=16, zmin=0, zmax=1600)

    output_html = f"results/html/Run{run_number}/ChannelMaps/FERSBoards.html"
    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"] +
        [output_name + "_Cer_quartz.png", output_name + "_Cer_plastic.png"],
        output_dir,
        2,
        output_html=output_html
    )
    return output_html


def DrawDRSBoards(run_number=316):
    """
    Draws the DRS boards for a given run_number in a TH2D
    """
    drs_boards = build_drs_boards(run_number)

    output_dir = f"results/plots/Run{run_number}/ChannelMaps/"
    output_name = f"DRS_Boards_Run{run_number}"
    [h2_DRS_Cer, h2_DRS_Cer_3mm], [h2_DRS_Sci, h2_DRS_Sci_3mm] = visualizeDRSBoards(
        drs_boards, suffix=f"Run{run_number}")
    [h2_DRS_Cer_quartz, h2_DRS_Cer_3mm_quartz], _ = visualizeDRSBoards(
        drs_boards, suffix=f"Run{run_number}_quartz", quartzOnly=1)
    [h2_DRS_Cer_plastic, h2_DRS_Cer_3mm_plastic], _ = visualizeDRSBoards(
        drs_boards, suffix=f"Run{run_number}_plastic", quartzOnly=2)

    DrawHistos([h2_DRS_Cer, h2_DRS_Cer_3mm], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Cer", run_number=run_number, zmax=1600, zmin=0, ncolors=16)
    DrawHistos([h2_DRS_Sci, h2_DRS_Sci_3mm], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Sci", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Sci", run_number=run_number, zmax=1600, zmin=0, ncolors=16)
    DrawHistos([h2_DRS_Cer_quartz, h2_DRS_Cer_3mm_quartz], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer_quartz", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Quartz", run_number=run_number, zmax=1600, zmin=0, ncolors=16)
    DrawHistos([h2_DRS_Cer_plastic, h2_DRS_Cer_3mm_plastic], "", FERS_XMIN_DISPLAY, FERS_XMAX_DISPLAY, "X [cm]", FERS_YMIN_DISPLAY,
               FERS_YMAX_DISPLAY, "Y [cm]", output_name + "_Cer_plastic", dology=False, drawoptions=["text,col", "text,col"],
               outdir=output_dir, doth2=True, W_ref=FERS_W_REF, H_ref=FERS_H_REF, extra_text="Plastic", run_number=run_number, zmax=1600, zmin=0, ncolors=16)

    output_html = f"results/html/Run{run_number}/ChannelMaps/DRSBoards.html"
    generate_html(
        [output_name + "_Cer.png", output_name + "_Sci.png"] +
        [output_name + "_Cer_quartz.png", output_name + "_Cer_plastic.png"],
        output_dir,
        2,
        output_html=output_html
    )
    return output_html


if __name__ == "__main__":
    # Example usage
    output_htmls = {}
    for run_number in [1175]:
        output_htmls[f"fers mapping {run_number}"] = DrawFERSBoards(run_number=run_number)
        output_htmls[f"drs mapping {run_number}"] = DrawDRSBoards(run_number=run_number)

    for key, value in output_htmls.items():
        print(f"{key}: {value}")

    print("Mapping plots generated for FERS and DRS boards.")
