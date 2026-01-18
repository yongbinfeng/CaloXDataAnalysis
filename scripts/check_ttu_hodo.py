import ROOT
from plotting.my_function import DrawHistos
from utils.html_generator import generate_html
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import get_ttu_hodo_ranges
from utils.parser import get_args
from utils.plot_helper import save_hists_to_file
from utils.root_setup import setup_root
from utils.timing import auto_timer

auto_timer("Total Execution Time")

setup_root(n_threads=10, batch_mode=True, load_functions=True)

args = get_args()
run_number = args.run

analysis = (CaloXAnalysisManager(args)
            .prepare()
            .apply_hole_veto(flag_only=True)
            )
paths = analysis.paths
rdf_org = analysis.get_rdf()

hodo_min, hodo_max, hodo_nbins = get_ttu_hodo_ranges()


def makeTTUHodoHits(rdf, suffix=""):
    h_nhit_x = rdf.Histo1D(
        (f"h_nhit_x_{suffix}", "Number of hits in X;N_{hits} (X);Entries",
         10, -0.5, 9.5),
        "TTU_Hodo_nHitX")
    h_nhit_y = rdf.Histo1D(
        (f"h_nhit_y_{suffix}", "Number of hits in Y;N_{hits} (Y);Entries",
         10, -0.5, 9.5),
        "TTU_Hodo_nHitY")
    # hit position
    h_x_pos = rdf.Histo1D(
        (f"h_x_pos_{suffix}", "Hodo X Position;X Channel;Entries",
         hodo_nbins, hodo_min, hodo_max),
        "TTU_Hodo_X")
    h_y_pos = rdf.Histo1D(
        (f"h_y_pos_{suffix}", "Hodo Y Position;Y Channel;Entries",
         hodo_nbins, hodo_min, hodo_max),
        "TTU_Hodo_Y")
    h_xy_pos = rdf.Histo2D(
        (f"h_xy_pos_{suffix}", "Hodo XY Position;X Channel;Y Channel;Entries",
         hodo_nbins, hodo_min, hodo_max,
         hodo_nbins, hodo_min, hodo_max),
        "TTU_Hodo_X",
        "TTU_Hodo_Y")

    return [h_nhit_x, h_nhit_y, h_x_pos, h_y_pos, h_xy_pos]


def plotTTUHodoHits():
    infilename = paths['root'] + "/ttuhodo_nhits.root"
    infile = ROOT.TFile.Open(infilename, "READ")
    plots = []
    outdir = paths['plots'] + "/ttuhodo"
    for suffix in ['inc', 'holeveto']:
        # number of hits
        h_nhit_x = infile.Get(f"h_nhit_x_{suffix}")
        h_nhit_y = infile.Get(f"h_nhit_y_{suffix}")
        # get efficiency
        eff_x = h_nhit_x.Integral(h_nhit_x.FindBin(1), h_nhit_x.GetNbinsX(
        )+1) / (h_nhit_x.Integral(0, h_nhit_x.GetNbinsX()+1) + 1e-6)
        eff_y = h_nhit_y.Integral(h_nhit_y.FindBin(1), h_nhit_y.GetNbinsX(
        )+1) / (h_nhit_y.Integral(0, h_nhit_y.GetNbinsX()+1) + 1e-6)

        extraToDraw = ROOT.TPaveText(0.20, 0.80, 0.90, 0.90, "NDC")
        extraToDraw.SetTextAlign(11)
        extraToDraw.SetFillColorAlpha(0, 0)
        extraToDraw.SetBorderSize(0)
        extraToDraw.SetTextFont(42)
        extraToDraw.SetTextSize(0.04)
        extraToDraw.AddText(f"Eff X : {eff_x:.3f}")
        extraToDraw.AddText(f"Eff Y : {eff_y:.3f}")

        DrawHistos(
            [h_nhit_x, h_nhit_y], [
                'X', 'Y'], -1, 6, f"nHits ({suffix})", 1, None, "Events",
            outputname=f"ttuhodo_nhits_{suffix}", outdir=outdir,
            mycolors=[1, 2], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True, extraToDraw=extraToDraw)
        plots.append(f'ttuhodo_nhits_{suffix}.png')

        # hit position
        h_x_pos = infile.Get(f"h_x_pos_{suffix}")
        h_y_pos = infile.Get(f"h_y_pos_{suffix}")
        DrawHistos(
            [h_x_pos, h_y_pos], [
                'X', 'Y'], hodo_min, hodo_max, f"Hit Position ({suffix}) [cm]", 1, None, "Events",
            outputname=f"ttuhodo_hitpos_{suffix}", outdir=outdir,
            mycolors=[1, 2], drawashist=True, run_number=run_number, addOverflow=True, addUnderflow=True)
        plots.append(f'ttuhodo_hitpos_{suffix}.png')

        # 2D hit position
        h_xy_pos = infile.Get(f"h_xy_pos_{suffix}")
        DrawHistos(
            [h_xy_pos], [], hodo_min, hodo_max, "Hodo X [cm]", hodo_min, hodo_max, "Hodo Y [cm]", outputname=f"ttuhodo_hitpos_xy_{suffix}", outdir=outdir,
            drawoptions="COLz", zmin=1, zmax=None, dologz=True,
            dology=False, run_number=run_number, addOverflow=True, doth2=True)
        plots.append(f'ttuhodo_hitpos_xy_{suffix}.png')

    output_html = paths['html'] + "/TTUHodo/ttuhodo_hits.html"
    output_html = generate_html(plots, outdir,
                                plots_per_row=3, output_html=output_html, title="TTU Hodoscope Hits",
                                intro_text="Plots showing the TTU Hodoscope hit distributions: number of hits and hit positions, both inclusive and with hole veto applied.\n If no hits are recorded, the corresponding entries will be underflow bins (-1).")
    print(f"Generated HTML: {output_html}")


def main():
    output_path = paths['root']
    # rdf = rdf_org.Filter("TTU_Hodo_nHitX > 0 && TTU_Hodo_nHitY > 0",
    #                     "At least one hit in both X and Y TTU Hodoscope")
    rdf = rdf_org
    hists = makeTTUHodoHits(rdf, "inc")
    # with hole veto
    rdf_hole_veto = rdf.Filter("is_HoleVeto_vetoed", "Apply Hole Veto")
    hists_veto = makeTTUHodoHits(rdf_hole_veto, "holeveto")
    # run the rdf
    ROOT.RDF.RunGraphs(hists + hists_veto)
    save_hists_to_file(
        hists + hists_veto,
        output_path + "/ttuhodo_nhits.root"
    )
    plotTTUHodoHits()


if __name__ == "__main__":
    main()
