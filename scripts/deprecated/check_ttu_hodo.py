"""
TTU Hodoscope Check Script

Plots TTU Hodoscope hit distributions and positions.
"""

import ROOT
from collections import OrderedDict
from core.analysis_manager import CaloXAnalysisManager
from configs.plot_config import get_ttu_hodo_ranges
from utils.parser import get_args
from utils.plot_helper import save_hists_to_file
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from plotting.calox_plot_helper import create_pave_text
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

# Particle-selected RDFs — inclusive first, then beam-type-dependent selections
rdfs = OrderedDict()
rdfs["inc"] = rdf_org
if analysis.beam_type in ["e+", "positron", "e-", "electron", "positrons", "electrons"]:
    rdfs["electron"] = analysis.get_particle_analysis("electron")
    rdfs["pion"] = analysis.get_particle_analysis("pion")
else:
    rdfs["proton"] = analysis.get_particle_analysis("proton")
    rdfs["pion"] = analysis.get_particle_analysis("pion")
    rdfs["muon"] = analysis.get_particle_analysis("muon")

hodo_min, hodo_max, hodo_nbins = get_ttu_hodo_ranges()

# Styles
STYLE_XY = PlotStyle(dology=False, drawoptions="HIST",
                     mycolors=[1, 2])
STYLE_2D_LOG = PlotStyle(dology=False, dologz=True,
                         drawoptions="COLZ", zmin=1, zmax=None)
STYLE_PROFILE = PlotStyle(dology=False, drawoptions="COLz",
                          zmin=0, zmax=4e4,
                          zlabel="Average PSD Sum")


def makeTTUHodoHits(rdf, suffix=""):
    """Book TTU Hodoscope histograms."""
    h_nhit_x = rdf.Histo1D(
        (f"h_nhit_x_{suffix}", "Number of hits in X;N_{hits} (X);Entries",
         10, -0.5, 9.5),
        "TTU_Hodo_nHitX")
    h_nhit_y = rdf.Histo1D(
        (f"h_nhit_y_{suffix}", "Number of hits in Y;N_{hits} (Y);Entries",
         10, -0.5, 9.5),
        "TTU_Hodo_nHitY")
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


def makePSDvsHodo(rdf, suffix=""):
    """Book 2D profile of PSD energy sum vs TTU hodoscope X and Y positions."""
    rdf = rdf.Define("PSD_Sum_neg", "-PSD_Sum")
    prof = rdf.Profile2D(
        (f"PSD_sum_profile_{suffix}",
         "Average PSD Sum vs Hodoscope X and Y"
         ";TTU Hodo X (cm);TTU Hodo Y (cm);Average PSD Sum",
         hodo_nbins, hodo_min, hodo_max,
         hodo_nbins, hodo_min, hodo_max),
        "TTU_Hodo_X", "TTU_Hodo_Y", "PSD_Sum_neg"
    )
    count = rdf.Histo2D(
        (f"PSD_hodo_count_{suffix}",
         "Event count vs Hodoscope X and Y"
         ";TTU Hodo X (cm);TTU Hodo Y (cm);Counts",
         hodo_nbins, hodo_min, hodo_max,
         hodo_nbins, hodo_min, hodo_max),
        "TTU_Hodo_X", "TTU_Hodo_Y"
    )
    return [prof, count]


def plotPSDvsHodo():
    """Plot 2D profile of PSD energy sum vs TTU hodoscope positions."""
    infilename = f"{paths['root']}/ttuhodo_nhits.root"

    with PlotManager(paths['root'], paths['plots'], paths['html'], run_number) as pm:
        pm.set_output_dir("ttuhodo")
        infile = ROOT.TFile.Open(infilename, "READ")

        for suffix in rdfs.keys():
            hist_prof = infile.Get(f"PSD_sum_profile_{suffix}")
            if hist_prof:
                pm.plot_2d(hist_prof,
                           f"PSD_sum_profile_{suffix}",
                           "TTU Hodo X [cm]", (hodo_min, hodo_max),
                           "TTU Hodo Y [cm]", (hodo_min, hodo_max),
                           style=STYLE_PROFILE)

            hist_count = infile.Get(f"PSD_hodo_count_{suffix}")
            if hist_count:
                pm.plot_2d(hist_count,
                           f"PSD_hodo_count_{suffix}",
                           "TTU Hodo X [cm]", (hodo_min, hodo_max),
                           "TTU Hodo Y [cm]", (hodo_min, hodo_max),
                           style=STYLE_2D_LOG)

        infile.Close()

        intro_text = """2D profile of the average PSD energy sum as a function of TTU hodoscope X and Y
hit positions [cm], shown for inclusive and each particle selection."""

        output_html = pm.generate_html(
            "TTUHodo/psd_vs_hodo.html",
            plots_per_row=2,
            title="PSD Energy vs Hodoscope Position",
            intro_text=intro_text
        )
        print(f"Generated HTML: {output_html}")
        return output_html


def plotTTUHodoHits():
    """Plot TTU Hodoscope distributions using PlotManager."""
    infilename = f"{paths['root']}/ttuhodo_nhits.root"

    with PlotManager(paths['root'], paths['plots'], paths['html'], run_number) as pm:
        pm.set_output_dir("ttuhodo")
        infile = ROOT.TFile.Open(infilename, "READ")

        for suffix in ['inc', 'holeveto']:
            # Number of hits
            h_nhit_x = infile.Get(f"h_nhit_x_{suffix}")
            h_nhit_y = infile.Get(f"h_nhit_y_{suffix}")

            if h_nhit_x and h_nhit_y:
                # Calculate efficiency
                eff_x = h_nhit_x.Integral(h_nhit_x.FindBin(1), h_nhit_x.GetNbinsX() + 1) / \
                    (h_nhit_x.Integral(0, h_nhit_x.GetNbinsX() + 1) + 1e-6)
                eff_y = h_nhit_y.Integral(h_nhit_y.FindBin(1), h_nhit_y.GetNbinsX() + 1) / \
                    (h_nhit_y.Integral(0, h_nhit_y.GetNbinsX() + 1) + 1e-6)

                pave = create_pave_text(0.20, 0.80, 0.90, 0.90)
                pave.AddText(f"Eff X : {eff_x:.3f}")
                pave.AddText(f"Eff Y : {eff_y:.3f}")

                pm.plot_1d(
                    [h_nhit_x, h_nhit_y],
                    f"ttuhodo_nhits_{suffix}",
                    f"nHits ({suffix})",
                    (-1, 6),
                    yrange=(1, None),
                    legends=['X', 'Y'],
                    style=STYLE_XY,
                    extraToDraw=pave
                )

            # Hit positions
            h_x_pos = infile.Get(f"h_x_pos_{suffix}")
            h_y_pos = infile.Get(f"h_y_pos_{suffix}")

            if h_x_pos and h_y_pos:
                pm.plot_1d(
                    [h_x_pos, h_y_pos],
                    f"ttuhodo_hitpos_{suffix}",
                    f"Hit Position ({suffix}) [cm]",
                    (hodo_min, hodo_max),
                    yrange=(1, None),
                    legends=['X', 'Y'],
                    style=STYLE_XY
                )

            # 2D hit position
            h_xy_pos = infile.Get(f"h_xy_pos_{suffix}")
            if h_xy_pos:
                pm.plot_2d(
                    h_xy_pos,
                    f"ttuhodo_hitpos_xy_{suffix}",
                    "Hodo X [cm]",
                    (hodo_min, hodo_max),
                    "Hodo Y [cm]",
                    (hodo_min, hodo_max),
                    style=STYLE_2D_LOG
                )

        infile.Close()

        intro_text = """Plots showing the TTU Hodoscope hit distributions: number of hits and hit positions, 
both inclusive and with hole veto applied.
If no hits are recorded, the corresponding entries will be underflow bins (-1)."""

        output_html = pm.generate_html(
            "TTUHodo/ttuhodo_hits.html",
            plots_per_row=3,
            title="TTU Hodoscope Hits",
            intro_text=intro_text
        )
        print(f"Generated HTML: {output_html}")
        return output_html


def main():
    output_path = paths['root']

    # Hodoscope occupancy: inclusive + hole-veto only
    hists_hodo = makeTTUHodoHits(rdf_org, "inc")
    rdf_hole_veto = rdf_org.Filter("is_HoleVeto_vetoed", "Apply Hole Veto")
    hists_hodo += makeTTUHodoHits(rdf_hole_veto, "holeveto")

    # PSD vs hodo profile: inclusive + all particle selections
    hists_psd = []
    for suffix, rdf in rdfs.items():
        hists_psd += makePSDvsHodo(rdf, suffix)
    # hists_psd += makePSDvsHodo(rdf_org, "inc")

    # Run the RDF graph and save
    all_hists = hists_hodo + hists_psd
    ROOT.RDF.RunGraphs(all_hists)
    save_hists_to_file(all_hists, f"{output_path}/ttuhodo_nhits.root")

    plotTTUHodoHits()
    plotPSDvsHodo()

    PlotManager.print_html_summary()


if __name__ == "__main__":
    main()
