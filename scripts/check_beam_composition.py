import ROOT
import os
import json
import matplotlib.pyplot as plt
from core.analysis_manager import CaloXAnalysisManager
from utils.parser import get_args
from utils.root_setup import setup_root
from utils.html_generator import generate_html


def check_composition_with_table():
    args = get_args()
    setup_root(n_threads=10, batch_mode=True)

    # 1. Initialize and branch
    analysis = (CaloXAnalysisManager(args)
                .prepare()
                .apply_hole_veto(flag_only=True))

    particle_types = ["electron", "pion", "muon", "proton"]
    initial_rdf = analysis.get_rdf()

    # Booking Initial vs Final
    initial_results = {p: analysis.get_particle_analysis(
        p, flag_only=False).Count() for p in particle_types}
    initial_total = initial_rdf.Count()

    vetoed_rdf = initial_rdf.Filter("isHoleVetoFired == 0")
    final_results = {}
    for p in particle_types:
        from core.selection_manager import SelectionManager
        p_mgr = SelectionManager(vetoed_rdf, analysis.run_number)
        p_mgr.apply_particle_selection(p, flag_only=False)
        final_results[p] = p_mgr.get_rdf().Count()

    final_total = vetoed_rdf.Count()

    all_proxies = list(initial_results.values()) + \
        list(final_results.values()) + [initial_total, final_total]
    ROOT.RDF.RunGraphs(all_proxies)

    i_tot = initial_total.GetValue()
    f_tot = final_total.GetValue()

    # Calculate "others" (Unidentified)
    i_others = i_tot - sum(r.GetValue() for r in initial_results.values())
    f_others = f_tot - sum(r.GetValue() for r in final_results.values())

    output_data = {
        "pre_veto": {
            "total": i_tot,
            "electron": initial_results["electron"].GetValue(),
            "pion": initial_results["pion"].GetValue(),
            "muon": initial_results["muon"].GetValue(),
            "proton": initial_results["proton"].GetValue(),
            "others": i_others
        },
        "post_veto": {
            "total": f_tot,
            "electron": final_results["electron"].GetValue(),
            "pion": final_results["pion"].GetValue(),
            "muon": final_results["muon"].GetValue(),
            "proton": final_results["proton"].GetValue(),
            "others": f_others
        }
    }

    # 4. Save to JSON
    json_path = os.path.join(
        analysis.paths['root'], f"beam_composition.json")
    if not os.path.exists(os.path.dirname(json_path)):
        os.makedirs(os.path.dirname(json_path))
    with open(json_path, "w") as f:
        json.dump(output_data, f, indent=4)

    print(f"✅ Simplified composition results saved to: {json_path}")

    def save_tdr_pie(data_dict, total_val, others_val, title, filename):
        counts = [data_dict[p].GetValue()
                  for p in particle_types] + [others_val]
        labels = [p.capitalize() for p in particle_types] + ["Others"]

        fig, ax = plt.subplots(figsize=(8, 7))
        # HEP-inspired color palette
        colors = ['#3f90da', '#ffa90e', '#bd1f01', '#94a4a2', '#832db6']

        ax.pie(counts, labels=labels, autopct='%1.1f%%', startangle=140,
               colors=colors, textprops={'fontsize': 14})
        ax.set_title(f"{title}\n(Total Events: {total_val})",
                     fontsize=16, pad=20)

        full_path = os.path.join(
            analysis.paths['plots'] + "/ServiceDRS", filename)
        if not os.path.exists(os.path.dirname(full_path)):
            os.makedirs(os.path.dirname(full_path))
        plt.savefig(full_path, bbox_inches='tight')
        plt.close()
        return filename

    plot_before = save_tdr_pie(initial_results, i_tot, i_others,
                               "Composition Pre-Veto", f"comp_before_run{analysis.run_number}.png")
    plot_after = save_tdr_pie(final_results, f_tot, f_others,
                              "Composition Post-Veto", f"comp_after_run{analysis.run_number}.png")

    output_html = os.path.join(
        analysis.paths['html'], f"ServiceDRS/beam_composition.html")
    generate_html(
        png_files=[plot_before,
                   plot_after],
        png_dir=analysis.paths['plots'] + "/drs_service",
        plots_per_row=4,
        output_html=output_html,
        intro_text=f"Beam Composition Summary before and after hole veto"
    )


def main():
    check_composition_with_table()


if __name__ == "__main__":
    main()
