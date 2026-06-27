import os
import json
import ROOT

# Import modules from your codebase
from core.plot_manager import PlotManager
from configs.plot_style import PlotStyle
from plotting.calox_plot_helper import create_pave_text

ROOT.gROOT.SetBatch(True)  # Run in batch mode to avoid opening canvases


def parse_energy(energy_str):
    """Extract numeric energy from strings like '80 GeV', '120GeV', etc."""
    try:
        clean_str = energy_str.lower().replace("gev", "").strip()
        return float(clean_str)
    except Exception:
        return None


def make_energy_scan_plots(run_list=[1422]):
    """
    Process energy scan for a specific list of runs.

    Args:
        run_list (list): A list of run numbers (strings or ints) to process. 
                         If None, processes all runs in Runlist.json.
    """
    rootdir = "results/root"
    plotdir = "results/plots"
    htmldir = "results/html"

    # Load Runlist
    runlist_path = "data/Runlist.json"
    if not os.path.exists(runlist_path):
        print(f"Error: Cannot find {runlist_path}")
        return

    with open(runlist_path, "r") as f:
        full_runlist = json.load(f)

    detectors = ["Sci", "Cer"]

    # Dictionary to hold our TGraphs and point tracking
    graphs = {det: {
        "resp_raw": ROOT.TGraphErrors(),
        "resp_corr": ROOT.TGraphErrors(),
        "sig_raw": ROOT.TGraphErrors(),
        "sig_corr": ROOT.TGraphErrors(),
        "res_raw": ROOT.TGraphErrors(),   # For Resolution: Sigma / E
        "res_corr": ROOT.TGraphErrors(),  # For Resolution: Sigma / E
        "pt_idx": 0
    } for det in detectors}

    # Ensure histograms detach from files so they don't get deleted on close
    ROOT.TH1.AddDirectory(False)

    # Use your PlotManager to save outputs & build the HTML page
    with PlotManager(rootdir, plotdir, htmldir, run_number="EnergyScan") as pm:
        pm.set_output_dir("EnergyScan_Fits")

        # 1. Loop through runs and perform individual fits
        for run_str in run_list:
            pm.add_newline()
            energy_str = full_runlist[str(run_str)].get("beam energy", "")
            energy = parse_energy(energy_str)

            # Require valid, positive energy
            if energy is None or energy <= 0:
                continue

            filepath = f"{rootdir}/Run{run_str}/fers_energy_sum_electron.root"
            if not os.path.exists(filepath):
                print(f"Skipping Run {run_str} - ROOT file not found.")
                continue

            root_file = ROOT.TFile.Open(filepath, "READ")
            if not root_file or root_file.IsZombie():
                print(f"Skipping Run {run_str} - Could not open ROOT file.")
                continue

            # Loop over Scintillator and Cherenkov
            for det in detectors:
                hist_raw_name = f"hist_FERS_energyMix_{det}_calib_sum_electron"
                hist_corr_name = f"hist_FERS_energyMix_{det}_calib_sum_PSDCorr_electron"

                h_raw = root_file.Get(hist_raw_name)
                h_corr = root_file.Get(hist_corr_name)

                # Skip if histograms are missing or empty
                if not h_raw or not h_corr or h_raw.GetEntries() < 10:
                    continue

                h_raw_clone = h_raw.Clone(f"raw_{det}_{run_str}")
                h_corr_clone = h_corr.Clone(f"corr_{det}_{run_str}")

                # PaveText for the fit results on the plot
                pave = create_pave_text(0.15, 0.70, 0.55, 0.90)
                tf1s = []

                # Define the fit range strictly around the beam energy (e.g. +/- 30%)
                fit_min = energy * 0.7
                fit_max = energy * 1.3

                idx = graphs[det]["pt_idx"]

                # --- Fit Raw ---
                fit_raw = h_raw_clone.Fit("gaus", "SQ0", "", fit_min, fit_max)
                if int(fit_raw) == 0:
                    func_raw = h_raw_clone.GetFunction("gaus")
                    func_raw.SetLineColor(ROOT.kBlue)
                    tf1s.append(func_raw)

                    mean = func_raw.GetParameter(1)
                    sigma = func_raw.GetParameter(2)
                    mean_err = func_raw.GetParError(1)
                    sigma_err = func_raw.GetParError(2)
                    chi2 = func_raw.GetChisquare()
                    ndf = func_raw.GetNDF()

                    # Calculate Response (<fit> / E)
                    response = mean / energy
                    response_err = mean_err / energy

                    # Calculate Resolution (Sigma / E)
                    resolution = sigma / mean
                    resolution_err = sigma_err / mean

                    pave.AddText(
                        f"Raw: #mu={mean:.1f}, #sigma={sigma:.1f}, #chi^{{2}}/ndf={chi2:.1f}/{ndf}")

                    graphs[det]["resp_raw"].SetPoint(idx, energy, response)
                    graphs[det]["resp_raw"].SetPointError(idx, 0, response_err)

                    graphs[det]["sig_raw"].SetPoint(idx, energy, sigma)
                    graphs[det]["sig_raw"].SetPointError(idx, 0, sigma_err)

                    graphs[det]["res_raw"].SetPoint(idx, energy, resolution)
                    graphs[det]["res_raw"].SetPointError(
                        idx, 0, resolution_err)

                # --- Fit PSD Corrected ---
                fit_corr = h_corr_clone.Fit(
                    "gaus", "SQ0", "", fit_min, fit_max)
                if int(fit_corr) == 0:
                    func_corr = h_corr_clone.GetFunction("gaus")
                    func_corr.SetLineColor(ROOT.kRed)
                    tf1s.append(func_corr)

                    mean_c = func_corr.GetParameter(1)
                    sigma_c = func_corr.GetParameter(2)
                    mean_err_c = func_corr.GetParError(1)
                    sigma_err_c = func_corr.GetParError(2)
                    chi2_c = func_corr.GetChisquare()
                    ndf_c = func_corr.GetNDF()

                    # Calculate Response (<fit> / E)
                    response_c = mean_c / energy
                    response_err_c = mean_err_c / energy

                    # Calculate Resolution (Sigma / E)
                    resolution_c = sigma_c / mean_c
                    resolution_err_c = sigma_err_c / mean_c

                    pave.AddText(
                        f"PSD Corr: #mu={mean_c:.1f}, #sigma={sigma_c:.1f}, #chi^{{2}}/ndf={chi2_c:.1f}/{ndf_c}")

                    graphs[det]["resp_corr"].SetPoint(idx, energy, response_c)
                    graphs[det]["resp_corr"].SetPointError(
                        idx, 0, response_err_c)

                    graphs[det]["sig_corr"].SetPoint(idx, energy, sigma_c)
                    graphs[det]["sig_corr"].SetPointError(idx, 0, sigma_err_c)

                    graphs[det]["res_corr"].SetPoint(idx, energy, resolution_c)
                    graphs[det]["res_corr"].SetPointError(
                        idx, 0, resolution_err_c)

                # Draw individual run histograms with fits overlaid
                xmax = energy * 2
                pm.plot_1d(
                    [h_raw_clone, h_corr_clone],
                    f"Fit_{det}_Run{run_str}",
                    f"{det} Energy [GeV]",
                    (0, xmax),
                    legends=["Raw", "PSD Corrected"],
                    style=PlotStyle(dology=False, drawoptions="HIST", mycolors=[
                                    ROOT.kBlue, ROOT.kRed]),
                    run_number=run_str,
                    extraToDraw=[pave] + tf1s,
                    legendPos=[0.6, 0.6, 0.9, 0.7]
                )

                graphs[det]["pt_idx"] += 1

            root_file.Close()

        # 2. Draw the Summary Plots
        total_fits = sum(graphs[det]["pt_idx"] for det in detectors)
        if total_fits == 0:
            print("No valid fits found across the requested runs.")
            return

        # Setup styling for the TGraphs
        style_map = {
            "Sci": {"raw_color": ROOT.kBlue, "corr_color": ROOT.kAzure+1, "raw_marker": 20, "corr_marker": 24},
            "Cer": {"raw_color": ROOT.kRed,  "corr_color": ROOT.kPink+1,  "raw_marker": 21, "corr_marker": 25}
        }

        # Legends
        leg_resp = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
        leg_resp.SetBorderSize(0)

        leg_sigma = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
        leg_sigma.SetBorderSize(0)

        leg_res = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
        leg_res.SetBorderSize(0)

        extra_draw_resp = []
        graphs_resp = []
        extra_draw_sigma = []
        extra_draw_res = []
        graphs_res = []

        for det in detectors:
            c_raw = style_map[det]["raw_color"]
            c_corr = style_map[det]["corr_color"]
            m_raw = style_map[det]["raw_marker"]
            m_corr = style_map[det]["corr_marker"]

            # Format Response Graphs
            for key, color, marker, label in [
                ("resp_raw", c_raw, m_raw, f"{det} Raw"),
                ("resp_corr", c_corr, m_corr, f"{det} PSD Corr")
            ]:
                gr = graphs[det][key]
                if gr.GetN() > 0:
                    gr.SetMarkerStyle(marker)
                    gr.SetMarkerColor(color)
                    gr.SetLineColor(color)
                    gr.SetDrawOption("PE")
                    leg_resp.AddEntry(gr, label, "pe")
                    graphs_resp.append(gr)
                    # extra_draw_resp.append(gr)

            # Format Sigma Graphs
            for key, color, marker, label in [
                ("sig_raw", c_raw, m_raw, f"{det} Raw"),
                ("sig_corr", c_corr, m_corr, f"{det} PSD Corr")
            ]:
                gr = graphs[det][key]
                if gr.GetN() > 0:
                    gr.SetMarkerStyle(marker)
                    gr.SetMarkerColor(color)
                    gr.SetLineColor(color)
                    gr.SetDrawOption("PE")
                    leg_sigma.AddEntry(gr, label, "pe")
                    extra_draw_sigma.append(gr)

            # Format Resolution (Sigma/E) Graphs
            for key, color, marker, label in [
                ("res_raw", c_raw, m_raw, f"{det} Raw"),
                ("res_corr", c_corr, m_corr, f"{det} PSD Corr")
            ]:
                gr = graphs[det][key]
                if gr.GetN() > 0:
                    gr.SetMarkerStyle(marker)
                    gr.SetMarkerColor(color)
                    gr.SetLineColor(color)
                    gr.SetDrawOption("PE")
                    leg_res.AddEntry(gr, label, "pe")
                    # extra_draw_res.append(gr)
                    graphs_res.append(gr)

        # Ideal y=1 line for the Response summary
        f_ideal_response = ROOT.TF1("ideal_response", "1", 0, 200)
        f_ideal_response.SetLineStyle(2)
        f_ideal_response.SetLineColor(ROOT.kBlack)
        # leg_resp.AddEntry(f_ideal_response, "Ideal (y=1)", "l")

        extra_draw_resp.append(f_ideal_response)
        extra_draw_resp.append(leg_resp)
        extra_draw_sigma.append(leg_sigma)
        extra_draw_res.append(leg_res)

        # Draw Response Summary (Prepend pushes it to the top of the HTML output)
        pm.plot_1d(
            graphs_resp,
            "Summary_Response_vs_Energy",
            "Beam Energy [GeV]",
            (0, 120),
            ylabel="Response (#mu / E)",
            yrange=(0.8, 1.4),  # Window clearly showing deviations around 1.0
            style=PlotStyle(drawoptions="PL"),
            extraToDraw=extra_draw_resp,
            prepend=True,
            includeRunNumber=False
        )

        # Draw Resolution (Sigma/E) Summary
        pm.plot_1d(
            graphs_res,
            "Summary_Resolution_vs_Energy",
            "Beam Energy [GeV]",
            (0, 120),
            ylabel="Resolution (#sigma / #mu)",
            yrange=(0, 0.5),
            style=PlotStyle(drawoptions="PL"),
            extraToDraw=extra_draw_res,
            prepend=True,
            includeRunNumber=False
        )

        # 3. Generate HTML output
        html_out = pm.generate_html(
            "EnergyScan/index.html", plots_per_row=4, title="Energy Scan Fit Results (Sci & Cer)")
        print(f"\nAll done! You can view the full report at: {html_out}")


if __name__ == "__main__":
    run_list = [1424, 1423, 1416, 1514, 1526, 1355, 1527]
    make_energy_scan_plots(run_list=run_list)
