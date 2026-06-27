"""
Overlay simulation and data DRS timing distributions for pions and electrons.

Simulation (40 GeV):  {DRSOutput,DRSOutput1stPeak}_meridional_{Pla,Qua}__{pi,ele}
    from hitsZ_plots_40GeV.root  (x-axis already in ns)

Data:
    pion     (Run 1469) — drs_vs_ts_calibrated_combined_pion_run1469.root
    electron (Run 1514) — drs_vs_ts_calibrated_combined_electron_1514.root
    (x-axis in time slices, 1 TS = 0.2 ns)

A single x-axis shift is derived from the DRSOutput Quartz-pion peak comparison
and applied consistently to all simulation histograms (both DRS types, both
particles, both fibers).

Colors:  Quartz = red,  Plastic = pink
Styles:  Data = solid,  MC = dashed
Y-scale: MC peak scaled to match corresponding data peak.
"""

import os
import sys
import ROOT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, ".."))   # CaloXDataAnalysis root
from plotting.my_function import DrawHistos

ROOT.TH1.AddDirectory(False)
ROOT.gROOT.SetBatch(True)

# ── constants ─────────────────────────────────────────────────────────────────
NS_PER_TS     = 0.2    # 1 time slice = 200 ps
RUN_PION      = 1469
RUN_ELECTRON  = 1514

COLOR_QUA  = ROOT.kRed
COLOR_PLA  = ROOT.kPink + 1
STYLE_DATA = 1           # solid
STYLE_MC   = 2           # dashed

# ── input files ──────────────────────────────────────────────────────────────
SIM_FILE      = "/Users/yongbinfeng/Desktop/plots/root/hitsZ_plots_40GeV.root"
DATA_FILE_PI  = "/Users/yongbinfeng/Desktop/plots/root/drs_vs_ts_calibrated_combined_pion_run1469.root"
DATA_FILE_ELE = "/Users/yongbinfeng/Desktop/plots/root/drs_vs_ts_calibrated_combined_electron_1514.root"

# ── output ───────────────────────────────────────────────────────────────────
OUTDIR = "results/plots/SimVsData/DRS_Timing/"
os.makedirs(OUTDIR, exist_ok=True)


# ── helpers ───────────────────────────────────────────────────────────────────
def peak_position(h):
    return h.GetXaxis().GetBinCenter(h.GetMaximumBin())


def rescale_xaxis(h, scale, new_name):
    nbins = h.GetNbinsX()
    xmin  = h.GetXaxis().GetXmin() * scale
    xmax  = h.GetXaxis().GetXmax() * scale
    h_new = ROOT.TH1D(new_name, h.GetTitle(), nbins, xmin, xmax)
    for i in range(1, nbins + 1):
        h_new.SetBinContent(i, h.GetBinContent(i))
        h_new.SetBinError(i, h.GetBinError(i))
    return h_new


def shift_xaxis(h, shift, new_name):
    nbins = h.GetNbinsX()
    xmin  = h.GetXaxis().GetXmin() + shift
    xmax  = h.GetXaxis().GetXmax() + shift
    h_new = ROOT.TH1D(new_name, h.GetTitle(), nbins, xmin, xmax)
    for i in range(1, nbins + 1):
        h_new.SetBinContent(i, h.GetBinContent(i))
        h_new.SetBinError(i, h.GetBinError(i))
    return h_new


def scale_mc_to_data_peak(h_mc, h_data):
    peak_mc = h_mc.GetMaximum()
    if peak_mc > 0:
        h_mc.Scale(h_data.GetMaximum() / peak_mc)


def draw(hists, labels, output_name, run_number, colors, linestyles, legend_pos=None):
    DrawHistos(
        hists, labels,
        xmin, xmax,
        "DRS output time (ns)",
        0, None,
        "ADC",
        output_name,
        dology=False,
        drawoptions=["hist,C"] * len(hists),
        mycolors=colors,
        linestyles=linestyles,
        addOverflow=False,
        addUnderflow=False,
        outdir=OUTDIR,
        legendPos=legend_pos or [0.30, 0.75, 0.40, 0.90],
        legendoptions=["L"] * len(hists),
        noCMS=False,
        noLumi=False,
        run_number=run_number,
        usePDF=True,
        usePNG=False,
    )


# ── load and convert data histograms (TS → ns) ────────────────────────────────
f_pi = ROOT.TFile(DATA_FILE_PI, "READ")
h_data_pi_qua = rescale_xaxis(f_pi.Get("prof_DRS_vs_TS_Cer_Quartz_Combined_pion"),  NS_PER_TS, "h_data_pi_qua")
h_data_pi_pla = rescale_xaxis(f_pi.Get("prof_DRS_vs_TS_Cer_Plastic_Combined_pion"), NS_PER_TS, "h_data_pi_pla")
f_pi.Close()

f_ele = ROOT.TFile(DATA_FILE_ELE, "READ")
h_data_ele_qua = rescale_xaxis(f_ele.Get("prof_DRS_vs_TS_Cer_Quartz_Combined_electron"),  NS_PER_TS, "h_data_ele_qua")
h_data_ele_pla = rescale_xaxis(f_ele.Get("prof_DRS_vs_TS_Cer_Plastic_Combined_electron"), NS_PER_TS, "h_data_ele_pla")
f_ele.Close()

# ── compute shift once from DRSOutput Quartz pion ────────────────────────────
f_sim = ROOT.TFile(SIM_FILE, "READ")
h_sim_pi_qua_ref = f_sim.Get("DRSOutput_meridional_Qua__pi").Clone("h_sim_pi_qua_ref")
f_sim.Close()
shift = peak_position(h_data_pi_qua) - peak_position(h_sim_pi_qua_ref)
print(f"Shift (from DRSOutput Quartz pion): {shift:.3f} ns")

# ── loop over DRS output types ────────────────────────────────────────────────
DRS_TYPES = ["DRSOutput", "DRSOutput1stPeak"]

for drs_type in DRS_TYPES:
    tag = "" if drs_type == "DRSOutput" else "_1stPeak"

    f_sim = ROOT.TFile(SIM_FILE, "READ")
    h_sim_pi_qua_raw  = f_sim.Get(f"{drs_type}_meridional_Qua__pi" ).Clone(f"h_sim_pi_qua_raw{tag}")
    h_sim_pi_pla_raw  = f_sim.Get(f"{drs_type}_meridional_Pla__pi" ).Clone(f"h_sim_pi_pla_raw{tag}")
    h_sim_ele_qua_raw = f_sim.Get(f"{drs_type}_meridional_Qua__ele").Clone(f"h_sim_ele_qua_raw{tag}")
    h_sim_ele_pla_raw = f_sim.Get(f"{drs_type}_meridional_Pla__ele").Clone(f"h_sim_ele_pla_raw{tag}")
    f_sim.Close()

    print(f"\n{drs_type}:")
    print(f"  pi  Qua: sim {peak_position(h_sim_pi_qua_raw)  + shift:.3f} ns  data {peak_position(h_data_pi_qua):.3f} ns")
    print(f"  pi  Pla: sim {peak_position(h_sim_pi_pla_raw)  + shift:.3f} ns  data {peak_position(h_data_pi_pla):.3f} ns")
    print(f"  ele Qua: sim {peak_position(h_sim_ele_qua_raw) + shift:.3f} ns  data {peak_position(h_data_ele_qua):.3f} ns")
    print(f"  ele Pla: sim {peak_position(h_sim_ele_pla_raw) + shift:.3f} ns  data {peak_position(h_data_ele_pla):.3f} ns")

    h_sim_pi_qua  = shift_xaxis(h_sim_pi_qua_raw,  shift, f"h_sim_pi_qua{tag}")
    h_sim_pi_pla  = shift_xaxis(h_sim_pi_pla_raw,  shift, f"h_sim_pi_pla{tag}")
    h_sim_ele_qua = shift_xaxis(h_sim_ele_qua_raw, shift, f"h_sim_ele_qua{tag}")
    h_sim_ele_pla = shift_xaxis(h_sim_ele_pla_raw, shift, f"h_sim_ele_pla{tag}")

    scale_mc_to_data_peak(h_sim_pi_qua,  h_data_pi_qua)
    scale_mc_to_data_peak(h_sim_pi_pla,  h_data_pi_pla)
    scale_mc_to_data_peak(h_sim_ele_qua, h_data_ele_qua)
    scale_mc_to_data_peak(h_sim_ele_pla, h_data_ele_pla)

    xmin, xmax = -17.0, -7.0

    # pion
    draw([h_sim_pi_qua, h_data_pi_qua],
         ["Simulation", "Data"],
         f"DRS_Timing_Quartz_SimVsData_pion{tag}",
         RUN_PION, [COLOR_QUA, COLOR_QUA], [STYLE_MC, STYLE_DATA])

    draw([h_sim_pi_pla, h_data_pi_pla],
         ["Simulation", "Data"],
         f"DRS_Timing_Plastic_SimVsData_pion{tag}",
         RUN_PION, [COLOR_PLA, COLOR_PLA], [STYLE_MC, STYLE_DATA])

    draw([h_sim_pi_qua, h_data_pi_qua, h_sim_pi_pla, h_data_pi_pla],
         ["Quartz, Simulation", "Quartz, Data", "Plastic, Simulation", "Plastic, Data"],
         f"DRS_Timing_Combined_SimVsData_pion{tag}",
         RUN_PION,
         [COLOR_QUA, COLOR_QUA, COLOR_PLA, COLOR_PLA],
         [STYLE_MC, STYLE_DATA, STYLE_MC, STYLE_DATA],
         legend_pos=[0.30, 0.65, 0.40, 0.90])

    # electron
    draw([h_sim_ele_qua, h_data_ele_qua],
         ["Simulation", "Data"],
         f"DRS_Timing_Quartz_SimVsData_electron{tag}",
         RUN_ELECTRON, [COLOR_QUA, COLOR_QUA], [STYLE_MC, STYLE_DATA])

    draw([h_sim_ele_pla, h_data_ele_pla],
         ["Simulation", "Data"],
         f"DRS_Timing_Plastic_SimVsData_electron{tag}",
         RUN_ELECTRON, [COLOR_PLA, COLOR_PLA], [STYLE_MC, STYLE_DATA])

    draw([h_sim_ele_qua, h_data_ele_qua, h_sim_ele_pla, h_data_ele_pla],
         ["Quartz, Simulation", "Quartz, Data", "Plastic, Simulation", "Plastic, Data"],
         f"DRS_Timing_Combined_SimVsData_electron{tag}",
         RUN_ELECTRON,
         [COLOR_QUA, COLOR_QUA, COLOR_PLA, COLOR_PLA],
         [STYLE_MC, STYLE_DATA, STYLE_MC, STYLE_DATA],
         legend_pos=[0.30, 0.65, 0.40, 0.90])

print(f"\nPlots saved to {OUTDIR}")
