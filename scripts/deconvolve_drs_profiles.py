"""
Deconvolve DRS vs TS profiles with the average SiPM pulse shape.

Reads already-saved per-channel DRS vs TS profiles from
  <rootdir>/drs_vs_ts_calibrated_{suffix}.root
and deconvolves each with a reference pulse shape using a Wiener filter.
The deconvolved histograms are written to
  <rootdir>/drs_vs_ts_deconvolved_{suffix}.root
Plots are saved as PDFs and an HTML summary is generated.

Run this after make_drs_timing_plots.py has produced the profile ROOT files.
"""

import numpy as np
import ROOT

from channels.channel_map import build_drs_boards
from configs.plot_style import PlotStyle
from core.plot_manager import PlotManager
from plotting.calox_plot_helper import create_pave_text
from plotting.my_function import LHistos2Hist
from utils.parser import get_args
from utils.plot_helper import get_run_paths
from utils.utils import number_to_string

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

args = get_args()
run_number = args.run

DRSBoards = build_drs_boards(run=run_number)
paths = get_run_paths(run_number)
rootdir = paths["root"]
plotdir = paths["plots"]
htmldir = paths["html"]


# ---------------------------------------------------------------------------
# Default pulse-shape file (path relative to the working directory)
# ---------------------------------------------------------------------------
DEFAULT_PULSE_FILE = (
    "data/AvePulse3mmSIPM_50ps_2200ts_max682ts_Inter4_R1592.npy"
)

def normalize_to_peak(h):
    """Scale histogram so its peak absolute value is 1."""
    peak = max(abs(h.GetMaximum()), abs(h.GetMinimum()))
    if peak > 0:
        h.Scale(1.0 / peak)


# Common plot styles (mirrored from make_drs_timing_plots.py)
STYLE_CER_QUARTZ_PLASTIC_SCI = PlotStyle(
    dology=False,
    drawoptions="HIST",
    mycolors=[2, 6, 4],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.25, 0.75, 0.40, 0.90]
)

STYLE_CER_QUARTZ_PLASTIC = PlotStyle(
    dology=False,
    drawoptions=["hist,C", "hist,C"],
    mycolors=[2, 6],
    addOverflow=False,
    addUnderflow=False,
    legendPos=[0.30, 0.80, 0.40, 0.90]
)


# ---------------------------------------------------------------------------
# Deconvolution
# ---------------------------------------------------------------------------
def deconvolveDRSvsTS(
    suffix: str = "",
    pulse_file: str = DEFAULT_PULSE_FILE,
    noise_level: float = 0.01,
    ts_per_pulse_ts: int = 4,
    save: bool = True,
):
    """Deconvolve DRS vs TS profiles with the SiPM average pulse shape.

    Parameters
    ----------
    suffix : str
        Category suffix used when the profiles were saved (e.g. ``'pion'``).
    pulse_file : str
        Path to the ``.npy`` file containing the average SiPM pulse shape.
        Expected: 2200 time slices at 50 ps/ts (peak at ts 682, inter x4).
    noise_level : float
        Wiener regularisation parameter — ratio of noise power to signal power.
        Increase to suppress ringing at the cost of resolution.
    ts_per_pulse_ts : int
        Number of pulse time slices per DRS TS (4 = 200 ps DRS / 50 ps pulse).
    save : bool
        If True, write deconvolved histograms to
        ``<rootdir>/drs_vs_ts_deconvolved_{suffix}.root``.

    Returns
    -------
    dict
        ``{"{var}_{sTowerX}_{sTowerY}": TH1D}`` of deconvolved profiles.
    """

    # ------------------------------------------------------------------
    # Load and prepare pulse shape
    # ------------------------------------------------------------------
    pulse_raw = np.load(pulse_file).astype(float)

    # Down-sample from 50 ps to DRS TS steps (200 ps) by averaging blocks
    n_full = (len(pulse_raw) // ts_per_pulse_ts) * ts_per_pulse_ts
    pulse_ds = pulse_raw[:n_full].reshape(-1, ts_per_pulse_ts).mean(axis=1)

    # Normalise to unit sum so deconvolved amplitude stays in ADC counts
    pulse_ds /= pulse_ds.sum()

    # Roll peak to index 0 so the FFT does not introduce a global time shift
    # (profiles are already time-aligned via AlignedTS)
    peak_idx = int(np.argmax(pulse_ds))
    pulse_ds = np.roll(pulse_ds, -peak_idx)

    # ------------------------------------------------------------------
    # Open input file
    # ------------------------------------------------------------------
    infile_name = f"{rootdir}/drs_vs_ts_calibrated_{suffix}.root"
    infile = ROOT.TFile(infile_name, "READ")
    if not infile or infile.IsZombie():
        print(f"Error: cannot open {infile_name}")
        return {}

    deconv_hists = {}

    for _, DRSBoard in DRSBoards.items():
        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            for var in ["Cer", "Sci"]:
                hprof_name = (
                    f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}"
                )
                hprof = infile.Get(hprof_name)
                if not hprof:
                    print(f"Warning: {hprof_name} not found, skipping")
                    continue

                hprof_1d = hprof.ProjectionX(f"proj_{hprof_name}", "E")

                n_bins = hprof_1d.GetNbinsX()
                data = np.array(
                    [hprof_1d.GetBinContent(i + 1) for i in range(n_bins)]
                )

                # Pad (or truncate) pulse template to data length
                pulse_padded = np.zeros(n_bins)
                n_copy = min(len(pulse_ds), n_bins)
                pulse_padded[:n_copy] = pulse_ds[:n_copy]

                # Wiener deconvolution:  X = Y * H* / (|H|^2 + lambda)
                DATA = np.fft.rfft(data)
                PULSE = np.fft.rfft(pulse_padded)
                PULSE_CONJ = np.conj(PULSE)
                deconv_freq = DATA * PULSE_CONJ / (PULSE * PULSE_CONJ + noise_level)
                deconv = np.fft.irfft(deconv_freq, n=n_bins)

                # Store result in a clone of the projection histogram
                h_out = hprof_1d.Clone(
                    f"deconv_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}"
                )
                h_out.SetTitle(
                    f"Deconvolved DRS vs TS {var}, "
                    f"Tower({sTowerX},{sTowerY}), {suffix};"
                    f"TS;Deconvolved ADC"
                )
                h_out.Reset()
                for i in range(n_bins):
                    h_out.SetBinContent(i + 1, deconv[i])
                    h_out.SetBinError(i + 1, 0.0)

                key = f"{var}_{sTowerX}_{sTowerY}"
                deconv_hists[key] = h_out

    infile.Close()

    if save and deconv_hists:
        outfile_name = f"{rootdir}/drs_vs_ts_deconvolved_{suffix}.root"
        outfile = ROOT.TFile(outfile_name, "RECREATE")
        for h in deconv_hists.values():
            h.SetDirectory(outfile)
            h.Write()
        outfile.Close()
        print(f"Saved {len(deconv_hists)} deconvolved histograms -> {outfile_name}")

    return deconv_hists


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def makeDRSDeconvPlots(suffix: str = ""):
    """Plot deconvolved DRS vs TS overlaid with the measured profile.

    Only 3mm boards are included.  For each channel the deconvolved histogram
    is drawn with a solid line and the measured profile with a dashed line,
    both normalised to unit area so they can be compared directly.
    """
    pm = PlotManager(rootdir, plotdir, htmldir, run_number)
    pm.set_output_dir(f"DRS_vs_ts_deconvolved_{suffix}")

    infile_deconv_name = f"{rootdir}/drs_vs_ts_deconvolved_{suffix}.root"
    infile_meas_name = f"{rootdir}/drs_vs_ts_calibrated_{suffix}.root"

    infile_deconv = ROOT.TFile(infile_deconv_name, "READ")
    if not infile_deconv or infile_deconv.IsZombie():
        print(f"Error: cannot open {infile_deconv_name}")
        return None

    infile_meas = ROOT.TFile(infile_meas_name, "READ")
    if not infile_meas or infile_meas.IsZombie():
        print(f"Error: cannot open {infile_meas_name}")
        return None

    # Collect combined-plot inputs for deconvolved and measured, 3mm channels
    hists_Cer_Quartz = []
    hists_Cer_Plastic = []
    hists_Cer = []
    hists_Sci = []
    meas_Cer_Quartz = []
    meas_Cer_Plastic = []
    meas_Cer = []
    meas_Sci = []
    xtitle = "TS"
    ytitle = "Normalised ADC"

    for _, DRSBoard in DRSBoards.items():
        board_no = DRSBoard.board_no
        pm.add_newline()

        for i_tower_x, i_tower_y in DRSBoard.get_list_of_towers():
            sTowerX = number_to_string(i_tower_x)
            sTowerY = number_to_string(i_tower_y)

            hists_deconv = {}
            hists_meas = {}
            isQuartz = False
            channelNos = {}

            for var in ["Cer", "Sci"]:
                chan = DRSBoard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))

                # 3mm channels only
                if chan is None or chan.is6mm:
                    continue

                if var == "Cer" and chan.isQuartz:
                    isQuartz = True

                # Deconvolved histogram
                h_deconv = infile_deconv.Get(
                    f"deconv_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}")
                hists_deconv[var] = h_deconv if h_deconv else None

                # Measured profile (project to 1D)
                hprof = infile_meas.Get(
                    f"prof_DRS_vs_TS_{var}_{sTowerX}_{sTowerY}_{suffix}")
                if hprof:
                    h_meas = hprof.ProjectionX(
                        f"meas_proj_{var}_{sTowerX}_{sTowerY}_{suffix}", "E")
                    hists_meas[var] = h_meas
                else:
                    hists_meas[var] = None

                if hists_deconv[var] or hists_meas[var]:
                    channelNos[var] = (
                        f"({chan.board_no},{chan.group_no},{chan.channel_no})"
                    )

            if not any(hists_deconv.values()) and not any(hists_meas.values()):
                print(
                    f"Warning: No histograms for Board {board_no}, "
                    f"Tower ({i_tower_x}, {i_tower_y})")
                continue

            pave = create_pave_text(0.20, 0.75, 0.70, 0.90)
            pave.AddText(f"Tower: ({i_tower_x}, {i_tower_y})")

            # Build overlay: [deconv_Cer, meas_Cer, deconv_Sci, meas_Sci]
            # Colors repeat per var; linestyles alternate solid/dashed per pair
            hists_to_draw = []
            colors = []
            linestyles = []
            labels = []

            for var, color in [("Cer", 2 if isQuartz else 6), ("Sci", 4)]:
                h_d = hists_deconv.get(var)
                h_m = hists_meas.get(var)
                if not h_d and not h_m:
                    continue
                var_label = ("Quartz Cer" if isQuartz else "Plastic Cer") if var == "Cer" else "Sci"
                if var in channelNos:
                    pave.AddText(f"{var}: {channelNos[var]}")
                # Clone for drawing so originals stay unscaled for accumulation
                if h_d:
                    h_d_draw = h_d.Clone(h_d.GetName() + "_draw")
                    normalize_to_peak(h_d_draw)
                    hists_to_draw.append(h_d_draw)
                    colors.append(color)
                    linestyles.append(1)   # solid
                    labels.append(f"{var_label} (deconv)")
                if h_m:
                    h_m_draw = h_m.Clone(h_m.GetName() + "_draw")
                    normalize_to_peak(h_m_draw)
                    hists_to_draw.append(h_m_draw)
                    colors.append(color)
                    linestyles.append(2)   # dashed
                    labels.append(f"{var_label} (meas)")

            # Accumulate originals (unscaled) for combined summary plots
            if hists_deconv.get("Cer") and hists_deconv.get("Sci"):
                if isQuartz:
                    hists_Cer_Quartz.append(hists_deconv["Cer"])
                else:
                    hists_Cer_Plastic.append(hists_deconv["Cer"])
                hists_Cer.append(hists_deconv["Cer"])
                hists_Sci.append(hists_deconv["Sci"])
                if hists_meas.get("Cer"):
                    if isQuartz:
                        meas_Cer_Quartz.append(hists_meas["Cer"])
                    else:
                        meas_Cer_Plastic.append(hists_meas["Cer"])
                    meas_Cer.append(hists_meas["Cer"])
                if hists_meas.get("Sci"):
                    meas_Sci.append(hists_meas["Sci"])

            pm.plot_1d(
                hists_to_draw,
                f"deconv_DRS_vs_TS_{sTowerX}_{sTowerY}_{suffix}",
                xtitle,
                (-200, 100),
                ylabel=ytitle,
                yrange=(None, None),
                legends=labels,
                style=PlotStyle(
                    dology=False,
                    drawoptions="HIST",
                    mycolors=colors,
                    linestyles=linestyles,

                    addOverflow=False,
                    addUnderflow=False
                ),
                extraToDraw=pave
            )

    # Summary / combined plots: deconvolved (solid) + measured (dashed)
    if hists_Cer:
        hCer_Combined        = LHistos2Hist(hists_Cer,        f"deconv_DRS_vs_TS_Cer_Combined_{suffix}")
        hSci_Combined        = LHistos2Hist(hists_Sci,        f"deconv_DRS_vs_TS_Sci_Combined_{suffix}")
        hCer_Quartz_Combined = LHistos2Hist(hists_Cer_Quartz, f"deconv_DRS_vs_TS_Cer_Quartz_Combined_{suffix}")
        hCer_Plastic_Combined= LHistos2Hist(hists_Cer_Plastic,f"deconv_DRS_vs_TS_Cer_Plastic_Combined_{suffix}")

        mSci_Combined        = LHistos2Hist(meas_Sci,         f"meas_DRS_vs_TS_Sci_Combined_{suffix}")
        mCer_Quartz_Combined = LHistos2Hist(meas_Cer_Quartz,  f"meas_DRS_vs_TS_Cer_Quartz_Combined_{suffix}")
        mCer_Plastic_Combined= LHistos2Hist(meas_Cer_Plastic, f"meas_DRS_vs_TS_Cer_Plastic_Combined_{suffix}")

        def scale_to_deconv(h_meas, h_deconv):
            """Scale measured histogram so its peak matches the deconvolved peak."""
            if not h_meas or not h_deconv:
                return
            peak_d = max(abs(h_deconv.GetMaximum()), abs(h_deconv.GetMinimum()))
            peak_m = max(abs(h_meas.GetMaximum()), abs(h_meas.GetMinimum()))
            if peak_m > 0 and peak_d > 0:
                h_meas.Scale(peak_d / peak_m)

        scale_to_deconv(mCer_Quartz_Combined,  hCer_Quartz_Combined)
        scale_to_deconv(mCer_Plastic_Combined, hCer_Plastic_Combined)
        scale_to_deconv(mSci_Combined,         hSci_Combined)

        # Cer + Sci combined
        deconv_items = [(hCer_Quartz_Combined, mCer_Quartz_Combined, "Cer Quartz", 2),
                        (hCer_Plastic_Combined, mCer_Plastic_Combined, "Cer Plastic", 6),
                        (hSci_Combined, mSci_Combined, "Sci", 4)]
        cer_sci_hists, cer_sci_labels, cer_sci_colors, cer_sci_ls = [], [], [], []
        for h_d, h_m, label, color in deconv_items:
            if h_d:
                cer_sci_hists.append(h_d);  cer_sci_labels.append(f"{label} (deconv)")
                cer_sci_colors.append(color); cer_sci_ls.append(1)
            if h_m:
                cer_sci_hists.append(h_m);  cer_sci_labels.append(f"{label} (meas)")
                cer_sci_colors.append(color); cer_sci_ls.append(2)

        pm.plot_1d(
            cer_sci_hists,
            f"deconv_DRS_vs_TS_Cer_Sci_Combined_{suffix}",
            xtitle, (-75, 0),
            ylabel=ytitle, yrange=(None, None),
            legends=cer_sci_labels,
            style=PlotStyle(dology=False, drawoptions="HIST",
                            mycolors=cer_sci_colors, linestyles=cer_sci_ls,
                            addOverflow=False, addUnderflow=False,
                            legendPos=[0.25, 0.70, 0.40, 0.90]),
            prepend=True
        )

        # Cer-only combined
        cer_items = [(hCer_Quartz_Combined, mCer_Quartz_Combined, "Cer Quartz", 2),
                     (hCer_Plastic_Combined, mCer_Plastic_Combined, "Cer Plastic", 6)]
        cer_hists, cer_labels, cer_colors, cer_ls = [], [], [], []
        for h_d, h_m, label, color in cer_items:
            if h_d:
                cer_hists.append(h_d);  cer_labels.append(f"{label} (deconv)")
                cer_colors.append(color); cer_ls.append(1)
            if h_m:
                cer_hists.append(h_m);  cer_labels.append(f"{label} (meas)")
                cer_colors.append(color); cer_ls.append(2)

        if cer_hists:
            pm.plot_1d(
                cer_hists,
                f"deconv_DRS_vs_TS_Cer_Combined_{suffix}",
                xtitle, (-80, -50),
                ylabel=ytitle, yrange=(None, None),
                legends=cer_labels,
                style=PlotStyle(dology=False, drawoptions="HIST",
                                mycolors=cer_colors, linestyles=cer_ls,
                                addOverflow=False, addUnderflow=False,
                                legendPos=[0.30, 0.75, 0.40, 0.90]),
                prepend=True
            )

        # Save combined histograms
        outfile_name = f"{rootdir}/drs_vs_ts_deconvolved_combined_{suffix}.root"
        outfile = ROOT.TFile(outfile_name, "RECREATE")
        for h in [hCer_Combined, hSci_Combined,
                  hCer_Quartz_Combined, hCer_Plastic_Combined]:
            if h:
                h.SetDirectory(outfile)
                h.Write()
        outfile.Close()

    pm.add_newline()
    infile_deconv.Close()
    infile_meas.Close()

    return pm.generate_html(
        f"DRS/{suffix}/DRS_VS_TS_deconvolved.html",
        plots_per_row=4,
        title=f"Deconvolved DRS ADC vs TS for {suffix} (3mm boards)",
        intro_text=(
            "Per-channel DRS pulse profiles for 3mm SiPM boards. "
            "Solid: Wiener-deconvolved with the average pulse shape. "
            "Dashed: measured profile. Both normalised to unit area. "
            f"Pulse template: {DEFAULT_PULSE_FILE}."
        )
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    suffixes = ["inclusive", "pion"]

    for suffix in suffixes:
        print(f"\n--- Deconvolving suffix: {suffix} ---")
        deconvolveDRSvsTS(suffix=suffix)

        print(f"--- Plotting suffix: {suffix} ---")
        html = makeDRSDeconvPlots(suffix=suffix)
        print(f"  HTML: {html}")


if __name__ == "__main__":
    main()
