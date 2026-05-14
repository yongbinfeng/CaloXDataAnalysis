"""
Histogram-booking functions for each CaloX sequence.

Each public function has the signature:
    book_<name>(ctx: CaloXAnalysisManager) -> (lazy_list, save_fn)

lazy_list  — flat list of lazy ROOT objects for ROOT.RDF.RunGraphs
save_fn()  — writes ROOT files / JSONs after the graph has been triggered
"""

from configs.selection_config import get_service_drs_cut
from configs.plot_config import get_service_drs_processed_info_ranges
import json
from configs.plot_config import get_drs_plot_ranges, get_fers_saturation_value
from utils.utils import number_to_string
from variables.drs import get_ts_arr_name
from channels.channel_map import get_mcp_channels, get_pid_channels

_N_EVENTS = 60000
_N_BINS_EVENT = 500


# ---------------------------------------------------------------------------
# FERS: monitor conditions (HV, current, temperature vs event)
# ---------------------------------------------------------------------------

def book_monitor_conditions(ctx):
    hprofs = []
    for fersboard in ctx.fersboards.values():
        hprofs.append(ctx.rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_hv_name()}_VS_Event",
            "FERS Board - SipmHV VS Event;Event;SipmHV (V)",
            _N_BINS_EVENT, 0, _N_EVENTS),
            "event_n", fersboard.get_sipm_hv_name()))
        hprofs.append(ctx.rdf.Profile1D((
            f"hprof_{fersboard.get_sipm_i_name()}_VS_Event",
            "FERS Board - SipmI VS Event;Event;SipmI (mA)",
            _N_BINS_EVENT, 0, _N_EVENTS),
            "event_n", fersboard.get_sipm_i_name()))
        hprofs.append(ctx.rdf.Profile1D((
            f"hprof_{fersboard.get_temp_det_name()}_VS_Event",
            "FERS Board - TempDET VS Event;Event;TempDET (C)",
            _N_BINS_EVENT, 0, _N_EVENTS),
            "event_n", fersboard.get_temp_det_name()))
        hprofs.append(ctx.rdf.Profile1D((
            f"hprof_{fersboard.get_temp_fpga_name()}_VS_Event",
            "FERS Board - TempFPGA VS Event;Event;TempFPGA (C)",
            _N_BINS_EVENT, 0, _N_EVENTS),
            "event_n", fersboard.get_temp_fpga_name()))

    ctx.hbook.add("conditions_vs_event.root", hprofs)


# ---------------------------------------------------------------------------
# FERS: energy sum vs event
# ---------------------------------------------------------------------------

def book_fers_energy_sum(ctx):
    hprofs = []
    for fersboard in ctx.fersboards.values():
        for var in ["Cer", "Sci"]:
            for gain in ["HG", "LG"]:
                hprofs.append(ctx.rdf.Profile1D((
                    f"hprof_{fersboard.get_energy_sum_name(gain=gain, isCer=(var == 'Cer'))}_VS_Event",
                    f"FERS Board - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                    _N_BINS_EVENT, 0, _N_EVENTS),
                    "event_n",
                    fersboard.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))))

    for var in ["Cer", "Sci"]:
        for gain in ["HG", "LG"]:
            hprofs.append(ctx.rdf.Profile1D((
                f"hprof_{ctx.fersboards.get_energy_sum_name(gain=gain, isCer=(var == 'Cer'))}_VS_Event",
                f"FERS - {var} Energy {gain} sum VS Event;Event;{var} Energy {gain} sum",
                _N_BINS_EVENT, 0, _N_EVENTS),
                "event_n",
                ctx.fersboards.get_energy_sum_name(gain=gain, isCer=(var == "Cer"))))

    ctx.hbook.add("fers_energysum_vs_event.root", hprofs)


# ---------------------------------------------------------------------------
# FERS: per-channel 1D distributions + pedestals
# ---------------------------------------------------------------------------

def _make_fers_1d_hists(ctx, gain):
    hists = []
    for fersboard in ctx.fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hists.append(ctx.rdf.Histo1D((
                    f"hist_FERS_Board{board_no}_{var}_{s_x}_{s_y}",
                    f"FERS Board {board_no} - {var} iTowerX {s_x} iTowerY {s_y};"
                    f"{var} Energy;Counts",
                    3000, 0, 9000),
                    chan.get_channel_name(gain=gain)))
    return hists


def _collect_pedestals(ctx, hists, gain):
    """Compute per-channel pedestal from the peak of the low-range bin content."""
    pedestals = {}
    hist_by_name = {}
    for h in hists:
        actual = h.GetPtr() if hasattr(h, "GetPtr") else h
        hist_by_name[actual.GetName()] = actual

    for fersboard in ctx.fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                channel_name = chan.get_channel_name(gain=gain)
                hname = f"hist_FERS_Board{board_no}_{var}_{s_x}_{s_y}"
                hist = hist_by_name.get(hname)
                if hist is None:
                    print(
                        f"Warning: {hname} not found for pedestal extraction")
                    pedestals[channel_name] = None
                    continue
                ibinmin = hist.GetXaxis().FindBin(100.0)
                ibinmax = hist.GetXaxis().FindBin(500.0)
                max_val, ibin_ped = -1, -1
                for ibin in range(ibinmin, ibinmax):
                    if hist.GetBinContent(ibin) > max_val:
                        max_val = hist.GetBinContent(ibin)
                        ibin_ped = ibin
                pedestals[channel_name] = hist.GetXaxis(
                ).GetBinCenter(ibin_ped)
    return pedestals


def book_fers_channels(ctx):
    hists_hg = _make_fers_1d_hists(ctx, gain="HG")
    hists_lg = _make_fers_1d_hists(ctx, gain="LG")

    def post_save():
        pedestals_hg = _collect_pedestals(ctx, hists_hg, gain="HG")
        pedestals_lg = _collect_pedestals(ctx, hists_lg, gain="LG")
        with open(f"{ctx.paths['root']}/fers_pedestals_hg.json", "w") as f:
            json.dump(pedestals_hg, f, indent=4)
        with open(f"{ctx.paths['root']}/fers_pedestals_lg.json", "w") as f:
            json.dump(pedestals_lg, f, indent=4)

    ctx.hbook.add("fers_all_channels_1d.root", hists_hg + hists_lg, post_save)


# ---------------------------------------------------------------------------
# FERS: channel statistics (mean, max, saturation frequency)
# ---------------------------------------------------------------------------

def book_fers_stats(ctx):
    saturation_value = get_fers_saturation_value()
    stats_lazy = {}
    for fersboard in ctx.fersboards.values():
        for chan in fersboard:
            for gain in ["HG", "LG"]:
                ch = chan.get_channel_name(gain=gain)
                stats_lazy[ch] = (
                    ctx.rdf.Mean(ch),
                    ctx.rdf.Max(ch),
                    ctx.rdf.Filter(f"{ch} >= {saturation_value}").Count(),
                )
    all_lazy = [item for triple in stats_lazy.values() for item in triple]

    def post_save():
        import os
        results = {}
        for ch, (mean, vmax, sat_count) in stats_lazy.items():
            results[ch] = (
                mean.GetValue(),
                vmax.GetValue(),
                float(sat_count.GetValue()) / _N_EVENTS,
            )
        os.makedirs(ctx.paths["root"], exist_ok=True)
        with open(f"{ctx.paths['root']}/fers_stats.json", "w") as f:
            json.dump(results, f, indent=4)

    ctx.hbook.add(None, all_lazy, post_save)


# ---------------------------------------------------------------------------
# FERS: max readout value per board per event
# ---------------------------------------------------------------------------

def book_fers_max(ctx):
    hists = []
    nbins, xmin, xmax = 100, 7500, 8500
    for fersboard in ctx.fersboards.values():
        for gain in ["HG", "LG"]:
            for isCer, var in [(True, "CER"), (False, "SCI")]:
                name = fersboard.get_energy_max_name(gain=gain, isCer=isCer)
                hists.append(ctx.rdf.Histo1D((
                    f"hist_{name}",
                    f"FERS Board - {var} Energy {gain} max;{var} Energy {gain} max per board;Counts",
                    nbins, xmin, xmax),
                    name))
    for gain in ["HG", "LG"]:
        for isCer, var in [(True, "CER"), (False, "SCI")]:
            name = ctx.fersboards.get_energy_max_name(gain=gain, isCer=isCer)
            hists.append(ctx.rdf.Histo1D((
                f"hist_{name}",
                f"FERS - {var} Energy {gain} max;{var} Energy {gain} max;Counts",
                nbins, xmin, xmax),
                name))

    ctx.hbook.add("fers_max_values.root", hists)


# ---------------------------------------------------------------------------
# FERS: Cer vs Sci and HG vs LG 2D correlations
# ---------------------------------------------------------------------------

def book_fers_2d(ctx):
    hists = []
    for fersboard in ctx.fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            chan_cer = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=True)
            chan_sci = fersboard.get_channel_by_tower(
                i_tower_x, i_tower_y, isCer=False)
            i_cer = chan_cer.channel_no
            i_sci = chan_sci.channel_no

            hists.append(ctx.rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{s_x}_{s_y}",
                f"CER {i_cer} VS SCI {i_sci} iTowerX {s_x} iTowerY {s_y};"
                f"Sci Energy HG;Cer Energy HG",
                300, 0, 9000, 300, 0, 9000),
                chan_sci.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="HG")))
            hists.append(ctx.rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_VS_Sci_{s_x}_{s_y}_zoom",
                f"CER {i_cer} VS SCI {i_sci} iTowerX {s_x} iTowerY {s_y} (zoomed);"
                f"SCI Energy HG;CER Energy HG",
                300, 0, 1000, 200, 0, 2000),
                chan_sci.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="HG")))
            hists.append(ctx.rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Sci_{s_x}_{s_y}_hg_VS_lg",
                f"SCI {i_sci} LG VS HG;SCI Energy HG;SCI Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_sci.get_channel_name(gain="HG"),
                chan_sci.get_channel_name(gain="LG")))
            hists.append(ctx.rdf.Histo2D((
                f"hist_FERS_Board{board_no}_Cer_{s_x}_{s_y}_hg_VS_lg",
                f"CER {i_cer} LG VS HG;CER Energy HG;CER Energy LG",
                300, 0, 9000, 300, 0, 3000),
                chan_cer.get_channel_name(gain="HG"),
                chan_cer.get_channel_name(gain="LG")))

    ctx.hbook.add("fers_all_channels_2d.root", hists)


# ---------------------------------------------------------------------------
# FERS: per-channel vs event 2D tracking
# ---------------------------------------------------------------------------

def book_fers_track(ctx):
    hists = []
    for fersboard in ctx.fersboards.values():
        board_no = fersboard.board_no
        for i_tower_x, i_tower_y in fersboard.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan = fersboard.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                hists.append(ctx.rdf.Histo2D((
                    f"hist_FERS_Board{board_no}_{var}_VS_Event_{s_x}_{s_y}",
                    f"FERS Board {board_no} - Event VS {var} {chan.channel_no} "
                    f"iTowerX {s_x} iTowerY {s_y};Event;{var} Energy HG",
                    _N_BINS_EVENT, 0, _N_EVENTS, 1000, 0, 9000),
                    "event_n", chan.get_channel_name(gain="HG")))

    ctx.hbook.add("fers_all_channels_2D_VS_event.root", hists)


# ---------------------------------------------------------------------------
# DRS: waveform 2D histograms (vs TS, raw / ref / mcp)
# ---------------------------------------------------------------------------

def _book_drs_channel_waveforms(ctx, channel_name, ymin, ymax):
    hists = []
    ch_blsub = f"{channel_name}_blsub"
    hists.append(ctx.rdf.Histo2D((
        f"hist_{ch_blsub}_VS_ts",
        "DRS values vs raw TS;TS;DRS values",
        1024, 0, 1024, 50, ymin, ymax),
        "ts", ch_blsub))
    hists.append(ctx.rdf.Profile1D((
        f"prof_{ch_blsub}_VS_ts",
        "Mean DRS values vs raw TS;TS;Mean DRS values",
        1024, 0, 1024),
        "ts", ch_blsub))
    hists.append(ctx.rdf.Histo2D((
        f"hist_{ch_blsub}_VS_ts_ref",
        "DRS values vs ref-corrected TS;ts_ref;DRS values",
        1024, 0, 1024, 50, ymin, ymax),
        get_ts_arr_name(channel_name), ch_blsub))
    hists.append(ctx.rdf.Profile1D((
        f"prof_{ch_blsub}_VS_ts_ref",
        "Mean DRS values vs ref-corrected TS;ts_ref;Mean DRS values",
        1024, 0, 1024),
        get_ts_arr_name(channel_name), ch_blsub))
    hists.append(ctx.rdf.Histo2D((
        f"hist_{ch_blsub}_VS_ts_mcp",
        "DRS values vs ref+MCP-corrected TS;ts_mcp;DRS values",
        1024, 0, 1024, 50, ymin, ymax),
        get_ts_arr_name(channel_name, use_mcp=True), ch_blsub))
    hists.append(ctx.rdf.Profile1D((
        f"prof_{ch_blsub}_VS_ts_mcp",
        "Mean DRS values vs ref+MCP-corrected TS;ts_mcp;Mean DRS values",
        1024, 0, 1024),
        get_ts_arr_name(channel_name, use_mcp=True), ch_blsub))
    return hists


def book_drs_waveforms(ctx):
    hists = []
    for _, board in ctx.drsboards.items():
        for chan in board:
            if chan is None:
                continue
            ymin, ymax = get_drs_plot_ranges(
                subtractMedian=True,
                is_amplified=chan.is_amplified,
                is6mm=chan.is6mm,
                is_reference=chan.is_reference)
            hists.extend(_book_drs_channel_waveforms(
                ctx, chan.get_channel_name(blsub=False), ymin, ymax))

    for _, mcp_channel in get_mcp_channels(ctx.run_number).items():
        ymin, ymax = get_drs_plot_ranges(subtractMedian=True, is_mcp=True)
        hists.extend(_book_drs_channel_waveforms(ctx, mcp_channel, ymin, ymax))

    ctx.hbook.add("drs_vs_ts.root", hists)


# ---------------------------------------------------------------------------
# DRS: pulse statistics (peak value, energy, timing)
# ---------------------------------------------------------------------------

def book_drs_stats(ctx):
    hists = []
    for _, board in ctx.drsboards.items():
        for chan in board:
            if chan.is_reference:
                continue
            channel_name = chan.get_channel_name(blsub=False)
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_peak_value",
                "DRS peak value;Peak value;Counts",
                200, 0, 800),
                f"{channel_name}_peak_value"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_energy",
                "DRS CFD energy integral;CFD energy;Counts",
                650, -500, 6000),
                f"{channel_name}_energy"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_TS_peak_ref",
                "DRS peak TS (ref-corrected);TS_peak_ref;Counts",
                1024, 0, 1024),
                f"{channel_name}_TS_peak_ref"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_TS_cfd_ref",
                "DRS CFD TS (ref-corrected);TS_cfd_ref;Counts",
                1024, 0, 1024),
                f"{channel_name}_TS_cfd_ref"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_TS_peak_mcp",
                "DRS peak TS (ref+MCP-corrected);TS_peak_mcp;Counts",
                1024, 0, 1024),
                f"{channel_name}_TS_peak_mcp"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_TS_cfd_mcp",
                "DRS CFD TS (ref+MCP-corrected);TS_cfd_mcp;Counts",
                1024, 0, 1024),
                f"{channel_name}_TS_cfd_mcp"))
            hists.append(ctx.rdf.Histo1D((
                f"hist_{channel_name}_TS_cfd_mcp_finebins",
                "DRS CFD TS (ref+MCP-corrected);TS_cfd_mcp;Counts",
                500, 420, 520),
                f"{channel_name}_TS_cfd_mcp"))

    ctx.hbook.add("drs_stats.root", hists)


# ---------------------------------------------------------------------------
# DRS: peak time-slice distributions
# ---------------------------------------------------------------------------

def book_drs_peak_ts(ctx):
    hists_cer, hists_sci, hists_2d = [], [], []
    for _, board in ctx.drsboards.items():
        board_no = board.board_no
        for i_tower_x, i_tower_y in board.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            channel_names = {}
            for var in ["Cer", "Sci"]:
                chan = board.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan is None:
                    if var == "Cer":
                        print(f"Warning: DRS channel not found for "
                              f"Board{board_no}, Tower({s_x},{s_y}), {var}")
                    continue
                ch_ts = chan.get_channel_peak_ts_name()
                channel_names[var] = ch_ts
                h1 = ctx.rdf.Histo1D((
                    f"hist_DRSPeakTS_{var}_{s_x}_{s_y}",
                    f"DRS Peak TS Board{board_no} Tower({s_x},{s_y}) {var};"
                    f"Peak TS;Counts",
                    1000, 0, 1000),
                    f"{ch_ts}_good")
                (hists_cer if var == "Cer" else hists_sci).append(h1)

            if len(channel_names) < 2:
                continue
            hists_2d.append(ctx.rdf.Histo2D((
                f"hist_DRSPeakTS_Cer_VS_Sci_{s_x}_{s_y}",
                f"DRS Peak TS CER VS SCI Board{board_no} Tower({s_x},{s_y});"
                f"SCI Peak TS;CER Peak TS",
                1000, 0, 1000, 1000, 0, 1000),
                f'{channel_names["Sci"]}_good',
                f'{channel_names["Cer"]}_good'))

    all_hists = hists_cer + hists_sci + hists_2d

    ctx.hbook.add("drspeakts.root", all_hists)


# ---------------------------------------------------------------------------
# DRS ↔ FERS correlations
# ---------------------------------------------------------------------------

def book_drs_sum_vs_fers(ctx):
    xymax = {"Cer": (20000, 8500), "Sci": (30000, 8500)}
    xymax_LG = {"Cer": (20000, 2000), "Sci": (30000, 4000)}
    hists = []
    for _, board in ctx.drsboards.items():
        board_no = board.board_no
        for i_tower_x, i_tower_y in board.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan_drs = board.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_drs is None:
                    if var == "Cer":
                        print(f"Warning: DRS channel not found for "
                              f"Board{board_no}, Tower({s_x},{s_y}), {var}")
                    continue
                chan_fers = None
                for fb in ctx.fersboards.values():
                    chan_fers = fb.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_fers is not None:
                        break
                if chan_fers is None:
                    print(f"Warning: FERS channel not found for "
                          f"Board{board_no}, Tower({s_x},{s_y}), {var}")
                    continue
                hists.append(ctx.rdf.Histo2D((
                    f"hist_DRSSum_VS_FERS_{var}_{s_x}_{s_y}",
                    f"DRS sum VS FERS HG Board{board_no} Tower({s_x},{s_y}) {var};"
                    f"FERS Energy HG;DRS Sum",
                    100, 0, xymax[var][1], 100, 0, xymax[var][0]),
                    chan_fers.get_channel_name(gain="HG"),
                    chan_drs.get_channel_sum_name()))
                hists.append(ctx.rdf.Histo2D((
                    f"hist_DRSSum_VS_FERSLG_{var}_{s_x}_{s_y}",
                    f"DRS sum VS FERS LG Board{board_no} Tower({s_x},{s_y}) {var};"
                    f"FERS Energy LG;DRS Sum",
                    100, 0, xymax_LG[var][1], 100, 0, xymax_LG[var][0]),
                    chan_fers.get_channel_name(gain="LG"),
                    chan_drs.get_channel_sum_name()))

    ctx.hbook.add("drssum_vs_fers.root", hists)


def book_drs_peak_vs_fers(ctx):
    hists = []
    for _, board in ctx.drsboards.items():
        board_no = board.board_no
        for i_tower_x, i_tower_y in board.get_list_of_towers():
            s_x = number_to_string(i_tower_x)
            s_y = number_to_string(i_tower_y)
            for var in ["Cer", "Sci"]:
                chan_drs = board.get_channel_by_tower(
                    i_tower_x, i_tower_y, isCer=(var == "Cer"))
                if chan_drs is None:
                    if var == "Cer":
                        print(f"Warning: DRS channel not found for "
                              f"Board{board_no}, Tower({s_x},{s_y}), {var}")
                    continue
                chan_fers = None
                for fb in ctx.fersboards.values():
                    chan_fers = fb.get_channel_by_tower(
                        i_tower_x, i_tower_y, isCer=(var == "Cer"))
                    if chan_fers is not None:
                        break
                if chan_fers is None:
                    print(f"Warning: FERS channel not found for "
                          f"Board{board_no}, Tower({s_x},{s_y}), {var}")
                    continue
                _, ymax = get_drs_plot_ranges(
                    subtractMedian=True,
                    is_amplified=chan_drs.is_amplified,
                    is6mm=chan_drs.is6mm)
                hists.append(ctx.rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERS_{var}_{s_x}_{s_y}",
                    f"DRS Peak VS FERS HG Board{board_no} Tower({s_x},{s_y}) {var};"
                    f"FERS Energy HG;DRS Peak",
                    100, 0, 9000, 100, 0, ymax),
                    chan_fers.get_channel_name(gain="HG"),
                    chan_drs.get_channel_peak_name()))
                hists.append(ctx.rdf.Histo2D((
                    f"hist_DRSPeak_VS_FERSLG_{var}_{s_x}_{s_y}",
                    f"DRS Peak VS FERS LG Board{board_no} Tower({s_x},{s_y}) {var};"
                    f"FERS Energy LG;DRS Peak",
                    100, 0, 3000, 100, 0, ymax),
                    chan_fers.get_channel_name(gain="LG"),
                    chan_drs.get_channel_peak_name()))

    ctx.hbook.add("drspeak_vs_fers.root", hists)


# ===========================================================================
# Service DRS sequences  (ctx: CaloXAnalysisManager)
# ===========================================================================


def _analyze_pulse(rdf, channels):
    """Book all histogram lazy objects for one set of service-DRS channels.

    Returns a flat list of lazy ROOT objects.
    """
    hists = {}

    for det, channel in channels.items():
        hists[channel] = {}
        wf_ymin, wf_ymax = get_service_drs_processed_info_ranges(
            det, "waveform")
        hists[channel]["ADC_VS_TS"] = rdf.Histo2D(
            (f"{det}_ADC_vs_TS",
             f"ADC vs Time Slice {channel};Time Slice;ADC Counts;Counts",
             1024, 0, 1024, 500, wf_ymin, wf_ymax),
            "ts", f"{channel}_blsub")
        hists[channel]["ADC_VS_TS_prof"] = rdf.Profile1D(
            (f"{det}_ADC_vs_TS_prof",
             f"Mean ADC vs Time Slice {channel};Time Slice;Mean ADC Counts",
             1024, 0, 1024, wf_ymin, wf_ymax),
            "ts", f"{channel}_blsub")
        hists[channel]["cfd_ts"] = rdf.Histo1D(
            (f"{det}_cfd_ts",
             f"CFD TS {channel};CFD TS;Counts", 128, 0, 1024),
            f"{channel}_TS_cfd")
        xmin, xmax = get_service_drs_processed_info_ranges(det, "peak_value")
        hists[channel]["peak_value_vs_cfd_ts"] = rdf.Histo2D(
            (f"{det}_peak_value_vs_cfd_ts",
             f"Peak Value vs CFD TS {channel};CFD TS;Peak Value;Counts",
             128, 0, 1024, 100, xmin, xmax),
            f"{channel}_TS_cfd", f"{channel}_peak_value")
        hists[channel]["peak_value"] = rdf.Histo1D(
            (f"{det}_peak_value",
             f"Peak Value {channel};ADC Counts;Counts", 50, xmin, xmax),
            f"{channel}_peak_value")
        xmin, xmax = get_service_drs_processed_info_ranges(det, "sum")
        hists[channel]["energy"] = rdf.Histo1D(
            (f"{det}_energy",
             f"Energy {channel};Energy (ADC);Counts", 500, xmin, xmax),
            f"{channel}_energy")
        hists[channel]["integral_to_peak"] = rdf.Histo1D(
            (f"{det}_integral_to_peak",
             f"Integral/Peak {channel};Integral/Peak [TS];Counts", 120, -5, 20),
            f"{channel}_integral_to_peak")

    return [h for hmap in hists.values() for h in hmap.values()]


def _analyze_detector_pair_correlations(rdf, channels):
    """Book 2D detector-pair correlation histograms for all valid channel pairs."""
    hists = {}
    det_list = list(channels.keys())
    for idx1, det1 in enumerate(det_list):
        for idx2, det2 in enumerate(det_list):
            if idx2 <= idx1:
                continue
            special_dets = {"KT1", "KT2", "T3", "T4"}
            if (det1 in special_dets) ^ (det2 in special_dets):
                continue
            channel1, channel2 = channels[det1], channels[det2]
            _, _, _, method1 = get_service_drs_cut(det1)
            _, _, _, method2 = get_service_drs_cut(det2)
            var1 = "peak_value" if method1 == "PeakValue" else "energy"
            var2 = "peak_value" if method2 == "PeakValue" else "energy"
            xmin, xmax = get_service_drs_processed_info_ranges(
                det1, "peak_value" if method1 == "PeakValue" else "sum")
            ymin, ymax = get_service_drs_processed_info_ranges(
                det2, "peak_value" if method2 == "PeakValue" else "sum")
            hists[f"{det1}_vs_{det2}"] = rdf.Histo2D(
                (f"{det1}_{var1}_vs_{det2}_{var2}",
                 f"{det1} vs {det2};{det1} {var1};{det2} {var2};Counts",
                 500, xmin, xmax, 500, ymin, ymax),
                f"{channel1}_{var1}", f"{channel2}_{var2}")

    return list(hists.values())


def _analyze_mcp_timing_diff(rdf, channels_mcp):
    """Book MCP CFD timing-difference histograms for all unique pairs."""
    hists = []
    dets = list(channels_mcp.keys())

    for i, det1 in enumerate(dets):
        col1 = channels_mcp[det1]
        for det2 in dets[i + 1:]:
            col2 = channels_mcp[det2]
            diff_col = f"{col1}_minus_{col2}_cfd_diff"
            rdf = rdf.Define(
                diff_col, f"{col1}_TS_cfd_ref - {col2}_TS_cfd_ref")
            hists.append(rdf.Histo1D(
                (f"{det1}_cfd_diff_vs_{det2}",
                 f"{det1} - {det2};#Delta t_{{CFD,ref}} [TS];Counts",
                 800, -10, 10),
                diff_col))

    return hists


def _analyze_hodo_peak(rdf, run_number):
    """Book hodoscope peak timing histograms."""
    from channels.channel_map import build_hodo_pos_channels
    hodo_pos_channels = build_hodo_pos_channels(run=run_number)

    histos = {key: {} for key in [
        "diff", "diff_relative", "sum", "left", "right",
        "left_vs_right", "left_peak", "right_peak", "LR_vs_UD"
    ]}
    conditions = ["is_HoleVeto_vetoed", "passNone"]
    rdf_filtered = rdf

    for group, channels in hodo_pos_channels.items():
        for channel in channels:
            rdf_filtered = rdf_filtered.Define(
                f"{channel}_subtracted_sum",
                f"ROOT::VecOps::Sum({channel}_blsub) + 1e-6")
            rdf_filtered = rdf_filtered.Define(
                f"{channel}_subtracted_norm",
                f"{channel}_blsub / {channel}_subtracted_sum")

        rdf_filtered = rdf_filtered.Define(
            f"{group}_delta_peak",
            f"(int){channels[1]}_peak_position - (int){channels[0]}_peak_position")
        rdf_filtered = rdf_filtered.Define(
            f"{group}_sum_peak",
            f"(int){channels[0]}_peak_position + (int){channels[1]}_peak_position")
        rdf_filtered = rdf_filtered.Define(
            f"{group}_delta_peak_relative",
            f"{group}_delta_peak / ({group}_sum_peak + 1e-6)")

        for sel in conditions:
            cat = f"{group}_{sel}"
            histos["diff"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak", f"Delta Peak {cat}", 256, -1024, 1024),
                f"{group}_delta_peak", sel)
            histos["diff_relative"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_delta_peak_relative",
                 f"Delta Peak Relative {cat}", 256, -1, 1),
                f"{group}_delta_peak_relative", sel)
            histos["sum"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_sum_peak", f"Sum Peak {cat}", 256, 0, 2048),
                f"{group}_sum_peak", sel)
            histos["left"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak", f"Left Peak {cat}", 256, 0, 1024),
                f"{channels[0]}_peak_position", sel)
            histos["right"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak", f"Right Peak {cat}", 256, 0, 1024),
                f"{channels[1]}_peak_position", sel)
            histos["left_vs_right"][cat] = rdf_filtered.Histo2D(
                (f"{cat}_left_peak_vs_right_peak",
                 f"Left vs Right Peak {cat}",
                 256, 0, 1024, 256, 0, 1024),
                f"{channels[0]}_peak_position", f"{channels[1]}_peak_position", sel)
            histos["left_peak"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_left_peak_value",
                 f"Left Peak Value {cat}", 200, -2500, 1999),
                f"{channels[0]}_peak_value", sel)
            histos["right_peak"][cat] = rdf_filtered.Histo1D(
                (f"{cat}_right_peak_value",
                 f"Right Peak Value {cat}", 200, -2500, 1999),
                f"{channels[1]}_peak_value", sel)

    for sel in conditions:
        for dwc in ["LR1_vs_UD1", "LR2_vs_UD2"]:
            if dwc[:3] not in hodo_pos_channels:
                continue
            cat = f"{dwc}_{sel}"
            histos["LR_vs_UD"][cat] = rdf_filtered.Histo2D(
                (cat, "LR vs UD", 256, -1, 1, 256, -1, 1),
                f"{dwc.split('_vs_')[0]}_delta_peak_relative",
                f"{dwc.split('_vs_')[1]}_delta_peak_relative", sel)

    return [h for hmap in histos.values() for h in hmap.values()]


# ---------------------------------------------------------------------------
# Public booking functions (ctx: CaloXAnalysisManager)
# ---------------------------------------------------------------------------

def book_service_drs_pid(ctx):
    pid_channels = get_pid_channels(ctx.run_number)
    hists = _analyze_pulse(ctx.rdf, pid_channels)
    hists += _analyze_detector_pair_correlations(ctx.rdf, pid_channels)
    ctx.hbook.add("drs_services.root", hists)


def book_service_drs_mcp(ctx):
    hists = _analyze_pulse(ctx.rdf, get_mcp_channels(ctx.run_number))
    ctx.hbook.add("drs_mcp.root", hists)


def book_service_drs_mcp_timing(ctx):
    # ctx.rdf is already the MCP-filtered view (set by define_selection in registry)
    channels_mcp = get_mcp_channels(ctx.run_number)
    hists = _analyze_mcp_timing_diff(ctx.rdf, channels_mcp)
    ctx.hbook.add("drs_mcp_timing_diff.root", hists)


def book_service_drs_hodo(ctx):
    hists = _analyze_hodo_peak(ctx.rdf, ctx.run_number)
    ctx.hbook.add("hodoscope_peaks.root", hists)
