"""
Sequence registry.

There are three sequence lists:

ALL_SEQUENCES         — main DQM (FERS + DRS readout), used by
                        make_dqm_hists.py and make_dqm_plots.py
SERVICE_DRS_SEQUENCES — service DRS (PID detectors, MCP, hodoscope), used by
                        check_service_drs.py
DRS_MCP_SEQUENCES     — DRS-only sequences run under MCP timing selection,
                        used by check_drs_mcp.py

To add a new monitoring variable to any list:
  1. Write book_<name>(ctx) in analysis/hist_functions.py
  2. Write plot_<name>(ctx) in analysis/plot_functions.py
  3. Add a CaloXSequence entry in the appropriate list below.

Set enabled_by_default=False to keep a sequence opt-in (--sequences on CLI).
"""

from core.sequence import CaloXSequence
from analysis.hist_functions import (
    book_monitor_conditions,
    book_fers_energy_sum,
    book_fers_channels,
    book_fers_stats,
    book_fers_max,
    book_fers_2d,
    book_fers_track,
    book_drs_waveforms,
    book_drs_stats,
    book_drs_peak_ts,
    book_drs_sum_vs_fers,
    book_drs_peak_vs_fers,
    book_service_drs_pid,
    book_service_drs_mcp,
    book_service_drs_mcp_timing,
    book_service_drs_hodo,
)
from analysis.plot_functions import (
    plot_monitor_conditions,
    plot_fers_energy_sum,
    plot_fers_mapping,
    plot_fers_channels,
    plot_fers_stats,
    plot_fers_max,
    plot_fers_2d,
    plot_fers_track,
    plot_drs_waveforms,
    plot_drs_stats,
    plot_drs_peak_ts,
    plot_drs_sum_vs_fers,
    plot_drs_peak_vs_fers,
    plot_service_drs_pid,
    plot_service_drs_mcp,
    plot_service_drs_mcp_timing,
    plot_service_drs_hodo,
)

# ---------------------------------------------------------------------------
# Canonical sequence list — order determines execution order.
# ---------------------------------------------------------------------------

ALL_SEQUENCES: list[CaloXSequence] = [
    # --- FERS monitoring ---
    CaloXSequence(
        name="monitor_conditions",
        book_hists=book_monitor_conditions,
        make_plots=plot_monitor_conditions,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="fers_energy_sum",
        book_hists=book_fers_energy_sum,
        make_plots=plot_fers_energy_sum,
        enabled_by_default=True,
    ),
    # plot-only: channel layout maps (no histogram step)
    CaloXSequence(
        name="fers_mapping",
        book_hists=None,
        make_plots=plot_fers_mapping,
        enabled_by_default=True,
    ),
    # per-channel 1D distributions + pedestal extraction
    CaloXSequence(
        name="fers_channels",
        book_hists=book_fers_channels,
        make_plots=plot_fers_channels,
        enabled_by_default=False,  # off by default; needed by fers_stats plots
    ),
    # board-level statistics (mean, max, saturation, pedestal maps)
    CaloXSequence(
        name="fers_stats",
        book_hists=book_fers_stats,
        make_plots=plot_fers_stats,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="fers_max",
        book_hists=book_fers_max,
        make_plots=plot_fers_max,
        enabled_by_default=False,
    ),
    CaloXSequence(
        name="fers_2d",
        book_hists=book_fers_2d,
        make_plots=plot_fers_2d,
        enabled_by_default=False,
    ),
    CaloXSequence(
        name="fers_track",
        book_hists=book_fers_track,
        make_plots=plot_fers_track,
        enabled_by_default=False,
    ),

    # --- DRS ---
    CaloXSequence(
        name="drs_waveforms",
        book_hists=book_drs_waveforms,
        make_plots=plot_drs_waveforms,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_stats",
        book_hists=book_drs_stats,
        make_plots=plot_drs_stats,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_peak_ts",
        book_hists=book_drs_peak_ts,
        make_plots=plot_drs_peak_ts,
        enabled_by_default=False,
    ),

    # --- DRS ↔ FERS correlations ---
    CaloXSequence(
        name="drs_sum_vs_fers",
        book_hists=book_drs_sum_vs_fers,
        make_plots=plot_drs_sum_vs_fers,
        enabled_by_default=False,
    ),
    CaloXSequence(
        name="drs_peak_vs_fers",
        book_hists=book_drs_peak_vs_fers,
        make_plots=plot_drs_peak_vs_fers,
        enabled_by_default=False,
    ),
]

# ---------------------------------------------------------------------------
# Service DRS sequences  (ctx: CaloXAnalysisManager)
# ---------------------------------------------------------------------------


def _select_mcp(ctx):
    """Filter to events where all MCPs have clean pulses."""
    from channels.channel_map import get_mcp_channels
    rdf = ctx.rdf
    for det in get_mcp_channels(ctx.run_number):
        rdf = rdf.Filter(
            f"{det}_integral_to_peak > 4 && {det}_peak_value > 10"
        )
    return rdf


def _select_mcp_diff(ctx):
    """Filter to events where all MCPs have clean pulses and small CFD timing differences."""
    from channels.channel_map import get_mcp_channels
    rdf = ctx.rdf
    rdf = _select_mcp(ctx)  # apply clean pulse filter first

    # require MCP timing difference between two MCPs
    dets = list(get_mcp_channels(ctx.run_number).keys())
    det1 = dets[0]
    det2 = dets[4]
    # for i in range(len(dets)):
    #    for j in range(i + 1, len(dets)):
    #        det1, det2 = dets[i], dets[j]
    rdf = rdf.Filter(
        f"abs({det1}_TS_cfd_ref - {det2}_TS_cfd_ref) > 2 && abs({det1}_TS_cfd_ref - {det2}_TS_cfd_ref) < 5")
    return rdf


SERVICE_DRS_SEQUENCES: list[CaloXSequence] = [
    CaloXSequence(
        name="service_drs_pid",
        book_hists=book_service_drs_pid,
        make_plots=plot_service_drs_pid,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="service_drs_mcp",
        book_hists=book_service_drs_mcp,
        make_plots=plot_service_drs_mcp,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="service_drs_mcp_timing",
        define_selection=_select_mcp,
        book_hists=book_service_drs_mcp_timing,
        make_plots=plot_service_drs_mcp_timing,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="service_drs_hodo",
        book_hists=book_service_drs_hodo,
        make_plots=plot_service_drs_hodo,
        enabled_by_default=False,
    ),
]

# ---------------------------------------------------------------------------
# DRS-only sequences under MCP timing selection
# ---------------------------------------------------------------------------

DRS_MCP_SEQUENCES: list[CaloXSequence] = [
    CaloXSequence(
        name="drs_waveforms",
        define_selection=_select_mcp_diff,
        book_hists=book_drs_waveforms,
        make_plots=plot_drs_waveforms,
    ),
    CaloXSequence(
        name="drs_stats",
        define_selection=_select_mcp_diff,
        book_hists=book_drs_stats,
        make_plots=plot_drs_stats,
    ),
]

# ---------------------------------------------------------------------------
# Lookup helpers — work for all lists
# ---------------------------------------------------------------------------

_ALL_BY_NAME: dict[str, CaloXSequence] = {
    s.name: s for s in ALL_SEQUENCES + SERVICE_DRS_SEQUENCES
}


def get_sequence(name: str) -> CaloXSequence:
    if name not in _ALL_BY_NAME:
        available = ", ".join(_ALL_BY_NAME)
        raise ValueError(f"Unknown sequence '{name}'. Available: {available}")
    return _ALL_BY_NAME[name]


def default_sequences() -> list[CaloXSequence]:
    """Main DQM sequences enabled by default."""
    return [s for s in ALL_SEQUENCES if s.enabled_by_default]


def default_service_drs_sequences() -> list[CaloXSequence]:
    """Service DRS sequences enabled by default."""
    return [s for s in SERVICE_DRS_SEQUENCES if s.enabled_by_default]


def default_drs_mcp_sequences() -> list[CaloXSequence]:
    """DRS-only sequences run under MCP timing selection."""
    return list(DRS_MCP_SEQUENCES)


def select_sequences(names: list[str],
                     pool: list[CaloXSequence] = None) -> list[CaloXSequence]:
    """Return sequences matching *names*, preserving pool order.

    pool defaults to ALL_SEQUENCES + SERVICE_DRS_SEQUENCES.
    """
    if pool is None:
        pool = ALL_SEQUENCES + SERVICE_DRS_SEQUENCES
    selected_names = {get_sequence(n).name for n in names}
    return [s for s in pool if s.name in selected_names]
