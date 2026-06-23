"""
Sequence registry.

Sequence lists:

FERS_DQM_SEQUENCES    — FERS monitoring only; runs without DRS data.
                        Used by make_fers_dqm_hists.py / make_fers_dqm_plots.py
DRS_DQM_SEQUENCES     — DRS monitoring; requires DRS data.
                        Used by make_dqm_hists.py / make_dqm_plots.py
ALL_SEQUENCES         — FERS_DQM + DRS_DQM combined (backward compat)
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
    book_fers_esum_vs_event,
    book_fers_channels,
    book_fers_stats,
    book_fers_max,
    book_fers_2d,
    book_fers_track,
    book_fers_energy_sum,
    book_fers_energy_sum_board,
    book_fers_lg_vs_mix,
    book_fers_cer_vs_sci,
    book_fers_dr,
    book_fers_ewc,
    book_fers_ewc_vs_hodo,
    book_fers_shower_shape,
    book_drs_waveforms,
    book_drs_profiles,
    book_drs_prof_combined,
    book_drs_stats,
    book_drs_finebins_combined,
    book_drs_finebins_corr_combined,
    book_drs_prof_corr_combined,
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
    plot_fers_esum_vs_event,
    plot_fers_mapping,
    plot_drs_mapping,
    plot_fers_channels,
    plot_fers_stats,
    plot_fers_max,
    plot_fers_2d,
    plot_fers_track,
    plot_fers_energy_sum,
    plot_fers_energy_sum_board,
    plot_fers_lg_vs_mix,
    plot_fers_cer_vs_sci,
    plot_fers_dr,
    plot_fers_ewc,
    plot_fers_ewc_vs_hodo,
    plot_fers_shower_shape,
    plot_drs_waveforms,
    plot_drs_profiles,
    plot_drs_stats,
    plot_drs_cfd_mpv,
    plot_drs_peak_ts,
    plot_drs_sum_vs_fers,
    plot_drs_peak_vs_fers,
    plot_service_drs_pid,
    plot_service_drs_mcp,
    plot_service_drs_mcp_timing,
    plot_service_drs_hodo,
)

# ---------------------------------------------------------------------------
# FERS-only DQM sequences — run without DRS data.
# ---------------------------------------------------------------------------

FERS_DQM_SEQUENCES: list[CaloXSequence] = [
    CaloXSequence(
        name="monitor_conditions",
        book_hists=book_monitor_conditions,
        make_plots=plot_monitor_conditions,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="fers_esum_vs_event",
        book_hists=book_fers_esum_vs_event,
        make_plots=plot_fers_esum_vs_event,
        enabled_by_default=True,
    ),
    # plot-only: channel layout map (no histogram step)
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
]

# ---------------------------------------------------------------------------
# DRS DQM sequences — require DRS data.
# ---------------------------------------------------------------------------

DRS_DQM_SEQUENCES: list[CaloXSequence] = [
    # plot-only: channel layout map (no histogram step)
    CaloXSequence(
        name="drs_mapping",
        book_hists=None,
        make_plots=plot_drs_mapping,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_waveforms",
        book_hists=book_drs_waveforms,
        make_plots=plot_drs_waveforms,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_profiles",
        book_hists=book_drs_profiles,
        make_plots=plot_drs_profiles,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_prof_combined",
        book_hists=book_drs_prof_combined,
        make_plots=None,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_stats",
        book_hists=book_drs_stats,
        make_plots=plot_drs_stats,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_finebins_combined",
        book_hists=book_drs_finebins_combined,
        make_plots=None,
        enabled_by_default=True,
    ),
    CaloXSequence(
        name="drs_peak_ts",
        book_hists=book_drs_peak_ts,
        make_plots=plot_drs_peak_ts,
        enabled_by_default=False,
    ),
]

# Combined list for backward compatibility.
ALL_SEQUENCES: list[CaloXSequence] = FERS_DQM_SEQUENCES + DRS_DQM_SEQUENCES

# ---------------------------------------------------------------------------
# Service DRS sequences  (ctx: CaloXAnalysisManager)
# ---------------------------------------------------------------------------


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
    #CaloXSequence(
    #    name="drs_waveforms",
    #    book_hists=book_drs_waveforms,
    #    make_plots=plot_drs_waveforms,
    #),
    CaloXSequence(
        name="drs_profiles",
        book_hists=book_drs_profiles,
        make_plots=plot_drs_profiles,
    ),
    CaloXSequence(
        name="drs_prof_combined",
        book_hists=book_drs_prof_combined,
        make_plots=None,
    ),
    CaloXSequence(
        name="drs_prof_corr_combined",
        book_hists=book_drs_prof_corr_combined,
        make_plots=None,
    ),
    CaloXSequence(
        name="drs_stats",
        book_hists=book_drs_stats,
        make_plots=plot_drs_stats,
    ),
    CaloXSequence(
        name="drs_finebins_combined",
        book_hists=book_drs_finebins_combined,
        make_plots=None,
    ),
    CaloXSequence(
        name="drs_finebins_corr_combined",
        book_hists=book_drs_finebins_corr_combined,
        make_plots=None,
    ),
    CaloXSequence(
        name="drs_cfd_mpv",
        book_hists=None,
        make_plots=plot_drs_cfd_mpv,
    ),
]

# ---------------------------------------------------------------------------
# DRS ↔ FERS correlation sequences  (require DRS integral variables defined)
# ---------------------------------------------------------------------------

DRS_FERS_SEQUENCES: list[CaloXSequence] = [
    CaloXSequence(
        name="drs_sum_vs_fers",
        book_hists=book_drs_sum_vs_fers,
        make_plots=plot_drs_sum_vs_fers,
    ),
    CaloXSequence(
        name="drs_peak_vs_fers",
        book_hists=book_drs_peak_vs_fers,
        make_plots=plot_drs_peak_vs_fers,
    ),
]

# ---------------------------------------------------------------------------
# FERS physics sequences  (require define_physics_variables on the manager)
# ---------------------------------------------------------------------------

FERS_SEQUENCES: list[CaloXSequence] = [
    CaloXSequence(
        name="fers_energy_sum",
        book_hists=book_fers_energy_sum,
        make_plots=plot_fers_energy_sum,
    ),
    CaloXSequence(
        name="fers_energy_sum_board",
        book_hists=book_fers_energy_sum_board,
        make_plots=plot_fers_energy_sum_board,
        enabled_by_default=False,  # opt-in via --sequences fers_energy_sum_board
    ),
    CaloXSequence(
        name="fers_lg_vs_mix",
        book_hists=book_fers_lg_vs_mix,
        make_plots=plot_fers_lg_vs_mix,
        enabled_by_default=False,  # opt-in via --sequences fers_lg_vs_mix
    ),
    #CaloXSequence(
    #    name="fers_stats",
    #    book_hists=book_fers_stats,
    #    make_plots=plot_fers_stats,
    #),
    CaloXSequence(
        name="fers_cer_vs_sci",
        book_hists=book_fers_cer_vs_sci,
        make_plots=plot_fers_cer_vs_sci,
    ),
    CaloXSequence(
        name="fers_dr",
        book_hists=book_fers_dr,
        make_plots=plot_fers_dr,
    ),
    CaloXSequence(
        name="fers_ewc",
        book_hists=book_fers_ewc,
        make_plots=plot_fers_ewc,
    ),
    CaloXSequence(
        name="fers_ewc_vs_hodo",
        book_hists=book_fers_ewc_vs_hodo,
        make_plots=plot_fers_ewc_vs_hodo,
    ),
    CaloXSequence(
        name="fers_shower_shape",
        book_hists=book_fers_shower_shape,
        make_plots=plot_fers_shower_shape,
    ),
]

# ---------------------------------------------------------------------------
# Lookup helpers — work for all lists
# ---------------------------------------------------------------------------

_ALL_BY_NAME: dict[str, CaloXSequence] = {
    s.name: s
    for s in (
        ALL_SEQUENCES
        + SERVICE_DRS_SEQUENCES
        + DRS_MCP_SEQUENCES
        + DRS_FERS_SEQUENCES
        + FERS_SEQUENCES
    )
}


def get_sequence(name: str) -> CaloXSequence:
    if name not in _ALL_BY_NAME:
        available = ", ".join(_ALL_BY_NAME)
        raise ValueError(f"Unknown sequence '{name}'. Available: {available}")
    return _ALL_BY_NAME[name]


def default_sequences() -> list[CaloXSequence]:
    """All DQM sequences (FERS + DRS) enabled by default."""
    return [s for s in ALL_SEQUENCES if s.enabled_by_default]


def default_fers_dqm_sequences() -> list[CaloXSequence]:
    """FERS-only DQM sequences enabled by default (no DRS required)."""
    return [s for s in FERS_DQM_SEQUENCES if s.enabled_by_default]


def default_drs_dqm_sequences() -> list[CaloXSequence]:
    """DRS DQM sequences enabled by default."""
    return [s for s in DRS_DQM_SEQUENCES if s.enabled_by_default]


def default_service_drs_sequences() -> list[CaloXSequence]:
    """Service DRS sequences enabled by default."""
    return [s for s in SERVICE_DRS_SEQUENCES if s.enabled_by_default]


def default_drs_mcp_sequences() -> list[CaloXSequence]:
    """DRS-only sequences run under MCP timing selection."""
    return list(DRS_MCP_SEQUENCES)


def default_fers_sequences() -> list[CaloXSequence]:
    """FERS physics sequences (require define_physics_variables)."""
    return [s for s in FERS_SEQUENCES if s.enabled_by_default]


def select_sequences(names: list[str],
                     pool: list[CaloXSequence] = None) -> list[CaloXSequence]:
    """Return sequences matching *names*, preserving pool order.

    pool defaults to ALL_SEQUENCES + SERVICE_DRS_SEQUENCES.
    """
    if pool is None:
        pool = ALL_SEQUENCES + SERVICE_DRS_SEQUENCES
    selected_names = {get_sequence(n).name for n in names}
    return [s for s in pool if s.name in selected_names]
