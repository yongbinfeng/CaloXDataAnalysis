"""Run-number boundaries that separate the CaloX data-taking eras.

Every map in this package switches on these, so they are defined once here.
_for_run is the shared lookup over a table of (run_min, run_max) -> value with
half-open ranges; None on either side means open-ended.
"""


# Run number at which DRS branches gained the "Brg{N}_" prefix.
_DRS_BRG_RUN = 1700


# Run number at which DRS boards switched from 3mm to 6mm fibers.
_DRS_6MM_RUN = 1748


_DRS_FULL_RUN_TB2026 = 1828


_DRS_TESTFIBER_RUN_TB2026 = 1994


_DRS_PHASE2_RUN = 1896


def _for_run(table, run_number, default=None):
    """Look up the entry of a run-range table (see above)."""
    for (run_min, run_max), value in table:
        if (run_min is None or run_min <= run_number) and \
           (run_max is None or run_number < run_max):
            return value
    return default
