"""Timing-resolution study: paths, runs and constants shared by the scripts of this directory.
The z-scan code (studies/zscan: datasets, event selection, in-time hits) is reused for the dispersion part."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ZSCAN = os.path.join(HERE, "..", "zscan")
DATA = os.environ.get("TIMING_DATA", "/data1/yofeng/CaloX/timing")      # refs/ (dump_refs.py), *.json
PLOTS = "results/plots/Timing"
HTML = "results/html/Timing"
TS_PS = 200.0                         # DRS sample [ps]

# electronics runs (dump_refs.py): reference edges, MCP copies and channel CFD times, one or two runs per period
REF_RUNS = {"2025": [1507, 1511], "2026A": [1910, 1920], "2026B": [2005, 2030]}
REF_EVENTS = 20000
CELL_RUN = 2005                       # run of the DRS cell-pattern plot

# dispersion (dispersion.py): z-scan datasets
DISP_DATASETS = ("2025", "2026A", "2026B", "testfibres")
