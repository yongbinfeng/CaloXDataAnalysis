"""The z-scan datasets: 90-degree positron runs in which the table moves the shower along the fibres (table X; larger X = closer to
the readout). Run lists, groups and run sets are set here; everything else per run (channel map, particle ID, MCP channels and cuts)
comes from the framework configuration."""
import json, os

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
DATA = os.environ.get("ZSCAN_DATA", "/data1/yofeng/CaloX/zscan_dumps")     # dumps (drs/, fers/), results zscan_<dataset>.json, checks/
PLOTS = "results/plots/ZScan"                                                 # per dataset: results/plots/ZScan/<dataset>/
HTML = "results/html/ZScan"
RUNLIST = json.load(open(os.path.join(REPO, "data/Runlist.json")))

# Event selection on the MCPs: a pulse in MCP_DS_0 above its configured cut (configs/selection_config.py: 20 ADC in tb2025,
# 50 ADC from run 1829; in tb2026 this is the test-fibre mcp_clean), and a valid time reference (the DRS times are referenced to
# MCP_DS_1, channels.maps.services.get_mcp_reference). About 30% of the positrons cross the MCPs.
MCP_SELECT = "MCP_DS_0"

DATASETS = {
    "2025": dict(
        title="tb2025, table Y = 1064 mm",
        desc="tb2025, 40 GeV e<sup>+</sup>, table Y &asymp; 1064 mm, five runs at X = &minus;400 &hellip; &minus;55 mm. 1513 and 15130 are the two "
             "data files of one DAQ run. 1501 / 1507 are logged at 90&deg;, the others at 93.3&deg; (notes: \"90 degrees\").",
        runs=[15130, 1511, 1507, 1501, 1513],
        daq={15130: 1513},              # second data file of DAQ run 1513 (own runlist entry; dumped from the repository mirror)
        map_run=1507,
        groups=lambda r: "40 GeV",
        # 1501: the shower sits about one tower away from the other runs in the other transverse direction
        runsets={"main": [15130, 1511, 1507, 1513], "all": [15130, 1511, 1507, 1501, 1513]},
        fers=True,                      # calorimeter FERS recorded: whole-module sums and position matching
        offsets="global",               # run-to-run offsets common to all channels (five runs: too few for one per DRS chip)
        main_fit="plain",
        min_events=200,
    ),
    "2026A": dict(
        title="2026 scan A, table Y = &minus;2577 mm",
        desc="2026 scan A: runs 1897&ndash;1931, 40 GeV e<sup>+</sup>, table Y = &minus;2577 mm, X = &minus;296 &hellip; +36 mm. No calorimeter FERS "
             "(only the hodoscope boards were read out).",
        runs=list(range(1897, 1932)),
        map_run=1897,
        groups=lambda r: "40 GeV",
        fers=False,
        offsets="chip",                 # run-to-run offsets per DRS chip (the 8 channels of a chip vary together by ~2%)
        main_fit="offsets",
        min_events=1000,
    ),
    "2026B": dict(
        title="2026 scan B (test-fibre runs), table Y = &minus;233 mm",
        desc="2026 scan B: the test-fibre runs, 40 GeV (1999&ndash;2016) and 60 GeV (2026&ndash;2039) e<sup>+</sup>, table Y = &minus;233 mm, "
             "X = &minus;295 &hellip; +50 mm, with a lead block in front of the module. The DRS channels rewired to the test fibres and the "
             "MCP channels are left out. No calorimeter FERS.",
        runs=list(range(1998, 2017)) + [2026, 2027, 2029, 2030, 2031, 2032, 2034, 2035, 2036, 2037, 2038, 2039],
        map_run=1897,
        # one normalisation per beam energy: the runs after the brick move (2012-2016) and the 60 GeV revisit (2037-2039) line up
        groups=lambda r: "40 GeV" if r <= 2016 else "60 GeV",
        xmax=0.050,                     # beyond +50 mm the integrals get heavier tails (end of the copper)
        test_fibres=True,
        fers=False,
        offsets="chip",
        main_fit="offsets",
        min_events=1000,
    ),
    # the test fibres of the 2026B runs, each fibre on its own (testfibres.py): single small bundles behind the lead block, no sums
    "testfibres": dict(
        title="2026 test fibres",
        desc="2026 test fibres (channels/data/testingfibers.json) in the positron runs 1998&ndash;2016 (40 GeV) and 2026&ndash;2039 (60 GeV), "
             "table Y = &minus;233 mm, behind the lead block. From run 2012 on the lead block was moved by 48.2 cm (runlist): X along the "
             "test fibres = table X &minus; 482 mm, as in the framework scan. Left out: pion runs (2017&ndash;2025, 2028, 2033), muons "
             "(2041) and 2040, 2042 (beam mostly missing the table).",
        runs=list(range(1998, 2017)) + [2026, 2027, 2029, 2030, 2031, 2032, 2034, 2035, 2036, 2037, 2038, 2039],
        x_offset=lambda r: -0.482 if r >= 2012 else 0.0,
        # one normalisation per continuous stretch of the scan: the move of the lead block changes how the shower meets the fibres
        groups=lambda r: "40 GeV" if r < 2012 else ("40 GeV, block moved" if r <= 2016 else "60 GeV"),
        channel_map="channels/data/testingfibers.json",   # dumped and analysed: these channels only (lines with // are switched off)
        dump_dir="testfibres",
        offsets=None,
        min_events=1000,
    ),
}


def xpos(run):
    """table X [m] from the runlist"""
    return float(str(RUNLIST[str(run)]["Table position X"]).replace("mm", "").strip()) / 1000.0


def drs_dump(name, run):
    """per-event DRS dump of a run for a dataset (datasets with their own channel map have their own dump directory)"""
    return os.path.join(DATA, DATASETS[name].get("dump_dir", "drs"), f"drs_{run}.npz")
