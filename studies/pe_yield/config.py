"""Photoelectrons per GeV in the 2025 FERS (tb2025, 6 mm SiPMs): runs, gain reference and paths."""
import os

DATA = os.environ.get("PE_YIELD_DATA", "/data1/yofeng/CaloX/pe_yield")   # fers/ (zscan dump_fers.py), raw/ (dump_raw.py), json
PLOTS = "results/plots/PEYield"
HTML = "results/html/PEYield"
RESPONSE = "data/fers/FERS_response_Sep.json"   # GeV per Mix ADC, inverted to get the Mix ADC back from the calibrated dump

# positron runs at position #1 (table X ~ 0, Y = 1032 mm, 3.3 deg): beam energy [GeV]
RUNS = {1424: 10, 1423: 20, 1422: 30, 1416: 30, 1514: 40, 1526: 60, 1527: 100}
# energy scan of the night Sep 25-26 at the same position: (S+C)/E ~ 0.3-0.5 with a tail to zero, the central channels at pedestal in
# most events (FERS sampling off the pulse peak?); shown, not used
BAD_RUNS = {1409: 80, 1410: 100, 1411: 120, 1412: 60, 1415: 40}
MAP_RUN = 1526                   # light map per channel

# gain (HG ADC per p.e.): template fit of the beam spectra of a low-energy run, with the pedestal run of the same evening as the
# zero-p.e. shape (pedestal + dark counts). Only boards 9 and 10 (about 5x lower SiPM current) show the single-p.e. peaks; on the
# other boards the dark counts, sampled at random times by the FERS hold, smear them out.
PED_RUN = 1425
GAIN_RUN = 1424
GAIN_BOARDS = (9, 10)
N_PED, N_BEAM = 100000, 30000    # events for the raw dumps

SEL_MIN_FRACTION = 0.7           # events: calibrated (S+C)/E above this (S and C are each calibrated to E, so S+C ~ 2E)
FAR = 6                          # [towers] from the shower centre: light outside this radius is listed (offset check)
