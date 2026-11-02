# Z scans: light speed and attenuation along the fibres

90° positron runs in which the table moves the shower along the fibres (table X; larger X = closer to the readout).
Datasets (`datasets.py`):

| Dataset | Runs | Readout | Notes |
|---|---|---|---|
| `2025` | 15130, 1511, 1507, 1501, 1513 (tb2025, Y = 1064 mm, 40 GeV) | DRS + calorimeter FERS | 15130 = second data file of 1513; 1501 outside the main run set (shower one tower away) |
| `2026A` | 1897–1931 (Y = −2577 mm, 40 GeV) | DRS | no calorimeter FERS recorded |
| `2026B` | test-fibre runs 1999–2016 (40 GeV), 2026–2039 (60 GeV), Y = −233 mm | DRS | lead block in front; test-fibre and MCP channels left out |
| `testfibres` | the same runs from 1998, the test fibres only (`channels/data/testingfibers.json`) | DRS | each fibre on its own (`testfibres.py`); own dumps; X along the fibres = table X − 482 mm from run 2012 (lead block moved) |

## Results in short

- **Light speed** is the robust measurement: the arrival time depends only on the distance to the readout. The speeds agree
  between the datasets (scintillating ≈ 0.53–0.55 c, quartz ≈ 0.62 c, plastic ≈ 0.60 c).
- **Attenuation** needs a sum over the whole module: the fibres are tilted relative to the table travel, so the shower moves across
  the channels during a scan, and a single channel or a partial sum follows the shower. Only 2025 has the calorimeter FERS;
  its whole-module sums with position matching (the hodoscope window shifted by the measured drift) are the attenuation measurement.
  DRS sums cover part of the module and are shown for comparison.
- **Timing resolution** (breakdown: MCP, DRS groups and boards, fibre dispersion): `studies/timing`.
- **Test fibres**: the same MCP requirement and in-time hits, each fibre on its own, without the hodoscope cuts (positrons with an MCP
  pulse; the z-scan events are kept as a comparison): speed, and attenuation from the light per crossing
  (light per hit × (1 − e<sup>−μ</sup>)) with one normalisation per stretch of the scan and the runs with partial overlap left out.

## Pipeline (repository root)

```bash
studies/zscan/dump.sh 2026A                 # per-event dumps (dump_drs.py, dump_fers.py), skipped if present
python3 studies/zscan/analysis.py 2026A     # -> $ZSCAN_DATA/zscan_2026A.json
python3 studies/zscan/gallery.py 2026A      # per-channel galleries
python3 studies/zscan/page.py 2026A         # results/html/ZScan/ZScan_2026A.html
python3 studies/zscan/overview.py           # results/html/ZScan/ZScan_overview.html (all datasets)
```

Test fibres: `dump.sh testfibres`, `analysis.py testfibres`, `page.py testfibres` (no gallery: the per-fibre plots are on the page).

Data: `$ZSCAN_DATA` (default `/data1/yofeng/CaloX/zscan_dumps`): `drs/`, `testfibres/` (test-fibre dumps), `fers/`, `zscan_<dataset>.json`,
`checks/`, `logs/`.
Plots: `results/plots/ZScan/<dataset>/`.

| File | Content |
|---|---|
| `datasets.py` | run lists, groups, run sets, MCP selection |
| `common.py` | dataset setup (channel map from the framework, MCP / test-fibre channels excluded), event selection, in-time hits, fits, drift fit, position matching, FERS loader |
| `analysis.py` | per-channel DRS (light per hit / per event, arrival time), sums per category with variants, FERS, drift, position check |
| `testfibres.py` | test fibres: per-fibre speed and attenuation, page (called by `analysis.py` / `page.py testfibres`) |
| `gallery.py`, `page.py`, `overview.py` | per-channel galleries, one page per dataset, overview |
| `dump_drs.py`, `dump_fers.py`, `dump.sh` | per-event dumps |
| `checks/` | finished one-off studies, kept runnable: integral variants, pulse containment, baseline drift and shape, MCP check, integral shapes |

## Method

- **Events**: positrons (beam-PID flags of the period), 1–3 hodoscope fibres per axis, beam within 1 cm of the run median (2025:
  vertically of a common centre after the drift shift). DRS: a pulse in MCP_DS_0 above its configured cut (in 2026 the test-fibre
  `mcp_clean`) and a valid time reference (MCP_DS_1).
- **In-time hits**: CFD time within ±4 ns (scintillating) / ±1.2 ns (Cherenkov) of the arrival peak and peak above 40 ADC (3 mm SiPM)
  or 25 ADC (6 mm); randoms from sidebands subtracted. Long integral (10 ns) for scintillating, short (3 ns) for Cherenkov.
- **Light per event** (events without a hit count as zero) is proportional to the light; light per hit leaves those events out and
  gives too flat a slope for few-p.e. channels.
- **Saturation**: DRS channels with more than 1% of events above 1500 ADC peak are dropped.
- **Fits**: ln(signal) = a<sub>group</sub> + X/λ, one normalisation per beam-energy group, with run offsets (per DRS chip in 2026,
  common to all channels in 2025) kept free of a constant and a linear trend in X. Time: straight line, speed = 1/|slope|.
