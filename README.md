# CaloXDataAnalysis

Data analysis framework for CaloX, built on **ROOT RDataFrame**.

### Prerequisites
Ensure you have the following installed on your system:
* **ROOT** (6.24+) with Python bindings (`PyROOT`)
* **Python** >= 3.8
* **C++ Compiler** (GCC or Clang)

### Installation
The repository is packaged for easy installation. Use the provided setup script to install the package in editable mode and compile the C++ utilities:

```bash
git clone https://github.com/yongbin-feng/CaloXDataAnalysis.git
cd CaloXDataAnalysis
chmod +x setup.sh
./setup.sh
```

### Data Setup
Download the data files from the CaloX Data Portal.

Configure your local data paths in data/datafiles.json.

Update configs/runconfig.py to set your target run number and event range.

### Usage
The framework provides several command-line entry points for standard analysis workflows:

Generate histograms and DQM plots for a specific run:

```bash
# Generate histograms (ROOT files)
calox-dqm-hists --run 1350

# Generate plots and HTML reports
calox-dqm-plots --run 1350
```

Energy & Timing Analysis

```bash
# FERS Energy analysis (pedestal subtraction, calibration, and sums)
calox-energy-plots --run 1350

# DRS Timing analysis (MCP-relative timing and Peak TS)
calox-timing-plots --run 1350
```

Auxiliary detector analysis
```bash
# Check Service DRS channels
calox-servicedrs

# Validate MCP signals
calox-checkmcp
```

### Test-fibre analysis
See [studies/testing_fiber/README.md](studies/testing_fiber/README.md) for the per-run and position-scan workflow of the SPS 2026 test fibres.

### DRS integrals
Default (`variables/drs.py`): constant baseline (median of samples 0-200), largest pulse in the **whole** record, 20% CFD crossing,
`<ch>_energy` (50 samples) and `<ch>_energy_short` (15 samples).

`prepare(drs_intime=True)` adds `<ch>_energy_it`, `<ch>_energy_short_it`, `<ch>_peak_value_it`, `<ch>_TS_cfd_mcp_it`: the same for
the largest pulse in raw samples 380-620 only. Recommendations (`studies/zscan`):
- Few-photoelectron channels: never use the all-event mean (dark counts and noise dominate). Hits: in time (6 TS, 20 TS
  scintillating) **and** peak > ~40 ADC (half a p.e.); subtract randoms from sidebands (inside 380-620 for the in-time columns,
  where every event has a noise "pulse"). Quote light per hit and hit fraction.
- With these cuts the default and in-time columns give the same results; the in-time ones are a convenience.
- Long integral (10 ns) for scintillating fibres: the short one (3 ns) holds a single-p.e. pulse but only ~40% of a
  scintillation pulse and is not proportional to the light. Short only for prompt (Cherenkov, single-p.e.) pulses.
- Baseline drifts ~1-2 ADC along the record: harmless for ratios and slopes; for absolute light subtract the average baseline shape.
- Saturation: drop whole channels with >1% of events above ~1500 ADC, never single events.

### Repository Structure
- `scripts/`: Main execution logic for plots and histogram generation.
- `utils/`: Core utilities including the data loader, fitter, and `html_generator`.
- `channels/`: Hardware mapping and geometry definitions for FERS and DRS boards. The tables they read (DRS->FERS CSVs, frozen per-era DRS maps, trigger delay/map, channel selections such as `testingfibers.json`) are in `channels/data/`.
- `variables/`: Definitions for high-level analysis variables (e.g., energy sums, timing).
- `configs/`: Global configurations for plot ranges and run parameters.
- `data/`: JSON maps for run lists, pedestals, and dead channels.