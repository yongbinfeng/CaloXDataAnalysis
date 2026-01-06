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

### Repository Structure
- `scripts/`: Main execution logic for plots and histogram generation.
- `utils/`: Core utilities including the data loader, fitter, and `html_generator`.
- `channels/`: Hardware mapping and geometry definitions for FERS and DRS boards.
- `variables/`: Definitions for high-level analysis variables (e.g., energy sums, timing).
- `configs/`: Global configurations for plot ranges and run parameters.
- `data/`: JSON maps for run lists, pedestals, and dead channels.