# Test-fibre analysis (SPS 2026)

How to run the analysis of the 20 candidate fibre bundles read out by the DRS digitizer,
from raw runs to the position-scan summaries used in the paper draft.

The analysis has two steps:

1. **Per run:** histograms, timing and energy, plus the in-time / out-of-time analysis,
   with `scripts/check_drs_mcp.py`.
2. **Across runs:** the position scan (speed, light yield, attenuation, pulse shape,
   timing and position resolution) with `scripts/compare/plot_drs_mpv_scan_testfibers.py`.

Run every command from the repository root.

## Setup

- ROOT with PyROOT and a C++ compiler. Minuit is **not** required: all fits are hand-written
  least-squares fits.
- `utils/functions.cc` (CFD, integral, pulse edges) is compiled on first use with ACLiC.
  After editing it, the next run recompiles it automatically.
- Data file locations come from the JSON set by `jsonFile` in `configs/run_config.py`
  (currently `data/datafiles_local.json`), or pass `--json-file`.
- Beam energy and table position of each run come from `data/Runlist.json`.
- The fibre label to DRS channel map is `channels/data/testingfibers.json`.
  The channel labelled `Quartz_1` actually carries a split copy of the upstream MCP,
  and the two `GradedIndex` bundles have no signal.

## Step 1: per-run analysis

```bash
python3 scripts/check_drs_mcp.py --run 2030 \
    --channels channels/data/testingfibers.json --mcp-clean --intime
```

| Option | What it does |
|---|---|
| `--channels` | Only the listed test-fibre channels are processed (MCP references are always kept). Required by `--intime`. |
| `--mcp-clean` | Keeps only events where all MCPs fire cleanly. Always use it: the timing is relative to the MCP. |
| `--intime` | After the standard DRS histograms, runs a second event loop for the in-time vs out-of-time analysis (`analysis/intime_analysis.py`). |
| `--first-event`, `--last-event` | Restrict the event range, useful for quick tests. |

The in-time analysis needs the timing MPVs from the first pass of the same command, so
keep both in one call.

### Many runs

```bash
scripts/batch_check_drs_mcp.sh 1994-2042            # ranges and single runs, any mix
EXTRA_ARGS="--mcp-clean --intime" scripts/batch_check_drs_mcp.sh 1998-2039
DRY_RUN=1 scripts/batch_check_drs_mcp.sh 2030       # print the commands only
```

- Logs go to `results/logs/check_drs_mcp/run<N>.log`, with one line per run in `status.txt`.
- Runs are processed one at a time on purpose. The raw files sit on one spinning disk,
  and parallel jobs are limited by disk reads, not CPU.
- Remember `--intime` in `EXTRA_ARGS`: the default is `--mcp-clean` only, and the `pe` and
  `shape` scans below need the in-time output.

### Per-run outputs

| Path | Content |
|---|---|
| `results/root/Run<N>/drs_stats.root` | Histograms, including the fine-binned CFD time (0.1 TS bins, 400–600 TS) |
| `results/root/Run<N>/drs_energy_stats.json` | Integral means per channel |
| `results/root/Run<N>/drs_intime.root` | In-time analysis histograms |
| `results/plots/Run<N>/DRS_InTime/intime_summary.json` | Per-fibre summary read by the scan: hit fraction, light per hit, p.e., pulse shape, timing widths, time walk |
| `results/html/Run<N>/DRS/DRS_Time.html`, `DRS_Sum.html` | Timing and energy-sum pages |
| `results/html/Run<N>/DRS/DRS_InTime.html` | Per-fibre in-time page: amplitude, energy, profiles, edges, timing, σ and walk vs amplitude |
| `results/html/Run<N>/DRS/DRS_InTime_*_overlay.html` | Interactive overlays of in-time / out-of-time / all, per fibre |

### Main definitions (in `analysis/intime_analysis.py`)

- Time: CFD at 20 % of the peak, relative to the MCP reference, `TS_cfd_mcp`. One TS is 0.2 ns.
- In time: within 6 TS of the timing MPV for Cherenkov fibres, 20 TS for scintillating fibres
  (`IN_TIME_WINDOW_TS`, `IN_TIME_WINDOW_SCINT_TS`).
- Random coincidences: estimated from the sideband 100–300 TS away from the MPV
  (`SIDEBAND_TS`) and subtracted from the hit fraction and the light per hit.
- Photoelectrons: the single-p.e. integral is calibrated per run from out-of-time dark counts.
- Timing width: half the 16–84 % interquantile range, before and after the time-walk correction.

## Step 2: position scan over runs

```bash
python3 scripts/compare/plot_drs_mpv_scan_testfibers.py --mode all \
    --channels channels/data/testingfibers.json
```

| `--mode` | Needs | Produces |
|---|---|---|
| `time` | `drs_stats.root` | CFD MPV vs X per fibre, propagation speed and intercept |
| `energy` | `drs_energy_stats.json` | Mean integral vs X |
| `pe` | `intime_summary.json` | Light per hit, p.e. per hit, hit fraction vs X; exponential attenuation fits |
| `shape` | `intime_summary.json`, and the `time` output for the speeds | Rise, fall, FWHM, timing width vs X; dispersion in ps/m; position resolution σx = v σt |
| `both` | | `time` + `energy` (default) |
| `all` | | All four |

Run `time` before `shape`, or use `all`, which runs them in that order.
Runs whose outputs are missing are skipped with a message.

### Outputs

| Mode | Plots and tables | HTML |
|---|---|---|
| `time` | `results/plots/PositionScan/DRS_CFD_MPV_Scan/` | `results/html/PositionScan/DRS_CFD_MPV_Scan_channels.html` |
| `energy` | `results/plots/PositionScan/DRS_SumMean_Scan/` | `.../DRS_SumMean_Scan_channels.html` |
| `pe` | `results/plots/PositionScan/DRS_PE_Scan/` | `.../DRS_PE_Scan_channels.html` |
| `shape` | `results/plots/PositionScan/DRS_Shape_Scan/` | `.../DRS_Shape_Scan_channels.html` |

Each directory also holds the summary tables as `.json` and `.md`, for example
`drs_cfd_mpv_scan_channels.json` (speeds), `drs_pe_scan_channels.json` (plateaus and
attenuation lengths), `drs_shape_scan_channels.json` (shape, timing, dispersion) and
`drs_position_resolution_channels.json`.

### Run selection and position (top of `plot_drs_mpv_scan_testfibers.py`)

- `SCAN_RUNS` is 1998–2042 minus `SKIP_RUNS`. Runs 2028 and 2033 are pion runs, and
  2040–2042 have the beam missing the table. Edit `SKIP_RUNS` by hand to change this.
- X is the table position from `data/Runlist.json`, shifted by −482 mm for runs from 2012 on,
  when the holder was moved.
- 40 GeV runs (1998–2016) and 60 GeV runs (2026–2039) are drawn with different markers and
  get separate normalisations in the fits, with a common slope or attenuation length.
- Plateau values average the runs with −220 ≤ X ≤ 50 mm (`PE_FIT_XMIN`, `PE_FIT_XMAX`).
  The edge runs at −375 and −395 mm are left out of the fits.

## Single-p.e. pulse linearity (any DRS channels)

```bash
python3 scripts/check_drs_pe_linearity.py --run 2030 \
    --channels channels/data/testingfibers.json --mcp-clean
python3 scripts/check_drs_pe_linearity.py --combine 2026 2027 2029-2032 2034-2039 \
    --channels channels/data/testingfibers.json --classes charge   # from the per-run caches
```

Cache many runs first, one at a time (the disk is the bottleneck), with `--dump-only`.
Options: `--sample in|out|both` (default `in`: in time with the beam), `--nmax` (default 8),
`--classes resolved|charge`. `resolved` (default) requires amplitude and integral to match n
for the peaks resolved in every run, which keeps synchronous photoelectrons, and uses charge bins
above. `charge` uses charge bins, n ± 0.5 p.e. of integral, for every n, the same definition at all n.
Combined pages go to `results/html/PE_Linearity/Runs<a>-<b>_<sample>_<classes>/`.

- Two integrals are available per DRS channel. `<ch>_energy` sums 50 samples, from 5 before to 45 after
  the first sample above the 20 % CFD threshold, and measures the total charge, including scintillation tails.
  `<ch>_energy_short` sums a fixed 15 samples, from 4 before to 11 after that sample, which covers the positive
  lobe of the single-p.e. pulse. It has the same noise for every pulse and resolves the p.e. peaks far better.
  Use it to count photoelectrons in the Cherenkov fibres.
- Pulses are sorted into p.e. classes using the amplitude and the 15-sample integral.
  The averaged waveforms per class are compared with n times the 1 p.e. pulse, aligned at the
  20 % CFD time and at the peak. Out-of-time pulses (dark counts and crosstalk) and in-time pulses are kept apart.
- A linear-pulse test adds n times the 1 p.e. shape to real noise-only waveforms and runs the same chain,
  to separate method bias from real non-linearity.
- The first call reads the raw data once (about a minute) and caches the pulses in
  `results/root/Run<N>/drs_pe_pulses.npz`. Later calls reuse it; `--redo` re-reads.
- Page: `results/html/Run<N>/DRS/DRS_PE_Linearity.html` (combined runs: `results/html/PE_Linearity/Runs<a>-<b>/`).
  Tables: `pe_linearity_summary.{json,md}` in the plot directory. `--channels` accepts any label-to-channel map.

## Typical full refresh

```bash
EXTRA_ARGS="--mcp-clean --intime" scripts/batch_check_drs_mcp.sh 1998-2042
python3 scripts/compare/plot_drs_mpv_scan_testfibers.py --mode all \
    --channels channels/data/testingfibers.json
```

Then open `results/html/PositionScan/*.html`. Figures for the paper draft in
`/data1/yofeng/CaloX/Tex/testing_fibers_paper` are copied from these plot directories.
