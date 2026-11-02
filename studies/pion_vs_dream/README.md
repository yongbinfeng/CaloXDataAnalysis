# Pion response and resolution vs DREAM (tb2025)

Compares the CaloX pion data (SPS 2025, FERS readout) with the original DREAM paper on the same module
(N. Akchurin et al., NIM A 537 (2005) 537). Result page: `results/html/FERS_PionResponse/PionResponse_vs_DREAM.html`.

Standalone study scripts for now; meant to be folded into the analysis sequences later (see the end).

## Method in short

Per pion event:

1. **Depth**: displacement = shower centre (energy-weighted S) minus impact point (TTU fibre hodoscope, mapped to the detector
   frame with 100 GeV positrons, run 1527), relative to the positron runs taken just before the pion scan (1409, 1411, 1416, 1423).
   Depth z = displacement along the tilt direction (dY/dX = 0.71, measured) / tan(2.5 deg).
2. **Attenuation**: S, quartz C and plastic C scaled by exp(-z/lambda), with lambda from the 90-degree positron runs at Y ~ 1064 mm
   (1513, 1507, 1511, 15130; 1501 excluded): 3.07 m (scintillating), 6.7 m (quartz), 6.5 m (plastic). Impact-position terms fitted per run.
3. **Depth cut**: z < 70 cm (DREAM: light centre of gravity in the first metre).
4. **Electron scale**: divided by the positron response at the same energy (or kept on the calibrated energy scale).
5. **Q/S**: DREAM's method (leakage-corrected S = S + outer ring, f_em from C/S with e/h = 1.3 / 5, then S / (f + (1-f)/1.3)).

## Layout

| Where | What |
|---|---|
| `paths.py` | `STUDY_DIR`, `REPO_DIR`, `DATA_DIR` (default `/data1/yofeng/CaloX/pion_vs_dream_dumps`, override with `PION_DREAM_DATA`) |
| `*_dump.py`, `*_batch.sh` | event dumps (run from the repo through the framework), written to `DATA_DIR` as `.npz` |
| `edepth_test.py` -> `pion_final.py` -> `pion_final_res.py` -> `final_analysis.py` | the analysis chain (scripts load each other with `exec`) |
| `final_plots_only.py`, `make_final_page.py` | plots from the saved results, and the HTML page |
| `att_plot.py`, `scan90_rows.py` | attenuation from the 90-degree positron runs |
| `offset_plots.py` | uncorrected S vs displacement |
| `res_fit.py`, `res_fit_plot.py` | resolution fits a/sqrt(E) + b, with and without b fixed to DREAM's |
| `leak_opt.py`, `depth_cut.py`, `angle_scan.py`, `test12.py`, `test2w.py` | side studies: leakage correction, depth cut, tilt, electron scale (tests 1, 2) |
| `dream_fig11.json`, `dream_fig32.json` | DREAM response (Fig. 11, 35) and resolution (Fig. 32a) values read off the figures |
| `exploratory/` | earlier intermediate scripts, kept for reference; not needed for the results |

`DATA_DIR` also holds the intermediate results (`final_analysis.json`, `res_fit.json`, `att90_Y1064.json`, ...), the DRS outputs of
runs 15120 / 15130 (`mirror/results/`), and the DREAM paper (`dream_537.pdf`, `dream.txt`).

## Running

From the repository root (`setup.sh` sourced):

```bash
S=studies/pion_vs_dream
# dumps (one run at a time: /data2 is a single disk)
$S/escan_batch.sh          # positrons: central energy scan and off-centre 80 GeV runs (hodo_*.npz)
$S/hodo_batch.sh           # 60 GeV positrons and pions (hodo_*.npz); also run hodo_dump.py 1527 electron 25000, 1432 pion 25000
$S/hodoq_batch.sh          # pions with quartz / plastic C (hodoq_*.npz)
$S/hodoc_batch.sh          # pions with every scintillation channel (hodoc_*.npz, leakage study)
$S/scan90ch_batch.sh       # 90-degree positron runs, per-channel FERS (scan90ch_*.npz)
$S/scan90ch_batch2.sh      # run 1500
$S/scan90_2025_batch.sh    # DRS analysis of the 90-degree runs (results/root/Run15xx)
$S/split_runs_batch.sh     # runs 15120 / 15130 (second files of 1512 / 1513), DRS + FERS, via make_mirror.sh
$S/chan_batch.sh           # per-channel dumps 1432 (pion), 1409 (positron), for test 2
# analysis, plots, page
python3 $S/att_plot.py && python3 $S/final_analysis.py && python3 $S/res_fit_plot.py && python3 $S/offset_plots.py
python3 $S/make_final_page.py
```

## Notes

- **Runs 15120 / 15130** are the second data files of DAQ runs 1512 / 1513, with their own runlist entries (1512: Y = 1000 mm,
  15120: Y = 936 mm; 1513: X = -54.5 mm, 15130: X = -400.3 mm). Period thresholds in the code compare the raw run number
  (> 1700 = tb2026), so they are analysed as 1512 / 1513 from a mirror of the repository whose data-file map points to the second
  files (`make_mirror.sh`).
- **Hodoscope reference**: the hodoscope zero moved by 0.56 cm after the 13.3-93.3 degree table rotations (runs 1481-1513); the pions
  use the positron runs before them.
- **Tilt**: the logged 3.3 deg is not precise; 2.5 deg (total, along dY/dX = 0.71) gives the best Q/S resolution; 2-4 deg changes the
  response by up to +-4%.

## Towards the sequences

Possible split when merging:
- **Variables** (event level, in `variables/`): hodoscope centroid (`hodo_centroid`), quartz / plastic C sums, shower displacement and
  depth, attenuation-corrected S / C, Q/S.
- **Selection**: depth cut, hodoscope 1-3 fibres per axis.
- **Configuration**: tilt, lambda per fibre type, hodoscope map and positron-depth reference per period.
- **Post-processing**: response / resolution vs energy, fits, comparison with DREAM (from the per-run histograms instead of event dumps).
