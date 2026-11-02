# Timing resolution: breakdown

What makes up the DRS timing resolution (CFD time of a channel minus that of the reference MCP_DS_1): the MCP, the alignment of DRS
groups and boards, the dispersion in the fibres, and the rest. All terms add in quadrature.

## Results in short

- **MCP**: two MCPs on the same DRS chip differ by about 43 ps in 2025 (about 30 ps per MCP, used for all periods: the 2026 upstream MCP
  gives about half the amplitude of the reference, ~80 ADC, and dominates the pair).
- **DRS alignment**: channel times are aligned through the reference signal (Channel8) of each group. Its time used to be the integer sample
  (0.2 ns/√12 = 58 ps of quantization per group); it is now the interpolated 50 % crossing (`process_dynamic_led` in `utils/functions.cc`).
  The two copies of the same MCP signal on two boards measure what the alignment adds to real pulses: about 110 ps before, 65 ps now.
  The rest is mostly the non-uniform DRS4 cell widths (a per-cell time calibration from the StartIndexCell could remove it).
- **Fibre dispersion**: at fixed number of photoelectrons the width grows with the distance to the readout, N × dσ²/dX the same in every
  p.e. bin; strong for the high-NA plastic and Quartz700 test fibres, absent for Quartz and Quartz400, small for the calorimeter channels.

## Pipeline (repository root)

```bash
python3 studies/timing/dump_refs.py 2005            # reference / MCP / channel times of one run -> $TIMING_DATA/refs/ (runs: config.REF_RUNS)
python3 studies/timing/electronics.py               # -> $TIMING_DATA/electronics.json
python3 studies/timing/dispersion.py 2026A          # per z-scan dataset (studies/zscan dumps) -> $TIMING_DATA/dispersion_2026A.json
python3 studies/timing/page.py                      # results/html/Timing/Timing_resolution.html
```

Data: `$TIMING_DATA` (default `/data1/yofeng/CaloX/timing`). Plots: `results/plots/Timing/`.

| File | Content |
|---|---|
| `config.py` | paths, runs |
| `dump_refs.py` | per-event reference times (interpolated, integer, line fit), StartIndexCell, MCP copies, channel CFD times |
| `electronics.py` | MCP pairs and copies, group alignment, DRS cell pattern, channel widths with integer / interpolated references |
| `dispersion.py` | width vs position in p.e. bins for the z-scan datasets: per-photon slope, 1/N test, MCP control |
| `page.py` | breakdown and page |
