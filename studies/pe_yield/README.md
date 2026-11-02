# p.e. per GeV in the 2025 FERS

Photoelectrons per GeV of beam energy for the scintillating, quartz and plastic fibres, 2025 positrons at position #1 (3.3°):
the Mix ADC of each channel divided by the SiPM gain (HG ADC per p.e.), summed over the channels of a fibre type.

## Results in short

- About **57 p.e./GeV scintillating, 4.2 quartz, 6.5 plastic** (≈ 11 Cherenkov), flat from 20 to 100 GeV; no leakage correction.
- **Gain** from the single-p.e. structure of a 10 GeV run, with the pedestal run as the zero-p.e. shape: 76 ADC/p.e. on boards 9 and 10
  (peaks resolved, quieter SiPMs), 53–61 on the other boards (peaks smeared by dark counts; per-board medians of shallow channel fits).
- ADC / gain counts fired SiPM cells: optical crosstalk is included (primary p.e. lower by the crosstalk probability).
- Runs 1409–1415 (energy scan of the night Sep 25–26) are left out: (S+C)/E ≈ 0.15 instead of ≈ 2, the FERS signal is not right there.

## Pipeline (repository root)

```bash
studies/pe_yield/dump.sh                  # FERS dumps (studies/zscan/dump_fers.py) and raw HG for the gain fit, skipped if present
python3 studies/pe_yield/gain.py          # -> $PE_YIELD_DATA/gain.json
python3 studies/pe_yield/yields.py        # -> $PE_YIELD_DATA/yields.json
python3 studies/pe_yield/page.py          # results/html/PEYield/PEYield_2025.html
```

Data: `$PE_YIELD_DATA` (default `/data1/yofeng/CaloX/pe_yield`): `fers/`, `raw/`, `gain.json`, `yields.json`, `logs/`.

| File | Content |
|---|---|
| `config.py` | runs, gain reference (pedestal run, gain run, boards with resolved peaks), selection |
| `dump.sh`, `dump_raw.py` | dumps |
| `gain.py` | template fit of the HG spectra, gain per channel and per board |
| `yields.py` | p.e. per GeV per fibre type and run |
| `page.py` | page |
