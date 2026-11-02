#!/bin/bash
# Dumps for the p.e. yield study (config.py): calibrated FERS per event (studies/zscan/dump_fers.py) for the positron runs, raw HG for
# the gain fit (pedestal run and the low-energy gain run). Runs already dumped are skipped.
# Usage (any directory): studies/pe_yield/dump.sh
STUDY=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$STUDY/../.." && pwd)
export PE_YIELD_DATA=${PE_YIELD_DATA:-/data1/yofeng/CaloX/pe_yield}
mkdir -p $PE_YIELD_DATA/logs; cd $REPO
RUNS=$(cd $STUDY && python3 -c "from config import RUNS, BAD_RUNS; print(*RUNS, *BAD_RUNS)")
for r in $RUNS; do
  [ -f $PE_YIELD_DATA/fers/fers_$r.npz ] && continue
  ZSCAN_DATA=$PE_YIELD_DATA python3 $REPO/studies/zscan/dump_fers.py $r 30000 > $PE_YIELD_DATA/logs/fers_$r.log 2>&1; echo "FERS $r rc=$?"
done
read PED NPED GAIN NBEAM <<< $(cd $STUDY && python3 -c "from config import *; print(PED_RUN, N_PED, GAIN_RUN, N_BEAM)")
for rn in "$PED $NPED" "$GAIN $NBEAM"; do
  set -- $rn
  [ -f $PE_YIELD_DATA/raw/raw_$1.npz ] && continue
  python3 $STUDY/dump_raw.py $1 $2 > $PE_YIELD_DATA/logs/raw_$1.log 2>&1; echo "raw $1 rc=$?"
done
echo ALLDONE
