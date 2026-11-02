#!/bin/bash
# Per-event dumps of a z-scan dataset (datasets.py), one run at a time: DRS (dump_drs.py) and, for datasets with calorimeter FERS,
# FERS (dump_fers.py). Runs already dumped are skipped. A second data file with its own runlist entry (datasets.py "daq", e.g. 15130 =
# second file of 1513) is dumped as its DAQ run from the repository mirror whose data-file map points to that file
# (studies/pion_vs_dream/make_mirror.sh). A dataset with its own channel map (testfibres) is dumped to its own directory.
# Usage (any directory): studies/zscan/dump.sh DATASET
STUDY=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$STUDY/../.." && pwd)
export ZSCAN_DATA=${ZSCAN_DATA:-/data1/yofeng/CaloX/zscan_dumps}
MIRROR=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}/mirror
mkdir -p $ZSCAN_DATA/drs $ZSCAN_DATA/fers $ZSCAN_DATA/logs
# "label daq fers dump" per run
LIST=$(cd $STUDY && python3 -c "
from datasets import DATASETS, drs_dump; c = DATASETS['$1']
for r in c['runs']: print(r, c.get('daq', {}).get(r, r), int(c.get('fers', False)), drs_dump('$1', r))")
echo "$LIST" | while read r daq fers out; do
  dir=$REPO; [ "$r" != "$daq" ] && { [ -d $MIRROR/core ] || $REPO/studies/pion_vs_dream/make_mirror.sh; dir=$MIRROR; }
  cd $dir
  if [ ! -f $out ]; then
    start=$(date +%s); python3 $STUDY/dump_drs.py $daq 40000 $r $1 > $ZSCAN_DATA/logs/drs_$1_$r.log 2>&1
    echo "DRS $r rc=$? $(( $(date +%s) - start ))s"
  fi
  if [ "$fers" = 1 ] && [ ! -f $ZSCAN_DATA/fers/fers_$r.npz ]; then
    start=$(date +%s); python3 $STUDY/dump_fers.py $daq 40000 $r > $ZSCAN_DATA/logs/fers_$r.log 2>&1
    echo "FERS $r rc=$? $(( $(date +%s) - start ))s"
  fi
done
echo ALLDONE
