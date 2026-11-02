#!/bin/bash
# Runs 15120 / 15130 are the second data files of DAQ runs 1512 / 1513 (own runlist entries). The framework picks the period
# setup from the run number, so they are analysed as 1512 / 1513 from a mirror of the repository whose data-file map points
# these numbers to the second files (make_mirror.sh); outputs are renamed to Run15120 / Run15130.
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
[ -d $DATA/mirror/core ] || $STUDY/make_mirror.sh
cd $DATA/mirror
for pair in "1513 15130" "1512 15120"; do
  set -- $pair
  start=$(date +%s)
  python3 scripts/check_drs_mcp.py --run $1 --particle electron > $DATA/drs_$2.log 2>&1
  rc=$?; rm -rf results/root/Run$2; mv results/root/Run$1 results/root/Run$2
  echo "DRS $2 rc=$rc $(( $(date +%s) - start ))s"
  python3 $STUDY/scan90ch_dump.py $1 20000 $2 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-120
done
echo ALLDONE
