#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
LOG=results/logs/scan90
rm -f $LOG/status.txt
for r in $(seq 1897 1930); do
  if [ -f results/root/Run$r/drs_energy_stats.json ]; then echo "$r skip (exists)" >> $LOG/status.txt; continue; fi
  start=$(date +%s)
  python3 scripts/check_drs_mcp.py --run $r > $LOG/run$r.log 2>&1
  echo "$r rc=$? $(( $(date +%s) - start ))s" >> $LOG/status.txt
done
echo done >> $LOG/status.txt
