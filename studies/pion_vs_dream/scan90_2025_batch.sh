#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
LOG=results/logs/scan90_2025
mkdir -p $LOG; rm -f $LOG/status.txt
for r in 1499 1506 1501 1507 1508 1512 1500 1502 1503 1504 1505 1509 1510 1511 1513; do
  if [ -f results/root/Run$r/drs_energy_stats.json ]; then echo "$r skip (exists)" >> $LOG/status.txt; continue; fi
  start=$(date +%s)
  python3 scripts/check_drs_mcp.py --run $r --particle electron > $LOG/run$r.log 2>&1
  echo "$r rc=$? $(( $(date +%s) - start ))s" >> $LOG/status.txt
done
echo done >> $LOG/status.txt
