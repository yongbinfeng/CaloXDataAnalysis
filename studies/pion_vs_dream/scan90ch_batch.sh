#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for r in 1499 1501 1502 1504 1505 1506 1507 1508 1509 1510 1511 1512 1513; do
  python3 $STUDY/scan90ch_dump.py $r 20000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-120
done
echo ALLDONE
