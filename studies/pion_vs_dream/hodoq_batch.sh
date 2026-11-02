#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for r in 1432 1441 1439 1438 1437 1429 1433 1434; do
  python3 $STUDY/hodoq_dump.py $r pion 50000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-80
done
echo ALLDONE
