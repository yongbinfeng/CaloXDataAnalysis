#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for r in 1441 1439 1438 1437 1429 1432 1433 1434 1469 1471 1465 1470; do
  python3 $STUDY/resp_dump.py $r 25000 2>&1 | grep -E "^done|Error|Traceback" -A3
done
echo ALLDONE
