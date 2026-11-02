#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for spec in "1526 electron" "1441 pion" "1439 pion" "1438 pion" "1437 pion" "1429 pion" "1433 pion" "1434 pion"; do
  set -- $spec
  python3 $STUDY/hodo_dump.py $1 $2 25000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-80
done
echo ALLDONE
