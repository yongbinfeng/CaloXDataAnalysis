#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for r in 1409 1411 1423 1416 1514 1389 1390 1398 1399 1404 1405 1406 1408; do
  [ -f $DATA/hodo_$r.npz ] && continue
  python3 $STUDY/hodo_dump.py $r electron 25000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-80
done
echo ALLDONE
