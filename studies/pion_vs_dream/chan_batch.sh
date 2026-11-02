#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
cd "$STUDY/../.."
for spec in "1432 pion" "1409 electron"; do
  set -- $spec
  python3 $STUDY/chan_dump.py $1 $2 25000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-80
done
echo ALLDONE
