#!/bin/bash
STUDY=$(cd "$(dirname "$0")" && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
# FERS per-channel dump for run 1500 (same position as 1501), missing from scan90ch_batch.sh
cd "$STUDY/../.."
for r in 1500; do      # 1503 has no local data file
  python3 $STUDY/scan90ch_dump.py $r 20000 2>&1 | grep -E "^done|Error|Traceback" | head -1 | cut -c1-120
done
echo ALLDONE
