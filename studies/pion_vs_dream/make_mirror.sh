#!/bin/bash
# Mirror of the repository (symlinks) with its own results/ and a data-file map in which 1512 / 1513 point to the second files
# of those DAQ runs (15120 / 15130), so they get the tb2025 setup. Used by split_runs_batch.sh.
STUDY=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$STUDY/../.." && pwd)
DATA=${PION_DREAM_DATA:-/data1/yofeng/CaloX/pion_vs_dream_dumps}
M=$DATA/mirror
mkdir -p $M/results $M/data
for e in $REPO/*; do b=$(basename $e); [ "$b" = results ] || [ "$b" = data ] || ln -sfn $e $M/$b; done
for e in $REPO/data/*; do b=$(basename $e); [ "$b" = datafiles_local.json ] || ln -sfn $e $M/data/$b; done
python3 - "$REPO/data/datafiles_local.json" "$M/data/datafiles_local.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1])); d["1512"], d["1513"] = d["15120"], d["15130"]
json.dump(d, open(sys.argv[2], "w"), indent=4)
PY
echo "mirror ready: $M"
