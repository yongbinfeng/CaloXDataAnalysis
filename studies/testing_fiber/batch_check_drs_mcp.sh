#!/bin/bash
# Run scripts/check_drs_mcp.py over a set of runs, one job at a time.
#
#   studies/testing_fiber/batch_check_drs_mcp.sh 1994-2016            # inclusive range
#   studies/testing_fiber/batch_check_drs_mcp.sh 1994 1996 2001-2003  # single runs and ranges, any mix
#
# Knobs (environment variables):
#   CHANNELS    channel-map JSON passed as --channels   (default channels/data/testingfibers.json)
#   EXTRA_ARGS  further check_drs_mcp.py options        (default "--mcp-clean")
#   LOGDIR      one log per run goes here               (default results/logs/check_drs_mcp)
#   DRY_RUN=1   print the commands instead of running them
#
# Each run's stdout/stderr goes to $LOGDIR/run<N>.log; a one-line summary per
# run (exit code, wall time) is appended to $LOGDIR/status.txt. A failing run
# does not stop the others; the exit code is the number of failed runs.

set -u
cd "$(dirname "$0")/../.." || exit 1

SCRIPT="scripts/check_drs_mcp.py"
CHANNELS="${CHANNELS:-channels/data/testingfibers.json}"
EXTRA_ARGS="${EXTRA_ARGS:---mcp-clean}"
LOGDIR="${LOGDIR:-results/logs/check_drs_mcp}"
DRY_RUN="${DRY_RUN:-0}"

if [ $# -eq 0 ]; then
    sed -n '2,15p' "$0"
    exit 1
fi

# expand "a-b" ranges and bare run numbers into one list
RUNS=()
for spec in "$@"; do
    case "$spec" in
        *-*) lo=${spec%-*}; hi=${spec#*-}
             if ! [[ "$lo" =~ ^[0-9]+$ && "$hi" =~ ^[0-9]+$ && "$lo" -le "$hi" ]]; then
                 echo "bad range: $spec" >&2; exit 1
             fi
             for ((r = lo; r <= hi; r++)); do RUNS+=("$r"); done ;;
        *)   if ! [[ "$spec" =~ ^[0-9]+$ ]]; then echo "bad run number: $spec" >&2; exit 1; fi
             RUNS+=("$spec") ;;
    esac
done

mkdir -p "$LOGDIR"
echo "${#RUNS[@]} runs: ${RUNS[*]}"
echo "command: python3 $SCRIPT --run <N> --channels $CHANNELS $EXTRA_ARGS"
echo "logs:    $LOGDIR/run<N>.log, summary in $LOGDIR/status.txt"

nfail=0
for run in "${RUNS[@]}"; do
    cmd=(python3 "$SCRIPT" --run "$run" --channels "$CHANNELS" $EXTRA_ARGS)
    if [ "$DRY_RUN" != 0 ]; then
        echo "${cmd[*]}"
        continue
    fi
    echo "=== Run $run  ($(date '+%F %T')) ==="
    start=$(date +%s)
    "${cmd[@]}" > "$LOGDIR/run$run.log" 2>&1
    rc=$?
    line="$run rc=$rc $(( $(date +%s) - start ))s $(date '+%F %T')"
    echo "$line" >> "$LOGDIR/status.txt"
    if [ $rc -ne 0 ]; then
        nfail=$((nfail + 1))
        echo "    FAILED (rc=$rc), see $LOGDIR/run$run.log; last lines:"
        grep -v "Warning in\|Info in" "$LOGDIR/run$run.log" | tail -5 | sed 's/^/    /'
    else
        echo "    done in $(( $(date +%s) - start ))s"
    fi
done

[ "$DRY_RUN" != 0 ] && exit 0
echo "finished: $(( ${#RUNS[@]} - nfail )) ok, $nfail failed"
exit $nfail
