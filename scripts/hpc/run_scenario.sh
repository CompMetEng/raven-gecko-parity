#!/usr/bin/env bash
# Run one scenario on both implementations and compare.
#   run_scenario.sh <scenario-id>
# Writes $PARITY_OUT/{python,matlab}/<id>.json, logs/<id>.log and status/<id>.txt, whose
# content is RUNNING while the job lives, then MATCH, DIFFER or ERROR. SEEDS lists the
# PYTHONHASHSEED values for the Python side (default "0 1 2"); a result that changes with
# the seed is a finding by itself. MATLAB_TIMEOUT (seconds, default 5400) limits the MATLAB
# step; genome-scale scenarios need a larger value. Exit code: 0 for MATCH and DIFFER (both are answers),
# 1 for ERROR, 2 for a usage error.
set -uo pipefail
source "$(dirname "$0")/env.sh" || exit 1

id=${1:-}
case $id in ""|*[!A-Za-z0-9_]*) echo "usage: run_scenario.sh <scenario-id>" >&2; exit 2 ;; esac
[ -f "$PARITY_ROOT/scenarios/$id/scenario.yml" ] || { echo "no scenario $id" >&2; exit 2; }
set -f; seeds=( ${SEEDS:-0 1 2} ); set +f      # splits on any whitespace, newlines included
[ "${#seeds[@]}" -ge 1 ] || { echo "SEEDS is empty" >&2; exit 2; }
for s in "${seeds[@]}"; do [[ $s =~ ^[0-9]{1,9}$ ]] || { echo "SEEDS must be integers below 10^9: $s" >&2; exit 2; }; done

mkdir -p "$PARITY_OUT"/{python,matlab,logs,status} || { echo "$id ERROR (cannot create $PARITY_OUT)"; exit 1; }
status_file=$PARITY_OUT/status/$id.txt
log=$PARITY_OUT/logs/$id.log
sdir=$PARITY_ROOT/scenarios/$id

# The context and result files live in the scenario folder of the shared checkout, so two
# runs of the same scenario would overwrite each other's files and could return the other
# run's result. One run per scenario at a time, across all batches.
exec 9>"$sdir/.run.lock" || { echo "$id ERROR (cannot open $sdir/.run.lock)"; echo ERROR >"$status_file"; exit 1; }
flock -n -E 200 9; lockrc=$?
if [ "$lockrc" -eq 200 ]; then echo "$id ERROR (another run of this scenario is active)"; exit 1; fi
if [ "$lockrc" -ne 0 ]; then echo "$id ERROR (flock failed with exit $lockrc; does this filesystem support it? see check_licences.sh)"; echo ERROR >"$status_file"; exit 1; fi

template=$WORK/prefs-template/matlabprefs.mat
[ -f "$template" ] || { echo "$id ERROR (no $template; run make_prefs_template.sh once)"; echo ERROR >"$status_file"; exit 1; }
scratch=
cleanup() { rm -f "$sdir/.matlab_result.json" "$sdir/.context.json"; [ -z "$scratch" ] || rm -rf "$scratch"; }
trap cleanup EXIT
trap 'echo ERROR >"$status_file"; exit 143' TERM INT
echo RUNNING >"$status_file"
rm -f "$PARITY_OUT/python/$id".seed*.json "$PARITY_OUT/matlab/$id.json"   # no result of an earlier batch survives

scratch=$(mktemp -d "${TMPDIR:-/tmp}/parity-$id-XXXXXX") || { echo "$id ERROR (mktemp failed)"; echo ERROR >"$status_file"; exit 1; }
export MATLAB_PREFDIR=$scratch/prefs      # private: concurrent sessions can corrupt a shared one
mkdir -p "$MATLAB_PREFDIR" && cp "$template" "$MATLAB_PREFDIR/" || { echo "$id ERROR (cannot prepare $MATLAB_PREFDIR)"; echo ERROR >"$status_file"; exit 1; }
rm -f "$sdir/.matlab_result.json" "$sdir/.context.json"    # never read a result left by an earlier run
mkdir -p "$PARITY_ROOT/results"   # git-ignored, so a fresh clone has none; some scenarios write files there
cd "$scratch"

status=ERROR      # becomes MATCH only inside the run block below, so a block that never runs stays ERROR
raise() { case $1 in ERROR) status=ERROR ;; DIFFER) [ "$status" = ERROR ] || status=DIFFER ;; esac; }
compare() {   # <left> <right> <message if they differ>
  # Exit code 1 means DIFFER, but an uncaught Python exception also exits 1, so the
  # verdict line the tool prints decides.
  local out; out=$(parity compare "$1" "$2" --scenario "$id" 2>&1); local rc=$?
  printf '%s\n' "$out"
  if   [ "$rc" -eq 0 ] && grep -q '^result:   MATCH' <<<"$out"; then :
  elif [ "$rc" -eq 1 ] && grep -q '^result:   DIFFER' <<<"$out"; then echo "$3"; raise DIFFER
  else echo "parity compare failed without a verdict (exit $rc)"; raise ERROR
  fi
}

{
  status=MATCH
  echo "scenario: $id  host: $(hostname)  start: $(date -u +%FT%TZ)  seeds: ${seeds[*]}"
  ok_seeds=
  for seed in "${seeds[@]}"; do
    if PYTHONHASHSEED=$seed parity run "$id" --impl python -o "$PARITY_OUT/python/$id.seed$seed.json"; then
      ok_seeds="$ok_seeds $seed"
    else
      echo "PYTHON FAILED (seed $seed)"; raise ERROR
    fi
  done
  # A MATLAB session that hangs (one in ten runs of a 40 s scenario stalled at 0% CPU, with no
  # MATLAB stack and a loaded node) would otherwise hold its licence seat until Slurm's time
  # limit. Each attempt runs in its own process group so that the timeout reaches MATLAB
  # itself, not only the parity process that started it. A scenario that timed out is tried
  # once more; a second timeout is an ERROR. MATLAB_TIMEOUT is in seconds (default 1800).
  # fd 9 holds the scenario lock; the MATLAB step and the watchdog must not inherit it, or a
  # leftover process would keep the scenario locked after this script has ended.
  limit=${MATLAB_TIMEOUT:-1800}
  for attempt in 1 2; do
    rm -f "$scratch/timed_out" "$sdir/.matlab_result.json"
    setsid parity run "$id" --impl matlab --matlab "$MATLAB" -o "$PARITY_OUT/matlab/$id.json" 9>&- &
    mpid=$!
    ( exec 9>&-; sleep "$limit" || exit 0; touch "$scratch/timed_out"; kill -TERM -- "-$mpid" 2>/dev/null
      sleep 30; kill -KILL -- "-$mpid" 2>/dev/null ) &
    wpid=$!
    wait "$mpid"; mrc=$?
    pkill -P "$wpid" 2>/dev/null; kill "$wpid" 2>/dev/null; wait "$wpid" 2>/dev/null
    if [ -f "$scratch/timed_out" ]; then
      echo "MATLAB TIMED OUT after $limit s (attempt $attempt of 2)"; mrc=124
      continue
    fi
    break
  done
  if [ "$mrc" -eq 0 ]; then
    for seed in $ok_seeds; do
      echo "--- python (PYTHONHASHSEED=$seed) vs matlab"
      compare "$PARITY_OUT/python/$id.seed$seed.json" "$PARITY_OUT/matlab/$id.json" "DIFFERS FROM MATLAB (seed $seed)"
    done
  else
    echo "MATLAB FAILED"; raise ERROR
  fi
  [ -n "$ok_seeds" ] || raise ERROR
  set -- $ok_seeds; first=${1:-}; [ $# -eq 0 ] || shift
  for seed in "$@"; do
    echo "--- python seed $first vs seed $seed"
    compare "$PARITY_OUT/python/$id.seed$first.json" "$PARITY_OUT/python/$id.seed$seed.json" "PYTHON RESULT DEPENDS ON THE HASH SEED"
  done
  echo "end: $(date -u +%FT%TZ)  status: $status"
} >"$log" 2>&1

echo "$status" >"$status_file" || { echo "$id ERROR (cannot write $status_file)"; exit 1; }
echo "$id $status"
[ "$status" != ERROR ]
