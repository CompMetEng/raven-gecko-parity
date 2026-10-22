#!/usr/bin/env bash
# usage: hashseeds.sh <scenario-id> <matlab-result-batch>; Python side under PYTHONHASHSEED 0..19, twice each
set -uo pipefail
export PARITY_OUT=${PARITY_OUT:-$PROBE_OUT}; source "$(dirname "$0")/../hpc-scripts/env.sh" || exit 1
id=$1; ref=$2/matlab/$id.json
d=$PROBE_OUT/hashseeds/$id; mkdir -p "$d"
[ -f "$ref" ] || { echo "$id NO-MATLAB-RESULT ($ref)"; exit 1; }
match=0; stable=0; bad=""
for s in $(seq 0 19); do
  for rep in a b; do
    PYTHONHASHSEED=$s parity run "$id" --impl python -o "$d/seed$s$rep.json" > "$d/seed$s$rep.log" 2>&1 || echo "PYTHON FAILED seed $s$rep" >> "$d/failures.txt"
  done
  if parity compare "$d/seed${s}a.json" "$ref" --scenario "$id" > "$d/compare-seed$s.txt" 2>&1 && grep -q '^result:   MATCH' "$d/compare-seed$s.txt"; then match=$((match+1)); else bad="$bad $s"; fi
  cmp -s "$d/seed${s}a.json" "$d/seed${s}b.json" && stable=$((stable+1))
  rm -f "$d/seed${s}b.log"
done
echo "$id seeds_match=$match/20 repeat_identical=$stable/20 differ_seeds:${bad:- none}" | tee "$d/summary.txt"
