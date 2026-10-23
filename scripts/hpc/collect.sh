#!/usr/bin/env bash
# Bundle one batch for sending back: result files, logs, test tables, STATUS.tsv and
# MANIFEST.txt (revisions and versions). Does not start MATLAB.
#   collect.sh [batch-folder]        default $PARITY_OUT; writes <batch-folder>.tar.gz
set -uo pipefail
source "$(dirname "$0")/env.sh" || exit 1

out=${1:-$PARITY_OUT}
[ -d "$out" ] || { echo "no such folder: $out" >&2; exit 1; }
out=$(cd "$out" && pwd -P)

{
  echo "collected: $(date -u +%FT%TZ) on $(hostname)"
  # The revisions are read now. If a repository moved since the batch ran, they are wrong
  # for it; `dirty` counts modified files, which would make the results unreproducible.
  for r in RAVEN GECKO raven-toolbox geckopy raven-gecko-parity; do
    printf '%s %s %s dirty=%s\n' "$r" "$(git -C "$WORK/$r" rev-parse HEAD 2>&1)" "$(git -C "$WORK/$r" rev-parse --abbrev-ref HEAD 2>&1)" \
      "$(git -C "$WORK/$r" status --porcelain 2>/dev/null | wc -l)"
  done
  # Each result document records the runtime that wrote it.
  grep -h -o '"runtime": *"[^"]*"' "$out"/matlab/*.json "$out"/python/*.json 2>/dev/null | sort | uniq -c
  echo "gurobipy: $(python -c 'import gurobipy; print(gurobipy.gurobi.version())' 2>&1 | tail -1)"
  echo "cobra: $(python -c 'import cobra; print(cobra.__version__)' 2>&1 | tail -1)"
  echo "raven_toolbox: $(python -c 'import raven_toolbox; print(raven_toolbox.__file__)' 2>&1 | tail -1)"
} >"$out/MANIFEST.txt"

: >"$out/STATUS.tsv"
for f in "$out"/status/*.txt; do
  [ -f "$f" ] && printf '%s\t%s\n' "$(basename "$f" .txt)" "$(head -1 "$f")" >>"$out/STATUS.tsv"
done
sort -o "$out/STATUS.tsv" "$out/STATUS.tsv"
echo "status of the jobs that reported:"
cut -f2 "$out/STATUS.tsv" | awk '{print $1}' | sort | uniq -c

# Jobs that were submitted but never reported (still queued, or killed before they started).
if ls "$out"/lists/scenarios-*.txt >/dev/null 2>&1; then
  cat "$out"/lists/scenarios-*.txt | sort -u | while read -r id; do
    [ -f "$out/status/$id.txt" ] || echo "NO STATUS  $id"
  done | tee "$out/MISSING.txt"
fi
if grep -q RUNNING "$out/STATUS.tsv"; then
  echo "warning: jobs still marked RUNNING (running, or killed); the archive may be incomplete:" >&2
  grep RUNNING "$out/STATUS.tsv" >&2
fi

tmp=$out.tar.gz.tmp
tar -czf "$tmp" -C "$(dirname "$out")" "$(basename "$out")"; rc=$?
# tar exits 1 when a file changed while it was read (a job still writing); that archive is usable.
if [ "$rc" -le 1 ] && mv "$tmp" "$out.tar.gz"; then
  [ "$rc" -eq 0 ] || echo "warning: files changed while the archive was written" >&2
  echo "wrote $out.tar.gz"
else
  rm -f "$tmp"; echo "tar failed (exit $rc); nothing written" >&2; exit 1
fi
