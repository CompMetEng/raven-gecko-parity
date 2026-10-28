#!/usr/bin/env bash
# Submit work as Slurm job arrays. Other schedulers: the unit of work is one call to
# run_scenario.sh <id>, run_matlab_tests.sh raven <class> or run_matlab_tests.sh gecko.
#   submit_slurm.sh all [raven|gecko]     scenarios, then the RAVEN tests, then the GECKO tests,
#                                         each starting when the previous one has ended
#   submit_slurm.sh scenarios [raven|gecko]
#   submit_slurm.sh tests                 the RAVEN test classes
#   submit_slurm.sh gecko-tests           the GECKO unit tests; run after the scenarios ended
# PARITY_OUT must be set (it names the batch). MAX_PARALLEL (default 8) caps the concurrent
# tasks of one array. Every task holds a MATLAB licence seat and may hold a Gurobi session;
# `all` keeps the arrays apart, and the single genome-scale scenario runs next to the
# scenario array, so at most MAX_PARALLEL + 1 seats are in use.
# tInstallation is not submitted: it needs the network and checkRaven saves the MATLAB path.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/env.sh" || exit 1
case $here in *[[:space:]\'\"\`\$\\,%:\;=]*) echo "the scripts path must not contain spaces, quotes or , % : ; = : $here" >&2; exit 2 ;; esac
command -v sbatch >/dev/null || { echo "sbatch not found" >&2; exit 2; }

MAX_PARALLEL=${MAX_PARALLEL:-8}
HEAVY="ec_model_full_yeastgem"      # genome scale: more memory and time, submitted separately
batch=$(date -u +%Y%m%dT%H%M%S)-$$
lists=$PARITY_OUT/lists; mkdir -p "$lists" "$PARITY_OUT/slurm"
export_opt=ALL,PARITY_OUT=$PARITY_OUT,HPC_SCRIPTS=$here

submit() {   # <list-file> <script> <script args...> -- <sbatch args...>; prints the job id
  local list=$1 script=$2; shift 2
  local fixed=(); while [ "$1" != -- ]; do fixed+=("$1"); shift; done; shift
  local n; n=$(grep -c . "$list" || true)
  [ "$n" -gt 0 ] || return 0
  sbatch --parsable --array="1-$n%$MAX_PARALLEL" --job-name=parity --export="$export_opt" \
    --output="$PARITY_OUT/slurm/%x-%A_%a.out" "$@" "$here/array_task.sh" "$list" "$script" ${fixed[@]+"${fixed[@]}"}
}

submit_scenarios() {   # [pair]; prints the job ids, space separated
  local pair=${1:-}
  case $pair in ""|raven|gecko) ;; *) echo "pair must be raven or gecko: $pair" >&2; return 2 ;; esac
  local all=$lists/scenarios-$batch.txt
  # ids of this pair, read with a YAML parser (a quoted or commented `pair:` line is fine)
  python - "$PARITY_ROOT/scenarios" "$pair" >"$all" <<'PY'
import sys, pathlib, yaml
root, pair = pathlib.Path(sys.argv[1]), sys.argv[2]
for yml in sorted(root.glob("*/scenario.yml")):
    if not pair or yaml.safe_load(yml.read_text(encoding="utf-8")).get("pair") == pair:
        print(yml.parent.name)
PY
  [ -s "$all" ] || { echo "no scenarios found" >&2; return 1; }
  grep -vxF -f <(printf '%s\n' $HEAVY) "$all" >"$lists/normal-$batch.txt" || true
  grep -xF  -f <(printf '%s\n' $HEAVY) "$all" >"$lists/heavy-$batch.txt"  || true
  echo "$(wc -l <"$all") scenarios (list: $all)" >&2
  local a b
  a=$(submit "$lists/normal-$batch.txt" run_scenario.sh -- --cpus-per-task=2 --mem=16G --time=02:00:00)
  b=$(MATLAB_TIMEOUT=${MATLAB_TIMEOUT_HEAVY:-18000} submit "$lists/heavy-$batch.txt" run_scenario.sh -- --cpus-per-task=4 --mem=64G --time=12:00:00)
  echo $a $b
}

submit_tests() {   # [dependency]; prints the job id
  local dep=${1:-} classes=$lists/raven-tests-$batch.txt
  ls "$WORK"/RAVEN/tests/function_tests/t*.m | xargs -n1 basename | sed 's/\.m$//' \
    | grep -vxE 'RavenTestCase|tBinaries|tInstallation' >"$classes"
  echo "$(wc -l <"$classes") RAVEN test classes" >&2
  submit "$classes" run_matlab_tests.sh raven -- --cpus-per-task=2 --mem=16G --time=04:00:00 ${dep:+--dependency=afterany:$dep}
}

submit_gecko() {   # [dependency]; prints the job id
  local dep=${1:-}
  sbatch --parsable --job-name=parity-gecko --export="$export_opt" --output="$PARITY_OUT/slurm/%x-%j.out" \
    --cpus-per-task=2 --mem=16G --time=04:00:00 ${dep:+--dependency=afterany:$dep} \
    --wrap "$here/run_matlab_tests.sh gecko"
}

join() { local IFS=:; echo "$*"; }

case ${1:-} in
  scenarios)   ids=$(submit_scenarios "${2:-}"); echo "submitted: $ids" ;;
  tests)       echo "submitted: $(submit_tests)" ;;
  gecko-tests) echo "submitted: $(submit_gecko)" ;;
  all)
    ids=$(submit_scenarios "${2:-}")
    t=$(submit_tests "$(join $ids)")
    g=$(submit_gecko "$t")
    echo "submitted: scenarios $ids, then tests $t, then gecko tests $g"
    ;;
  *) echo "usage: submit_slurm.sh all [pair] | scenarios [pair] | tests | gecko-tests" >&2; exit 2 ;;
esac
echo "results in $PARITY_OUT"
