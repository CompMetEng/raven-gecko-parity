#!/usr/bin/env bash
# Run one MATLAB test class of RAVEN, or the GECKO core tests.
#   run_matlab_tests.sh raven <t-class>      e.g. run_matlab_tests.sh raven tINIT
#   run_matlab_tests.sh gecko
# Writes $PARITY_OUT/matlab-tests/<name>.csv (one row per test), logs/<name>.log and
# status/test-<name>.txt: PASS, FAIL, INCOMPLETE (no failure, but tests were filtered or
# skipped), or ERROR (MATLAB failed before the tests finished), followed by the counts.
# Exit code: 0 unless ERROR.
#
# Which solver the tests use is decided by the tests, not by this script: RavenTestCase
# sets glpk when a class starts, and tests that need a MILP solver switch to Gurobi only
# if the function gurobi() is on the MATLAB path; otherwise they run on GLPK or SCIP
# without reporting a skip. The job solves one small LP with gurobi() first: the status
# line ends in NO-GUROBI when that failed (not on the path, no licence, or no free
# session), and GUROBI_ERROR in the log gives the reason. A PASS with NO-GUROBI means the
# MILP tests ran on the fallback solver.
set -uo pipefail
source "$(dirname "$0")/env.sh" || exit 1

suite=${1:-}
mkdir -p "$PARITY_OUT"/{matlab-tests,logs,status}
case $suite in
  raven)
    cls=${2:-}
    case $cls in ""|*[!A-Za-z0-9_]*) echo "usage: run_matlab_tests.sh raven <class>" >&2; exit 2 ;; esac
    [ -f "$WORK/RAVEN/tests/function_tests/$cls.m" ] || { echo "no test class $cls" >&2; exit 2; }
    name=raven-$cls
    run_body="cd('$WORK/RAVEN/tests/function_tests'); r = runtests({'$cls'});"
    ;;
  gecko)
    # GECKO's installer drops every GECKO folder from the path when the checkout path
    # contains "doc", and it saves the MATLAB path; run it when no scenario job is active.
    case $WORK in *doc*) echo "WORK must not contain 'doc' for the GECKO tests: $WORK" >&2; exit 2 ;; esac
    name=gecko-core
    # CI runs only geckoCoreFunctionTests; the other *Tests.m files of the folder run here too.
    classes=$(cd "$WORK/GECKO/test/unit_tests" && ls *Tests.m | sed "s/\.m\$//; s/.*/'&'/" | paste -sd, -)
    [ -n "$classes" ] || { echo "no *Tests.m files in $WORK/GECKO/test/unit_tests" >&2; exit 2; }
    run_body="cd('$WORK/GECKO'); GECKOInstaller.install; setRavenSolver('glpk'); cd('test/unit_tests'); r = runtests({$classes});"
    ;;
  *) echo "usage: run_matlab_tests.sh raven <class> | gecko" >&2; exit 2 ;;
esac

status_file=$PARITY_OUT/status/test-$name.txt
log=$PARITY_OUT/logs/$name.log
csv=$PARITY_OUT/matlab-tests/$name.csv
template=$WORK/prefs-template/matlabprefs.mat
[ -f "$template" ] || { echo "$name ERROR (no $template; run make_prefs_template.sh once)"; echo ERROR >"$status_file"; exit 1; }
scratch=
trap '[ -z "$scratch" ] || rm -rf "$scratch"' EXIT
trap 'echo ERROR >"$status_file"; exit 143' TERM INT
echo RUNNING >"$status_file"; rm -f "$csv"

scratch=$(mktemp -d "${TMPDIR:-/tmp}/parity-tests-XXXXXX") || { echo "$name ERROR (mktemp failed)"; echo ERROR >"$status_file"; exit 1; }
export MATLAB_PREFDIR=$scratch/prefs
mkdir -p "$MATLAB_PREFDIR" && cp "$template" "$MATLAB_PREFDIR/" || { echo "$name ERROR (cannot prepare $MATLAB_PREFDIR)"; echo ERROR >"$status_file"; exit 1; }
cd "$scratch"

cat >"$scratch/job.m" <<MATLAB
try
    restoredefaultpath;
    addpath(genpath('$WORK/RAVEN'));
    if ~isempty('$GUROBI_MATLAB') && isfolder('$GUROBI_MATLAB'), addpath('$GUROBI_MATLAB'); end
    fprintf('GUROBI_ON_PATH=%s\n', which('gurobi'));
    % A Gurobi that is on the path but unlicensed (or over its session limit) would make
    % every MILP test fail and read as a RAVEN problem. Solve one tiny LP to find out.
    gurobiOk = 0;
    if ~isempty(which('gurobi'))
        try
            gm.A = sparse(1); gm.obj = 1; gm.rhs = 1; gm.sense = '>'; gm.modelsense = 'min';
            gres = gurobi(gm, struct('outputflag', 0));
            gurobiOk = double(strcmp(gres.status, 'OPTIMAL'));
        catch gerr
            fprintf('GUROBI_ERROR=%s\n', strtrim(gerr.message));
        end
    end
    fprintf('GUROBI_OK=%d\n', gurobiOk);
    $run_body
    t = table(r);
    if ismember('Details', t.Properties.VariableNames), t.Details = []; end
    writetable(t, '$csv');
    fprintf('SUMMARY passed=%d failed=%d incomplete=%d total=%d\n', nnz([r.Passed]), nnz([r.Failed]), nnz([r.Incomplete]), numel(r));
    exit(0);
catch err
    fprintf(2, '%s\n', getReport(err));
    exit(2);
end
MATLAB

"$MATLAB" -batch "run('$scratch/job.m')" >"$log" 2>&1
code=$?

summary=$(grep -m1 '^SUMMARY' "$log" || true)
gurobi_ok=$(sed -n 's/^GUROBI_OK=//p' "$log" | head -1)
if [ "$code" -ne 0 ] || [ -z "$summary" ]; then
  state=ERROR
else
  failed=$(sed -n 's/.*failed=\([0-9]*\).*/\1/p' <<<"$summary")
  incomplete=$(sed -n 's/.*incomplete=\([0-9]*\).*/\1/p' <<<"$summary")
  total=$(sed -n 's/.*total=\([0-9]*\).*/\1/p' <<<"$summary")
  if   [ "${total:-0}" -eq 0 ]; then state=ERROR
  elif [ "${failed:-0}" -gt 0 ]; then state=FAIL
  elif [ "${incomplete:-0}" -gt 0 ]; then state=INCOMPLETE
  else state=PASS
  fi
fi
note=
[ "${gurobi_ok:-0}" = 1 ] || note=" NO-GUROBI"
echo "$state ${summary#SUMMARY }$note" >"$status_file"
echo "$name $state ${summary#SUMMARY }$note"
[ "$state" != ERROR ]
