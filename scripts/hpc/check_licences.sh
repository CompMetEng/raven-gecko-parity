#!/usr/bin/env bash
# Check, before submitting anything, that MATLAB starts and that Gurobi solves a tiny LP
# from both MATLAB and Python, in the same environment the jobs will use.
#   check_licences.sh
set -uo pipefail
source "$(dirname "$0")/env.sh" || exit 1

scratch=$(mktemp -d "${TMPDIR:-/tmp}/parity-check-XXXXXX"); trap 'rm -rf "$scratch"' EXIT
export MATLAB_PREFDIR=$scratch/prefs; mkdir -p "$MATLAB_PREFDIR"
rc=0

echo "== flock on the scenario folders (run_scenario.sh needs it to work across nodes)"
lock=$PARITY_ROOT/scenarios/.flock_check.$$
exec 8>"$lock" && flock -n 8 && echo "lock acquired" || { echo "flock failed"; rc=1; }
( exec 7>"$lock"; flock -n -E 200 7; [ $? -eq 200 ] ) && echo "a second holder is refused" || { echo "a second holder was NOT refused: flock does not work here"; rc=1; }
exec 8>&-; rm -f "$lock"

echo "== MATLAB"
"$MATLAB" -batch "disp(['MATLAB ' version('-release')]); disp(['userpath: ' userpath]); disp(['startup.m: ' which('startup')]);" || rc=1

echo "== Gurobi from MATLAB"
"$MATLAB" -batch "m.A=sparse(1); m.obj=1; m.rhs=1; m.sense='>'; m.modelsense='min'; r=gurobi(m, struct('outputflag',0)); disp(['gurobi status: ' r.status]); assert(strcmp(r.status,'OPTIMAL'));" || rc=1

echo "== Gurobi from Python"
python - <<'PY' || rc=1
import gurobipy as gp
m = gp.Model(); x = m.addVar(lb=1.0); m.setObjective(x); m.setParam("OutputFlag", 0); m.optimize()
print("gurobipy status:", m.Status, "(2 = optimal)"); assert m.Status == 2
PY

[ "$rc" = 0 ] && echo "all checks passed" || echo "a check failed; fix it before submitting" >&2
exit $rc
