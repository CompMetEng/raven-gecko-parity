#!/usr/bin/env bash
# usage: run_probe.sh <probe script name without .m>; writes JSON to $PROBE_OUT
set -uo pipefail
export PARITY_OUT=${PARITY_OUT:-$PROBE_OUT}; source "$(dirname "$0")/../hpc-scripts/env.sh" || exit 1
here=$(cd "$(dirname "$0")" && pwd)
p=$1
export RAVEN_ROOT=$WORK/RAVEN
scratch=$(mktemp -d "${TMPDIR:-/tmp}/probe-$p-XXXXXX"); trap 'rm -rf "$scratch"' EXIT
export MATLAB_PREFDIR=$scratch/prefs; mkdir -p $MATLAB_PREFDIR; cp $WORK/prefs-template/matlabprefs.mat $MATLAB_PREFDIR/
mkdir -p "$PROBE_OUT"; cd "$scratch"
"$MATLAB" -batch "restoredefaultpath; addpath(genpath('$WORK/RAVEN')); addpath('$here'); setRavenSolver('glpk'); try, run('$here/$p.m'); catch err, disp(getReport(err)); pr_common('save', '${p}_FAILED', struct('error', err.message)); end" > "$PROBE_OUT/$p.log" 2>&1
echo "$p exit $?"
