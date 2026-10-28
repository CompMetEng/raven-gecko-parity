#!/usr/bin/env bash
# Run once per installation. Creates $WORK/prefs-template/matlabprefs.mat with RAVEN's
# solver set to glpk. Every job copies it into its own private MATLAB preference folder,
# because RAVEN stores its solver choice in the MATLAB preferences and a fresh folder has
# none (optimizeProb then stops with "RAVEN solver not defined"). A scenario that needs
# another solver sets it itself and restores this one afterwards.
set -uo pipefail
source "$(dirname "$0")/env.sh" || exit 1

export MATLAB_PREFDIR=$WORK/prefs-template
mkdir -p "$MATLAB_PREFDIR"
"$MATLAB" -batch "restoredefaultpath; addpath(genpath('$WORK/RAVEN')); setRavenSolver('glpk'); if ispref('RAVEN','ravenPath'), rmpref('RAVEN','ravenPath'); end; disp(['solver: ' getpref('RAVEN','solver')]);" || exit 1
[ -f "$MATLAB_PREFDIR/matlabprefs.mat" ] || { echo "no matlabprefs.mat was written to $MATLAB_PREFDIR" >&2; exit 1; }
echo "prefs template written to $MATLAB_PREFDIR"
