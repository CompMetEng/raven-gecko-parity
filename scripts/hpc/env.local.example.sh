#!/usr/bin/env bash
# Copy to env.local.sh (ignored by git) and edit. env.sh reads it first, so it may set any
# subset of these variables. The `: "${VAR:=value}"` form keeps a value given on the command
# line or in the job environment ahead of the one written here.

: "${WORK:=$HOME/parity-work}"                  # shared, fast filesystem; absolute, no spaces, quotes or : , % ; = ; no "doc" for the GECKO tests
: "${MATLAB:=matlab}"                           # full path, or load the site module before submitting
: "${GRB_LICENSE_FILE:=$HOME/gurobi.lic}"
: "${GUROBI_HOME:=/opt/gurobi1303/linux64}"     # must contain matlab/ (the MATLAB interface) and lib/
# : "${PARITY_PYTHON_ENV:=$WORK/venv}"          # virtualenv with raven-gecko-parity, raven-toolbox and gurobipy
# PARITY_OUT is deliberately not set here: name every batch, e.g. export PARITY_OUT=$WORK/out/2026-10-21
export WORK MATLAB GRB_LICENSE_FILE GUROBI_HOME
