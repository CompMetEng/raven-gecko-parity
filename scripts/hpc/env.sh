#!/usr/bin/env bash
# Settings shared by the scripts in this folder. Do not copy or run this file: put site
# values in env.local.sh (start from env.local.example.sh). Scripts call
# `source env.sh || exit 1`.
#
# Layout under $WORK: RAVEN GECKO raven-toolbox geckopy raven-gecko-parity (parity.toml
# finds the code repositories as siblings of raven-gecko-parity), plus venv/ and out/.
#
# PARITY_OUT, the folder of one batch, has no default: give every batch a name, for example
#   export PARITY_OUT=$WORK/out/2026-10-21
# and use the same value for submit_slurm.sh, the gecko tests, reruns and collect.sh. A
# default derived from the date sends a rerun on the next day to a different folder.

if [ "${BASH_SOURCE[0]}" = "$0" ]; then echo "env.sh is meant to be sourced by the other scripts" >&2; exit 2; fi
_hpc_here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
[ -f "$_hpc_here/env.local.sh" ] && source "$_hpc_here/env.local.sh"

_hpc_die() { echo "env.sh: $*" >&2; return 1; }

: "${WORK:=$HOME/parity-work}"
: "${PARITY_ROOT:=$WORK/raven-gecko-parity}"
: "${PARITY_PYTHON_ENV:=$WORK/venv}"
: "${MATLAB:=matlab}"
: "${GRB_LICENSE_FILE:=$HOME/gurobi.lic}"
: "${GUROBI_HOME:=/opt/gurobi1303/linux64}"
[ -n "${PARITY_OUT:-}" ] || { _hpc_die "PARITY_OUT is not set; name the batch, e.g. export PARITY_OUT=$WORK/out/2026-10-21"; return 1; }
export WORK PARITY_ROOT PARITY_PYTHON_ENV MATLAB GRB_LICENSE_FILE GUROBI_HOME PARITY_OUT

# The paths go into MATLAB statements, Slurm options and path lists unquoted.
for _v in WORK PARITY_ROOT PARITY_OUT PARITY_PYTHON_ENV GUROBI_HOME; do
  case ${!_v} in
    /*) ;;
    *) _hpc_die "$_v must be an absolute path: ${!_v}"; return 1 ;;
  esac
  case ${!_v} in
    *[[:space:]\'\"\`\$\\,%:\;=]*) _hpc_die "$_v must not contain spaces, quotes, \$, backslash or any of , % : ; = : ${!_v}"; return 1 ;;
  esac
done

[ -d "$PARITY_ROOT/scenarios" ] || { _hpc_die "PARITY_ROOT is not a raven-gecko-parity checkout: $PARITY_ROOT"; return 1; }
[ -f "$PARITY_PYTHON_ENV/bin/activate" ] || { _hpc_die "no virtualenv at PARITY_PYTHON_ENV=$PARITY_PYTHON_ENV"; return 1; }
source "$PARITY_PYTHON_ENV/bin/activate"
command -v parity >/dev/null || { _hpc_die "parity is not installed in $PARITY_PYTHON_ENV (pip install -e raven-gecko-parity)"; return 1; }
export PARITY_CONFIG=$PARITY_ROOT/parity.toml

# Gurobi's MATLAB interface. Prepended once, however often this file is sourced.
GUROBI_MATLAB=$GUROBI_HOME/matlab
if [ -d "$GUROBI_MATLAB" ]; then
  case ":${MATLABPATH:-}:" in *":$GUROBI_MATLAB:"*) ;; *) export MATLABPATH=$GUROBI_MATLAB${MATLABPATH:+:$MATLABPATH} ;; esac
  case ":${LD_LIBRARY_PATH:-}:" in *":$GUROBI_HOME/lib:"*) ;; *) export LD_LIBRARY_PATH=$GUROBI_HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH} ;; esac
else
  GUROBI_MATLAB=
  echo "env.sh: warning: $GUROBI_HOME/matlab does not exist; Gurobi scenarios and MILP tests will fail or fall back" >&2
fi
export GUROBI_MATLAB
