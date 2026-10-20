#!/usr/bin/env bash
# usage: run_probe.sh <script.m>  (runs MATLAB with a private preference folder, glpk)
source /etc/profile.d/modules.sh; module purge
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
export MATLAB_PREFDIR=$(mktemp -d); cp -r $WORK/prefs-template/. $MATLAB_PREFDIR/ 2>/dev/null
export MATLABPATH=/apps/Arch/software/Gurobi/13.0.1-GCCcore-14.3.0/matlab
export LD_LIBRARY_PATH=/apps/Arch/software/Gurobi/13.0.1-GCCcore-14.3.0/lib:$LD_LIBRARY_PATH
export GRB_LICENSE_FILE=/apps/Arch/software/Gurobi/13.0.1-GCCcore-14.3.0/gurobi.lic
cd $(dirname "$1")
/apps/Common/software/MATLAB/2024b/bin/matlab -batch "run('$1')"
