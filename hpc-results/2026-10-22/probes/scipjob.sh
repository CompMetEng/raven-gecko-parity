#!/usr/bin/env bash
#SBATCH -J scip
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 03:00:00
# usage: sbatch scipjob.sh <out>; the three test classes without Gurobi on the path (SCIP, else GLPK)
export PARITY_OUT=$1
cd /cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/probes2
for c in tINIT tGapfilling tTasks; do ./run_tests_scip.sh raven $c; done
