#!/usr/bin/env bash
#SBATCH -J exportscen
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
# usage: sbatch exportscen.sh <raven-commit> <toolbox-commit> <batch> <scenario>...; the commits are checked out in place and restored
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
R=$WORK/RAVEN; T=$WORK/raven-toolbox; rc=$1; tc=$2; b=$3; shift 3
trap 'git -C $R checkout -q develop; git -C $T checkout -q develop; git -C $R status -sb | head -2; git -C $T status -sb | head -2' EXIT
git -C $R checkout -q $rc && git -C $T checkout -q $tc || exit 1
export PARITY_ROOT=$WORK/raven-gecko-parity-t1 PARITY_OUT=$WORK/out/$b
cd $WORK/hpc-scripts
for s in "$@"; do ./run_scenario.sh $s; done
