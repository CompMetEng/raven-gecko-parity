#!/usr/bin/env bash
#SBATCH -J pr139
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
# Scenarios of PR 139 (9326b3f) with RAVEN at PR 11 and raven-toolbox at PR 181; both restored afterwards.
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
R=$WORK/RAVEN; T=$WORK/raven-toolbox
restore() { git -C $R checkout -q develop; git -C $T checkout -q develop; git -C $R status -sb | head -3; git -C $T status -sb | head -3; }
trap restore EXIT
git -C $R checkout -q 875e45b && git -C $T checkout -q 8fbaed4 || exit 1
git -C $R rev-parse --short HEAD; git -C $T rev-parse --short HEAD
export PARITY_ROOT=$WORK/raven-gecko-parity-139
git -C $PARITY_ROOT rev-parse --short HEAD
mkdir -p $PARITY_ROOT/results
cd $WORK/raven-gecko-parity/scripts/hpc
export PARITY_OUT=$WORK/out/pr139
for s in remove_low_score_genes_ties_smallyeast task_list_edge_cases yaml_float_digits_smallyeast task_list_parsing replace_metabolite_toy; do ./run_scenario.sh $s; done
