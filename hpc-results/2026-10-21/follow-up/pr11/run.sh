#!/usr/bin/env bash
#SBATCH -J pr11
#SBATCH -n 2
#SBATCH --mem=24G
#SBATCH -t 05:00:00
# RAVEN PR 11 head and raven-toolbox PR 181 head checked out in place, restored at the end.
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
R=$WORK/RAVEN; T=$WORK/raven-toolbox
restore() { git -C $R checkout -q develop; git -C $T checkout -q develop; git -C $R status -sb | head -3; git -C $T status -sb | head -3; }
trap restore EXIT
git -C $R checkout -q 875e45b && git -C $T checkout -q 8fbaed4 || exit 1
git -C $R rev-parse --short HEAD; git -C $T rev-parse --short HEAD
cd $WORK/raven-gecko-parity/scripts/hpc
export PARITY_OUT=$WORK/out/pr11
for c in tTasks tIO tManipulation tSyntax tINIT; do ./run_matlab_tests.sh raven $c; done
for s in remove_low_score_genes_smallyeast replace_metabolite_toy task_checking_smallyeast task_list_parsing yaml_roundtrip_smallyeast; do ./run_scenario.sh $s; done
