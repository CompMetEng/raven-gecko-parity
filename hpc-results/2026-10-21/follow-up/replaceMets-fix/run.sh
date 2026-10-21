#!/usr/bin/env bash
#SBATCH -J peer
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
P=$WORK/peer; R=$WORK/RAVEN
cd $WORK/raven-gecko-parity/scripts/hpc
(cd $R && git checkout -- . && git apply --recount --ignore-whitespace $P/test.patch)
export PARITY_OUT=$WORK/out/peer-A; ./run_matlab_tests.sh raven tManipulation
(cd $R && git apply --recount --ignore-whitespace $P/src.patch)
export PARITY_OUT=$WORK/out/peer-B; ./run_matlab_tests.sh raven tManipulation
export PARITY_OUT=$WORK/out/peer-C; SEEDS=0 ./run_scenario.sh replace_metabolite_toy
(cd $R && git checkout -- . ; git status -s | head; git rev-parse --short HEAD)
