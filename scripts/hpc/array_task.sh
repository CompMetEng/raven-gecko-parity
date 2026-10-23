#!/usr/bin/env bash
# Slurm array task: run line $SLURM_ARRAY_TASK_ID of a list file with one of the scripts.
#   array_task.sh <list-file> <script> [fixed args...]
# sbatch runs a copy of this file from Slurm's spool folder, so $0 does not locate the
# other scripts; submit_slurm.sh exports HPC_SCRIPTS for that. The list is read when the
# task starts, so every submission has its own.
set -uo pipefail
dir=${HPC_SCRIPTS:?HPC_SCRIPTS is not set; submit with submit_slurm.sh}
list=$1 script=$2; shift 2
line=$(sed -n "${SLURM_ARRAY_TASK_ID:?not a Slurm array task}p" "$list")
[ -n "$line" ] || { echo "no line ${SLURM_ARRAY_TASK_ID} in $list" >&2; exit 2; }
set -f                                  # the line is an id or a pair of words, never a pattern
exec "$dir/$script" "$@" $line
