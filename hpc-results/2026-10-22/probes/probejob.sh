#!/usr/bin/env bash
#SBATCH -J probe
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 01:00:00
# usage: sbatch --array=... probejob.sh <batch-folder>; task i runs the i-th line of the list
export PROBE_OUT=$1
here=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/probes2
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $here/list.txt)
case $line in p5) exec $here/run_p5.sh ;; *) exec $here/run_probe.sh $line ;; esac
