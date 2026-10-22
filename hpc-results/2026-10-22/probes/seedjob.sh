#!/usr/bin/env bash
#SBATCH -J seeds
#SBATCH -n 2
#SBATCH --mem=12G
#SBATCH -t 03:00:00
# usage: sbatch --array=1-N seedjob.sh <probe-out> <matlab-batch>
here=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/probes2
export PROBE_OUT=$1
id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $here/seedlist.txt)
exec $here/hashseeds.sh $id $2
