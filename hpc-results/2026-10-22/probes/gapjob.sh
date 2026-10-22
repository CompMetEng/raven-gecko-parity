#!/usr/bin/env bash
#SBATCH -J gapp1
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 01:30:00
# gap_report_toy with cobra processes=1 (sitecustomize), then the 20 hash seeds against its MATLAB result
here=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/probes2
export PARITY_PYTHON_ENV=$here/venv-p1
export PARITY_OUT=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/out/gap-p1 PROBE_OUT=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/out/probes-A
cd /cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/hpc-scripts
./run_scenario.sh gap_report_toy
$here/hashseeds.sh gap_report_toy $PARITY_OUT
