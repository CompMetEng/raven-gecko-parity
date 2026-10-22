#!/usr/bin/env bash
#SBATCH -J pyhang
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 01:30:00
# usage: sbatch pyhang.sh <out-dir> <scenario>; Python side 40 times; a stack dump after 150 s shows where a hang sits
export PARITY_OUT=$1/x PROBE_OUT=$1
source /cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/hpc-scripts/env.sh || exit 1
id=$2; d=$1/$id; mkdir -p $d
for n in $(seq 1 40); do
  s=$((n % 4))
  PYTHONHASHSEED=$s python -X faulthandler -c "
import faulthandler, sys
faulthandler.dump_traceback_later(150, exit=True, file=open('$d/stack-$n.txt','w'))
from parity.cli import main
sys.argv=['parity','run','$id','--impl','python','-o','$d/out-$n.json']
main()" > $d/run-$n.log 2>&1
  rc=$?
  echo "$n seed=$s rc=$rc $(wc -c < $d/stack-$n.txt)" >> $d/summary.txt
  [ -s $d/stack-$n.txt ] || rm -f $d/stack-$n.txt
done
