#!/usr/bin/env bash
#SBATCH -J pyhang2
#SBATCH -n 8
#SBATCH --mem=24G
#SBATCH -t 00:20:00
# usage: sbatch pyhang2.sh <out-dir> <scenario> <seed>...; one Python run per seed, all at once, stack dump after 150 s
export PARITY_OUT=$1/x PROBE_OUT=$1
source /cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work/hpc-scripts/env.sh || exit 1
id=$2; shift 2; d=$PROBE_OUT/$id; mkdir -p $d
i=0
for s in "$@"; do
  i=$((i+1))
  ( PYTHONHASHSEED=$s python -X faulthandler -c "
import faulthandler, sys
faulthandler.dump_traceback_later(150, exit=True, file=open('$d/stack2-$i-seed$s.txt','w'))
from parity.cli import main
sys.argv=['parity','run','$id','--impl','python','-o','$d/out2-$i-seed$s.json']
main()" > $d/run2-$i-seed$s.log 2>&1; echo "run $i seed=$s rc=$? stack_bytes=$(wc -c < $d/stack2-$i-seed$s.txt)" >> $d/summary2.txt ) &
done
wait
