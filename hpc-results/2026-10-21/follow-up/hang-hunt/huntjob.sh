#!/usr/bin/env bash
#SBATCH -J hunt
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 01:30:00
# usage: sbatch huntjob.sh <tag>; loops run_scenario with trace marks until one run exceeds 4 min
TAG=$1
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
export PARITY_ROOT=$WORK/raven-gecko-parity-hunt SEEDS=0
D=$WORK/probes/hunt-$TAG; mkdir -p $D
cd $WORK/raven-gecko-parity/scripts/hpc
for n in $(seq 1 40); do
  export PARITY_OUT=$WORK/out/hunt-$TAG-$n PARITY_TRACE=$TMPDIR/trace-$n.txt
  ./run_scenario.sh add_reactions_from_model_smallyeast > $D/run$n.txt 2>&1 &
  pid=$!
  for i in $(seq 1 48); do kill -0 $pid 2>/dev/null || break; sleep 5; done
  if kill -0 $pid 2>/dev/null; then
    echo "iteration $n HANG on $(hostname), load $(cut -d' ' -f1-3 /proc/loadavg)" > $D/HANG.txt
    cp $PARITY_TRACE $D/hang-trace.txt
    m=$(pgrep -u $USER -f "MATLAB -batch" | head -1); echo "matlab pid $m" >> $D/HANG.txt
    ps -u $USER -o pid,etime,pcpu,stat,wchan:20,args --sort=-pcpu | head -6 | cut -c1-180 >> $D/HANG.txt
    for t in /proc/$m/task/*; do echo "$(cat $t/comm) $(cat $t/wchan 2>/dev/null) $(awk '{print $3}' $t/stat)"; done | sort | uniq -c | sort -rn >> $D/HANG.txt
    ls -l /proc/$m/fd 2>/dev/null | grep -v -E "\.jar|/dev/null|socket|pipe|anon_inode" | cut -c40-220 >> $D/HANG.txt
    kill -SEGV $m; sleep 15; pkill -u $USER -f run_scenario.sh; wait $pid 2>/dev/null
    echo "hang found at iteration $n"; break
  fi
  wait $pid; cp $PARITY_TRACE $D/trace-$n.txt
  echo "$n $(tail -1 $D/run$n.txt)" >> $D/summary.txt
done
echo "done after $n iterations" >> $D/summary.txt
