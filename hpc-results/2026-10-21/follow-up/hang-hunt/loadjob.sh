#!/usr/bin/env bash
#SBATCH -J huntload
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 01:30:00
# usage: sbatch loadjob.sh <tag>; 8 concurrent workers, each with its own copy of the traced checkout
TAG=$1
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
D=$WORK/probes/hunt-$TAG; mkdir -p $D
cd $WORK/raven-gecko-parity/scripts/hpc
worker() {
  k=$1
  mkdir -p $TMPDIR/w$k; cp -r $WORK/raven-gecko-parity-hunt $TMPDIR/w$k/raven-gecko-parity
  for r in RAVEN GECKO raven-toolbox geckopy; do ln -s $WORK/$r $TMPDIR/w$k/$r; done
  export PARITY_ROOT=$TMPDIR/w$k/raven-gecko-parity SEEDS=0
  for n in $(seq 1 15); do
    export PARITY_OUT=$WORK/out/hunt-$TAG-w$k-$n PARITY_TRACE=$TMPDIR/trace-w$k-$n.txt
    ./run_scenario.sh add_reactions_from_model_smallyeast > $D/w$k-run$n.txt 2>&1 &
    pid=$!
    for i in $(seq 1 48); do kill -0 $pid 2>/dev/null || break; sleep 5; done
    if kill -0 $pid 2>/dev/null; then
      F=$D/HANG-w$k-$n.txt
      echo "worker $k iteration $n HANG on $(hostname), load $(cut -d' ' -f1-3 /proc/loadavg)" > $F
      cp $PARITY_TRACE $D/hang-trace-w$k-$n.txt
      for m in $(pgrep -u $USER -f "MATLAB -batch"); do
        tr '\0' '\n' < /proc/$m/environ 2>/dev/null | grep -qx "PARITY_TRACE=$PARITY_TRACE" || continue
        echo "matlab pid $m" >> $F
        for t in /proc/$m/task/*; do echo "$(cat $t/comm) $(cat $t/wchan 2>/dev/null) $(awk '{print $3}' $t/stat)"; done | sort | uniq -c | sort -rn >> $F
        kill -SEGV $m
      done
      sleep 15; wait $pid 2>/dev/null; break
    fi
    wait $pid; cp $PARITY_TRACE $D/trace-w$k-$n.txt
    echo "w$k $n $(tail -1 $D/w$k-run$n.txt)" >> $D/summary.txt
  done
}
for k in 1 2 3 4 5 6 7 8; do worker $k & done
wait
echo "done" >> $D/summary.txt
