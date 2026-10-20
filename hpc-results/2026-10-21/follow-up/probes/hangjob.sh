#!/usr/bin/env bash
#SBATCH -J hang
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 00:40:00
# usage: sbatch hangjob.sh <N>; one run of add_reactions_from_model_smallyeast, watchdog at 8 min
N=$1
export WORK=/cephyr/NOBACKUP/groups/compmeteng/eduardk/parity-work
export PARITY_OUT=$WORK/out/hang-$N
D=$WORK/probes/hang; mkdir -p $D
cd $WORK/raven-gecko-parity/scripts/hpc
echo "run $N on $(hostname) start $(date -u +%T)  node load: $(cut -d' ' -f1-3 /proc/loadavg)" > $D/run$N.txt
./run_scenario.sh add_reactions_from_model_smallyeast >> $D/run$N.txt 2>&1 &
pid=$!
for i in $(seq 1 96); do kill -0 $pid 2>/dev/null || break; sleep 5; done
if kill -0 $pid 2>/dev/null; then
  echo "HANG: still running after 8 min" >> $D/run$N.txt
  ps -u $USER -o pid,etime,pcpu,args --sort=-pcpu | head -8 | cut -c1-200 >> $D/run$N.txt
  m=$(pgrep -u $USER -f "MATLAB -batch" | head -1)
  if [ -n "$m" ]; then top -H -b -n1 -p $m | head -20 >> $D/run$N.txt; kill -SEGV $m; sleep 20; fi
  cp $PARITY_OUT/logs/add_reactions_from_model_smallyeast.log $D/run$N.log 2>/dev/null
  pkill -u $USER -f run_scenario.sh; sleep 2
fi
wait $pid 2>/dev/null
echo "end $(date -u +%T) status: $(cat $PARITY_OUT/status/add_reactions_from_model_smallyeast.txt 2>/dev/null)" >> $D/run$N.txt
