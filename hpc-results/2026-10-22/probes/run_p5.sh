#!/usr/bin/env bash
set -uo pipefail
export PARITY_OUT=${PARITY_OUT:-$PROBE_OUT}; source "$(dirname "$0")/../hpc-scripts/env.sh" || exit 1
here=$(cd "$(dirname "$0")" && pwd)
export P5_DIR=$PROBE_OUT/p5-files
python -I "$here/p5_make.py" "$P5_DIR" "$WORK/RAVEN/tutorials/data/models/smallYeast.yml" > "$PROBE_OUT/p5_python.log" 2>&1 || { echo python failed; cat "$PROBE_OUT/p5_python.log"; exit 1; }
exec "$here/run_probe.sh" p5_readback
