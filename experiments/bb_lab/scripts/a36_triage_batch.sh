#!/bin/zsh
# A36 triage assembly line — sequential per-group censuses+scans, both
# axes, ordered by q_cover.  Run from experiments/bb_lab/ AFTER the
# transport sign is validated (pass it as $1, e.g. +1 or -1).
set -e
SIGN=${1:?usage: a36_triage_batch.sh <sign>}
run() {
  echo "=== triage $1 axis $2 W=$3 ==="
  uv run python scripts/a36_light_census.py triage --sign "$SIGN" \
    --only-group "$1" --axis "$2" --w-override "$3" --threads 4
}
# q_cover 19.2 (d=12 -> W=18 pre-tier; W=22 exactness via T1.5/certify)
run Z5xZ12 y 18
run Z5xZ12 x 18
run Z6xZ10 x 18
run Z6xZ10 y 18
# q_cover 19.05 (k=12 -> W=14 fast tier)
run Z7xZ9 x 14
run Z7xZ9 y 14
# q_cover 17.78 / 17.63
run Z15xZ3 x 18
run Z15xZ3 y 18
run Z7xZ7 x 18
run Z7xZ7 y 18
# q_cover 14.81 (bb_108's parameter point - other codes than bb_108)
run Z9xZ6 x 18
run Z9xZ6 y 18
# q_cover 15.43 (d=12 -> W=18 pre-tier)
run Z7xZ8 x 18
run Z7xZ8 y 18
# q_cover 14.29
run Z6xZ7 x 18
run Z6xZ7 y 18
echo "=== triage line COMPLETE ==="
