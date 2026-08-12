#!/bin/zsh
# A36 orchestrator: transport verify (both axes on the cheap T6 point)
# -> full triage line with the validated sign.  Logs to data/a36/.
set -e
cd "$(dirname "$0")/.."
LOG=data/a36/orchestrator.log
echo "[$(date +%H:%M:%S)] verify T6 x" >> "$LOG"
uv run python scripts/a36_light_census.py verify --point T6 --axis x \
  --W 18 --threads 4 > data/a36/verify_T6x.log 2>&1
SIGNX=$(grep -o 'sign [+-]1' data/a36/verify_T6x.log | awk '{print $2}')
echo "[$(date +%H:%M:%S)] T6 x sign: $SIGNX" >> "$LOG"
echo "[$(date +%H:%M:%S)] verify T6 y" >> "$LOG"
uv run python scripts/a36_light_census.py verify --point T6 --axis y \
  --W 18 --threads 4 > data/a36/verify_T6y.log 2>&1
SIGNY=$(grep -o 'sign [+-]1' data/a36/verify_T6y.log | awk '{print $2}')
echo "[$(date +%H:%M:%S)] T6 y sign: $SIGNY" >> "$LOG"
if [ -z "$SIGNX" ] || [ "$SIGNX" != "$SIGNY" ]; then
  echo "[$(date +%H:%M:%S)] SIGN MISMATCH/EMPTY ($SIGNX vs $SIGNY) — abort" >> "$LOG"
  exit 1
fi
echo "[$(date +%H:%M:%S)] launching triage line sign=$SIGNX" >> "$LOG"
zsh scripts/a36_triage_batch.sh "$SIGNX" > data/a36/triage_line.log 2>&1
echo "[$(date +%H:%M:%S)] triage line done" >> "$LOG"
