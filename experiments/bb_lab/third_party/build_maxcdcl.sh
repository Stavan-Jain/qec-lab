#!/usr/bin/env bash
# Fetch MaxCDCL (MSE 2023 source artifact), apply the qec-lab patch,
# and build both binaries:
#   maxcdcl_stock   — pristine MSE 2023 solver (baseline)
#   maxcdcl_release — qec-lab fork (-cost-step, -prime-vars)
#
# Usage:  ./build_maxcdcl.sh [target-dir]     (default: ./maxcdcl)
# Needs:  curl, unzip, make, a C++ compiler, zlib headers.
set -euo pipefail

DIR="${1:-maxcdcl}"
URL="https://maxsat-evaluations.github.io/2023/mse23-solver-src/exact/MaxCDCL.zip"
HERE="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$DIR" && cd "$DIR"
[ -f MaxCDCL.zip ] || curl -sL -o MaxCDCL.zip "$URL"
unzip -q -o MaxCDCL.zip

# stock build first, preserved as the baseline binary
cd MaxCDCL/code/simp
MROOT=.. make r >/dev/null
cp maxcdcl_release maxcdcl_stock

# apply the qec-lab patch (idempotent: skip if already applied) and rebuild
cd ../..
if ! grep -q costStep code/core/Solver.h; then
    patch -p1 -d code < "$HERE/maxcdcl-qeclab.patch"
fi
cd code/simp
MROOT=.. make r >/dev/null

echo "built:"
echo "  $(pwd)/maxcdcl_stock    (pristine MSE 2023)"
echo "  $(pwd)/maxcdcl_release  (qec-lab fork)"
