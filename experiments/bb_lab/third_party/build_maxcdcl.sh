#!/usr/bin/env bash
# Fetch MaxCDCL (MSE 2023 source artifact), apply the qec-lab patch,
# and build both binaries:
#   maxcdcl_stock — pristine MSE 2023 solver (baseline)
#   tandem        — the qec-lab fork (-cost-step, -prime-vars)
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

# apply the qec-lab patch (idempotent: skip if already applied) and rebuild.
# The patch's paths are already `a/code/...`, so -p1 lands them relative to
# the MaxCDCL root — do NOT add `-d code` (that looks for code/code/...).
# --batch keeps a mismatch from stalling on patch's interactive prompt.
cd ../..
if ! grep -q costStep code/core/Solver.h; then
    patch -p1 --batch --forward < "$HERE/maxcdcl-qeclab.patch"
    grep -q costStep code/core/Solver.h || {
        echo "error: patch reported success but costStep is absent from" \
             "code/core/Solver.h — the fork would be a stock binary." >&2
        exit 1
    }
fi
cd code/simp
# clean first: the patch changes the Solver class layout (Solver.h) and
# the stock makefile's dependency tracking is not reliable across that —
# stale objects from the stock build segfault at startup.
MROOT=.. make clean >/dev/null 2>&1 || true
MROOT=.. make r >/dev/null
cp maxcdcl_release tandem

echo "built:"
echo "  $(pwd)/maxcdcl_stock  (pristine MSE 2023)"
echo "  $(pwd)/tandem         (qec-lab fork)"
