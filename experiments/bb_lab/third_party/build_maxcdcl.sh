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

# Stock build first, preserved as the baseline binary.
#
# `make clean` is load-bearing on a *rerun*. The previous run left release
# objects (*.or) compiled from patched sources, and `unzip -o` restores
# the pristine sources with their 2023 archive timestamps — older than
# those objects. make would judge the stale patched objects up to date and
# relink them straight into `maxcdcl_stock`, so the "pristine baseline"
# would silently BE the fork and every A/B would compare it to itself.
cd MaxCDCL/code/simp
MROOT=.. make clean >/dev/null 2>&1 || true
MROOT=.. make r >/dev/null
cp maxcdcl_release maxcdcl_stock

# Apply the qec-lab patch (idempotent: skip if already applied) and rebuild.
#
# The merged patch (v5, fiber-lb included) carries `a/core/...` paths, so
# it applies with -p1 from the MaxCDCL root WITHOUT `-d code` (one strip
# level from `a/code/...` under `-d code` was the historical mislabelled-
# binary bug — both lines found it independently). --batch keeps a
# mismatch from dropping into patch's interactive "File to patch:" prompt;
# the explicit failure branch refuses to ship an unpatched 'tandem'.
cd ../..
if ! grep -q costStep code/core/Solver.h; then
    if ! patch --batch --forward -p1 < "$HERE/maxcdcl-qeclab.patch"; then
        echo "ERROR: maxcdcl-qeclab.patch did not apply — refusing to" >&2
        echo "       build an unpatched binary named 'tandem'. The" >&2
        echo "       upstream zip layout may have changed." >&2
        exit 1
    fi
fi
grep -q costStep code/core/Solver.h || {
    echo "ERROR: patch reported success but costStep is absent" >&2
    exit 1
}

cd code/simp
# clean first: the patch changes the Solver class layout (Solver.h) and
# the stock makefile's dependency tracking is not reliable across that —
# stale objects from the stock build segfault at startup.
MROOT=.. make clean >/dev/null 2>&1 || true
MROOT=.. make r >/dev/null
cp maxcdcl_release tandem

# The fork's entire value is in its flags, so prove the binary carries
# them rather than trusting the build. This also distinguishes the two
# binaries: stock does not know -cost-step.
#
# The help text is captured once rather than piped into grep: under
# `set -o pipefail`, `solver --help | grep -q` fails the whole pipeline
# because grep exits at the first match and SIGPIPEs the solver.
tandem_help="$(./tandem --help 2>&1 || true)"
for flag in cost-step prime-vars init-lb phase-file fiber-lb; do
    case "$tandem_help" in
        *"-$flag"*) ;;
        *) echo "ERROR: built 'tandem' is missing -$flag" >&2; exit 1 ;;
    esac
done
stock_help="$(./maxcdcl_stock --help 2>&1 || true)"
case "$stock_help" in
    *-cost-step*)
        echo "ERROR: 'maxcdcl_stock' knows -cost-step — the baseline is" >&2
        echo "       not pristine; delete $DIR and rebuild." >&2
        exit 1 ;;
esac

echo "built:"
echo "  $(pwd)/maxcdcl_stock  (pristine MSE 2023)"
echo "  $(pwd)/tandem         (qec-lab fork; 4 flags verified present)"
