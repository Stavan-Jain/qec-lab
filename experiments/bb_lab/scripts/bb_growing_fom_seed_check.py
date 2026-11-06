#!/usr/bin/env python3
"""Exact finite checks for the weight-six separable fractal BB candidate.

The enumeration is over a q-dimensional CLASSICAL kernel, not over quantum
logical cosets.  The square-matrix hypergraph-product distance theorem transfers
the exact classical minimum to the BB code.  No external packages are needed.

Run from the repository root:
  uv run --project experiments/bb_lab python \
    experiments/bb_lab/scripts/bb_growing_fom_seed_check.py

The default enumerates s=1,...,4 completely and checks only the algebraic
identity and upper witness for s=5,6.  --native-q32 also completely enumerates
s=5 using an optional compiled ARM NEON loop (requires an ARM64 C++ compiler).
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


Q32_NATIVE_SOURCE = r"""
#include <arm_neon.h>
#include <cstdint>
#include <cstdio>
int main() {
  constexpr unsigned q=32, n=1023;
  uint64_t h[16]={};
  // Independent combinatorial construction of (p^(q-1)+1)/x.
  for (unsigned a=0;a<q;++a) for(unsigned b=0;b<q;++b)
    if (!(a&b) && (a||b)) {
      unsigned e=a+q*b-1;
      h[e/64]^=uint64_t(1)<<(e%64);
    }
  uint8x16_t basis[32][8];
  for(unsigned j=0;j<q;++j) {
    uint64_t row[16]={};
    for(unsigned e=0;e<n-j;++e) if((h[e/64]>>(e%64))&1)
      row[(e+j)/64]^=uint64_t(1)<<((e+j)%64);
    for(unsigned k=0;k<8;++k)
      basis[j][k]=vld1q_u8(reinterpret_cast<uint8_t*>(row)+16*k);
  }
  uint8x16_t word[8];
  for(auto &part:word) part=vdupq_n_u8(0);
  unsigned minimum=n;
  uint64_t minimum_count=0, attaining_message=0, total_words=0, total_weight=0;
  for(uint64_t i=1;i<(uint64_t(1)<<q);++i) {
    unsigned j=__builtin_ctzll(i);
    uint8x16_t counts=vdupq_n_u8(0);
    for(unsigned k=0;k<8;++k) {
      word[k]=veorq_u8(word[k],basis[j][k]);
      // Each byte sums eight popcounts, hence is <=64: no byte overflow.
      counts=vaddq_u8(counts,vcntq_u8(word[k]));
    }
    unsigned weight=vaddlvq_u8(counts);
    ++total_words;
    total_weight+=weight;
    if(weight<minimum) {
      minimum=weight; minimum_count=1; attaining_message=i^(i>>1);
    } else if(weight==minimum) ++minimum_count;
  }
  printf("{\"minimum\":%u,\"minimum_word_count\":%llu,"
         "\"attaining_message\":%llu,\"nonzero_messages_checked\":%llu,"
         "\"total_weight\":%llu}\n",minimum,
         static_cast<unsigned long long>(minimum_count),
         static_cast<unsigned long long>(attaining_message),
         static_cast<unsigned long long>(total_words),
         static_cast<unsigned long long>(total_weight));
}
"""


def poly_mul(a: int, b: int) -> int:
    """Carryless multiplication over F_2; bit j is the x^j coefficient."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        a <<= 1
        b >>= 1
    return result


def poly_divmod(a: int, b: int) -> tuple[int, int]:
    if b == 0:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = 0
    while a and a.bit_length() >= b.bit_length():
        shift = a.bit_length() - b.bit_length()
        quotient ^= 1 << shift
        a ^= b << shift
    return quotient, a


def poly_pow(a: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = poly_mul(result, a)
        exponent >>= 1
        if exponent:
            a = poly_mul(a, a)
    return result


def exhaustive_minimum(generator: int, dimension: int) -> dict[str, object]:
    """Visit all nonzero messages in reflected Gray-code order exactly once."""
    word = 0
    histogram: Counter[int] = Counter()
    minimum = generator.bit_count()
    attaining_message = 1
    for step in range(1, 1 << dimension):
        changed_bit = (step & -step).bit_length() - 1
        word ^= generator << changed_bit
        weight = word.bit_count()
        if weight == 0:
            raise AssertionError("nonzero message encoded as zero")
        histogram[weight] += 1
        if weight < minimum:
            minimum = weight
            attaining_message = step ^ (step >> 1)
    assert sum(histogram.values()) == (1 << dimension) - 1
    witness = poly_mul(generator, attaining_message)
    assert witness.bit_count() == minimum
    return {
        "minimum": minimum,
        "attaining_message_hex": hex(attaining_message),
        "attaining_word_hex": hex(witness),
        "nonzero_messages_checked": sum(histogram.values()),
        "weight_distribution": dict(sorted(histogram.items())),
    }


def exhaustive_q32_native(generator: int) -> dict[str, object]:
    """Run the complete 2^32-1 enumeration; validate totals and witness in Python."""
    if platform.machine().lower() not in {"arm64", "aarch64"}:
        raise RuntimeError("--native-q32 requires ARM64 NEON")
    compiler = shutil.which("clang++") or shutil.which("c++")
    if compiler is None:
        raise RuntimeError("--native-q32 requires a C++ compiler")
    with tempfile.TemporaryDirectory(prefix="bb-fractal-q32-") as directory:
        source = Path(directory) / "enumerate.cpp"
        executable = Path(directory) / "enumerate"
        source.write_text(Q32_NATIVE_SOURCE)
        subprocess.run(
            [compiler, "-O3", "-std=c++17", str(source), "-o", str(executable)],
            check=True,
        )
        completed = subprocess.run(
            [str(executable)], check=True, text=True, capture_output=True
        )
    result = json.loads(completed.stdout)
    assert result["nonzero_messages_checked"] == (1 << 32) - 1
    assert result["total_weight"] == 1023 * (1 << 31)
    message = result.pop("attaining_message")
    witness = poly_mul(generator, message)
    assert 0 < message < (1 << 32)
    assert witness.bit_count() == result["minimum"]
    assert witness.bit_length() <= 1023
    result["attaining_message_hex"] = hex(message)
    result["attaining_word_hex"] = hex(witness)
    result["method"] = "complete Gray enumeration, ARM NEON"
    return result


def check_seed(s: int, max_dimension: int, native_q32: bool) -> dict[str, object]:
    q = 1 << s
    period = q * q - 1
    check = (1 << q) | 3
    modulus = (1 << period) | 1
    generator, remainder = poly_divmod(modulus, check)
    assert remainder == 0
    assert generator.bit_length() - 1 == period - q
    assert poly_mul(check, generator) == modulus

    # Independent closed-form construction: x*g = p^(q-1)+1.
    fractal = poly_pow(check, q - 1)
    assert fractal & 1
    assert generator == (fractal ^ 1) >> 1
    assert fractal.bit_count() == 3**s
    assert generator.bit_count() == 3**s - 1

    # Each shifted generator is an actual kernel vector. Distinct degrees give
    # linear independence; polynomial cancellation proves these span the kernel.
    for j in range(q):
        assert poly_divmod(poly_mul(check, generator << j), modulus)[1] == 0

    result: dict[str, object] = {
        "s": s,
        "q": q,
        "classical_length": period,
        "classical_dimension": q,
        "check_support": [0, 1, q],
        "generator_hex": hex(generator),
        "proved_upper_witness_weight": generator.bit_count(),
        "distance_status": "UPPER_BOUND_ONLY",
        "quantum_n": 2 * period * period,
        "quantum_k": 2 * q * q,
        "quantum_check_weight": 6,
    }
    if q <= max_dimension or (q == 32 and native_q32):
        exact = exhaustive_q32_native(generator) if q == 32 else exhaustive_minimum(generator, q)
        result["enumeration"] = exact
        distance = int(exact["minimum"])
        ratio = Fraction(q * q * distance * distance, period * period)
        result.update(
            distance_status="EXACT_CLASSICAL_ENUMERATION_PLUS_HGP_THEOREM",
            classical_distance=distance,
            quantum_distance=distance,
            quantum_fom_fraction=str(ratio),
            quantum_fom_decimal=float(ratio),
        )
        assert distance == 3**s - 1
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s", nargs="+", type=int, default=list(range(1, 7)))
    parser.add_argument("--max-enumerated-dimension", type=int, default=16)
    parser.add_argument("--native-q32", action="store_true", help="also enumerate all 2^32-1 messages for s=5 using ARM NEON")
    parser.add_argument("--json", action="store_true", help="emit full certificates")
    args = parser.parse_args()
    if any(s < 1 or s > 10 for s in args.s):
        parser.error("require 1 <= s <= 10 to bound polynomial allocation")
    if not 0 <= args.max_enumerated_dimension <= 24:
        parser.error("enumeration dimension must be between 0 and 24")
    results = [check_seed(s, args.max_enumerated_dimension, args.native_q32) for s in args.s]
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print("s  classical [N,k,d]       BB [[n,k,d]]            k*d^2/n   status")
        for r in results:
            distance = r.get("classical_distance")
            distance_label = str(distance) if distance is not None else f"<={r['proved_upper_witness_weight']}"
            classical = f"[{r['classical_length']},{r['classical_dimension']},{distance_label}]"
            quantum = f"[[{r['quantum_n']},{r['quantum_k']},{distance_label}]]"
            ratio = f"{r['quantum_fom_decimal']:.6f}" if distance is not None else "unknown"
            status = "exact finite enumeration + HGP" if distance is not None else "upper witness only"
            print(f"{r['s']:<2} {classical:<23} {quantum:<24} {ratio:<10} {status}")


if __name__ == "__main__":
    main()
