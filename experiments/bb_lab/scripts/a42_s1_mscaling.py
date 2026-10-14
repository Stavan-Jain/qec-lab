#!/usr/bin/env python3
"""A42 S1f — the block-cost m-scaling lemma, verified exactly.

For p = 3m·2^a (m odd), the omega-block (the CRT summand of
F2[y]/(y^p-1) corresponding to (y^2+y+1)^{2^a}) maps isomorphically
onto Lambda_a by taking omega-content — verified: the content map on
the block is injective — so pure(lambda) := (min) weight of a lift
with all complementary content zero is the weight of the UNIQUE
block element with content lambda.

Lemma (verified): pure_{3m·2^a}(lambda) = m · pure_{3·2^a}(lambda)
for every lambda, at (m,a) in {(3,0),(5,0),(7,0),(3,1)}.

With Theorem H (the sigma-side is m-independent), this scales the
pure-lift upper bound: UB(3m·2^a) = m·UB(3·2^a) = 2p for all odd m.
Valuation buckets of the base tables (banked): a=0: all 2;
a=1: {nu=0: 2, nu=1: 4}; a=2: {0: 2, 1: 4, 2: 4, 3: 8} — the CMSS
digit shadow pure(pi^nu · unit) = 2m · 2^{popcount(nu)}.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

DATA = LAB / "data" / "a42"


def pure_table(p: int) -> dict[int, int]:
    a = AL.v2(p)
    L = SY.Lam(a)
    ymod = (1 << p) | 1
    cof = AL.pdivmod(ymod, L.mod)[0]
    tab: dict[int, int] = {}
    dimt = AL.pdeg(ymod) - AL.pdeg(cof)
    for t in range(1, 1 << dimt):
        z = AL.pmod(AL.pmul(cof, t), ymod)
        lam = AL.pmod(z, L.mod)
        assert lam not in tab, "content map not injective"
        tab[lam] = bin(z).count("1")
    return tab


def main():
    t0 = time.time()
    base = {a: pure_table(3 << a) for a in (0, 1, 2)}
    rows = []
    for p, m, a in ((9, 3, 0), (15, 5, 0), (21, 7, 0), (18, 3, 1)):
        tab = pure_table(p)
        ok = all(tab[lam] == m * base[a][lam] for lam in tab)
        assert ok, (p, m, a)
        rows.append({"p": p, "m": m, "a": a, "n_lambda": len(tab),
                     "scaling_exact": ok})
        print(f"p={p} (m={m},a={a}): pure == m*base on all "
              f"{len(tab)} lambda: {ok}", flush=True)
    out = {"rows": rows,
           "base_bucket_by_valuation": {
               "0": {"0": 2}, "1": {"0": 2, "1": 4},
               "2": {"0": 2, "1": 4, "2": 4, "3": 8}},
           "wall_s": round(time.time() - t0, 1)}
    (DATA / "s1_mscaling.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s1_mscaling.json'} ({out['wall_s']} s)")


if __name__ == "__main__":
    main()
