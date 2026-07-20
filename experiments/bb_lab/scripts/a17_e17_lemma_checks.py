#!/usr/bin/env python3
"""A17 E17 — verification battery for the (1,3)-vacuity and (1,5)
structure lemmas (confirmation of hand proofs, never an ingredient).

V1 (Lemma V bookkeeping): for Sidon 5-sets B (members' A and B sides
    AND random Sidon 5-sets — the lemma quantifies over ALL Sidon B)
    and exhaustive 3-sets T: verify n₂ + 3n₃ = a (= #unordered T-pairs
    with difference in dB), |B·T| = 15 − 2(n₂+n₃) ≥ 15 − 2a ≥ 9, and
    in particular |B·T| ≠ 5 everywhere.

V2 (equality analysis): exhaustive 5-sets T with |B·T| = 5: verify
    dT = dB exactly, multiplicity profile has n₃ = n₄ = n₅ = 0
    (all collisions simple pairs), image = {φ(t) + t} with φ(t) the
    unique non-tail, and record whether φ is injective and whether
    T is a translate of B (completeness census).

V3 (5-clique census, completeness probe): enumerate ALL 5-sets T with
    dT ⊆ dB (cliques in Cay(G, dB) with Sidon-distinct pair
    differences), classified as B-translates / (−B)-translates /
    EXOTIC. Exotics with all multiplicities ≤ 2 would be
    counterexamples to the completeness conjecture — expect none.

Usage:
    uv run python scripts/a17_e17_lemma_checks.py \
        --in data/a17/members_7x9.jsonl --rand 40 --out data/a17/e17_checks.jsonl
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from bb_lab.group import AbelianGroup
from bb_lab.poly import Poly


def diffset_counter(S, G) -> Counter:
    d: Counter = Counter()
    for a in S:
        for b in S:
            if a != b:
                d[G.sub(a, b)] += 1
    return d


def is_sidon(S, G) -> bool:
    d = diffset_counter(S, G)
    return all(v == 1 for v in d.values())


def conv_cells(S, T, G) -> Counter:
    c: Counter = Counter()
    for s in S:
        for t in T:
            c[G.add(s, t)] += 1
    return c


def is_translate(S: frozenset, T: frozenset, G) -> bool:
    t0 = next(iter(T))
    return any(S == frozenset(G.add(x, G.sub(s0, t0)) for x in T)
               for s0 in S)


def check_B(B: tuple, G: AbelianGroup, tag: str, do_v2: bool) -> dict:
    dB = diffset_counter(B, G)
    Bset = frozenset(B)
    negB = frozenset(G.neg(b) for b in B)
    elems = [g for g in G]
    zero = tuple(0 for _ in range(G.rank))
    nonzero = [g for g in elems if g != zero]
    out = {"tag": tag, "B": [list(b) for b in B]}

    # V1: all 3-sets
    v1_bad = 0
    min_img = 99
    for t1, t2 in itertools.combinations(nonzero, 2):
        T = (zero, t1, t2)
        a = sum(1 for x, y in itertools.combinations(T, 2)
                if G.sub(x, y) in dB)
        cells = conv_cells(B, T, G)
        prof = Counter(cells.values())
        n2, n3 = prof.get(2, 0), prof.get(3, 0)
        img = sum(1 for v in cells.values() if v % 2)
        ok = (n2 + 3 * n3 == a) and (img == 15 - 2 * (n2 + n3)) \
            and (img >= 15 - 2 * a) and (img >= 9) and (img != 5)
        if not ok:
            v1_bad += 1
        min_img = min(min_img, img)
    out["V1_bad"] = v1_bad
    out["V1_min_image"] = min_img

    # V3: 5-cliques in Cay(G, dB)
    nbrs = {g: {G.add(g, d) for d in dB} for g in elems}
    cliques = []
    # translation-normalized: 0 in T; extend greedily over dB-neighbors of 0
    cands = sorted(nbrs[zero])
    for combo in itertools.combinations(cands, 4):
        T = (zero, *combo)
        if all(G.sub(x, y) in dB for x, y in itertools.combinations(T, 2)):
            dT = diffset_counter(T, G)
            if all(v == 1 for v in dT.values()):  # Sidon-distinct pairs
                cliques.append(T)
    cls = Counter()
    exotic_m2 = []
    for T in cliques:
        Tset = frozenset(T)
        if is_translate(Tset, Bset, G):
            cls["B_translate"] += 1
        elif is_translate(Tset, negB, G):
            cls["negB_translate"] += 1
        else:
            cls["EXOTIC"] += 1
            cells = conv_cells(B, T, G)
            if max(cells.values()) <= 2:
                exotic_m2.append([list(t) for t in T])
    out["V3_cliques"] = dict(cls)
    out["V3_exotic_multleq2"] = exotic_m2

    # V2: exhaustive 5-sets with image 5 (skippable — slow)
    if do_v2:
        t0 = time.time()
        hits = []
        for rest in itertools.combinations(nonzero, 4):
            T = (zero, *rest)
            cells = conv_cells(B, T, G)
            img = frozenset(g for g, v in cells.items() if v % 2)
            if len(img) != 5:
                continue
            dT = diffset_counter(T, G)
            prof = Counter(cells.values())
            # equality-analysis invariants
            eq_ok = (dict(dT) == dict(dB)) and prof.get(3, 0) == 0 \
                and prof.get(4, 0) == 0 and prof.get(5, 0) == 0
            # phi map: for each t, tails are v(t,t') over t' != t
            phi = {}
            phi_ok = True
            for t in T:
                tails = set()
                for t2 in T:
                    if t2 == t:
                        continue
                    delta = G.sub(t, t2)
                    # unique (u,v) in B with u - v = delta
                    pair = [(u, v) for u in B for v in B
                            if G.sub(u, v) == delta]
                    if len(pair) != 1:
                        phi_ok = False
                        break
                    tails.add(pair[0][1])
                if not phi_ok:
                    break
                non_tails = Bset - tails
                if len(non_tails) != 1:
                    phi_ok = False
                    break
                phi[t] = next(iter(non_tails))
            img_pred = (frozenset(G.add(phi[t], t) for t in T)
                        if phi_ok else None)
            hits.append({
                "T_is_B_translate": is_translate(frozenset(T), Bset, G),
                "eq_invariants": eq_ok,
                "phi_wellformed": phi_ok,
                "phi_injective": phi_ok and len(set(phi.values())) == 5,
                "image_is_phi_formula": img_pred == img,
            })
        out["V2_hits"] = len(hits)
        out["V2_all_ok"] = all(
            h["eq_invariants"] and h["phi_wellformed"]
            and h["image_is_phi_formula"] for h in hits)
        out["V2_all_B_translates"] = all(
            h["T_is_B_translate"] for h in hits)
        out["V2_secs"] = round(time.time() - t0, 1)
    return out


def random_sidon_5set(G: AbelianGroup, rng: random.Random) -> tuple:
    elems = [g for g in G]
    while True:
        S = tuple(rng.sample(elems, 5))
        if is_sidon(S, G):
            return S


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default=None,
                    help="member jsonl (uses both A and B sides)")
    ap.add_argument("--rand", type=int, default=0,
                    help="additionally check N random Sidon 5-sets per frame")
    ap.add_argument("--frames", type=str, default="7x9,6x9",
                    help="frames for the random checks")
    ap.add_argument("--v2", action="store_true",
                    help="run the (slow) exhaustive V2 5-set pass")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20260718)
    args = ap.parse_args()

    jobs: list[tuple[tuple, AbelianGroup, str]] = []
    if args.inp:
        with open(args.inp) as f:
            for i, line in enumerate(f):
                r = json.loads(line)
                if "A" not in r:
                    continue
                G = AbelianGroup(tuple(r["frame"]))
                jobs.append((tuple(sorted(
                    Poly.from_string(r["B"], G).support)), G, f"member{i}B"))
                jobs.append((tuple(sorted(
                    Poly.from_string(r["A"], G).support)), G, f"member{i}A"))
    rng = random.Random(args.seed)
    for tok in args.frames.split(","):
        ell, m = (int(x) for x in tok.split("x"))
        G = AbelianGroup((ell, m))
        for i in range(args.rand):
            jobs.append((random_sidon_5set(G, rng), G, f"rand{tok}_{i}"))

    print(f"{len(jobs)} Sidon 5-sets to check "
          f"(V2={'on' if args.v2 else 'off'})", flush=True)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    tot_bad = 0
    tot_exotic = 0
    with outp.open("w") as out:
        for i, (B, G, tag) in enumerate(jobs, 1):
            res = check_B(B, G, tag, do_v2=args.v2)
            out.write(json.dumps(res) + "\n")
            out.flush()
            tot_bad += res["V1_bad"]
            tot_exotic += len(res["V3_exotic_multleq2"])
            if i % 10 == 0 or res["V1_bad"] or res["V3_exotic_multleq2"]:
                print(json.dumps({"done": i, "tag": tag,
                                  "V1_bad_total": tot_bad,
                                  "min_img": res["V1_min_image"],
                                  "cliques": res["V3_cliques"],
                                  "exotic_m2_total": tot_exotic}),
                      flush=True)
    print(f"\nV1 violations: {tot_bad} (expect 0); "
          f"m<=2 exotic cliques: {tot_exotic} (expect 0)")


if __name__ == "__main__":
    main()
