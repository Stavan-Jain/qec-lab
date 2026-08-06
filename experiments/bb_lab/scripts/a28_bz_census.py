"""A28 — Brouwer–Zimmermann-style certified census of light boundaries.

Method (two disjoint information sets; the census specialization of BZ):
let C = C(A,B) (dim kappa), and let I1, I2 be DISJOINT information sets
(full-rank column sets, |I_j| = kappa) with systematic generator matrices
G1, G2 (G_j restricted to I_j = identity).  Every codeword c is the sum of
the G_j-rows selected by its restriction c|_{I_j}.  If |c| <= W then, by
disjointness, min(|c|_{I1}, |c|_{I2}) <= floor(W/2) =: r, so enumerating
all row-combinations of size <= r in BOTH systematic matrices visits every
codeword of weight <= W at least once.  The enumeration is exhaustive
DFS over subsets (C kernel, pthreads); the certificate is (I1, I2 disjoint,
both full rank, systematic forms verified, node counts = sum_{s<=r} C(kappa,s)).

This is an independent, solver-free verification lane for the LSC census:
no SAT, no fibering — plain linear algebra + counting.

Usage:
  uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_bz_census.py <code> [--threads T] [--dry]

Outputs (experiments/bb_lab/data/a28/):
  census_<code>.json      — translation classes with weights + representatives
  bzcert_<code>.json      — the completeness certificate
"""

import argparse
import json
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import DATA, REGISTRY, load_f2a6_census, rref

BUILD = DATA / "a28" / "build"

C_SRC = r"""
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#define NW 3
static uint64_t R0[4096], R1[4096], R2[4096];
static int KAPPA, RMAX, WMAX, NTHREADS;
static const char *TAG;

typedef struct { int tid; long long nodes; long long hits; FILE *out; } targ_t;

static void dfs(targ_t *ta, int start, int depth,
                uint64_t a0, uint64_t a1, uint64_t a2) {
    for (int i = start; i < KAPPA; i++) {
        uint64_t b0 = a0 ^ R0[i], b1 = a1 ^ R1[i], b2 = a2 ^ R2[i];
        ta->nodes++;
        int w = __builtin_popcountll(b0) + __builtin_popcountll(b1)
              + __builtin_popcountll(b2);
        if (w <= WMAX) {
            fprintf(ta->out, "%016llx%016llx%016llx\n",
                    (unsigned long long)b2, (unsigned long long)b1,
                    (unsigned long long)b0);
            ta->hits++;
        }
        if (depth + 1 < RMAX)
            dfs(ta, i + 1, depth + 1, b0, b1, b2);
    }
}

static void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    for (int i0 = ta->tid; i0 < KAPPA; i0 += NTHREADS) {
        uint64_t b0 = R0[i0], b1 = R1[i0], b2 = R2[i0];
        ta->nodes++;
        int w = __builtin_popcountll(b0) + __builtin_popcountll(b1)
              + __builtin_popcountll(b2);
        if (w <= WMAX) {
            fprintf(ta->out, "%016llx%016llx%016llx\n",
                    (unsigned long long)b2, (unsigned long long)b1,
                    (unsigned long long)b0);
            ta->hits++;
        }
        if (RMAX > 1)
            dfs(ta, i0 + 1, 1, b0, b1, b2);
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: bzkern matfile r W nthreads tag\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("matfile"); return 2; }
    KAPPA = 0;
    unsigned long long w2, w1, w0;
    while (fscanf(f, "%llx %llx %llx", &w2, &w1, &w0) == 3) {
        R2[KAPPA] = w2; R1[KAPPA] = w1; R0[KAPPA] = w0; KAPPA++;
    }
    fclose(f);
    RMAX = atoi(argv[2]); WMAX = atoi(argv[3]); NTHREADS = atoi(argv[4]);
    TAG = argv[5];
    pthread_t th[64]; targ_t ta[64];
    for (int t = 0; t < NTHREADS; t++) {
        char name[512];
        snprintf(name, sizeof name, "%s_t%02d.hits", TAG, t);
        ta[t].tid = t; ta[t].nodes = 0; ta[t].hits = 0;
        ta[t].out = fopen(name, "w");
        if (!ta[t].out) { perror("out"); return 2; }
        pthread_create(&th[t], NULL, worker, &ta[t]);
    }
    long long nodes = 0, hits = 0;
    for (int t = 0; t < NTHREADS; t++) {
        pthread_join(th[t], NULL);
        fclose(ta[t].out);
        nodes += ta[t].nodes; hits += ta[t].hits;
    }
    printf("KAPPA=%d nodes=%lld hits=%lld\n", KAPPA, nodes, hits);
    return 0;
}
"""


def systematic_forms(inst, seed=0):
    """Two disjoint information sets with systematic bases.

    A greedy first set can starve the complement (e.g. docket37: taking all
    pivots in the A-block leaves the B-block rank-deficient by dim Ann(B)),
    so try an interleaved column order first, then seeded random shuffles.
    """
    import random

    rows = inst.generator_rows()
    N = inst.G.N
    ncols = 2 * N
    kappa = rank_of = len(rref(rows, ncols)[0])
    orders = [[c for i in range(N) for c in (i, N + i)]]  # u,v interleave
    rng = random.Random(seed)
    for _ in range(50):
        o = list(range(ncols))
        rng.shuffle(o)
        orders.append(o)
    basis1 = piv1 = basis2 = piv2 = None
    for order in orders:
        b1x, p1x = rref(rows, ncols, col_order=order)
        comp = [c for c in range(ncols) if c not in set(p1x)]
        b2x, p2x = rref(rows, ncols, col_order=comp)
        if len(b2x) == kappa:
            basis1, piv1, basis2, piv2 = b1x, p1x, b2x, p2x
            break
    assert basis2 is not None, "no disjoint information-set pair found in 51 tries"
    assert len(basis1) == kappa and len(basis2) == kappa
    assert not (set(piv1) & set(piv2))
    # verify systematicity: basis_j restricted to piv_j is the identity
    for basis, piv in ((basis1, piv1), (basis2, piv2)):
        mask = 0
        for c in piv:
            mask |= 1 << c
        for i, b in enumerate(basis):
            assert b & mask == (1 << piv[i]), "not systematic"
    return (basis1, piv1), (basis2, piv2), kappa


def write_mat(path, basis):
    with open(path, "w") as fh:
        for b in basis:
            w0 = b & ((1 << 64) - 1)
            w1 = (b >> 64) & ((1 << 64) - 1)
            w2 = (b >> 128) & ((1 << 64) - 1)
            fh.write(f"{w2:016x} {w1:016x} {w0:016x}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("code", choices=sorted(REGISTRY))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dry", action="store_true", help="print node counts and exit")
    args = ap.parse_args()

    inst = REGISTRY[args.code]
    ncols = 2 * inst.G.N
    assert ncols <= 192, "C kernel is 3x64-bit words"
    W = inst.W
    r = W // 2

    (b1, p1), (b2, p2), kappa = systematic_forms(inst)
    exp_nodes = sum(math.comb(kappa, s) for s in range(1, r + 1))
    print(f"{args.code}: kappa={kappa} W={W} r={r} "
          f"nodes/set={exp_nodes:,} total={2*exp_nodes:,}")
    if args.dry:
        return

    BUILD.mkdir(parents=True, exist_ok=True)
    src = BUILD / "bzkern.c"
    binp = BUILD / "bzkern"
    src.write_text(C_SRC)
    subprocess.run(["cc", "-O3", "-o", str(binp), str(src), "-lpthread"], check=True)

    t0 = time.time()
    certs = []
    all_hits = set()
    for j, (basis, piv) in enumerate(((b1, p1), (b2, p2)), start=1):
        mat = BUILD / f"{args.code}_G{j}.mat"
        write_mat(mat, basis)
        tag = str(BUILD / f"{args.code}_G{j}")
        res = subprocess.run(
            [str(binp), str(mat), str(r), str(W), str(args.threads), tag],
            capture_output=True, text=True, check=True)
        line = res.stdout.strip()
        print(f"  set {j}: {line}  [{time.time()-t0:.1f}s]")
        parts = dict(kv.split("=") for kv in line.split())
        assert int(parts["KAPPA"]) == kappa
        assert int(parts["nodes"]) == exp_nodes, "node count mismatch — enumeration hole!"
        certs.append({"info_set": piv, "nodes": int(parts["nodes"]),
                      "hits": int(parts["hits"])})
        for t in range(args.threads):
            for line in open(f"{tag}_t{t:02d}.hits"):
                h = int(line.strip(), 16)
                # words were written w2|w1|w0 big-to-small
                all_hits.add(h)

    # canonicalize
    classes = {}
    for h in sorted(all_hits):
        uv = inst.unpack(h)
        c = inst.canonical(uv)
        if c not in classes:
            classes[c] = uv
    by_w = Counter(inst.pair_weight(uv) for uv in classes.values())
    wall = time.time() - t0
    print(f"  raw hits {len(all_hits):,} -> {len(classes)} classes, "
          f"weights {dict(sorted(by_w.items()))}  [{wall:.1f}s]")

    out = {
        "code": args.code, "W": W, "r": r, "kappa": kappa,
        "n_classes": len(classes),
        "weight_histogram": dict(sorted(by_w.items())),
        "classes": [
            {"weight": inst.pair_weight(uv),
             "u_support": inst.G.support(uv[0]),
             "v_support": inst.G.support(uv[1])}
            for uv in sorted(classes.values(),
                             key=lambda uv: (inst.pair_weight(uv), inst.pack(uv)))
        ],
    }
    (DATA / "a28").mkdir(exist_ok=True)
    with open(DATA / "a28" / f"census_{args.code}.json", "w") as fh:
        json.dump(out, fh)
    cert = {
        "code": args.code, "W": W, "r": r, "kappa": kappa,
        "coverage": "min(|c|_I1,|c|_I2) <= r for every |c|<=W (disjoint sets)",
        "expected_nodes_per_set": exp_nodes,
        "sets": certs, "wall_secs": wall,
    }
    with open(DATA / "a28" / f"bzcert_{args.code}.json", "w") as fh:
        json.dump(cert, fh, indent=1)

    # ground-truth comparison where available
    if args.code == "f2a6":
        want = set()
        for row in load_f2a6_census():
            u = inst.G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 0])
            v = inst.G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 1])
            want.add(inst.canonical((u, v)))
        got = set(classes)
        print(f"  vs SAT census: ours {len(got)}, SAT {len(want)}, "
              f"agree={got == want}")
        assert got == want, "CENSUS MISMATCH vs a17 ground truth"
    if args.code == "grossbase":
        assert dict(by_w) == {6: 1, 10: 6}, by_w
        print("  vs (CLASS): 1 hexagon + 6 D-pair classes, nothing else ✓")


if __name__ == "__main__":
    main()
