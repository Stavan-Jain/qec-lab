"""Certified coset enumeration (two-window Brouwer–Zimmermann), as a library.

The A28/A30 machinery, extracted for reuse: enumerate every element of a
coset t0 + rowspace(M) whose weight is at most W, with completeness
guaranteed by a counting invariant rather than a solver.

Method (A30 §1).  Let I1, I2 be DISJOINT information sets for the code
C = rowspace(M), with systematic generators G1, G2.  A coset element c is
the unique element matching its own restriction c|_Ij; enumerating the
combinations of up to r_j rows of G_j (XOR-seeded at the coset element
vanishing on I_j) visits every c with |c|_{I_j} <= r_j.  With the
asymmetric complete pair r1 + r2 + 2 > W, every |c| <= W is visited.  The
certificate is: windows disjoint + systematic identity blocks + exact node
counts sum_{s<=r_j} C(kappa, s) (asserted per run) + the pair inequality.

The inner walk is a small pthreads C kernel (3x64-bit words, so n <= 192),
built on demand under data/build/cosetbz/.  Session provenance:
scripts/a28_bz_census.py (censuses) and scripts/a30_coset_bz.py (coset
floors); this module is the library-grade copy — the session scripts stay
frozen as artifacts.
"""

from __future__ import annotations

import math
import subprocess
import time
from pathlib import Path
from typing import Callable, Optional

import numpy as np

LAB = Path(__file__).resolve().parents[2]
BUILD = LAB / "data" / "build" / "cosetbz"

NMAX = 192          # 3x64-bit words in the C kernel
NOFF_MAX = 256      # base-word slots in the C kernel

C_SRC = r"""
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>

/* coset-BZ kernel: enumerate XOR-combinations of <= RMAX generator rows
   (3x64-bit words); for each, weight-check against NOFF coset base words.
   matfile: first line "KAPPA RMAX WMAX NOFF"; then KAPPA row triples
   (hex w2 w1 w0); then NOFF base triples. */

static uint64_t R0[4096], R1[4096], R2[4096];
static uint64_t B0[256], B1[256], B2[256];
static int KAPPA, RMAX, WMAX, NOFF, NTHREADS;

typedef struct { int tid; long long nodes; long long hits; FILE *out; } targ_t;

static inline void check(targ_t *ta, uint64_t a0, uint64_t a1, uint64_t a2) {
    for (int j = 0; j < NOFF; j++) {
        uint64_t c0 = a0 ^ B0[j], c1 = a1 ^ B1[j], c2 = a2 ^ B2[j];
        int w = __builtin_popcountll(c0) + __builtin_popcountll(c1)
              + __builtin_popcountll(c2);
        if (w <= WMAX) {
            fprintf(ta->out, "%d %016llx%016llx%016llx\n", j,
                    (unsigned long long)c2, (unsigned long long)c1,
                    (unsigned long long)c0);
            ta->hits++;
        }
    }
}

static void dfs(targ_t *ta, int start, int depth,
                uint64_t a0, uint64_t a1, uint64_t a2) {
    for (int i = start; i < KAPPA; i++) {
        uint64_t b0 = a0 ^ R0[i], b1 = a1 ^ R1[i], b2 = a2 ^ R2[i];
        ta->nodes++;
        check(ta, b0, b1, b2);
        if (depth + 1 < RMAX)
            dfs(ta, i + 1, depth + 1, b0, b1, b2);
    }
}

static void *worker(void *arg) {
    targ_t *ta = (targ_t *)arg;
    for (int i0 = ta->tid; i0 < KAPPA; i0 += NTHREADS) {
        uint64_t b0 = R0[i0], b1 = R1[i0], b2 = R2[i0];
        ta->nodes++;
        check(ta, b0, b1, b2);
        if (RMAX > 1)
            dfs(ta, i0 + 1, 1, b0, b1, b2);
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: cosetbz matfile nthreads tag\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("matfile"); return 2; }
    if (fscanf(f, "%d %d %d %d", &KAPPA, &RMAX, &WMAX, &NOFF) != 4) return 2;
    unsigned long long w2, w1, w0;
    for (int i = 0; i < KAPPA; i++) {
        if (fscanf(f, "%llx %llx %llx", &w2, &w1, &w0) != 3) return 2;
        R2[i] = w2; R1[i] = w1; R0[i] = w0;
    }
    for (int j = 0; j < NOFF; j++) {
        if (fscanf(f, "%llx %llx %llx", &w2, &w1, &w0) != 3) return 2;
        B2[j] = w2; B1[j] = w1; B0[j] = w0;
    }
    fclose(f);
    NTHREADS = atoi(argv[2]);
    pthread_t th[64]; targ_t ta[64];
    for (int t = 0; t < NTHREADS; t++) {
        char name[512];
        snprintf(name, sizeof name, "%s_t%02d.hits", argv[3], t);
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


# ---------------------------------------------------------------- GF2 utils
def rref(M: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Row-reduce over GF(2); returns (nonzero reduced rows, pivot cols)."""
    M = M.copy() % 2
    piv: list[int] = []
    r = 0
    for c in range(M.shape[1]):
        rows = np.nonzero(M[r:, c])[0]
        if len(rows) == 0:
            continue
        M[[r, r + rows[0]]] = M[[r + rows[0], r]]
        for i in np.nonzero(M[:, c])[0]:
            if i != r:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
        if r == M.shape[0]:
            break
    return M[:r], piv


def in_rowspace(R: np.ndarray, piv: list[int], v: np.ndarray) -> bool:
    v = v.copy() % 2
    for i, c in enumerate(piv):
        if v[c]:
            v ^= R[i]
    return not v.any()


def kernel_of(M: np.ndarray) -> np.ndarray:
    """Basis (rows) of the right kernel of M over GF(2)."""
    Mr, piv = rref(M)
    n = M.shape[1]
    free = [c for c in range(n) if c not in set(piv)]
    out = []
    for c in free:
        v = np.zeros(n, dtype=np.uint8)
        v[c] = 1
        for i, pc in enumerate(piv):
            if Mr[i, c]:
                v[pc] ^= 1
        out.append(v)
    return (np.array(out, dtype=np.uint8) if out
            else np.zeros((0, n), dtype=np.uint8))


def pack3(v: np.ndarray) -> tuple[int, int, int]:
    x = 0
    for i in np.nonzero(v)[0]:
        x |= 1 << int(i)
    return x & ((1 << 64) - 1), (x >> 64) & ((1 << 64) - 1), x >> 128


def unpack3(hexstr: str, n: int) -> np.ndarray:
    x = int(hexstr, 16)
    return np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)


# --------------------------------------------------------------- the kernel
def build_kernel() -> Path:
    BUILD.mkdir(parents=True, exist_ok=True)
    src = BUILD / "cosetbz.c"
    binp = BUILD / "cosetbz"
    stale = (not binp.exists() or not src.exists()
             or src.read_text() != C_SRC)
    if stale:
        src.write_text(C_SRC)
        subprocess.run(
            ["cc", "-O3", "-march=native", "-o", str(binp), str(src),
             "-lpthread"],
            check=True,
        )
    return binp


def disjoint_info_sets(
    Mfull: np.ndarray, tries: int = 100, seed: int = 0
) -> tuple[list[int], np.ndarray, list[int], np.ndarray, int]:
    """Two disjoint information sets + systematic generators for rowspace.

    Greedy-leftmost can strand the complement below full rank; retry with
    random column orders until both windows are full rank."""
    R1n, piv0 = rref(Mfull)
    kappa = len(piv0)
    n = Mfull.shape[1]
    if 2 * kappa > n:
        raise ValueError(f"2*kappa = {2*kappa} > n = {n}: no disjoint pair")
    rng = np.random.default_rng(seed)
    piv1: list[int] = []
    piv2: list[int] = []
    for t in range(tries):
        order = (list(range(n)) if t == 0
                 else [int(c) for c in rng.permutation(n)])
        _, p1i = rref(R1n[:, order])
        cand1 = sorted(order[c] for c in p1i)
        rest = [c for c in range(n) if c not in set(cand1)]
        _, p2i = rref(R1n[:, rest])
        if len(p2i) == kappa:
            piv1 = cand1
            piv2 = sorted(rest[c] for c in p2i)
            break
    if not piv2:
        raise RuntimeError(
            f"no disjoint information-set pair found in {tries} tries"
        )

    def systematic(window: list[int]) -> np.ndarray:
        S, p = rref(np.concatenate([R1n[:, window], R1n], axis=1))
        assert p[: len(window)] == list(range(len(window))), "not systematic"
        return S[:, len(window):]

    return piv1, systematic(piv1), piv2, systematic(piv2), kappa


def coset_base(
    Gsys: np.ndarray, window: list[int], t0: np.ndarray
) -> np.ndarray:
    """The unique element of t0 + rowspace vanishing on the window."""
    c = t0.copy()
    for i, col in enumerate(window):
        if c[col]:
            c = (c + Gsys[i]) % 2
    assert not c[window].any()
    return c


def run_window(
    binp: Path,
    tag: str,
    Gsys: np.ndarray,
    bases: list[np.ndarray],
    r: int,
    W: int,
    deadline: float,
    threads: int = 8,
    workdir: Optional[Path] = None,
) -> dict:
    """One window pass; returns nodes/hits + raw hit rows.

    Enforces `deadline` (time.monotonic) via subprocess timeout; asserts
    the exact node count sum_{s<=r} C(kappa, s)."""
    assert len(bases) <= NOFF_MAX
    kappa, n = Gsys.shape[0], Gsys.shape[1]
    assert n <= NMAX
    wd = workdir or BUILD
    wd.mkdir(parents=True, exist_ok=True)
    mat = wd / f"{tag}.mat"
    with mat.open("w") as f:
        f.write(f"{kappa} {r} {W} {len(bases)}\n")
        for row in Gsys:
            w0, w1, w2 = pack3(row)
            f.write(f"{w2:x} {w1:x} {w0:x}\n")
        for b in bases:
            w0, w1, w2 = pack3(b)
            f.write(f"{w2:x} {w1:x} {w0:x}\n")
    left = deadline - time.monotonic()
    if left <= 0:
        raise TimeoutError("budget exhausted before window run")
    t0 = time.monotonic()
    out = subprocess.run(
        [str(binp), str(mat), str(threads), str(wd / tag)],
        capture_output=True, text=True, timeout=left,
    )
    assert out.returncode == 0, out.stderr
    parts = dict(p.split("=") for p in out.stdout.strip().split())
    nodes, hits = int(parts["nodes"]), int(parts["hits"])
    expect = sum(math.comb(kappa, s) for s in range(1, r + 1))
    assert nodes == expect, f"node count {nodes} != {expect}"
    raw = []
    for t in range(threads):
        p = wd / f"{tag}_t{t:02d}.hits"
        if p.exists():
            raw.extend(p.read_text().split())
            p.unlink()
    hit_rows = [(int(raw[i]), raw[i + 1]) for i in range(0, len(raw), 2)]
    assert len(hit_rows) == hits
    return {
        "nodes": nodes, "hits": hits, "expect": expect,
        "wall_s": round(time.monotonic() - t0, 2), "hit_rows": hit_rows,
    }


def pair_radii(W: int) -> tuple[int, int]:
    """The asymmetric complete pair (r1, r2): r1 + r2 + 2 > W."""
    r1 = W // 2
    return r1, max(W - r1 - 1, 0)


def coset_sweep(
    binp: Path,
    tag: str,
    Gsys_pair: tuple[list[int], np.ndarray, list[int], np.ndarray],
    bases_for: Callable[[list[int], np.ndarray], list[np.ndarray]],
    W: int,
    deadline: float,
    threads: int = 8,
    workdir: Optional[Path] = None,
    on_hit: Optional[Callable[[int, np.ndarray], None]] = None,
    n: Optional[int] = None,
) -> dict:
    """Run the complete (r1, r2) pair over both windows.

    `bases_for(window, Gsys)` supplies the per-window coset base words
    (called once per window; the S = empty elements are also reported to
    `on_hit`).  Every hit row (offset index j, vector) goes to `on_hit`.
    Returns per-window stats."""
    I1, G1, I2, G2 = Gsys_pair
    r1, r2 = pair_radii(W)
    stats = []
    for wi, (window, Gs, r) in enumerate([(I1, G1, r1), (I2, G2, r2)]):
        bases = bases_for(window, Gs)
        nn = n if n is not None else Gs.shape[1]
        if on_hit:
            for j, b in enumerate(bases):
                if int(b.sum()) <= W:
                    on_hit(j, b)
        res = run_window(binp, f"{tag}_w{wi}", Gs, bases, r, W, deadline,
                         threads=threads, workdir=workdir)
        rows = res.pop("hit_rows")
        if on_hit:
            for j, hx in rows:
                on_hit(j, unpack3(hx, nn))
        res["window"] = wi
        res["r"] = r
        stats.append(res)
    return {"windows": stats, "r_pair": [r1, r2], "W": W}
