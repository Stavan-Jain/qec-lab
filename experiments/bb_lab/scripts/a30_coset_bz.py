"""A30 — coset-BZ: certified safe-floor decisions for the d=10 doubling cells.

The A28 two-window certified enumeration (a28_bz_census.py), extended from
linear-code censuses to COSET floors — the I4 species of A27 §2 and exactly
the safe-floor obligation of the doubling template.

Method.  Safe class coset = t0 + C_AB, where C_AB = {(A*f, B*f) : f in
F2[G]} (dim kappa = |G| - dim K) and t0 = (seamC_u, seamC_v) per Prop
A14.1(3) (one kernel G-orbit rep per cell; coset minima are constant on
G-orbits).  Let I1, I2 be DISJOINT information sets for C_AB with
systematic generators G1, G2.  Any coset element c with |c| <= W has
min(|c|_I1, |c|_I2) <= floor(W/2) =: r, and c is the UNIQUE coset element
matching its own restriction c|_Ij; that element is

    c = e_j(S) := c_empty_j  ^  XOR_{i in S} G_j[i],   S = supp(c|_Ij),

with c_empty_j = t0 ^ XOR_{i: t0|_Ij[i]=1} G_j[i] (the coset element
vanishing on I_j).  Enumerating all |S| <= r in both windows therefore
visits every coset element of weight <= W; the certificate is (windows
disjoint + systematic identity blocks + exact node counts
sum_{s<=r} C(kappa, s) + the parity of the coset).  Same species as the
A28 census certificates: completeness is a counting invariant, no UNSAT.

Parity: |A|, |B| odd => |Af|+|Bf| even => coset weights == |t0| (mod 2).
All three floor-20 cells have even t0, so target 20 decides at W = 18.

Stages per code (budget: 15 min of compute per code, enforced):
  basefloor  d_base >= 10 both CSS sides (r = 4, seconds)
  r8         safe floor >= 18 anytime pre-stage (~1 min)
  r9         the decision at 20: certify (0 hits) or refute (witnesses,
             which at r = 9 pin d_safe exactly if <= 18)

Usage:
  uv run python scripts/a30_coset_bz.py validate          # ground-truth battery
  uv run python scripts/a30_coset_bz.py decide            # the 3 floor-20 cells
  uv run python scripts/a30_coset_bz.py decide --only 37a70e02
Outputs: data/a30/*.json (certificates + witnesses), build in data/a30/build.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))

from bb_lab.group import AbelianGroup  # noqa: E402
from bb_lab.poly import Poly  # noqa: E402
from bb_lab.checks import circulant  # noqa: E402
from bb_lab.fibering import kernel_basis, seam_offsets  # noqa: E402

DATA = LAB / "data" / "a30"
BUILD = DATA / "build"
BUDGET_S = 900.0  # 15 min of compute per code (hard)
NTHREADS = 8

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
    """Row-reduce a GF(2) matrix; returns (reduced rows without zeros, pivots)."""
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


def pack3(v: np.ndarray) -> tuple[int, int, int]:
    x = 0
    for i in np.nonzero(v)[0]:
        x |= 1 << int(i)
    return x & ((1 << 64) - 1), (x >> 64) & ((1 << 64) - 1), x >> 128


def unpack3(hexstr: str, n: int) -> np.ndarray:
    x = int(hexstr, 16)
    return np.array([(x >> i) & 1 for i in range(n)], dtype=np.uint8)


# ------------------------------------------------------------- code objects
def code_rows(A: Poly, B: Poly) -> np.ndarray:
    """Generators of C_AB = {(A*f, B*f)}: rows (A*e_g, B*e_g), g in G."""
    MA = circulant(A).astype(np.uint8)
    MB = circulant(B).astype(np.uint8)
    return np.concatenate([MA.T, MB.T], axis=1) % 2  # row g = (A*e_g, B*e_g)


def disjoint_info_sets(
    Mfull: np.ndarray, tries: int = 100, seed: int = 0
) -> tuple[list[int], np.ndarray, list[int], np.ndarray, int]:
    """Two disjoint information sets + systematic generators for rowspace.

    Greedy-leftmost can strand the complement below full rank; retry with
    random column orders until both windows are full rank (disjointness is
    by construction)."""
    R1n, piv0 = rref(Mfull)
    kappa = len(piv0)
    n = Mfull.shape[1]
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
    # systematic on a window: reduce so submatrix at window cols = identity
    def systematic(window: list[int]) -> np.ndarray:
        S, p = rref(np.concatenate([R1n[:, window], R1n], axis=1))
        assert p[: len(window)] == list(range(len(window))), "not systematic"
        return S[:, len(window):]

    return piv1, systematic(piv1), piv2, systematic(piv2), kappa


def coset_base(
    Gsys: np.ndarray, window: list[int], t0: np.ndarray
) -> np.ndarray:
    """The unique coset element t0 + C vanishing on the window."""
    c = t0.copy()
    for i, col in enumerate(window):
        if c[col]:
            c = (c + Gsys[i]) % 2
    assert not c[window].any()
    return c


# ---------------------------------------------------------------- C harness
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


def run_window(
    binp: Path,
    tag: str,
    Gsys: np.ndarray,
    bases: list[np.ndarray],
    r: int,
    W: int,
    deadline: float,
) -> dict:
    """One window pass; returns nodes/hits + raw hit rows; enforces deadline."""
    kappa = Gsys.shape[0]
    mat = BUILD / f"{tag}.mat"
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
        [str(binp), str(mat), str(NTHREADS), str(BUILD / tag)],
        capture_output=True, text=True, timeout=left,
    )
    assert out.returncode == 0, out.stderr
    line = out.stdout.strip()
    parts = dict(p.split("=") for p in line.split())
    nodes, hits = int(parts["nodes"]), int(parts["hits"])
    expect = sum(math.comb(kappa, s) for s in range(1, r + 1))
    assert nodes == expect, f"node count {nodes} != {expect}"
    raw = []
    for t in range(NTHREADS):
        p = BUILD / f"{tag}_t{t:02d}.hits"
        if p.exists():
            raw.extend(p.read_text().split())
            p.unlink()
    # raw alternates "j hex"; regroup
    hit_rows = [(int(raw[i]), raw[i + 1]) for i in range(0, len(raw), 2)]
    assert len(hit_rows) == hits
    return {
        "nodes": nodes, "hits": hits, "expect": expect,
        "wall_s": round(time.monotonic() - t0, 2), "hit_rows": hit_rows,
    }


# ------------------------------------------------------------------- floors
def coset_floor(
    binp: Path,
    tag: str,
    M: np.ndarray,
    offsets: list[tuple[str, np.ndarray]],
    target: int,
    deadline: float,
    prestage: bool = True,
) -> dict:
    """Decide min weight >= target for each labeled coset t0 + rowspace(M)."""
    n = M.shape[1]
    Rn, piv = rref(M)
    I1, G1, I2, G2, kappa = disjoint_info_sets(M)
    assert not (set(I1) & set(I2))
    pars = set()
    labels = [lab for lab, _ in offsets]
    t0s = [t for _, t in offsets]
    for t in t0s:
        assert not in_rowspace(Rn, piv, t), "offset IS in the rowspace — class 0?"
        pars.add(int(t.sum()) % 2)
    W = target - 1
    if len(pars) == 1:
        par = pars.pop()
        if (W - par) % 2 == 1:
            W -= 1  # coset weights == par (mod 2)
    else:
        par = -1  # mixed parities: keep W = target - 1 (sound, no refinement)
    r_full = W // 2
    report: dict = {
        "tag": tag, "n": n, "kappa": kappa, "target": target, "W": W,
        "r": r_full, "parity": par, "I1": I1, "I2": I2,
        "labels": labels, "t0_weights": [int(t.sum()) for t in t0s],
        "stages": [],
    }
    best: dict[str, tuple[int, np.ndarray]] = {}
    empty_w = []
    for wi, (window, Gs) in enumerate([(I1, G1), (I2, G2)]):
        for j, t in enumerate(t0s):
            ce = coset_base(Gs, window, t)
            w = int(ce.sum())
            empty_w.append(w)
            if w <= W and (labels[j] not in best or w < best[labels[j]][0]):
                best[labels[j]] = (w, ce)
    report["empty_weights"] = empty_w
    # Asymmetric complete pair: any |c| <= W has |c|_I1 <= r1 or |c|_I2 <= r2
    # (else |c| >= r1 + r2 + 2 = W + 1).  Schedule cheap window first; an
    # optional anytime insurance stage runs window 0 at depth 8 too.
    r1 = W // 2
    r2 = max(W - r1 - 1, 0)
    report["r_pair"] = [r1, r2]
    schedule: list[tuple[int, int]] = []
    if prestage and r1 > 8:
        schedule.append((1, min(8, r2)))
        schedule.append((0, 8))
        if r2 > 8:
            schedule.append((1, r2))
    else:
        schedule.append((1, r2))
    schedule.append((0, r1))
    windows = [(I1, G1), (I2, G2)]
    done = {0: -1, 1: -1}
    for wi, r in schedule:
        if r <= done[wi]:
            continue
        if all(lab in best for lab in labels):
            break  # every class already refuted; nothing left to decide
        window, Gs = windows[wi]
        bases = [coset_base(Gs, window, t) for t in t0s]
        try:
            res = run_window(
                binp, f"{tag}_r{r}_w{wi}", Gs, bases, r, W, deadline
            )
        except (TimeoutError, subprocess.TimeoutExpired) as e:
            report["aborted"] = f"window {wi} at r={r}: {e.__class__.__name__}"
            break
        for j, hx in res.pop("hit_rows"):
            c = unpack3(hx, n)
            w = int(c.sum())
            assert w <= W
            assert in_rowspace(Rn, piv, (c + t0s[j]) % 2), "hit not in coset"
            if labels[j] not in best or w < best[labels[j]][0]:
                best[labels[j]] = (w, c)
        res["window"] = wi
        res["r"] = r
        report["stages"].append(res)
        done[wi] = max(done[wi], r)
    decided = done[0] >= r1 and done[1] >= r2
    proved = (done[0] + done[1] + 2) if min(done.values()) >= 0 else 0
    report["min_found"] = {lab: best[lab][0] for lab in labels if lab in best}
    per_class = {}
    for lab in labels:
        if lab in best:
            w, c = best[lab]
            per_class[lab] = {
                "verdict": "REFUTED", "min_weight_found": w,
                "witness_hex": f"{int(''.join(map(str, c[::-1])), 2):x}",
                "exact": decided,
            }
        elif decided or proved > W:
            per_class[lab] = {"verdict": "CERTIFIED", "floor": target}
        else:
            per_class[lab] = {
                "verdict": "PARTIAL", "floor": min(target, proved),
            }
    report["per_class"] = per_class
    return report


def base_floor(
    binp: Path, tag: str, A: Poly, B: Poly, d: int, deadline: float
) -> dict:
    """Certify d_base >= d: every logical class coset has min weight >= d.

    X side: logicals = ker D1 / im D2, with D2 = [[MA],[MB]], D1 = [MB, MA]
    (D1 D2 = MB MA + MA MB = 0); Z side: ker D2^T / im D1^T.  Each side is
    2^k - 1 cosets of a rank-(|G| - dim K) stabilizer code — the same
    coset-BZ lane as the safe floors, with r = floor((d-1)/2).
    """
    n = 2 * A.group.cardinality
    out: dict = {"tag": tag, "d": d, "sides": {}}
    MA = circulant(A).astype(np.uint8) % 2
    MB = circulant(B).astype(np.uint8) % 2
    for side in ("X", "Z"):
        if side == "X":  # ker D1 / im D2
            H = np.concatenate([MB, MA], axis=1) % 2
            Mrows = np.concatenate([MA.T, MB.T], axis=1) % 2  # rows span im D2
        else:  # ker D2^T / im D1^T
            H = np.concatenate([MA.T, MB.T], axis=1) % 2
            Mrows = np.concatenate([MB, MA], axis=1) % 2  # rows span im D1^T
        # logical representatives: kernel(H) vectors independent mod Mrows
        Hn, hp = rref(H)
        free = [c for c in range(n) if c not in hp]
        Rcomb, pcomb = rref(Mrows)
        Rcomb = list(Rcomb)
        pcomb = list(pcomb)
        reps: list[np.ndarray] = []
        for c in free:
            v = np.zeros(n, dtype=np.uint8)
            v[c] = 1
            for i, pc in enumerate(hp):
                if Hn[i, c]:
                    v[pc] ^= 1
            assert not ((H @ v) % 2).any()
            w = v.copy()
            for i, pc2 in enumerate(pcomb):
                if w[pc2]:
                    w ^= Rcomb[i]
            if w.any():  # new logical direction; absorb into the RREF
                reps.append(v)
                lead = int(np.nonzero(w)[0][0])
                Rcomb.append(w)
                pcomb.append(lead)
        k_side = len(reps)
        combos = []
        for t in range(1, 1 << k_side):
            L = np.zeros(n, dtype=np.uint8)
            for j in range(k_side):
                if (t >> j) & 1:
                    L ^= reps[j]
            combos.append((f"L{t}", L))
        rep = coset_floor(
            binp, f"{tag}_{side}", Mrows, combos, d, deadline, prestage=False
        )
        mins = [pc["min_weight_found"]
                for pc in rep["per_class"].values()
                if pc["verdict"] == "REFUTED"]
        out["sides"][side] = {
            "k": k_side, "kappa": rep["kappa"], "r": rep["r"], "W": rep["W"],
            "n_classes": len(combos),
            "certified_ge_d": not mins,
            "light_logicals": sorted(mins)[:5],
        }
    out["certified"] = all(s["certified_ge_d"] for s in out["sides"].values())
    return out


# ------------------------------------------------------------------ targets
def poly_pair(group: tuple[int, ...], As: str, Bs: str) -> tuple[Poly, Poly]:
    G = AbelianGroup(group)
    return Poly.from_string(As, G), Poly.from_string(Bs, G)


DECIDE = {
    "37a70e02": {
        "group": (15, 6), "A": "1 + y + x", "B": "y^4 + x + x^11*y^2",
        "axes": [0], "d_base": 10, "floor": 20,
    },
    "5e50a976": {
        "group": (15, 6), "A": "1 + y + x", "B": "y^4 + x^8*y^2 + x^13",
        "axes": [0, 1], "d_base": 10, "floor": 20,
    },
}

VALIDATE = [
    # (tag, group, A, B, axis, target, expect) — recorded ground truths
    ("pair72_x@8", (3, 6), "x^2 + y + y^3", "1 + x + y^2", 0, 8,
     "CERTIFIED"),                                               # A14: minima 8
    ("pair72_x@10", (3, 6), "x^2 + y + y^3", "1 + x + y^2", 0, 10,
     "REFUTED@8"),
    ("f2a6_y@16", (5, 15), "1 + y + x", "x*y^6 + x*y^10 + x^2*y^12", 1, 16,
     "CERTIFIED"),                                               # A23/A29
    ("f2a6_y@18", (5, 15), "1 + y + x", "x*y^6 + x*y^10 + x^2*y^12", 1, 18,
     "REFUTED@16"),                                              # d_safe = 16
    ("ac46bbea_y@16", (5, 15), "1 + y + x", "1 + y^11 + x*y^2", 1, 16,
     "CERTIFIED"),                                               # A29 NEW
    ("ac46bbea_y@18", (5, 15), "1 + y + x", "1 + y^11 + x*y^2", 1, 18,
     "REFUTED@16"),                                              # exact 16
    ("38d3c884_x@16", (5, 15), "1 + y + x", "y^9 + y^12 + x^2*y^4", 0, 16,
     "CERTIFIED"),                                               # A29 NEW
    ("38d3c884_y@18", (5, 15), "1 + y + x", "y^9 + y^12 + x^2*y^4", 1, 18,
     "CERTIFIED"),                                               # A29: >= 18
    ("a8base_x@12", (6, 14), "1 + y + x^3*y^3", "1 + x + x^2*y^7", 0, 12,
     "CERTIFIED"),                                               # A29 §5.2
    ("a8base_x@14", (6, 14), "1 + y + x^3*y^3", "1 + x + x^2*y^7", 0, 14,
     "REFUTED@12"),                                              # min 12
]


def main() -> None:
    global NTHREADS
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["validate", "decide"])
    ap.add_argument("--only", default=None)
    ap.add_argument("--threads", type=int, default=NTHREADS)
    args = ap.parse_args()
    NTHREADS = args.threads
    DATA.mkdir(parents=True, exist_ok=True)
    binp = build_kernel()

    if args.stage == "validate":
        results = []
        for row in VALIDATE:
            tag, group, As, Bs, axis, target, expect = row
            if As is None or Bs is None:
                print(f"{tag}: SKIP (polynomials TBD)")
                continue
            A, B = poly_pair(group, As, Bs)
            offs = seam_offsets(A, B, axis)
            labeled = [(f"rep{i}", np.concatenate([su, sv]))
                       for i, (z, su, sv) in enumerate(offs)]
            deadline = time.monotonic() + 300
            rep = coset_floor(binp, tag, code_rows(A, B), labeled, target,
                              deadline, prestage=False)
            verd = {lab: pc["verdict"] +
                    (f"@{pc.get('min_weight_found')}"
                     if pc["verdict"] == "REFUTED" else "")
                    for lab, pc in rep["per_class"].items()}
            if expect.startswith("CERTIFIED"):
                ok = all(v == "CERTIFIED" for v in verd.values())
            else:  # REFUTED@w: the floor is min-over-classes — one dip refutes
                want = int(expect.split("@")[1])
                mins = [pc["min_weight_found"]
                        for pc in rep["per_class"].values()
                        if pc["verdict"] == "REFUTED"]
                ok = bool(mins) and min(mins) == want
            print(f"{tag}: {verd}  expect={expect}  "
                  f"{'OK' if ok else '** MISMATCH **'}")
            rep["expect"] = expect
            rep["ok"] = ok
            results.append(rep)
        (DATA / "validate.json").write_text(json.dumps(results, indent=1))
        bad = [r["tag"] for r in results if not r["ok"]]
        print(f"validate: {len(results) - len(bad)}/{len(results)} OK"
              + (f"  FAILURES: {bad}" if bad else ""))
        sys.exit(1 if bad else 0)

    # decide
    for cid, spec in DECIDE.items():
        if args.only and not cid.startswith(args.only):
            continue
        A, B = poly_pair(spec["group"], spec["A"], spec["B"])
        code_t0 = time.monotonic()
        deadline = code_t0 + BUDGET_S
        out: dict = {"id": cid, **{k: spec[k] for k in
                                   ("group", "A", "B", "axes", "d_base",
                                    "floor")}}
        out["base_floor"] = base_floor(
            binp, f"{cid}_base", A, B, spec["d_base"], deadline
        )
        labeled = []
        for ax in spec["axes"]:
            for i, (z, su, sv) in enumerate(seam_offsets(A, B, ax)):
                labeled.append(
                    (f"axis{'xy'[ax]}_rep{i}", np.concatenate([su, sv]))
                )
        try:
            out["safe_floor"] = coset_floor(
                binp, f"{cid}_safe", code_rows(A, B), labeled, spec["floor"],
                deadline, prestage=(len(labeled) == 1)
            )
        except (TimeoutError, subprocess.TimeoutExpired) as e:
            out["safe_floor"] = {"aborted": str(e)}
        out["compute_s"] = round(time.monotonic() - code_t0, 1)
        path = DATA / f"decide_{cid}.json"
        path.write_text(json.dumps(out, indent=1))
        print(f"{cid}: base_floor certified={out['base_floor']['certified']} "
              f"safe={json.dumps(out['safe_floor'].get('per_class', out['safe_floor']))} "
              f"compute={out['compute_s']}s -> {path.name}")


if __name__ == "__main__":
    main()
