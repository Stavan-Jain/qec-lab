#!/usr/bin/env python
"""A26: cover-detect-and-descend on the mitten codes (arXiv:2607.28795).

The mitten codes are non-abelian lifted-product codes LP(A,B) with 1x2 base
matrices over F2[G].  This tool ports the lab's descent question (descend.py's
"is there a base?" move, previously bivariate-only) to general finite groups:

  1. rebuild every Table XIII mitten code from GAP group tables
     (data/mitten_groups/*.txt, exported by a26_export_groups.g in GAP's
     Elements(G) ordering — the paper's element-index convention);
  2. census translation decks: central involutions (guaranteed decks) plus
     any non-central left/right conjugation-fixing involutions;
  3. for each central involution, build the quotient mitten code on G/<iota>,
     verify the deck and the pushforward intertwining exactly, and estimate
     the quotient's distance by random information-set decoding (ISD upper
     bounds with hit counters — never floors);
  4. classify the pushforward of min-weight cover logicals (dangerous sector:
     projects to a base logical / deck-odd sector: projects to zero);
  5. optionally (--yarn PATH, needs pynauty) cross-check against the paper's
     shipped processor_codes matrices: labeling identification, Tanner-graph
     isomorphism of my builds vs shipped, and quotient-vs-[[150,30,10]].

Findings of record (2026-08-03): see notes/A26_mitten_descent.md.
Run:  uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a26_mitten_descent.py
      [--yarn /path/to/yarn]  [--trials 20000]  [--cover-trials 2500]  [--skip-540]
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np

try:
    import pynauty  # optional: Tanner-graph iso / certificates

    HAVE_NAUTY = True
except ImportError:
    HAVE_NAUTY = False

REPO = Path(__file__).resolve().parents[3]
DATA = REPO / "experiments" / "bb_lab" / "instances" / "mitten_groups"


# ----------------------------------------------------------------- groups


class Group:
    """Finite group given by a multiplication table (0-based indices)."""

    def __init__(self, mt: np.ndarray, structure: str = "?"):
        self.mt = np.asarray(mt, dtype=np.int64)
        self.n = self.mt.shape[0]
        self.structure = structure
        eye = np.arange(self.n)
        e_cands = [i for i in range(self.n) if np.array_equal(self.mt[i], eye)]
        assert len(e_cands) == 1, f"identity not unique: {e_cands}"
        self.e = e_cands[0]
        self.inv = np.empty(self.n, dtype=np.int64)
        for g in range(self.n):
            w = np.nonzero(self.mt[g] == self.e)[0]
            assert len(w) == 1
            self.inv[g] = w[0]
        self.orders = np.empty(self.n, dtype=np.int64)
        for g in range(self.n):
            k, x = 1, g
            while x != self.e:
                x = int(self.mt[x, g])
                k += 1
            self.orders[g] = k

    @classmethod
    def from_file(cls, path) -> "Group":
        # GAP wraps long lines with backslash continuations; undo them.
        lines = Path(path).read_text().replace("\\\n", "").splitlines()
        n = int(lines[0])
        structure = lines[1].strip()
        mt = np.array([list(map(int, l.split())) for l in lines[4 : 4 + n]])
        g = cls(mt, structure)
        orders = np.array(list(map(int, lines[2].split())))
        inv = np.array(list(map(int, lines[3].split())))
        assert np.array_equal(orders, g.orders), "order table mismatch vs GAP"
        assert np.array_equal(inv, g.inv), "inverse table mismatch vs GAP"
        assert g.e == 0, "GAP Elements(G) does not start with identity"
        return g

    def mul(self, a: int, b: int) -> int:
        return int(self.mt[a, b])

    def conj_set(self, t: int, supp) -> tuple:
        ti = int(self.inv[t])
        return tuple(sorted(int(self.mt[self.mt[t, s], ti]) for s in supp))

    def center(self) -> list[int]:
        return [g for g in range(self.n) if np.array_equal(self.mt[g], self.mt[:, g])]

    def involutions(self) -> list[int]:
        return [g for g in range(self.n) if self.orders[g] == 2]

    def central_involutions(self) -> list[int]:
        c = set(self.center())
        return [g for g in self.involutions() if g in c]


# ----------------------------------------------------- ring-rep matrices


def Lmat(G: Group, supp) -> np.ndarray:
    """L(a)[y,x] = 1 iff y = g*x for some g in supp  (L(g)b(h)=b(gh))."""
    M = np.zeros((G.n, G.n), dtype=np.uint8)
    cols = np.arange(G.n)
    for g in supp:
        M[G.mt[g, cols], cols] ^= 1
    return M


def Rmat(G: Group, supp) -> np.ndarray:
    """R(a)[y,x] = 1 iff y = x*g^-1 for some g in supp  (R(g)b(h)=b(hg^-1))."""
    M = np.zeros((G.n, G.n), dtype=np.uint8)
    cols = np.arange(G.n)
    for g in supp:
        gi = int(G.inv[g])
        M[G.mt[cols, gi], cols] ^= 1
    return M


def star(G: Group, supp) -> tuple:
    return tuple(sorted(int(G.inv[g]) for g in supp))


def mitten_code(G: Group, a0, a1, b0, b1):
    """H_X, H_Z of LP([a0 a1],[b0 b1]) per arXiv:2607.28795 Eq. (J1)."""
    Z = np.zeros((G.n, G.n), dtype=np.uint8)
    HX = np.block(
        [
            [Lmat(G, a0), Z, Lmat(G, a1), Z, Rmat(G, star(G, b0))],
            [Z, Lmat(G, a0), Z, Lmat(G, a1), Rmat(G, star(G, b1))],
        ]
    )
    HZ = np.block(
        [
            [Rmat(G, b0), Rmat(G, b1), Z, Z, Lmat(G, star(G, a0))],
            [Z, Z, Rmat(G, b0), Rmat(G, b1), Lmat(G, star(G, a1))],
        ]
    )
    return HX.astype(np.uint8), HZ.astype(np.uint8)


# ----------------------------------------------------------- F2 algebra


def rref_f2(M: np.ndarray):
    M = M.copy().astype(np.uint8)
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        p = r + nz[0]
        if p != r:
            M[[r, p]] = M[[p, r]]
        mask = M[:, c].astype(bool).copy()
        mask[r] = False
        M[mask] ^= M[r]
        piv.append(c)
        r += 1
    return M[:r], piv


def rank_f2(M) -> int:
    return len(rref_f2(M)[1])


def nullspace_f2(M) -> np.ndarray:
    R, piv = rref_f2(M)
    cols = M.shape[1]
    pivset = set(piv)
    free = [c for c in range(cols) if c not in pivset]
    B = np.zeros((len(free), cols), dtype=np.uint8)
    for i, f in enumerate(free):
        B[i, f] = 1
        for r_i, p in enumerate(piv):
            B[i, p] = R[r_i, f]
    if len(B):
        assert not ((M.astype(np.uint8) @ B.T) % 2).any()
    return B


def in_rowspace_f2(rows_rref, piv, v) -> bool:
    v = v.copy().astype(np.uint8)
    for r_i, p in enumerate(piv):
        if v[p]:
            v ^= rows_rref[r_i]
    return not v.any()


# ----------------------------------------------------------------- ISD


class IsdEstimator:
    """Random information-set search; upper bounds only (sQetch Eq. H4 test)."""

    def __init__(self, HX, HZ):
        self.n = HX.shape[1]
        self.Nx = nullspace_f2(HX)
        self.Nz = nullspace_f2(HZ)

    def run(self, trials: int, rng: np.random.Generator):
        n, Nx, Nz = self.n, self.Nx, self.Nz
        best_w, best_v, hits, last = np.inf, None, 0, 0
        for t in range(trials):
            perm = rng.permutation(n)
            R, _ = rref_f2(Nx[:, perm])
            w = R.sum(axis=1)
            for i in np.argsort(w):
                wi = int(w[i])
                if wi > best_w:
                    break
                v = np.zeros(n, dtype=np.uint8)
                v[perm] = R[i]
                if ((Nz @ v) % 2).any():
                    if wi < best_w:
                        best_w, best_v, hits, last = wi, v, 1, t
                    elif wi == best_w and not np.array_equal(v, best_v):
                        hits += 1
                    break
        return int(best_w), best_v, hits, last


def estimate_distance(HX, HZ, trials, seed=0):
    rng = np.random.default_rng(seed)
    ez = IsdEstimator(HX, HZ).run(trials, rng)
    ex = IsdEstimator(HZ, HX).run(trials, rng)
    return {
        "dz_ub": ez[0], "dz_hits": ez[2], "dz_wit": ez[1],
        "dx_ub": ex[0], "dx_hits": ex[2], "dx_wit": ex[1],
        "d_ub": min(ez[0], ex[0]), "trials": trials,
    }


# ------------------------------------------------------ decks/quotients


def translation_deck_census(G: Group, a0, a1, b0, b1):
    A0, A1 = tuple(sorted(a0)), tuple(sorted(a1))
    B0, B1 = tuple(sorted(b0)), tuple(sorted(b1))
    invs = G.involutions()
    left = [t for t in invs if G.conj_set(t, A0) == A0 and G.conj_set(t, A1) == A1]
    right = [s for s in invs
             if G.conj_set(int(G.inv[s]), B0) == B0 and G.conj_set(int(G.inv[s]), B1) == B1]
    return {"central": G.central_involutions(), "left": left, "right": right,
            "n_involutions": len(invs)}


def central_quotient(G: Group, iota: int):
    assert G.orders[iota] == 2 and iota in G.center()
    n = G.n
    partner = G.mt[iota]
    rep = np.minimum(np.arange(n), partner)
    reps = sorted(set(int(r) for r in rep))
    idx = {r: i for i, r in enumerate(reps)}
    proj = np.array([idx[int(rep[x])] for x in range(n)], dtype=np.int64)
    m = n // 2
    qmt = np.zeros((m, m), dtype=np.int64)
    for i, x in enumerate(reps):
        qmt[i] = proj[G.mt[x, reps]]
    return Group(qmt, structure=f"{G.structure}/<iota#{iota}>"), proj


def push_support(proj: np.ndarray, supp):
    c = Counter(int(proj[g]) for g in supp)
    img = tuple(sorted(k for k, v in c.items() if v % 2))
    return img, any(v > 1 for v in c.values())


def pushforward_matrix(proj: np.ndarray, nblocks: int) -> np.ndarray:
    n = len(proj)
    m = int(proj.max()) + 1
    P1 = np.zeros((m, n), dtype=np.uint8)
    P1[proj, np.arange(n)] = 1
    P = np.zeros((nblocks * m, nblocks * n), dtype=np.uint8)
    for b in range(nblocks):
        P[b * m : (b + 1) * m, b * n : (b + 1) * n] = P1
    return P


# ------------------------------------------------------------- nauty


def tanner_certificate(HX, HZ) -> bytes:
    n, mx, mz = HX.shape[1], HX.shape[0], HZ.shape[0]
    adj = {}
    for r in range(mx):
        adj[n + r] = [int(c) for c in np.nonzero(HX[r])[0]]
    for r in range(mz):
        adj[n + mx + r] = [int(c) for c in np.nonzero(HZ[r])[0]]
    g = pynauty.Graph(
        n + mx + mz, directed=False, adjacency_dict=adj,
        vertex_coloring=[set(range(n)), set(range(n, n + mx)),
                         set(range(n + mx, n + mx + mz))],
    )
    return pynauty.certificate(g)


# --------------------------------------------------------------- data

# arXiv:2607.28795 Table XIII (0-based GAP Elements(G) indices).
MITTEN = {
    "150,30,10": dict(gid=(30, 1), a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29)),
    "200,40,12": dict(gid=(40, 5), a0=(10, 21, 29), a1=(0, 17, 18), b0=(2, 27, 38), b1=(0, 19, 21)),
    "300,60,14": dict(gid=(60, 11), a0=(38, 51, 54), a1=(0, 6, 45), b0=(25, 33, 48), b1=(0, 16, 58)),
    "500,100,16": dict(gid=(100, 9), a0=(19, 84, 87), a1=(0, 75, 78), b0=(39, 45, 71), b1=(0, 7, 77)),
    "540,108,18": dict(gid=(108, 9), a0=(20, 35, 52), a1=(0, 36, 39), b0=(38, 63, 104), b1=(0, 35, 94)),
    "630,126,20": dict(gid=(126, 1), a0=(50, 117, 123), a1=(0, 62, 104), b0=(4, 39, 82), b1=(0, 67, 87)),
    "780,156,22": dict(gid=(156, 13), a0=(38, 46, 88), a1=(0, 8, 59), b0=(13, 40, 131), b1=(0, 38, 133)),
    "975,195,24": dict(gid=(195, 1), a0=(112, 123, 135), a1=(0, 104, 185), b0=(52, 56, 132), b1=(0, 62, 75)),
}

# The shipped processor_codes [[300,60,14]] is NOT Table XIII's code: it lives
# in DirectProduct(C10,S3) Elements order (== row-major Kronecker order) with
# these sets, and is Tanner-INEQUIVALENT to the Table XIII build (whose R(b1)
# is singular and whose distance is 6 — Table XIII erratum, see the A26 note).
SHIPPED_300 = dict(group_file="group_dp_c10_s3.txt",
                   sets=((4, 9, 49), (15, 21, 58), (3, 10, 52), (15, 30, 38)))


# --------------------------------------------------------------- driver


def load_group(gid) -> Group:
    return Group.from_file(DATA / f"group_{gid[0]}_{gid[1]}.txt")


def descend(tag, G, sets, trials, cover_trials, published_d=None, seed=1):
    a0, a1, b0, b1 = sets
    HX, HZ = mitten_code(G, a0, a1, b0, b1)
    n = 5 * G.n
    k = n - rank_f2(HX) - rank_f2(HZ)
    print(f"\n=== {tag}: n={n} k={k} group={G.structure!r}")

    cents = G.central_involutions()
    census = translation_deck_census(G, a0, a1, b0, b1)
    extra_l = sorted(set(census["left"]) - set(cents))
    extra_r = sorted(set(census["right"]) - set(cents))
    print(f"    central involutions: {cents}; extra left/right translation decks: {extra_l}/{extra_r}")
    if not cents:
        print("    no central involution — chain stops here.")
        return None

    iota = cents[0]
    pi = G.mt[iota]
    data_perm = np.concatenate([b * G.n + pi for b in range(5)])
    chk_perm = np.concatenate([b * G.n + pi for b in range(2)])
    assert np.array_equal(HX[np.ix_(chk_perm, data_perm)], HX)
    assert np.array_equal(HZ[np.ix_(chk_perm, data_perm)], HZ)
    assert not (pi == np.arange(G.n)).any()
    print(f"    deck verified: iota=#{iota} free, preserves HX/HZ exactly")

    Gq, proj = central_quotient(G, iota)
    imgs, coll = {}, {}
    for nm, s in zip(("a0", "a1", "b0", "b1"), (a0, a1, b0, b1)):
        imgs[nm], coll[nm] = push_support(proj, s)
    print(f"    quotient order {Gq.n}; collisions {coll}")
    print(f"    quotient sets: a0={imgs['a0']} a1={imgs['a1']} b0={imgs['b0']} b1={imgs['b1']}")

    HXq, HZq = mitten_code(Gq, imgs["a0"], imgs["a1"], imgs["b0"], imgs["b1"])
    nq = 5 * Gq.n
    kq = nq - rank_f2(HXq) - rank_f2(HZq)
    p5 = pushforward_matrix(proj, 5)
    p2 = pushforward_matrix(proj, 2)
    assert not ((HXq @ p5) % 2 ^ (p2 @ HX) % 2).any()
    assert not ((HZq @ p5) % 2 ^ (p2 @ HZ) % 2).any()
    print(f"    pushforward intertwining verified; quotient: n={nq} k={kq} (cover k={k})")

    eq = estimate_distance(HXq, HZq, trials, seed=seed)
    print(f"    quotient d_ub={eq['d_ub']} (dx {eq['dx_ub']}/{eq['dx_hits']} hits, "
          f"dz {eq['dz_ub']}/{eq['dz_hits']} hits; {trials} trials/side)")
    ec = estimate_distance(HX, HZ, cover_trials, seed=seed + 1)
    print(f"    cover d_ub={ec['d_ub']} (published {published_d})")

    RZq, pivZq = rref_f2(HZq)
    RXq, pivXq = rref_f2(HXq)
    for lane, wit, Hq, stab in (("Z", ec["dz_wit"], HXq, (RZq, pivZq)),
                                ("X", ec["dx_wit"], HZq, (RXq, pivXq))):
        if wit is None:
            continue
        pw = (p5 @ wit) % 2
        assert not ((Hq @ pw) % 2).any()
        cls = ("ZERO (deck-odd)" if not pw.any()
               else f"base STABILIZER (wt {int(pw.sum())})" if in_rowspace_f2(*stab, pw)
               else f"base LOGICAL (wt {int(pw.sum())})")
        print(f"    min-wt cover {lane}-logical (wt {int(wit.sum())}) pushes to: {cls}")
    return dict(Gq=Gq, imgs=imgs, HXq=HXq, HZq=HZq, kq=kq, est=eq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--yarn", type=Path, default=None,
                    help="path to a checkout of github.com/a7b/yarn (optional cross-checks)")
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--cover-trials", type=int, default=2500)
    ap.add_argument("--skip-540", action="store_true")
    args = ap.parse_args()

    print("== Phase 1: rebuild Table XIII codes, census decks ==")
    built = {}
    for name, spec in MITTEN.items():
        G = load_group(spec["gid"])
        HX, HZ = mitten_code(G, spec["a0"], spec["a1"], spec["b0"], spec["b1"])
        n = HX.shape[1]
        assert not ((HX @ HZ.T) % 2).any()
        k = n - rank_f2(HX) - rank_f2(HZ)
        census = translation_deck_census(G, spec["a0"], spec["a1"], spec["b0"], spec["b1"])
        built[name] = (G, HX, HZ)
        print(f"[[{name}]] {G.structure}: n={n} k={k} central-inv decks {census['central']} "
              f"(involutions {census['n_involutions']})")
        if args.yarn and HAVE_NAUTY and (args.yarn / "processor_codes" / "mitten" / f"[[{name}]]").exists():
            d = args.yarn / "processor_codes" / "mitten" / f"[[{name}]]"
            HXs = np.load(d / "Hx.npy").astype(np.uint8)
            HZs = np.load(d / "Hz.npy").astype(np.uint8)
            exact = np.array_equal(HX, HXs) and np.array_equal(HZ, HZs)
            iso = exact or tanner_certificate(HX, HZ) == tanner_certificate(HXs, HZs)
            print(f"    vs shipped: exact={exact} Tanner-iso={iso}")

    print("\n== Phase 2: descent chains ==")
    spec = MITTEN["200,40,12"]
    r = descend("[[200,40,12]]", load_group(spec["gid"]),
                sets=(spec["a0"], spec["a1"], spec["b0"], spec["b1"]),
                trials=args.trials, cover_trials=args.cover_trials, published_d=12)
    lvl = 2
    while r is not None:
        im = r["imgs"]
        r = descend(f"[[200]] chain level {lvl}", r["Gq"],
                    sets=(im["a0"], im["a1"], im["b0"], im["b1"]),
                    trials=args.trials, cover_trials=args.cover_trials)
        lvl += 1

    spec = MITTEN["300,60,14"]
    r_tab = descend("[[300,60,14]] Table XIII variant (ERRATUM: d=6)", load_group(spec["gid"]),
                    sets=(spec["a0"], spec["a1"], spec["b0"], spec["b1"]),
                    trials=args.trials, cover_trials=args.cover_trials, published_d=14)
    r_ship = descend("[[300,60,14]] shipped variant", Group.from_file(DATA / SHIPPED_300["group_file"]),
                     sets=SHIPPED_300["sets"],
                     trials=args.trials, cover_trials=args.cover_trials, published_d=14)

    if args.yarn and HAVE_NAUTY:
        d150 = args.yarn / "processor_codes" / "mitten" / "[[150,30,10]]"
        c150 = tanner_certificate(np.load(d150 / "Hx.npy").astype(np.uint8),
                                  np.load(d150 / "Hz.npy").astype(np.uint8))
        for nm, rr in (("TableXIII", r_tab), ("shipped", r_ship)):
            if rr is not None:
                same = tanner_certificate(rr["HXq"], rr["HZq"]) == c150
                print(f"[[300]]/iota ({nm}) == published [[150,30,10]]?  {same}")

    if not args.skip_540:
        spec = MITTEN["540,108,18"]
        descend("[[540,108,18]]", load_group(spec["gid"]),
                sets=(spec["a0"], spec["a1"], spec["b0"], spec["b1"]),
                trials=args.trials, cover_trials=max(args.cover_trials, 6000), published_d=18)


if __name__ == "__main__":
    main()
