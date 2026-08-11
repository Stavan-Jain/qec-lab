#!/usr/bin/env python3
"""A32 M4 — exhaustive light-triple collector (vectorized split sweeps).

For each of the four triple instances (X0, X1, Z0, Z1) enumerate ALL
solutions with |u|+|w|+|t| <= 9, by sweeping (u, w) supports per split
and deriving t (t is uniquely determined by (u, w) on both sides since
b_beta / a1~ act invertibly ... on the X side t = (a0 u + a1 w) b^-1;
on the Z side D-alpha:  u b0~ + w b1~ = a_alpha~ t  with a0~ SINGULAR:
for alpha = 0 the map t -> a0~ t is 4-to-1 onto its range, so t is
derived as a 4-coset (particular + Ann(a0~)); alpha = 1 unique).

Everything is exact; the result is the complete classification list the
Lean sweeps will certify against.  Families (row fragments) are checked
to be contained in the lists; the JOIN of the two instances of each side
along shared t is simulated and must produce exactly the 60 generator
rows of that side — the M4 floor statement, verified end-to-end offline.

Output: data-a32/m4_triples.json
Run: cd experiments/bb_lab && uv run python scripts/a32_m4_collect.py
"""

from __future__ import annotations

import importlib.util
import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GROUPS = HERE.parent / "instances" / "mitten_groups"
OUT = HERE.parent / "data-a32"
SETS = dict(a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29))
N = 30


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def popcnt(a: np.ndarray) -> np.ndarray:
    if hasattr(np, "bitwise_count"):
        return np.bitwise_count(a)
    # 16-bit table fallback
    tbl = np.array([bin(i).count("1") for i in range(1 << 16)], dtype=np.uint8)
    return (tbl[a & 0xFFFF] + tbl[(a >> 16) & 0x3FFF]).astype(np.uint8)


def col_masks(M: np.ndarray) -> np.ndarray:
    """uint64 mask per column."""
    out = np.zeros(M.shape[1], dtype=np.uint64)
    for j in range(M.shape[1]):
        out[j] = sum(1 << i for i in np.flatnonzero(M[:, j] % 2))
    return out

_SUPP_CACHE: dict = {}


def supp_masks(k: int) -> tuple[np.ndarray, list[tuple]]:
    if k not in _SUPP_CACHE:
        sups = list(combinations(range(N), k))
        _SUPP_CACHE[k] = (np.array([sum(1 << i for i in s) for s in sups],
                                   dtype=np.uint64), sups)
    return _SUPP_CACHE[k]


def fold_cols(cols: np.ndarray, sups: list[tuple]) -> np.ndarray:
    """XOR of cols over each support tuple (vectorized per position)."""
    if not sups:
        return np.zeros(0, dtype=np.uint64)
    k = len(sups[0])
    arr = np.zeros(len(sups), dtype=np.uint64)
    if k == 0:
        return arr
    idx = np.array(sups, dtype=np.int64)
    for pos in range(k):
        arr ^= cols[idx[:, pos]]
    return arr


def sweep_instance(name: str, colU: np.ndarray, colW: np.ndarray,
                   ann_t: list[int] | None) -> list[tuple[int, int, int]]:
    """All (u,w,t) masks with the instance equation and total weight <= 9.

    colU[i] / colW[i]: t-mask contribution of u/w bit i (through the
    already-inverted map); ann_t: additive t-ambiguity (Z0 only) — the
    derived t is particular; solutions are t = part ^ ann for ann in the
    list whose OWN equation-check passes (we verify all 4 and keep those
    that solve; for unique instances ann_t = [0])."""
    sols = []
    anns = np.array(ann_t if ann_t else [0], dtype=np.uint64)
    # in-triple parity (augmentation hom, all four base sets odd):
    # |u|+|w| ≡ |t| (mod 2), so p+q = 9 forces t = 0 with odd parity —
    # impossible; skipping it removes the C(60,9)-scale bulk.
    for p in range(0, 9):
        mu, su = supp_masks(p)
        if mu.size == 0:
            continue
        fu = fold_cols(colU, su)
        for q in range(0, 9 - p):
            mw, sw = supp_masks(q)
            if mw.size == 0:
                continue
            fw = fold_cols(colW, sw)
            budget = 9 - p - q
            # chunk the outer to bound memory
            step = max(1, (1 << 22) // max(1, mw.size))
            for lo in range(0, mu.size, step):
                hi = min(mu.size, lo + step)
                t_part = fu[lo:hi, None] ^ fw[None, :]
                for ann in anns:
                    t = t_part ^ ann
                    w = popcnt(t)
                    ii, jj = np.nonzero(w <= budget)
                    for a, b in zip(ii, jj):
                        sols.append((int(mu[lo + a]), int(mw[b]),
                                     int(t[a, b])))
    return sorted(set(sols))


def main() -> None:
    pg = _load("a32_m4_probe_gen")
    a26 = _load("a26_mitten_descent")
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    mul, inv = G.mul, G.inv
    s = {k: set(v) for k, v in SETS.items()}

    def Lm(a):
        return np.array([[int(mul(h, int(inv[x])) in a) for x in range(N)]
                         for h in range(N)], dtype=np.uint8)

    def Rm(b):
        return np.array([[int(mul(int(inv[x]), h) in b) for x in range(N)]
                         for h in range(N)], dtype=np.uint8)

    def Lms(a):
        return np.array([[int(mul(x, int(inv[y])) in a) for x in range(N)]
                         for y in range(N)], dtype=np.uint8)

    def Rms(b):
        return np.array([[int(mul(int(inv[y]), x) in b) for x in range(N)]
                         for y in range(N)], dtype=np.uint8)

    A0, A1 = Lm(s["a0"]), Lm(s["a1"])
    Rb = {0: Rm(s["b0"]), 1: Rm(s["b1"])}
    B0s, B1s = Rms(s["b0"]), Rms(s["b1"])
    As = {0: Lms(s["a0"]), 1: Lms(s["a1"])}

    inv2 = pg.inv2
    out: dict = {}
    t_all0 = time.perf_counter()

    # ---- X side: t = Rb^-1 (A0 u + A1 w), unique
    for beta in (0, 1):
        Ri = inv2(Rb[beta])
        colU = col_masks(Ri @ A0 % 2)
        colW = col_masks(Ri @ A1 % 2)
        t0 = time.perf_counter()
        sols = sweep_instance(f"X{beta}", colU, colW, None)
        # verify each against the raw equation; collect split histogram
        splits: dict = {}
        for (mu, mw, mt) in sols:
            vu = np.array([mu >> i & 1 for i in range(N)], dtype=np.uint8)
            vw = np.array([mw >> i & 1 for i in range(N)], dtype=np.uint8)
            vt = np.array([mt >> i & 1 for i in range(N)], dtype=np.uint8)
            assert not ((A0 @ vu + A1 @ vw + Rb[beta] @ vt) % 2).any()
            k = f"({int(vu.sum())},{int(vw.sum())},{int(vt.sum())})"
            splits[k] = splits.get(k, 0) + 1
        # families contained?
        fam = set()
        for y in range(N):
            u = sum(1 << mul(y, t) for t in s[f"b{beta}"])
            fam.add((u, 0, sum(1 << mul(t, y) for t in s["a0"])))
            fam.add((0, u, sum(1 << mul(t, y) for t in s["a1"])))
        assert fam <= set(map(tuple, sols)), f"X{beta}: family missing!"
        print(f"[X{beta}] {len(sols)} light triples "
              f"({time.perf_counter() - t0:.0f}s), splits {splits}")
        out[f"X{beta}"] = dict(sols=sols, splits=splits, n=len(sols))

    # ---- Z side: a_alpha~ t = u b0~ + w b1~
    annZ = {1: [0]}
    P0, ln0, ann0 = pg.pseudo_inv_and_certs(As[0])
    annZ[0] = [int(sum(1 << i for i in np.flatnonzero(a))) for a in ann0]
    for alpha in (0, 1):
        if alpha == 1:
            Ai = inv2(As[1])
            colU = col_masks(Ai @ B0s % 2)
            colW = col_masks(Ai @ B1s % 2)
            sols = sweep_instance("Z1", colU, colW, None)
        else:
            # t_part = P0 (u b0~ + w b1~) when solvable; enumerate coset.
            colU = col_masks(P0 @ B0s % 2)
            colW = col_masks(P0 @ B1s % 2)
            raw = sweep_instance("Z0", colU, colW, annZ[0])
            # keep only genuine solutions (solvability filter via check)
            sols = []
            for (mu, mw, mt) in raw:
                vu = np.array([mu >> i & 1 for i in range(N)], dtype=np.uint8)
                vw = np.array([mw >> i & 1 for i in range(N)], dtype=np.uint8)
                vt = np.array([mt >> i & 1 for i in range(N)], dtype=np.uint8)
                if not ((B0s @ vu + B1s @ vw + As[0] @ vt) % 2).any():
                    sols.append((mu, mw, mt))
            sols = sorted(set(sols))
        splits = {}
        for (mu, mw, mt) in sols:
            vu = np.array([mu >> i & 1 for i in range(N)], dtype=np.uint8)
            vw = np.array([mw >> i & 1 for i in range(N)], dtype=np.uint8)
            vt = np.array([mt >> i & 1 for i in range(N)], dtype=np.uint8)
            assert not ((B0s @ vu + B1s @ vw + As[alpha] @ vt) % 2).any()
            k = f"({int(vu.sum())},{int(vw.sum())},{int(vt.sum())})"
            splits[k] = splits.get(k, 0) + 1
        fam = set()
        a_i = {int(inv[t]) for t in s[f"a{alpha}"]}
        b0i = {int(inv[t]) for t in s["b0"]}
        b1i = {int(inv[t]) for t in s["b1"]}
        for h in range(N):
            u = sum(1 << mul(t, h) for t in a_i)
            fam.add((u, 0, sum(1 << mul(h, t) for t in b0i)))
            fam.add((0, u, sum(1 << mul(h, t) for t in b1i)))
        assert fam <= set(map(tuple, sols)), f"Z{alpha}: family missing!"
        print(f"[Z{alpha}] {len(sols)} light triples, splits {splits}")
        out[f"Z{alpha}"] = dict(sols=sols, splits=splits, n=len(sols))

    # ---- JOIN simulation per side: must produce exactly the 60 rows
    for side, (n0, n1) in (("X", ("X0", "X1")), ("Z", ("Z0", "Z1"))):
        rows = set()
        l0, l1 = out[n0]["sols"], out[n1]["sols"]
        by_t: dict = {}
        for (mu, mw, mt) in l1:
            by_t.setdefault(mt, []).append((mu, mw))
        joined = []
        for (mu0, mw0, mt) in l0:
            for (mu1, mw1) in by_t.get(mt, []):
                wtot = (popcnt(np.array([mu0, mw0, mu1, mw1],
                                        dtype=np.uint64)).sum()
                        + popcnt(np.array([mt], dtype=np.uint64)).sum())
                if 0 < wtot <= 9:
                    joined.append((mu0, mw0, mu1, mw1, mt))
        # expected rows
        exp = set()
        if side == "X":
            for y in range(N):
                exp.add((sum(1 << mul(y, t) for t in s["b0"]),
                         0,
                         sum(1 << mul(y, t) for t in s["b1"]),
                         0,
                         sum(1 << mul(t, y) for t in s["a0"])))
                exp.add((0,
                         sum(1 << mul(y, t) for t in s["b0"]),
                         0,
                         sum(1 << mul(y, t) for t in s["b1"]),
                         sum(1 << mul(t, y) for t in s["a1"])))
        else:
            for h in range(N):
                a0i = {int(inv[t]) for t in s["a0"]}
                a1i = {int(inv[t]) for t in s["a1"]}
                b0i = {int(inv[t]) for t in s["b0"]}
                b1i = {int(inv[t]) for t in s["b1"]}
                exp.add((sum(1 << mul(t, h) for t in a0i), 0,
                         sum(1 << mul(t, h) for t in a1i), 0,
                         sum(1 << mul(h, t) for t in b0i)))
                exp.add((0, sum(1 << mul(t, h) for t in a0i),
                         0, sum(1 << mul(t, h) for t in a1i),
                         sum(1 << mul(h, t) for t in b1i)))
        got = set(joined)
        assert got == exp, (f"{side}: join mismatch: "
                            f"{len(got - exp)} extra, {len(exp - got)} missing")
        print(f"[{side}] JOIN = exactly the 60 generator rows  ✓")

    OUT.mkdir(exist_ok=True)
    (OUT / "m4_triples.json").write_text(json.dumps(
        {k: dict(n=v["n"], splits=v["splits"], sols=v["sols"])
         for k, v in out.items()}, indent=1))
    print(f"[done] {time.perf_counter() - t_all0:.0f}s total; "
          f"wrote {OUT / 'm4_triples.json'}")


if __name__ == "__main__":
    main()
