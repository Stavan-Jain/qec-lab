#!/usr/bin/env python3
"""A32 M4 sweep-driver rehearsal probe (X-side instance beta=0).

Emits a self-contained, import-free Lean file that runs the FULL
even-split sweep for the triple  a0*u + a1*w = t*b0  with real tables,
as ONE native_decide.  Timing THIS is the M4 budget rehearsal; it also
de-risks correctness (kill splits must actually kill, the two classify
splits must fire exactly at the 60 family triples).

Tables (little-endian bit i = GAP element index i):
  tUW[i]  = column i of A1inv @ A0      (w-contribution of u bit i)
  tTW[i]  = column i of A1inv @ Rb0     (w-contribution of t bit i)
  tUT[i]  = column i of Rb0inv @ A0     (t-contribution of u bit i)
  tWT[i]  = column i of Rb0inv @ A1     (t-contribution of w bit i)
  tWC[i]  = column i of A1              (c-contribution of w bit i)
  tTC[i]  = column i of Rb0             (c-contribution of t bit i)
  tP[i]   = column i of P (A0-pseudo-inverse: A0 @ P @ c = c on range)
  ln[2]   = left-null rows of A0 (solvability filter)
  ann[4]  = Ann(A0) = {u : A0 u = 0}  (all four elements, 0 included)
  fam[60] = the two row-fragment families (u,w,t) as mask triples
Everything is hard-asserted in numpy before emission.
Run: cd experiments/bb_lab && uv run python scripts/a32_m4_probe_gen.py <outdir>
"""

from __future__ import annotations

import importlib.util
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GROUPS = HERE.parent / "instances" / "mitten_groups"
SETS = dict(a0=(0, 14, 23), a1=(0, 2, 11), b0=(7, 20, 24), b1=(0, 2, 29))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rank2(M):
    A = M.copy().astype(np.uint8)
    r = 0
    for c in range(A.shape[1]):
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        p = r + nz[0]
        A[[r, p]] = A[[p, r]]
        for q in np.flatnonzero(A[:, c]):
            if q != r:
                A[q] ^= A[r]
        r += 1
    return r


def inv2(M):
    n = M.shape[0]
    A = np.concatenate([M.astype(np.uint8), np.eye(n, dtype=np.uint8)], 1)
    for c in range(n):
        nz = [i for i in range(c, n) if A[i, c]]
        assert nz, "singular"
        A[[c, nz[0]]] = A[[nz[0], c]]
        for q in range(n):
            if q != c and A[q, c]:
                A[q] ^= A[c]
    W = A[:, n:] % 2
    assert ((M @ W) % 2 == np.eye(n, dtype=int)).all()
    assert ((W @ M) % 2 == np.eye(n, dtype=int)).all()
    return W


def pseudo_inv_and_certs(A0):
    """P with A0 P c = c for all c in range(A0); left-null rows; Ann list."""
    n = 30
    # RREF of A0 to get pivot rows/cols and a particular-solution map
    A = A0.copy().astype(np.uint8)
    E = np.eye(n, dtype=np.uint8)
    piv = []
    r = 0
    for c in range(n):
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0:
            continue
        p = r + nz[0]
        A[[r, p]], E[[r, p]] = A[[p, r]], E[[p, r]]
        for q in np.flatnonzero(A[:, c]):
            if q != r:
                A[q] ^= A[r]
                E[q] ^= E[r]
        piv.append(c)
        r += 1
    assert r == 28 and len(piv) == 28
    # E @ A0 = A (RREF).  Particular solution of A0 u = c: u[piv[k]] = (E c)[k].
    P = np.zeros((n, n), dtype=np.uint8)
    for k, pc in enumerate(piv):
        P[pc] = E[k]
    ok = (A0 @ P) % 2
    # left-null rows = last 2 rows of E (A rows 28,29 are zero)
    assert not A[28:].any()
    ln = E[28:] % 2
    # verify: A0 P c = c for c in range (test on all columns of A0)
    for j in range(n):
        c = A0[:, j] % 2
        assert not (ln @ c % 2).any()
        assert ((A0 @ (P @ c % 2)) % 2 == c).all()
    # Ann(A0)
    ker = []
    free = [c for c in range(n) if c not in piv]
    assert len(free) == 2
    for f in free:
        u = np.zeros(n, dtype=np.uint8)
        u[f] = 1
        u ^= P @ (A0[:, f] % 2) % 2
        assert not (A0 @ u % 2).any() and u.any()
        ker.append(u)
    ann = [np.zeros(n, dtype=np.uint8), ker[0], ker[1], ker[0] ^ ker[1]]
    assert len({tuple(a) for a in ann}) == 4
    return P, ln, ann


def mask(v) -> int:
    return int(sum(1 << i for i in np.flatnonzero(v)))


def cols(M) -> list[int]:
    return [mask(M[:, j] % 2) for j in range(M.shape[1])]


def fmt(name: str, vals: list[int], per=6) -> str:
    lines = [f"def {name} : List Nat := ["]
    for i in range(0, len(vals), per):
        chunk = ", ".join(str(v) for v in vals[i:i + per])
        sep = "," if i + per < len(vals) else "]"
        lines.append(f"  {chunk}{sep}")
    return "\n".join(lines)


def main() -> None:
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parent / "data-a32"
    a26 = _load("a26_mitten_descent")
    G = a26.Group.from_file(GROUPS / "group_30_1.txt")
    mul, inv = G.mul, G.inv
    n = 30
    a0s, a1s = set(SETS["a0"]), set(SETS["a1"])
    b0s = set(SETS["b0"])

    A0 = np.array([[int(mul(h, int(inv[x])) in a0s) for x in range(n)]
                   for h in range(n)], dtype=np.uint8)
    A1 = np.array([[int(mul(h, int(inv[x])) in a1s) for x in range(n)]
                   for h in range(n)], dtype=np.uint8)
    Rb0 = np.array([[int(mul(int(inv[x]), h) in b0s) for x in range(n)]
                    for h in range(n)], dtype=np.uint8)
    assert rank2(A0) == 28 and rank2(A1) == 30 and rank2(Rb0) == 30

    A1i, Rb0i = inv2(A1), inv2(Rb0)
    P, ln, ann = pseudo_inv_and_certs(A0)

    tUW = cols(A1i @ A0 % 2)
    tTW = cols(A1i @ Rb0 % 2)
    tUT = cols(Rb0i @ A0 % 2)
    tWT = cols(Rb0i @ A1 % 2)
    tWC = cols(A1)
    tTC = cols(Rb0)
    tP = cols(P)
    lnm = [mask(r) for r in ln]
    annm = [mask(a) for a in ann]

    # classification list: the complete X0 light-triple census from the
    # vectorized collector (a32_m4_collect.py) when available, else the
    # 60 row fragments alone (in which case the simulation below will
    # fire at the first stray split — falsify-first guard).
    triples_json = HERE.parent / "data-a32" / "m4_triples.json"
    if triples_json.exists():
        import json
        fam = [tuple(t) for t in json.loads(triples_json.read_text())["X0"]["sols"]]
        fam = [t for t in fam if t != (0, 0, 0)]
        print(f"[probe-gen] using collector census: {len(fam)} X0 triples")
    else:
        fam = []
        for y in range(n):
            u = sum(1 << mul(y, t) for t in b0s)
            t1 = sum(1 << mul(s, y) for s in a0s)
            t2 = sum(1 << mul(s, y) for s in a1s)
            fam.append((u, 0, t1))
            fam.append((0, u, t2))
    # numpy check: every family member solves the triple
    for (mu, mw, mt) in fam:
        vu = np.array([mu >> i & 1 for i in range(n)], dtype=np.uint8)
        vw = np.array([mw >> i & 1 for i in range(n)], dtype=np.uint8)
        vt = np.array([mt >> i & 1 for i in range(n)], dtype=np.uint8)
        assert not ((A0 @ vu + A1 @ vw + Rb0 @ vt) % 2).any()

    # simulate the whole sweep in numpy-free python ints (bit-for-bit what
    # Lean will do), so a Lean 'false' could only be a Lean-side bug.
    def xor_fold(tbl, m):
        acc = 0
        i = 0
        while m:
            if m & 1:
                acc ^= tbl[i]
            m >>= 1
            i += 1
        return acc

    def pop(m):
        return bin(m).count("1")

    def par(m):
        return pop(m) & 1

    famset = set(fam)
    kmasks = {k: [sum(1 << i for i in c) for c in combinations(range(n), k)]
              for k in range(0, 7)}

    splits = []
    for s in range(0, 9, 2):
        for p in range(s + 1):
            for q in range(s + 1 - p):
                r = s - p - q
                c_ut = len(kmasks.get(p, [])) * len(kmasks.get(r, [])) \
                    if p <= 6 and r <= 6 else 10**18
                c_uw = len(kmasks.get(p, [])) * len(kmasks.get(q, [])) \
                    if p <= 6 and q <= 6 else 10**18
                c_wt = (len(kmasks.get(q, [])) * len(kmasks.get(r, [])) * 4
                        if q <= 6 and r <= 6 else 10**18)
                # weight >6 sides must be derived, never enumerated
                mode, cost = min((("ut", c_ut), ("uw", c_uw), ("wt", c_wt)),
                                 key=lambda x: x[1])
                assert cost < 10**18, (p, q, r)
                splits.append((p, q, r, mode, cost))

    def check(p, q, r, mode):
        if mode == "ut":
            for mu in kmasks[p]:
                wu = xor_fold(tUW, mu)
                for mt in kmasks[r]:
                    mw = wu ^ xor_fold(tTW, mt)
                    if pop(mw) == q and (mu, mw, mt) not in famset \
                            and (mu | mw | mt) != 0:
                        return False
        elif mode == "uw":
            for mu in kmasks[p]:
                tu = xor_fold(tUT, mu)
                for mw in kmasks[q]:
                    mt = tu ^ xor_fold(tWT, mw)
                    if pop(mt) == r and (mu, mw, mt) not in famset \
                            and (mu | mw | mt) != 0:
                        return False
        else:
            for mw in kmasks[q]:
                cw = xor_fold(tWC, mw)
                for mt in kmasks[r]:
                    c = cw ^ xor_fold(tTC, mt)
                    if par(lnm[0] & c) or par(lnm[1] & c):
                        continue
                    up = xor_fold(tP, c)
                    for a in annm:
                        mu = up ^ a
                        if pop(mu) == p and (mu, mw, mt) not in famset \
                                and (mu | mw | mt) != 0:
                            return False
        return True

    import time
    t0 = time.perf_counter()
    total = 0
    for (p, q, r, mode, cost) in splits:
        total += cost
        assert check(p, q, r, mode), f"SWEEP FIRES OFF-FAMILY at {(p,q,r)}!"
    print(f"[probe-gen] python simulation of all {len(splits)} splits PASSED "
          f"({total:,} classes, {time.perf_counter()-t0:.0f}s python)")

    split_lits = ", ".join(
        f"({p}, {q}, {r}, {dict(ut=0, uw=1, wt=2)[m]})"
        for (p, q, r, m, _) in splits)
    packed = sorted((u | (w << 30) | (t << 60)) for (u, w, t) in fam)
    fam_lits = ",\n  ".join(str(v) for v in packed)

    lean = f"""/- A32 M4 rehearsal probe (GENERATED by a32_m4_probe_gen.py; scratch only,
NOT a library module).  Full even-split sweep for the X-side beta=0 triple
a0*u + a1*w = t*b0 over F2[C5xS3], real tables, one native_decide. -/

set_option maxRecDepth 65536

namespace A32M4Probe

{fmt("tUW", tUW)}
{fmt("tTW", tTW)}
{fmt("tUT", tUT)}
{fmt("tWT", tWT)}
{fmt("tWC", tWC)}
{fmt("tTC", tTC)}
{fmt("tP", tP)}
def lnm : List Nat := [{lnm[0]}, {lnm[1]}]
def annm : List Nat := [{annm[0]}, {annm[1]}, {annm[2]}, {annm[3]}]
def fam : List Nat := [
  {fam_lits}]

def splits : List (Nat × Nat × Nat × Nat) := [{split_lits}]

def xorFold : List Nat → Nat → Nat
  | _, 0 => 0
  | [], _ => 0
  | c :: rest, m => (if m &&& 1 = 1 then c else 0) ^^^ xorFold rest (m >>> 1)

/-- Fueled popcount (proof-facing form; masks are < 2^30). -/
def popCntGo : Nat -> Nat -> Nat
  | 0, _ => 0
  | f + 1, m => (m &&& 1) + popCntGo f (m >>> 1)

/-- 15-bit popcount table, correct by construction. -/
def popTbl : Array Nat := (Array.range 32768).map (popCntGo 15)

/-- Driver popcount: two table lookups (= popCntGo 30 on masks < 2^30). -/
def popCnt (m : Nat) : Nat :=
  popTbl.getD (m &&& 32767) 0 + popTbl.getD (m >>> 15) 0

def parity (m : Nat) : Bool := popCnt m % 2 = 1

/-- Pair each mask with its folded table image (hoists the per-pair fold
out of the quadratic loop: the fold is per-ELEMENT, the pair body is
xor + popcount only). -/
def foldedPairs (tbl : List Nat) (ms : List Nat) : List (Nat × Nat) :=
  ms.map fun m => (m, xorFold tbl m)

/-- All n-bit masks of popcount k. -/
def masksOfWt : Nat → Nat → List Nat
  | _, 0 => [0]
  | 0, _ + 1 => []
  | n + 1, k + 1 =>
      masksOfWt n (k + 1) ++ (masksOfWt n k).map (fun m => m ||| (1 <<< n))

def okTriple (mu mw mt : Nat) : Bool :=
  (mu ||| mw ||| mt) == 0
    || fam.contains (mu ||| (mw <<< 30) ||| (mt <<< 60))

def checkSplit : Nat × Nat × Nat × Nat → Bool
  | (p, q, r, 0) =>
      let ups := foldedPairs tUW (masksOfWt 30 p)
      let tps := foldedPairs tTW (masksOfWt 30 r)
      ups.all fun (mu, wu) =>
        tps.all fun (mt, wt) =>
          let mw := wu ^^^ wt
          popCnt mw != q || okTriple mu mw mt
  | (p, q, r, 1) =>
      let ups := foldedPairs tUT (masksOfWt 30 p)
      let wps := foldedPairs tWT (masksOfWt 30 q)
      ups.all fun (mu, tu) =>
        wps.all fun (mw, tw) =>
          let mt := tu ^^^ tw
          popCnt mt != r || okTriple mu mw mt
  | (p, q, r, _) =>
      let wps := foldedPairs tWC (masksOfWt 30 q)
      let tps := foldedPairs tTC (masksOfWt 30 r)
      wps.all fun (mw, cw) =>
        tps.all fun (mt, ct) =>
          let c := cw ^^^ ct
          parity (lnm.getD 0 0 &&& c) || parity (lnm.getD 1 0 &&& c) ||
            (let up := xorFold tP c
             annm.all fun a =>
               let mu := up ^^^ a
               popCnt mu != p || okTriple mu mw mt)

theorem sweep_X0 : splits.all checkSplit = true := by native_decide

end A32M4Probe
"""
    outdir.mkdir(exist_ok=True)
    target = outdir / "A32M4Probe.lean"
    target.write_text(lean)
    print(f"[probe-gen] wrote {target}")


if __name__ == "__main__":
    main()
