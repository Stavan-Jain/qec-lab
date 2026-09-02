#!/usr/bin/env python3
"""A42 S5 Stage 3 — machine check of the tower statement at a = 3
(p = 24, the r = 4 member period):

  (i)  pi^* : H(2) -> H(3) is an ISOMORPHISM (A42 §2.13: for a >= 2 the
       pullback along the y-fold Z_{3q} -> Z_{3q/2} is the inclusion
       M[pi^{q/2}] <= M[pi^q] = M) — so every nontrivial class at p = 24
       is the pullback of a class at p = 12, and every nontrivial p = 24
       cycle is  pi^*(v_12) + boundary  with |pi^* v_12| = 2|v_12| >= 48;
  (ii) pi_* : H(3) -> H(2) is ZERO (the fold of every p = 24 cycle is a
       p = 12 boundary);
  (iii) the same at a = 2 -> 1 (iso / zero) and the a = 1 -> 0 exception
       (pi^* injective with image of dim 4, pi_* of rank 2).

Everything is exact window linear algebra on the omega-window engine of
a42_s1_syzygy (Lambda_a = F2[y]/((y^2+y+1)^{2^a}), the class functionals
of §2.2); the pullback R_{a-1} -> R_a is c(y) -> (1 + y^{3q/2}) c(y),
which on the omega-factor is multiplication by the image of 1 + y^{3q/2}
in Lambda_a (= pi^{q/2} x unit); the fold is the ring map y -> y with
y^{3q/2} = 1, i.e. reduction Lambda_a -> Lambda_{a-1} composed with the
sum over the two sheets (on the omega-factor: the quotient map).

Output: data/a42/s5_pullback_a3.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).parent / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SZ = _load("a42_s1_syzygy")
AL = SZ.AL
DATA = LAB / "data" / "a42"


def rank_f2(M):
    return SZ.rank_f2(np.array(M, dtype=np.uint8)) if len(M) else 0


def window_data(a: int, Wx: int):
    ow = SZ.OmegaWindow(a, Wx)
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    return ow, H, T, Z, Lrows, dimZ, ncls


def content(ow, v, blk, c):
    """Lambda-element (int) of a window vector at slot (blk, c)."""
    z = 0
    for i in range(ow.L.dim):
        if v[ow.bit(blk, c, i)]:
            z |= 1 << i
    return z


def vec_from_contents(ow, cont):
    v = np.zeros(ow.nbits, dtype=np.uint8)
    for (blk, c), z in cont.items():
        for i in range(ow.L.dim):
            if (z >> i) & 1:
                v[ow.bit(blk, c, i)] = 1
    return v


def pullback_vec(ow_lo, ow_hi, v):
    """pi^*: multiply every slot content by (1 + y^{3q/2}) in Lambda_hi."""
    q_lo = 1 << ow_lo.L.a          # period 3 q_lo below, 6 q_lo above
    eps = ow_hi.L.red(1 | (1 << (3 * q_lo)))
    cont = {}
    for blk in (0, 1):
        for c in range(ow_lo.Wx):
            z = content(ow_lo, v, blk, c)
            if z:
                cont[(blk, c)] = ow_hi.L.mul(ow_hi.L.red(z), eps)
    return vec_from_contents(ow_hi, cont)


def fold_vec(ow_hi, ow_lo, v):
    """pi_*: reduce every slot content modulo Lambda_lo (the quotient
    F2[y]/(y^{3q}+1) -> F2[y]/(y^{3q/2}+1) restricted to the omega
    factor is reduction mod (y^2+y+1)^{q/2})."""
    cont = {}
    for blk in (0, 1):
        for c in range(ow_hi.Wx):
            z = content(ow_hi, v, blk, c)
            if z:
                cont[(blk, c)] = ow_lo.L.red(z)
    return vec_from_contents(ow_lo, cont)


def classes_of(Lrows, vecs):
    return [(Lrows @ v) % 2 for v in vecs]


def main():
    Wx = 10
    out = {}
    win = {a: window_data(a, Wx) for a in (0, 1, 2, 3)}
    for a in (0, 1, 2, 3):
        ow, H, T, Z, Lrows, dimZ, ncls = win[a]
        print(f"a={a}: Lambda dim {ow.L.dim}, window bits {ow.nbits}, "
              f"dim Z {dimZ}, classes {ncls}", flush=True)
        out[f"a{a}"] = {"dimZ": int(dimZ), "classes": int(ncls)}
    for a in (1, 2, 3):
        ow_lo, H_lo, T_lo, Z_lo, L_lo, _, n_lo = win[a - 1]
        ow_hi, H_hi, T_hi, Z_hi, L_hi, _, n_hi = win[a]
        # class reps below: pick cycles spanning the classes
        reps, sigs = [], []
        for z in Z_lo:
            s = (L_lo @ z) % 2
            if s.any() and rank_f2(sigs + [s]) > len(sigs):
                reps.append(z)
                sigs.append(s)
        assert len(reps) == n_lo
        # (i) pullback: cycles above, classes independent?
        pb = [pullback_vec(ow_lo, ow_hi, z) for z in reps]
        for v in pb:
            assert not ((H_hi @ v) % 2).any(), "pullback not a cycle"
        pb_sigs = [(L_hi @ v) % 2 for v in pb]
        r_pb = rank_f2(pb_sigs)
        # (ii) fold of the classes above: reps above
        reps_hi, sigs_hi = [], []
        for z in Z_hi:
            s = (L_hi @ z) % 2
            if s.any() and rank_f2(sigs_hi + [s]) > len(sigs_hi):
                reps_hi.append(z)
                sigs_hi.append(s)
        assert len(reps_hi) == n_hi
        fd = [fold_vec(ow_hi, ow_lo, z) for z in reps_hi]
        for v in fd:
            assert not ((H_lo @ v) % 2).any(), "fold not a cycle"
        fd_sigs = [(L_lo @ v) % 2 for v in fd]
        r_fd = rank_f2(fd_sigs)
        # pi_* pi^* = 0 and pi^* pi_* = 1 + sigma (sanity)
        comp = [(L_lo @ fold_vec(ow_hi, ow_lo, v)) % 2 for v in pb]
        assert not any(c.any() for c in comp), "pi_* pi^* != 0 on H"
        print(f"a={a-1}->{a}: rank pi^* = {r_pb} (of {n_lo} classes below, "
              f"{n_hi} above); rank pi_* = {r_fd}; pi_* pi^* = 0 OK",
              flush=True)
        out[f"a{a-1}_to_a{a}"] = {"rank_pullback": int(r_pb),
                                  "classes_below": int(n_lo),
                                  "classes_above": int(n_hi),
                                  "rank_pushforward": int(r_fd),
                                  "pullback_iso": bool(r_pb == n_lo == n_hi),
                                  "pushforward_zero": bool(r_fd == 0)}
        if a >= 2:
            assert r_pb == n_lo == n_hi == 6 and r_fd == 0
        else:
            assert r_pb == 4 and n_lo == 4 and n_hi == 6 and r_fd == 2
    (DATA / "s5_pullback_a3.json").write_text(json.dumps(out, indent=1))
    print("->", DATA / "s5_pullback_a3.json")


if __name__ == "__main__":
    main()
