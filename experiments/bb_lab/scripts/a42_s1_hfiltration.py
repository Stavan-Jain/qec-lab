#!/usr/bin/env python3
"""A42 S1e — the exact pi-filtration of the cylinder omega-homology
H(a) = Syz/T over Lambda_a[x^pm], a = 0..3.

Computes, per depth a:
  * dim_F2 H(a) (window-stable: checked at two window widths),
  * the pi-action on H and the filtration dims
    F_l = im(pi^l : H -> H)   (NOT the same as classes representable
    with level-l contents, but the module-theoretic filtration),
  * the annihilator exponent (smallest l with pi^l H = 0),
  * elementary-divisor structure of H as an F4[pi]-module
    (dims of ker(pi^l) give the partition).

This decides the stabilization question (is H(a) = H(1) for a >= 1?)
and gives the exact class inventory the forall-a floor must price.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bb_lab.linalg import nullspace_f2  # noqa: E402
import a42_lib as AL  # noqa: E402
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "a42_s1_syzygy", Path(__file__).parent / "a42_s1_syzygy.py")
SY = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SY)

DATA = LAB / "data" / "a42"


def rank_f2(M):
    return SY.rank_f2(M)


def h_structure(a: int, Wx: int, log=print):
    ow = SY.OmegaWindow(a, Wx)
    L = ow.L
    H, T, Z, Lrows, dimZ, ncls = ow.functionals()
    log(f"a={a} Wx={Wx}: dim Syz_W={dimZ} classes={ncls}")
    # class coordinates: cv(v) = Lrows @ v  (Lrows chosen so that on
    # Syz_W the map v -> cv factors through H and is injective on a
    # complement — verify: rank of Lrows@Z^T equals ncls)
    PZ = (Lrows @ Z.T) % 2
    assert rank_f2(PZ) == ncls
    # pi action on window syzygies: multiply each column content by pi
    pi = L.pi

    def pi_vec(v):
        wv = np.zeros_like(v)
        for blk in (0, 1):
            for c in range(ow.Wx):
                z = 0
                for i in range(L.dim):
                    if v[ow.bit(blk, c, i)]:
                        z |= 1 << i
                if z:
                    zz = L.mul(pi, z)
                    for i in range(L.dim):
                        if (zz >> i) & 1:
                            wv[ow.bit(blk, c, i)] = 1
        return wv

    # basis of H: pick Z-rows independent in class coordinates
    picked = []
    basis = []
    for i in range(Z.shape[0]):
        cv = (Lrows @ Z[i]) % 2
        vv = cv.copy()
        for (bv, bp) in basis:
            if vv[bp]:
                vv ^= bv
        nz = np.flatnonzero(vv)
        if nz.size:
            basis.append((vv, nz[0]))
            picked.append(i)
        if len(picked) == ncls:
            break
    Hbasis = Z[picked]
    # pi-action matrix on class coordinates
    cls_of = lambda v: (Lrows @ v) % 2  # noqa: E731
    # express pi*b in the class-coordinate space directly
    Pcols = []
    for b in Hbasis:
        Pcols.append(cls_of(pi_vec(b)))
    Pmat_cc = np.array(Pcols, dtype=np.uint8).T  # (ncls-coord dim) x ncls
    # to get a square action we need coordinates: cls_of has image of
    # rank ncls inside F2^{len(Lrows)}; build a coordinate map via the
    # pivot rows of PZ
    # coordinate extraction: solve cls_of(b_j) basis
    B = np.array([cls_of(b) for b in Hbasis], dtype=np.uint8).T
    # B: len(Lrows) x ncls, full column rank; coords(x) = solve B c = x
    # via rref bookkeeping
    m, n = B.shape
    Aug = np.concatenate([B, np.eye(m, dtype=np.uint8)], axis=1)
    r = 0
    pivots = []
    for cc in range(n):
        piv = None
        for i in range(r, m):
            if Aug[i, cc]:
                piv = i
                break
        if piv is None:
            continue
        Aug[[r, piv]] = Aug[[piv, r]]
        for i in range(m):
            if i != r and Aug[i, cc]:
                Aug[i] ^= Aug[r]
        pivots.append(cc)
        r += 1
    assert r == n, "basis not independent"
    Sol = Aug[:, n:]

    def coords(x):
        y = (Sol @ x) % 2
        return y[:n]

    P = np.array([coords(cls_of(pi_vec(b))) for b in Hbasis],
                 dtype=np.uint8).T  # n x n: action of pi
    # filtration dims
    dims_im = []
    dims_ker = []
    M = np.eye(ncls, dtype=np.uint8)
    for lev in range(0, (1 << a) + 1):
        dims_im.append(int(rank_f2(M.T)))
        K = nullspace_f2(M)
        dims_ker.append(K.shape[0] if K.size else 0)
        M = (P @ M) % 2
    log(f"  dim H = {ncls}; im(pi^l) dims {dims_im}; "
        f"ker(pi^l) dims {dims_ker}")
    return {"a": a, "Wx": Wx, "dimH": ncls, "im_dims": dims_im,
            "ker_dims": dims_ker}


def main():
    t0 = time.time()
    out = {"rows": []}
    for a, Wx in ((0, 8), (0, 10), (1, 8), (1, 11), (2, 9), (2, 12),
                  (3, 10), (3, 13)):
        row = h_structure(a, Wx, log=lambda s: print(s, flush=True))
        out["rows"].append(row)
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s1_hfiltration.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA/'s1_hfiltration.json'} ({out['wall_s']} s)",
          flush=True)


if __name__ == "__main__":
    main()
