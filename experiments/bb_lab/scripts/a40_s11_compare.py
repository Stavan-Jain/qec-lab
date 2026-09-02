#!/usr/bin/env python3
"""A40 S11 — the comparison theorem engine (note §16).

Question: at a b = 1 member (6r+6, 6r), can a doubly-spanning
(toroidal) nontrivial X-logical weigh less than the windowed-sector
minimum?  §16.0 pins the statement through u-COMPACTNESS: v is
u-compact (u in Lambda = <(l,0),(0,m)>) if it has an injective compact
lift to the cylinder C_u = Z^2/<u>; the windowed sector W is
(0,m)-compact ∪ (l,0)-compact, D its complement (helical = u-compact
for a twisted u; 2D = no u), and Lemma U gives |v| >= N(u) for every
u-compact nontrivial v (N(u) = the compact floor of C_u).

Lanes
  classify  every banked/recomputed population (bb72 w6, gross <= 12,
            the (12,12) witness family, (18,12) witnesses), each object
            classified by gap sector AND by exact u-compactness over a
            box of directions (windowed lift search on C_u — positive
            verdicts are verified compact cycles, negative verdicts are
            solver UNSAT over a window that provably contains every
            connected compact lift); per-frame sector minima; the
            class-wise probe (H1 = Wx ⊕ Wy?) at every member r <= 3 and
            the b = 0 frames; class minima at gross.
  hunt      SAT witness hunts (existence only — witnesses are upper
            bounds, UNSAT is an observation, never a certificate) for
            gap-dense nontrivial logicals: (18,12) at 24 (the equality
            question), controls (18,12) <= 23 / (12,12) <= 18 / <= 17 /
            gross <= 12, and the conjecture-falsification probe at
            (24,18) <= 35 (time-capped).
  norm      the wrapping norm N(u) for small u and the member helix
            directions through the S5 twisted-atlas reduction (compact
            triviality by the deterministic generator march).

House rules: validate_banked before every lane; every vector
re-verified end to end (cycle, class, weight); RSS guard (current
RSS via ps, <= 2.5 GB); no /tmp; outputs data/a40/s11_*.json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

from bb_lab import cosetbz  # noqa: E402
from bb_lab.tower import (  # noqa: E402
    TowerCode, i2v, v2i, rep_for, validate_banked, rref_ints, kernel_basis,
)
from a38_c37xx_freeze import census_pass  # noqa: E402
from a40_s4_phase_triage import snf2  # noqa: E402

DATA = LAB / "data" / "a40"
A_L = [(0, 0), (0, 1), (3, -1)]
B_L = [(0, 0), (1, 0), (-1, -3)]
# Z-check at g reads block-1 cells g + s (s in supp B) and block-2 cells
# g + t (t in supp A): H_Z = [M_B^T | M_A^T], circulant M[g,h] = P(g-h).
READ = [(0, s) for s in B_L] + [(1, t) for t in A_L]
# X-stabilizer generator at g has block-1 cells g - a (a in supp A) and
# block-2 cells g - b (b in supp B): H_X = [M_A | M_B].
BND = [(0, (-a[0], -a[1])) for a in A_L] + [(1, (-b[0], -b[1])) for b in B_L]
SPAN = 4
RSS_CAP_GB = 2.5


def rss_gb() -> float:
    try:
        out = subprocess.run(["ps", "-o", "rss=", "-p", str(os.getpid())],
                             capture_output=True, text=True).stdout
        return int(out.strip()) / 1e6
    except Exception:
        return 0.0


def rss_guard(tag=""):
    r = rss_gb()
    if r > RSS_CAP_GB:
        raise RuntimeError(f"RSS {r:.2f} GB > cap {RSS_CAP_GB} GB at {tag}")
    return r


def red(supp, lm):
    return frozenset((e[0] % lm[0], e[1] % lm[1]) for e in supp)


_CODES: dict = {}


def member_code(l, m) -> TowerCode:
    key = (l, m)
    if key not in _CODES:
        _CODES[key] = TowerCode(f"tdg({l},{m})", key, red(A_L, key),
                                red(B_L, key))
    return _CODES[key]


def assert_conventions(code: TowerCode):
    """The READ / BND offset tables against the code's own matrices."""
    l, m = code.G.orders
    ng = code.ng
    got = {(int(i) // ng, code.G.from_index(int(i) % ng))
           for i in np.nonzero(code.HZ[0])[0]}
    want = {(blk, ((s[0]) % l, (s[1]) % m)) for blk, s in READ}
    assert got == want, ("READ convention", got, want)
    got = {(int(i) // ng, code.G.from_index(int(i) % ng))
           for i in np.nonzero(code.HX[0])[0]}
    want = {(blk, ((s[0]) % l, (s[1]) % m)) for blk, s in BND}
    assert got == want, ("BND convention", got, want)


# ------------------------------------------------------------ objects
def cells_of(code: TowerCode, v: np.ndarray):
    ng = code.ng
    return frozenset((int(i) // ng,) + tuple(code.G.from_index(int(i) % ng))
                     for i in np.nonzero(v)[0])


def vec_of(code: TowerCode, cells) -> np.ndarray:
    v = np.zeros(code.n, dtype=np.uint8)
    l, m = code.G.orders
    for blk, x, y in cells:
        v[blk * code.ng + code.G.index((x % l, y % m))] ^= 1
    return v


def gap_structure(vals, order):
    vs = sorted(set(vals))
    if not vs:
        return 0, order
    if len(vs) == 1:
        return 1, order - 1
    gaps = [(vs[(i + 1) % len(vs)] - vs[i]) % order - 1
            for i in range(len(vs))]
    return len(vs), max(gaps)


def components(code: TowerCode, cells):
    """Footprint-connected components (two cells adjacent iff some
    Z-check reads both) — the S6 connectivity lemma's adjacency."""
    l, m = code.G.orders
    idx = {c: i for i, c in enumerate(cells)}
    parent = list(range(len(idx)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    diffs = set()
    for b1, s1 in READ:
        for b2, s2 in READ:
            diffs.add((b1, b2, s2[0] - s1[0], s2[1] - s1[1]))
    for (blk, x, y), i in idx.items():
        for b1, b2, dx, dy in diffs:
            if b1 != blk:
                continue
            c2 = (b2, (x + dx) % l, (y + dy) % m)
            j = idx.get(c2)
            if j is not None:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[ra] = rb
    comps: dict = {}
    for c, i in idx.items():
        comps.setdefault(find(i), []).append(c)
    return list(comps.values())


# --------------------------------------------------------- the cover C_u
class Cover:
    """C_u = Z^2/<u> ~= Z x Z_g in coordinates X = w1.e (free),
    Y = w2.e mod g, with W = [[w1],[w2]] unimodular, W u = (0, g)."""

    def __init__(self, u):
        u0, u1 = int(u[0]), int(u[1])
        g = math.gcd(abs(u0), abs(u1))
        assert g > 0
        # bezout a u0 + b u1 = g
        a, b = _bezout(u0, u1)
        assert a * u0 + b * u1 == g
        w1 = (u1 // g, -u0 // g)
        w2 = (a, b)
        det = w1[0] * w2[1] - w1[1] * w2[0]
        assert det in (1, -1)
        if det == -1:
            w2 = (-a, -b)
        assert w1[0] * u0 + w1[1] * u1 == 0
        assert w2[0] * u0 + w2[1] * u1 == g
        self.u, self.g, self.w1, self.w2 = (u0, u1), g, w1, w2
        # inverse of W (det 1): [[w2[1], -w1[1]], [-w2[0], w1[0]]]
        self.Wi = ((w2[1], -w1[1]), (-w2[0], w1[0]))

    def XY(self, e):
        return (self.w1[0] * e[0] + self.w1[1] * e[1],
                (self.w2[0] * e[0] + self.w2[1] * e[1]) % self.g)

    def e_of(self, X, Y):
        """One plane representative of the cover cell (X, Y)."""
        return (self.Wi[0][0] * X + self.Wi[0][1] * Y,
                self.Wi[1][0] * X + self.Wi[1][1] * Y)

    def torus_cell(self, blk, X, Y, lm):
        e = self.e_of(X, Y)
        return (blk, e[0] % lm[0], e[1] % lm[1])


def _bezout(a, b):
    r0, r1, a0, a1, b0, b1 = a, b, 1, 0, 0, 1
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        a0, a1 = a1, a0 - q * a1
        b0, b1 = b1, b0 - q * b1
    if r0 < 0:
        return -a0, -b0
    return a0, b0


def complement_in_lambda(ab, lm):
    """For u = (a l, b m) with gcd(a,b) = 1, a w = (c l, d m) with
    a d - b c = 1, so Lambda = <u, w>."""
    a, b = ab
    if a == 0:
        assert abs(b) == 1
        return (lm[0] * (-b), 0) if b == 1 else (lm[0], 0)
    if b == 0:
        assert abs(a) == 1
        return (0, lm[1] * a)
    # find c, d: a d - b c = 1
    x, y = _bezout(a, -b)   # x a + y (-b) = gcd = 1  -> d = x, c = y
    assert x * a - y * b == 1, (a, b, x, y)
    return (y * lm[0], x * lm[1])


def lift_search(cells, lm, ab, want_lift=True, solver_time=60.0):
    """Exact u-compactness test for the torus object `cells` on
    (l, m), u = (a l, b m).  Returns (lift or None, info).

    The lift is searched among the preimages inside a window of X
    wide enough to contain every connected compact lift (X-step per
    footprint adjacency times the cell count, per component), so
    UNSAT is exhaustive for that window.  A found lift is re-verified
    as a cycle on C_u (every check touching it) and as bijective."""
    from pycryptosat import Solver
    l, m = lm
    a, b = ab
    u = (a * l, b * m)
    cv = Cover(u)
    w = complement_in_lambda(ab, lm)
    Xw = cv.w1[0] * w[0] + cv.w1[1] * w[1]
    assert Xw != 0
    if Xw < 0:
        w = (-w[0], -w[1])
        Xw = -Xw
    # footprint X-steps
    steps = []
    for b1, s1 in READ:
        for b2, s2 in READ:
            d = (s2[0] - s1[0], s2[1] - s1[1])
            steps.append(abs(cv.w1[0] * d[0] + cv.w1[1] * d[1]))
    maxstep = max(steps)
    ncell = len(cells)
    R = maxstep * ncell + 2 * maxstep + 2
    Xlo, Xhi = -R, Xw + R
    # candidate lifts per torus cell
    cells = sorted(cells)
    var = {}
    per_cell = []
    for ci, (blk, x, y) in enumerate(cells):
        opts = []
        # e = (x, y) + j w; X = X(x,y) + j Xw
        X0, Y0 = cv.XY((x, y))
        jlo = math.floor((Xlo - X0) / Xw) - 1
        jhi = math.ceil((Xhi - X0) / Xw) + 1
        for j in range(jlo, jhi + 1):
            e = (x + j * w[0], y + j * w[1])
            X, Y = cv.XY(e)
            if Xlo <= X <= Xhi:
                vid = len(var) + 1
                var[vid] = (ci, blk, X, Y)
                opts.append(vid)
        per_cell.append(opts)
    s = Solver()
    for opts in per_cell:
        s.add_clause(opts)
        for i in range(len(opts)):
            for j in range(i + 1, len(opts)):
                s.add_clause([-opts[i], -opts[j]])
    # checks: group vars by the checks reading them
    by_check: dict = {}
    for vid, (ci, blk, X, Y) in var.items():
        for rb, sft in READ:
            if rb != blk:
                continue
            # a check at (X', Y') reads blk cell (X' + X(sft), Y' + Y(sft))
            dX = cv.w1[0] * sft[0] + cv.w1[1] * sft[1]
            dY = cv.w2[0] * sft[0] + cv.w2[1] * sft[1]
            key = (X - dX, (Y - dY) % cv.g)
            by_check.setdefault(key, []).append(vid)
    for key, vids in by_check.items():
        s.add_xor_clause(vids, False)
    sat, sol = s.solve()
    info = dict(u=list(u), g=cv.g, nvars=len(var), nchecks=len(by_check),
                window=[Xlo, Xhi], Xw=Xw)
    if not sat:
        return None, info
    lift = [var[vid] for vid in var if sol[vid]]
    # verify: bijective onto cells, and a cycle on C_u
    assert sorted(ci for ci, *_ in lift) == list(range(ncell))
    lifted = {(blk, X, Y) for _, blk, X, Y in lift}
    assert len(lifted) == ncell
    touched: dict = {}
    for blk, X, Y in lifted:
        for rb, sft in READ:
            if rb != blk:
                continue
            dX = cv.w1[0] * sft[0] + cv.w1[1] * sft[1]
            dY = cv.w2[0] * sft[0] + cv.w2[1] * sft[1]
            key = (X - dX, (Y - dY) % cv.g)
            touched[key] = touched.get(key, 0) ^ 1
    assert not any(touched.values()), "lift is not a cycle on C_u"
    # every lifted cell projects to its own torus cell
    for ci, blk, X, Y in lift:
        assert cv.torus_cell(blk, X, Y, lm) == cells[ci]
    Xs = [X for _, _, X, _ in lift]
    info["X_extent"] = max(Xs) - min(Xs) + 1
    return sorted(lifted), info


DIRS = [(0, 1), (1, 0), (1, 1), (1, -1), (1, 2), (1, -2), (2, 1), (2, -1),
        (1, 3), (1, -3), (3, 1), (3, -1), (2, 3), (2, -3), (3, 2), (3, -2)]


def classify_object(code: TowerCode, v: np.ndarray, dirs=DIRS):
    lm = tuple(code.G.orders)
    cells = cells_of(code, v)
    xs = [c[1] for c in cells]
    ys = [c[2] for c in cells]
    nx, gx = gap_structure(xs, lm[0])
    ny, gy = gap_structure(ys, lm[1])
    comps = components(code, cells)
    rec = dict(w=int(v.sum()), nx=nx, ny=ny, gap_x=gx, gap_y=gy,
               ncomp=len(comps),
               gap_sector=("both-gap" if gx >= SPAN and gy >= SPAN else
                           "x-gap" if gx >= SPAN else
                           "y-gap" if gy >= SPAN else "gap-dense"))
    compact = []
    for ab in dirs:
        lift, info = lift_search(cells, lm, ab)
        if lift is not None:
            compact.append(dict(ab=list(ab), u=info["u"], g=info["g"],
                                X_extent=info["X_extent"]))
    rec["compact_dirs"] = compact
    if any(c["ab"] == [0, 1] for c in compact):
        sec = "W_x"
    elif any(c["ab"] == [1, 0] for c in compact):
        sec = "W_y"
    elif compact:
        sec = "helical"
    else:
        sec = "2D"
    rec["sector"] = sec
    return rec


# ------------------------------------------------------- populations
def census_all_classes(binp, code: TowerCode, W: int, tag: str):
    allb = list(range(1, 1 << code.k))
    CH = 51
    out = []
    for lo in range(0, len(allb), CH):
        chunk = allb[lo:lo + CH]
        hits = census_pass(binp, code,
                           [(f"C{c}", rep_for(code, c)) for c in chunk],
                           W, f"{tag}_{lo}")
        for c in chunk:
            for h in sorted(hits[f"C{c}"]):
                v = i2v(h, code.n)
                assert code.is_cycle(v) and not code.is_stab(v)
                assert v2i(code.sig(v)) == c
                out.append((c, v))
    return out


def a36_witness(tg: TowerCode) -> np.ndarray:
    wit = json.loads((LAB / "data" / "a36" /
                      "w18_witness_banked.json").read_text())
    tg_s = TowerCode("tg/stored", (12, 12), "x^3 + y^2 + y^7",
                     "y^3 + x + x^2")
    v_s = np.zeros(tg_s.n, dtype=np.uint8)
    v_s[wit["v_support"]] = 1
    ngt = tg.ng
    v_p = np.zeros(tg.n, dtype=np.uint8)
    for i in np.nonzero(v_s)[0]:
        blk, gi = divmod(int(i), ngt)
        h = tg_s.G.from_index(gi)
        s = (0, 7) if blk == 0 else (1, 0)
        uu = ((h[0] + s[0]) % 12, (7 * (h[1] + s[1])) % 12)
        v_p[blk * ngt + tg.G.index(uu)] = 1
    assert tg.is_cycle(v_p) and not tg.is_stab(v_p) and v_p.sum() == 18
    return v_p


def shear_pullbacks(binp, l, p, d, W, target: TowerCode):
    """All nontrivial cycles of weight <= W on the shear frame
    Z^2/<(l,0),(d,p)>, pulled back to the target torus (which must
    cover the frame); each pullback re-verified."""
    from a40_s4_phase_triage import quotient_code
    code, o = quotient_code(l, p, d)
    M = [[l, 0], [d, p]]
    D, U, V = snf2(M)
    o1, o2 = abs(D[0][0]), abs(D[1][1])
    assert (o1, o2) == tuple(o)
    L, Mm = target.G.orders

    def phi(e):
        return ((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                (e[0] * V[0][1] + e[1] * V[1][1]) % o2)
    # covering condition: (L,0),(0,Mm) in <(l,0),(d,p)>
    assert phi((L, 0)) == (0, 0) and phi((0, Mm)) == (0, 0), \
        "target does not cover the shear frame"
    objs = census_all_classes(binp, code, W, f"s11_shear_{l}_{p}_{d}")
    out = []
    for c, v in objs:
        vt = np.zeros(target.n, dtype=np.uint8)
        for blk in range(2):
            for x in range(L):
                for y in range(Mm):
                    f = phi((x, y))
                    if v[blk * code.ng + code.G.index(f)]:
                        vt[blk * target.ng + target.G.index((x, y))] = 1
        assert target.is_cycle(vt)
        if target.is_stab(vt):
            continue
        out.append((c, vt))
    return out, code.k


def l12_stack(t21: TowerCode) -> np.ndarray:
    L12 = [tuple(t) for t in
           json.loads((DATA / "s2_ub_bands.json").read_text())["L12"]]
    v = np.zeros(t21.n, dtype=np.uint8)
    for j in range(2):
        for blk, gx, gy in L12:
            v[blk * t21.ng + t21.G.index((gx, (gy + 6 * j) % 12))] ^= 1
    assert t21.is_cycle(v) and not t21.is_stab(v) and v.sum() == 24
    return v


# ------------------------------------------------ class decomposition
def windowed_class_span(code: TowerCode, axis: int):
    """Span (as sig-ints) of the classes of cycles supported in a
    window of width order - 4 along `axis`."""
    l, m = code.G.orders
    order = code.G.orders[axis]
    keep = []
    for blk in range(2):
        for x in range(l):
            for y in range(m):
                c = (x, y)[axis]
                if c <= order - 5:
                    keep.append(blk * code.ng + code.G.index((x, y)))
    keep = np.array(keep)
    Hsub = code.HZ[:, keep]
    kb = kernel_basis(Hsub)
    sigs = []
    for kv in kb:
        v = np.zeros(code.n, dtype=np.uint8)
        v[keep] = kv
        assert code.is_cycle(v)
        sigs.append(v2i(code.sig(v)))
    basis, _ = rref_ints(sigs)
    return basis


def decomposition(code: TowerCode) -> dict:
    Wx = windowed_class_span(code, 0)
    Wy = windowed_class_span(code, 1)
    both, _ = rref_ints(list(Wx) + list(Wy))
    dim_sum = len(both)
    dim_int = len(Wx) + len(Wy) - dim_sum
    return dict(k=code.k, dim_Wx=len(Wx), dim_Wy=len(Wy), dim_sum=dim_sum,
                dim_int=dim_int, direct_sum=bool(dim_int == 0 and
                                                dim_sum == code.k),
                Wx_basis=[int(x) for x in Wx], Wy_basis=[int(x) for x in Wy])


def class_split(sig_int: int, Wx_basis, Wy_basis):
    """Write a class as cx + cy (unique when the sum is direct)."""
    rows = list(Wx_basis) + list(Wy_basis)
    nb = len(Wx_basis)
    # solve over F2: sig = sum eps_i rows_i, by elimination with tracking
    basis, piv, tags = [], [], []
    for i, r in enumerate(rows):
        cur, tag = r, 1 << i
        for bb, pp, tt in zip(basis, piv, tags):
            if (cur >> pp) & 1:
                cur ^= bb
                tag ^= tt
        if cur:
            basis.append(cur)
            piv.append((cur & -cur).bit_length() - 1)
            tags.append(tag)
    cur, tag = sig_int, 0
    for bb, pp, tt in zip(basis, piv, tags):
        if (cur >> pp) & 1:
            cur ^= bb
            tag ^= tt
    assert cur == 0, "class outside Wx + Wy"
    cx = 0
    cy = 0
    for i in range(len(rows)):
        if (tag >> i) & 1:
            if i < nb:
                cx ^= rows[i]
            else:
                cy ^= rows[i]
    assert cx ^ cy == sig_int
    return cx, cy


# ---------------------------------------------------------------- lanes
def lane_classify(args):
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    binp = cosetbz.build_kernel()
    out: dict = {"frames": {}, "wall": {}}

    def frame_report(name, code, objs, extra=None):
        t1 = time.time()
        rows = []
        for c, v in objs:
            rec = classify_object(code, v)
            rec["class"] = int(c)
            rows.append(rec)
        # per-sector minima and counts
        sec = {}
        for rec in rows:
            s = rec["sector"]
            d = sec.setdefault(s, dict(n=0, min_w=None, gap_sectors={}))
            d["n"] += 1
            d["min_w"] = rec["w"] if d["min_w"] is None else \
                min(d["min_w"], rec["w"])
            d["gap_sectors"][rec["gap_sector"]] = \
                d["gap_sectors"].get(rec["gap_sector"], 0) + 1
        dirs = {}
        for rec in rows:
            for cd in rec["compact_dirs"]:
                key = str(tuple(cd["ab"]))
                dirs[key] = dirs.get(key, 0) + 1
        dec = decomposition(code)
        rep = dict(lm=list(code.G.orders), k=code.k, n_objects=len(rows),
                   sectors=sec, compact_dir_counts=dirs,
                   decomposition={k: v for k, v in dec.items()
                                  if not k.endswith("basis")},
                   wall_s=round(time.time() - t1, 1))
        if extra:
            rep.update(extra)
        # class-wise: split each object's class
        if dec["direct_sum"]:
            split_hist = {}
            for rec in rows:
                cx, cy = class_split(rec["class"], dec["Wx_basis"],
                                     dec["Wy_basis"])
                kind = ("x-only" if cy == 0 else "y-only" if cx == 0
                        else "mixed")
                rec["class_kind"] = kind
                key = f"{kind}|{rec['sector']}|w{rec['w']}"
                split_hist[key] = split_hist.get(key, 0) + 1
            rep["class_kind_x_sector_x_weight"] = dict(sorted(
                split_hist.items()))
        rep["objects"] = rows if args.keep else \
            [r for r in rows if r["sector"] != "W_x"][:200]
        out["frames"][name] = rep
        print(f"[{name}] {len(rows)} objects, k={code.k}, sectors="
              f"{ {k: (v['n'], v['min_w']) for k, v in sec.items()} }, "
              f"dirs={dirs}, decomposition={rep['decomposition']}, "
              f"{rep['wall_s']} s", flush=True)
        return rows

    # --- b = 0 control frames
    bb72 = member_code(6, 6)
    assert_conventions(member_code(12, 12))
    objs = census_all_classes(binp, bb72, 6, "s11_bb72")
    frame_report("bb72(6,6)", bb72, objs)
    rss_guard("bb72")

    tg = member_code(12, 12)
    wit = a36_witness(tg)
    objs = [(v2i(tg.sig(wit)), wit)]
    pb, k44 = shear_pullbacks(binp, 12, 4, 4, 6, tg)
    # dedupe against the witness
    seen = {v2i(wit)}
    for c, v in pb:
        if v2i(v) not in seen:
            seen.add(v2i(v))
            objs.append((v2i(tg.sig(v)), v))
    frame_report("two-gross(12,12)", tg, objs,
                 extra=dict(source="a36 witness + (12;4,4)-w6 shear "
                                   "pullbacks", shear_k=k44,
                            n_pullbacks=len(pb)))
    rss_guard("tg")

    # --- b = 1 members
    gross = member_code(12, 6)
    objs = census_all_classes(binp, gross, 12, "s11_gross")
    rows = frame_report("gross(12,6)", gross, objs)
    # class minima at gross: per class, min weight and the sectors at it
    per_class: dict = {}
    for rec in rows:
        d = per_class.setdefault(rec["class"], dict(min_w=None, secs=set()))
        if d["min_w"] is None or rec["w"] < d["min_w"]:
            d["min_w"], d["secs"] = rec["w"], set()
        if rec["w"] == d["min_w"]:
            d["secs"].add(rec["sector"])
    hist = {}
    for c, d in per_class.items():
        key = f"w{d['min_w']}|" + "+".join(sorted(d["secs"]))
        hist[key] = hist.get(key, 0) + 1
    out["frames"]["gross(12,6)"]["class_minima_hist"] = dict(sorted(
        hist.items()))
    out["frames"]["gross(12,6)"]["n_classes_with_objects"] = len(per_class)
    print(f"[gross] class minima: {hist}", flush=True)
    rss_guard("gross")

    t21 = member_code(18, 12)
    objs = [(v2i(t21.sig(l12_stack(t21))), l12_stack(t21))]
    hf = DATA / "s11_hunt.json"
    if hf.exists():
        H = json.loads(hf.read_text())
        for run in H.get("runs", []):
            if tuple(run["lm"]) == (18, 12):
                for sol in run.get("solutions", []):
                    v = vec_of(t21, [tuple(c) for c in sol["cells"]])
                    assert t21.is_cycle(v) and not t21.is_stab(v)
                    objs.append((v2i(t21.sig(v)), v))
    frame_report("tdg432(18,12)", t21, objs,
                 extra=dict(source="L12x2 witness + hunt solutions"))

    # decomposition-only at the larger members (no populations)
    for lm in [(24, 18), (18, 18), (30, 24)]:
        code = member_code(*lm)
        dec = decomposition(code)
        out["frames"][f"decomp{lm}"] = dict(
            lm=list(lm), k=code.k,
            decomposition={k: v for k, v in dec.items()
                           if not k.endswith("basis")})
        print(f"[decomp {lm}] k={code.k} {out['frames'][f'decomp{lm}']['decomposition']}",
              flush=True)
        rss_guard(f"decomp{lm}")
    out["wall"]["total_s"] = round(time.time() - t0, 1)
    (DATA / "s11_classify.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA / 's11_classify.json'} ({out['wall']['total_s']} s)")


# ------------------------------------------------------------- hunts
def hunt(code: TowerCode, Wmax: int, dense_x=True, dense_y=True,
         nsol=20, time_limit=600.0, extra_block=None, tag=""):
    """SAT hunt for nontrivial X-logicals of weight <= Wmax that are
    gap-dense (every cyclic 4-window along each dense axis occupied).
    Existence lane only.  Returns the list of verified solutions
    (translation-orbit distinct) and the status."""
    from pycryptosat import Solver
    from pysat.card import CardEnc, EncType
    l, m = code.G.orders
    n = code.n
    top = n
    clauses = []
    xors = []
    for row in code.HZ:
        lits = [int(i) + 1 for i in np.nonzero(row)[0]]
        xors.append((lits, False))
    # nontriviality: OR_i (pairing_i = 1) via auxiliaries p_i <-> XOR
    # Encode: p_i is a fresh var, XOR(z_i-lits + [p_i]) = 0, and OR p_i.
    pv = []
    for zr in code.zreps:
        top += 1
        pv.append(top)
        lits = [int(i) + 1 for i in np.nonzero(zr)[0]] + [top]
        xors.append((lits, False))
    clauses.append(pv)
    if dense_x:
        for x0 in range(l):
            lits = []
            for blk in range(2):
                for dx in range(SPAN):
                    for y in range(m):
                        lits.append(blk * code.ng +
                                    code.G.index(((x0 + dx) % l, y)) + 1)
            clauses.append(lits)
    if dense_y:
        for y0 in range(m):
            lits = []
            for blk in range(2):
                for dy in range(SPAN):
                    for x in range(l):
                        lits.append(blk * code.ng +
                                    code.G.index((x, (y0 + dy) % m)) + 1)
            clauses.append(lits)
    card = CardEnc.atmost(lits=list(range(1, n + 1)), bound=Wmax,
                          top_id=top, encoding=EncType.seqcounter)
    clauses.extend(card.clauses)
    top = max(top, card.nv)
    sols = []
    status = "unknown"
    t0 = time.time()
    perms = None
    blocked = 0
    while len(sols) < nsol:
        left = time_limit - (time.time() - t0)
        if left <= 0:
            status = "timeout"
            break
        s = Solver(time_limit=max(1.0, left))
        for cl in clauses:
            s.add_clause(cl)
        for lits, rhs in xors:
            s.add_xor_clause(lits, rhs)
        sat, sol = s.solve()
        if sat is None:
            status = "timeout"
            break
        if not sat:
            status = "unsat" if not sols else "exhausted"
            break
        v = np.array([1 if sol[i + 1] else 0 for i in range(n)],
                     dtype=np.uint8)
        assert code.is_cycle(v) and not code.is_stab(v)
        assert v.sum() <= Wmax
        cells = cells_of(code, v)
        if dense_x:
            assert gap_structure([c[1] for c in cells], l)[1] < SPAN
        if dense_y:
            assert gap_structure([c[2] for c in cells], m)[1] < SPAN
        sols.append(dict(w=int(v.sum()), cells=sorted(cells),
                         sig=int(v2i(code.sig(v)))))
        print(f"  [{tag}] solution {len(sols)}: w={int(v.sum())} "
              f"({time.time() - t0:.1f} s)", flush=True)
        # block the whole translation orbit of this support
        for dx in range(l):
            for dy in range(m):
                sh = [blk * code.ng + code.G.index(((x + dx) % l,
                                                    (y + dy) % m)) + 1
                      for blk, x, y in cells]
                clauses.append([-lit for lit in sh])
                blocked += 1
        status = "found"
    return sols, status, round(time.time() - t0, 1)


def lane_hunt(args):
    t0 = time.time()
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    out = {"runs": []}
    plan = [
        # (l, m, Wmax, nsol, time_limit, label)
        (12, 12, 17, 1, 300, "control: (12,12) <= 17 must be UNSAT (d = 18)"),
        (12, 12, 18, 6, 300, "control: (12,12) <= 18 must find the witness "
                             "family"),
        (12, 6, 12, 12, 300, "gross gap-dense <= 12 (consistency with the "
                             "census classification)"),
        (18, 12, 23, 1, 600, "control: (18,12) <= 23 must be UNSAT (d = 24)"),
        (18, 12, 24, 12, args.t1812, "THE EQUALITY QUESTION: gap-dense "
                                     "nontrivial at 24 at (18,12)"),
    ]
    if args.r3:
        plan.append((24, 18, 35, 1, args.t2418,
                     "falsification probe: (24,18) gap-dense <= 35"))
        plan.append((24, 18, 36, 3, args.t2418,
                     "(24,18) gap-dense at 36 (equality probe)"))
    for l, m, W, nsol, tl, label in plan:
        code = member_code(l, m)
        print(f"== {label}: (l,m)=({l},{m}) W<={W}", flush=True)
        sols, status, wall = hunt(code, W, nsol=nsol, time_limit=tl,
                                  tag=f"({l},{m})<={W}")
        run = dict(lm=[l, m], Wmax=W, label=label, status=status,
                   n_solutions=len(sols), wall_s=wall, solutions=sols)
        # classify every solution
        for sol in sols:
            v = vec_of(code, [tuple(c) for c in sol["cells"]])
            rec = classify_object(code, v)
            sol["classification"] = {k: rec[k] for k in
                                     ("sector", "gap_sector", "gap_x",
                                      "gap_y", "ncomp", "compact_dirs")}
            print(f"   -> w={sol['w']} sector={rec['sector']} "
                  f"gap={rec['gap_sector']} dirs="
                  f"{[c['ab'] for c in rec['compact_dirs']]}", flush=True)
        print(f"   status={status} n={len(sols)} wall={wall} s", flush=True)
        out["runs"].append(run)
        (DATA / "s11_hunt.json").write_text(json.dumps(out, indent=1))
        rss_guard("hunt")
    out["wall_s"] = round(time.time() - t0, 1)
    (DATA / "s11_hunt.json").write_text(json.dumps(out, indent=1))
    print(f"wrote {DATA / 's11_hunt.json'} ({out['wall_s']} s)")


# ---------------------------------------------------------- the norm
def lane_norm(args):
    """N(u) for a list of directions: twisted-atlas exhaustion at
    Wcap on C_u ~= Z x Z_g with the transported pair; compact
    triviality of every find by the deterministic generator march."""
    validate_banked(LAB / "data")
    print("validate_banked: PASS", flush=True)
    _argv = sys.argv
    sys.argv = [_argv[0]]
    from a40_s5_twisted_atlas import (  # noqa: E402
        transform_class, tr_supp, atlas_run,
    )
    from a40_s4_compact_triviality import march_generator  # noqa: E402
    from a40_s4_phase_atlas import bar  # noqa: E402
    sys.argv = _argv
    out = {"rows": []}
    t0 = time.time()
    dirs = [tuple(int(x) for x in d.split(",")) for d in args.dirs]
    for (t, p) in dirs:
        Wcap = args.wcap
        w1, w2, g = transform_class(p, t)
        P = tr_supp(A_L, w1, w2, g)
        Q = tr_supp(B_L, w1, w2, g)
        row = dict(u=[t, p], g=g, w1=list(w1), w2=list(w2), P=sorted(P),
                   Q=sorted(Q), Wcap=Wcap)
        t1 = time.time()
        try:
            rows, npop = atlas_run(P, Q, g, Wcap, max_states=args.max_states,
                                   max_paths=args.max_paths)
        except RuntimeError as e:
            row["status"] = f"cap: {e}"
            print(f"u=({t},{p}) g={g}: {row['status']}", flush=True)
            out["rows"].append(row)
            continue
        rss_guard(f"norm {t},{p}")
        # compact triviality of every find
        Pb, Qb = bar(P), bar(Q)
        spec = {}
        nontriv = []
        for r in rows:
            pts = r["pts"]
            v1 = {}
            v2 = {}
            for (c, y, blk) in pts:
                dd = v1 if blk == 0 else v2
                dd[c] = dd.get(c, 0) ^ (1 << y)
            s = march_generator(Pb, Qb, v1, v2, g)
            key = f"w{r['weight']}"
            spec.setdefault(key, dict(n=0, compact_trivial=0,
                                      torus_nontrivial=0))
            spec[key]["n"] += 1
            if s is not None:
                spec[key]["compact_trivial"] += 1
            if r["nontrivial"]:
                spec[key]["torus_nontrivial"] += 1
            if s is None:
                nontriv.append(r["weight"])
        row.update(status="complete", npop=npop, n_cycles=len(rows),
                   spectrum=dict(sorted(spec.items())),
                   N_lower=(min(nontriv) if nontriv else Wcap + 1),
                   N_exact=(min(nontriv) if nontriv else None),
                   wall_s=round(time.time() - t1, 1))
        print(f"u=({t},{p}) g={g}: {row['n_cycles']} compact cycles "
              f"<= {Wcap}, spectrum {row['spectrum']}, N(u) "
              f"{'= ' + str(row['N_exact']) if row['N_exact'] else '>= ' + str(Wcap + 1)}"
              f" ({row['wall_s']} s, {npop} pops)", flush=True)
        out["rows"].append(row)
        out["wall_s"] = round(time.time() - t0, 1)
        (DATA / f"s11_norm_{args.tag}.json").write_text(
            json.dumps(out, indent=1))
    print(f"wrote {DATA / f's11_norm_{args.tag}.json'}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="lane", required=True)
    c = sub.add_parser("classify")
    c.add_argument("--keep", action="store_true",
                   help="keep every object record in the json")
    h = sub.add_parser("hunt")
    h.add_argument("--r3", action="store_true")
    h.add_argument("--t1812", type=float, default=1800.0)
    h.add_argument("--t2418", type=float, default=1800.0)
    nrm = sub.add_parser("norm")
    nrm.add_argument("--dirs", nargs="+", required=True,
                     help="directions as t,p (u = (t, p))")
    nrm.add_argument("--wcap", type=int, required=True)
    nrm.add_argument("--tag", required=True)
    nrm.add_argument("--max-states", type=int, default=4_000_000)
    nrm.add_argument("--max-paths", type=int, default=400_000)
    args = ap.parse_args()
    # every lane tees its own log (shell redirection is unavailable in
    # the headless environment — S10 incident ii)
    logname = f"s11_{args.lane}" + (f"_{args.tag}" if args.lane == "norm"
                                    else "")
    logf = (DATA / f"{logname}.log").open("a")

    class Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, s):
            for st in self.streams:
                st.write(s)
                st.flush()

        def flush(self):
            for st in self.streams:
                st.flush()
    sys.stdout = Tee(sys.__stdout__, logf)
    sys.stderr = Tee(sys.__stderr__, logf)
    print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')} argv={sys.argv[1:]}")
    {"classify": lane_classify, "hunt": lane_hunt,
     "norm": lane_norm}[args.lane](args)


if __name__ == "__main__":
    main()
