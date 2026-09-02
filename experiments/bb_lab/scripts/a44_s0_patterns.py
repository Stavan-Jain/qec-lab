#!/usr/bin/env python3
"""A44 S0 — the banked D-minimizer families as Laurent patterns.

For every (kind, shape) family of A43 §5.1 at (6,6), (12,12) and gross
(objects re-read from data/a43/s1_l3_shape.json, every one kernel
re-verified), write a representative as (f, g) in F2[x^±, y^±]^2 and
TEST the structural readings:

  x-band   g = 0, rows c + 3i (i < m/3), p_{i+1} = (x^-1 + x^-2) p_i,
           closure ((1+x)^s + x^{2s}) p_0 = 0 in F2[x]/(x^l + 1), s = m/3;
  y-band   the mirror (f = 0, columns c + 3j, q_{j-1} = (y^-1 + y^-2) q_j);
  helix    a period lattice Pi of rank >= 1 inside Lambda (pullback from a
           shear frame) and, on the (6,6) strand, f = x^-3 y g.

Then the two falsify-first computations the charter asks for:
  * the exact minimum cost of an x-band at (6r, 6r), r = 1..6 (min over
    the annihilator ideal; the band is closed-form, so this is exact);
  * the single-helix stack g = x sum_j (x^2 y^2)^j, f = x^-3 y g at
    (6r, 6r): a cycle of weight 6r for every r, and WHY it is trivial for
    r >= 2 (the trivializer h is exhibited).

Output data/a44/s0_patterns.{json,log}.  No SAT.  ~1 min.
"""
from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a44_common import (  # noqa: E402
    LAB, DATA, Logger, member_code, pair_of, vec_of_pair, alpha_beta, pmul,
    padd, pshift, poly_str, is_cycle_pair, assert_laurent_conventions,
    picture, rss_guard, C,
)
from bb_lab.tower import validate_banked, kernel_basis  # noqa: E402

log = Logger(DATA / "s0_patterns.log")


# ------------------------------------------------------------ 1-var helpers
def x1(p: frozenset, l: int):
    """Row polynomial in x only (as a set of exponents mod l)."""
    return frozenset(e[0] % l for e in p)


def mul1(p, q, l):
    cnt = {}
    for a in p:
        for b in q:
            k = (a + b) % l
            cnt[k] = cnt.get(k, 0) ^ 1
    return frozenset(k for k, v in cnt.items() if v)


def pow1(p, e, l):
    r = frozenset({0})
    for _ in range(e):
        r = mul1(r, p, l)
    return r


def binom_pow(s, l):
    """(1+x)^s + x^{2s} mod x^l + 1."""
    return frozenset(pow1(frozenset({0, 1}), s, l) ^ {(2 * s) % l})


def ann_ideal_basis(P, l):
    """F2-basis of ann(P) in F2[x]/(x^l+1) (kernel of multiplication)."""
    M = np.zeros((l, l), dtype=np.uint8)
    for j in range(l):
        for e in mul1(frozenset({j}), P, l):
            M[e, j] ^= 1
    return [frozenset(int(i) for i in np.nonzero(k)[0]) for k in kernel_basis(M)]


def band_rows(p0, s, l):
    """p_i = (x^-1 + x^-2)^i p_0, i < s."""
    step = frozenset({(-1) % l, (-2) % l})
    rows = [p0]
    for _ in range(s - 1):
        rows.append(mul1(rows[-1], step, l))
    return rows


def band_min_cost(l, m, dim_cap=20):
    """Exact minimum weight of an x-band at (l, m) over the closure ideal."""
    s = m // 3
    P = binom_pow(s, l)
    basis = ann_ideal_basis(P, l)
    dim = len(basis)
    if dim > dim_cap:
        return None, dim, None
    best, arg = None, None
    for mask in range(1, 1 << dim):
        p0 = frozenset()
        for i in range(dim):
            if (mask >> i) & 1:
                p0 = frozenset(p0 ^ basis[i])
        rows = band_rows(p0, s, l)
        # closure check: (x^-1+x^-2)^s p0 == p0
        assert mul1(rows[-1], frozenset({(-1) % l, (-2) % l}), l) == p0
        w = sum(len(r) for r in rows)
        if best is None or w < best:
            best, arg = w, (p0, rows)
    return best, dim, arg


def band_vector(code, rows, c=1):
    """rows[i] = (x^-1+x^-2)^i p_0 sits at row c - 3i (downward recursion)."""
    l, m = code.G.orders
    f = set()
    for i, r in enumerate(rows):
        for e in r:
            f.add((e, (c - 3 * i) % m))
    return vec_of_pair(code, frozenset(f), frozenset())


# ------------------------------------------------------------ structure tests
def period_lattice(code, cells):
    l, m = code.G.orders
    cs = frozenset(cells)
    per = []
    for dx in range(l):
        for dy in range(m):
            if (dx, dy) == (0, 0):
                continue
            if frozenset((b, (x + dx) % l, (y + dy) % m) for b, x, y in cs) == cs:
                per.append((dx, dy))
    return per


def test_xband(f, g, lm):
    l, m = lm
    if g or not f:
        return None
    ys = sorted({e[1] for e in f})
    c = ys[0]
    if any((y - c) % 3 for y in ys):
        return dict(ok=False, why="rows not 3-spaced")
    s = m // 3
    rows = [x1(frozenset(e for e in f if e[1] == (c + 3 * i) % m), l)
            for i in range(s)]
    step = frozenset({(-1) % l, (-2) % l})
    # beta f = 0 row by row: (1 + x^-1) p_i = x p_{i-1}, i.e. the recursion
    # runs DOWNWARD, p_{i-1} = (x^-1 + x^-2) p_i
    rec = all(mul1(rows[i], step, l) == rows[(i - 1) % s] for i in range(s))
    P = binom_pow(s, l)
    closure = (mul1(P, rows[0], l) == frozenset())
    return dict(ok=bool(rec and closure), recursion=bool(rec), closure=bool(closure),
                c=c, rows=[sorted(r) for r in rows], row_weights=[len(r) for r in rows])


def test_yband(f, g, lm):
    l, m = lm
    if f or not g:
        return None
    xs = sorted({e[0] for e in g})
    c = xs[0]
    if any((x - c) % 3 for x in xs):
        return dict(ok=False, why="cols not 3-spaced")
    s = l // 3
    cols = [frozenset(e[1] for e in g if e[0] == (c + 3 * j) % l) for j in range(s)]
    step = frozenset({(-1) % m, (-2) % m})
    # alpha g = 0: (1+y^-1) q_j + y q_{j-1} = 0  =>  q_{j-1} = (y^-1+y^-2) q_j
    rec = all(mul1(cols[j], step, m) == cols[(j - 1) % s] for j in range(s))
    return dict(ok=bool(rec), recursion=bool(rec), c=c,
                cols=[sorted(q) for q in cols], col_weights=[len(q) for q in cols])


def test_helix_relation(f, g, lm):
    """f == x^-3 y g (the (6,6) strand relation) up to a translation?"""
    l, m = lm
    if not f or not g:
        return None
    tgt = pshift(g, ((-3) % l, 1), lm)
    for dx in range(l):
        for dy in range(m):
            if pshift(tgt, (dx, dy), lm) == f:
                return dict(ok=True, extra_shift=(dx, dy))
    return dict(ok=False)


def trivializer(code, v):
    """h with (f, g) = h (alpha, beta), i.e. v = H_X^T h; None if not stab."""
    if not code.is_stab(v):
        return None
    # solve H_X^T h = v over F2 (n x ng system)
    M = np.concatenate([code.HX.T, v.reshape(-1, 1)], axis=1).astype(np.uint8)
    rows, cols = M.shape
    piv = []
    r = 0
    M = M.copy()
    for c in range(cols - 1):
        pr = None
        for i in range(r, rows):
            if M[i, c]:
                pr = i
                break
        if pr is None:
            continue
        M[[r, pr]] = M[[pr, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
    h = np.zeros(cols - 1, dtype=np.uint8)
    for i, c in enumerate(piv):
        h[c] = M[i, -1]
    assert not (((code.HX.T @ h) % 2) ^ v).any()
    return frozenset(code.G.from_index(int(i)) for i in np.nonzero(h)[0])


# ------------------------------------------------------------ main
def main():
    t0 = time.time()
    validate_banked(LAB / "data")
    log("validate_banked: PASS")
    for lm in [(6, 6), (12, 12), (12, 6), (18, 18)]:
        assert_laurent_conventions(member_code(*lm))
    log("Laurent conventions (beta f = alpha g; stab = h(alpha,beta)): PASS at "
        "(6,6) (12,12) (12,6) (18,18)")

    shape = json.loads((LAB / "data" / "a43" / "s1_l3_shape.json").read_text())
    out = {"families": {}, "band_costs": {}, "helix_stack": {}, "stack_of_66": {}}

    # ---- families -------------------------------------------------------
    for fname, fr in shape["frames"].items():
        lm = tuple(fr["lm"])
        code = member_code(*lm)
        fams: dict = {}
        for o in fr["objects"]:
            cells = [tuple(c) for c in o["cells"]]
            colc = Counter(c[1] for c in cells)
            rowc = Counter(c[2] for c in cells)
            key = (f"{o['sector']}|{o['kind']}|w{o['w']}|nx{o['n_x']}|ny{o['n_y']}"
                   f"|gx{o['gap_x']}|gy{o['gap_y']}|blk{tuple(o['blocks'])}"
                   f"|cmax{max(colc.values())}|rmax{max(rowc.values())}")
            fams.setdefault(key, []).append(cells)
        log(f"=== {fname} {lm}: {len(fr['objects'])} objects, {len(fams)} families")
        recs = {}
        for key, members in sorted(fams.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            if fname.startswith("gross") and "pure-y" not in key and "W_x" not in key:
                # blobs: record only the count + one Laurent form
                pass
            cells = members[0]
            v = C.vec_of(code, cells)
            assert code.is_cycle(v) and not code.is_stab(v)
            f, g = pair_of(code, v)
            assert is_cycle_pair(f, g, lm)
            per = period_lattice(code, cells)
            xb = test_xband(f, g, lm)
            yb = test_yband(f, g, lm)
            hx = test_helix_relation(f, g, lm) if (f and g and len(per) > 0) else None
            reading = ("x-band" if xb and xb["ok"] else "y-band" if yb and yb["ok"]
                       else "helix/pullback" if per else "blob")
            centre = min(f | g) if (f | g) else (0, 0)
            rec = dict(n=len(members), reading=reading, periods=per,
                       f=poly_str(f), g=poly_str(g),
                       f_centred=poly_str(f, lm, centre), g_centred=poly_str(g, lm, centre),
                       xband=xb, yband=yb, helix_rel=hx,
                       cells=sorted(cells))
            recs[key] = rec
            log(f"  [{len(members):4d}] {key}\n"
                f"      reading={reading} periods={per}\n"
                f"      f = {rec['f']}\n      g = {rec['g']}"
                + (f"\n      x-band rows(weights)={xb['row_weights']} rec={xb['recursion']} "
                   f"closure={xb['closure']}" if xb else "")
                + (f"\n      y-band cols(weights)={yb['col_weights']} rec={yb['recursion']}"
                   if yb else "")
                + (f"\n      helix relation f = x^-3 y g: {hx}" if hx else ""))
            if reading != "blob" or len(members) <= 72:
                log(picture(code, cells))
        out["families"][fname] = dict(lm=lm, families=recs)
        rss_guard(fname)

    # ---- x-band exact minimum at (6r,6r) ---------------------------------
    log("=== x-band exact minimum cost at (6r,6r) (min over the closure ideal)")
    for r in range(1, 7):
        l = m = 6 * r
        best, dim, arg = band_min_cost(l, m)
        rec = dict(r=r, lm=[l, m], target=12 * r - 6, ideal_dim=dim, min_cost=best)
        if arg is not None:
            p0, rows = arg
            rec["p0"] = sorted(p0)
            rec["row_weights"] = [len(x) for x in rows]
            if 6 * r <= 24:
                code = member_code(l, m)
                v = band_vector(code, rows)
                rec["cycle"] = bool(code.is_cycle(v))
                rec["nontrivial"] = bool(not code.is_stab(v))
        out["band_costs"][r] = rec
        log(f"  r={r} ({l},{m}): target {12*r-6}, ideal dim {dim}, min x-band cost "
            f"{best} rows={rec.get('row_weights')} p0={rec.get('p0')} "
            f"cycle={rec.get('cycle')} nontrivial={rec.get('nontrivial')}")
        rss_guard(f"band r={r}")

    # ---- the single-helix stack ------------------------------------------
    log("=== single helix g = x sum_j (x^2y^2)^j, f = x^-3 y g at (6r,6r)")
    for r in range(1, 6):
        l = m = 6 * r
        lm = (l, m)
        code = member_code(l, m)
        g = frozenset(((1 + 2 * j) % l, (2 * j) % m) for j in range(3 * r))
        f = pshift(g, ((-3) % l, 1), lm)
        v = vec_of_pair(code, f, g)
        cyc = code.is_cycle(v)
        stab = code.is_stab(v) if cyc else None
        rec = dict(r=r, weight=int(v.sum()), cycle=bool(cyc), stab=bool(stab) if cyc else None)
        if cyc and stab:
            h = trivializer(code, v)
            rec["h"] = poly_str(h)
            rec["h_weight"] = len(h)
            hper = period_lattice(code, [(0,) + e for e in h])
            rec["h_periods"] = hper[:6]
        out["helix_stack"][r] = rec
        log(f"  r={r}: weight {rec['weight']} cycle={cyc} stab={rec['stab']}"
            + (f" trivializer h (weight {rec['h_weight']}) = {rec['h']}" if cyc and stab else ""))
        rss_guard(f"helix r={r}")

    # ---- stacking the (6,6) x-band carrier r times -----------------------
    log("=== stacking the (6,6) x-band carrier: sum_j y^{6j} (y p + y^4 q), "
        "p = 1+x^2, q = (x^-1+x^-2)p, at (6r,6r)")
    for r in range(1, 5):
        l = m = 6 * r
        code = member_code(l, m)
        p = frozenset({0, 2})
        q = mul1(p, frozenset({(-1) % l, (-2) % l}), l)
        f = set()
        for j in range(r):
            for e in p:
                f.add((e, (1 + 6 * j) % m))
            for e in q:
                f.add((e, (4 + 6 * j) % m))
        v = vec_of_pair(code, frozenset(f), frozenset())
        cyc = code.is_cycle(v)
        # the failing identity: (1+x^2+x^4) p = 1 + x^6 must vanish mod x^l+1
        fail = mul1(frozenset({0, 2, 4}), p, l)
        rec = dict(r=r, weight=int(v.sum()), cycle=bool(cyc),
                   nontrivial=bool(not code.is_stab(v)) if cyc else None,
                   closure_residual=sorted(fail))
        out["stack_of_66"][r] = rec
        log(f"  r={r}: weight {rec['weight']} cycle={cyc} nontrivial={rec['nontrivial']} "
            f"closure residual (1+x^2+x^4)p mod x^{l}+1 = {sorted(fail)}")

    out["wall_s"] = round(time.time() - t0, 1)
    out["rss_gb"] = round(rss_guard("end"), 2)
    (DATA / "s0_patterns.json").write_text(json.dumps(out, indent=1, default=str))
    log(f"wrote {DATA / 's0_patterns.json'} ({out['wall_s']} s, RSS {out['rss_gb']} GB)")


if __name__ == "__main__":
    main()
