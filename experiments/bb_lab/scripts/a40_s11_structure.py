#!/usr/bin/env python3
"""A40 S11 structure checks (note §16.3).

(1) PARITY LAW of the cylinder images: im(pi_u*) over 22 directions
    u = (a l, b m) at five frames depends only on (a, b) mod 2 — three
    pairwise-complementary 6-dim subspaces Wx, Wy, Wd, with Wd the graph
    of an isomorphism Wx -> Wy (every basis vector of Wd has both
    components nonzero).  Output data/a40/s11_parity_law.json.
(2) ker(rho^*) = W_u for the three Z2-covers T_u = Z^2/<u, 2 Lambda>
    (x-cover (2l, m), y-cover (l, 2m), diagonal cover <(l,m),(2l,0)>),
    rank rho^* = 6 each; the (12,12) witness pulls back NONTRIVIALLY to
    all three (a mixed class).  Output data/a40/s11_pullback.json.

Run: uv run python scripts/a40_s11_structure.py [parity|pullback|both]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(LAB / "src"))
sys.path.insert(0, str(LAB / "scripts"))

import a40_s11_compare as C  # noqa: E402
from bb_lab.tower import TowerCode, v2i, rref_ints, in_span  # noqa: E402
from a40_s4_phase_triage import snf2  # noqa: E402

DATA = LAB / "data" / "a40"


def parity():
    out = {}
    for lm in [(6, 6), (12, 6), (12, 12), (18, 12), (24, 18)]:
        code = C.member_code(*lm)
        bases = {}
        for ab in C.DIRS + [(4, 1), (1, 4), (4, 3), (3, 4), (5, 2), (2, 5)]:
            b, info = C.image_classes(code, ab, K=4)
            bases[ab] = sorted(b)
        groups = {}
        for ab, b in bases.items():
            key = (ab[0] % 2, ab[1] % 2)
            groups.setdefault(key, []).append(ab)
            ref = bases[groups[key][0]]
            rb, rp = rref_ints(list(ref))
            same = all(in_span(x, rb, rp) for x in b) and len(b) == len(ref)
            assert same, (lm, ab, key)
        Wx, Wy, Wd = bases[(0, 1)], bases[(1, 0)], bases[(1, 1)]

        def dim(*bs):
            r, _ = rref_ints([x for b in bs for x in b])
            return len(r)
        rec = dict(dims={str(k): len(bases[v[0]]) for k, v in groups.items()},
                   groups={str(k): [list(a) for a in v]
                           for k, v in groups.items()},
                   dim_xy=dim(Wx, Wy), dim_xd=dim(Wx, Wd), dim_yd=dim(Wy, Wd),
                   dim_all=dim(Wx, Wy, Wd))
        splits = []
        for d in Wd:
            cx, cy = C.class_split(d, Wx, Wy)
            splits.append((cx != 0, cy != 0))
        rec["Wd_split_both_nonzero"] = all(a and b for a, b in splits)
        out[str(lm)] = rec
        print(lm, rec["dims"], "xy", rec["dim_xy"], "xd", rec["dim_xd"],
              "yd", rec["dim_yd"], "all", rec["dim_all"], "Wd graph:",
              rec["Wd_split_both_nonzero"], flush=True)
    (DATA / "s11_parity_law.json").write_text(json.dumps(out, indent=1))


def cover_code(lm, gens):
    """BB code on Z^2/<gens> with the fixed pair, plus torus cell ->
    cover cells (index-2 covers of the member torus)."""
    M = [list(gens[0]), list(gens[1])]
    D, U, V = snf2(M)
    o1, o2 = abs(D[0][0]), abs(D[1][1])

    def phi(e):
        return ((e[0] * V[0][0] + e[1] * V[1][0]) % o1,
                (e[0] * V[0][1] + e[1] * V[1][1]) % o2)

    def tr(supp):
        cnt = {}
        for e in supp:
            k = phi(e)
            cnt[k] = cnt.get(k, 0) ^ 1
        return frozenset(k for k, v in cnt.items() if v)
    code = TowerCode(f"cover{gens}", (o1, o2), tr(C.A_L), tr(C.B_L))
    l, m = lm

    def lifts(x, y):
        out = set()
        for i in range(2):
            for j in range(2):
                out.add(phi((x + i * l, y + j * m)))
        return out
    return code, lifts


def pullback(code_T, code_C, lifts, v):
    vc = np.zeros(code_C.n, dtype=np.uint8)
    for blk, x, y in C.cells_of(code_T, v):
        for f in lifts(x, y):
            vc[blk * code_C.ng + code_C.G.index(f)] ^= 1
    assert code_C.is_cycle(vc)
    return vc


def pullback_checks():
    out = {}
    for lm in [(12, 12), (12, 6), (18, 12)]:
        l, m = lm
        T = C.member_code(*lm)
        covers = {"x": ((2 * l, 0), (0, m)), "y": ((l, 0), (0, 2 * m)),
                  "d": ((l, m), (2 * l, 0))}
        W = {ab: C.image_classes(T, ab, K=4)[0]
             for ab in [(0, 1), (1, 0), (1, 1)]}
        names = {(0, 1): "x", (1, 0): "y", (1, 1): "d"}
        rec = {}
        for cname, gens in covers.items():
            code_C, lifts = cover_code(lm, gens)
            imgs = []
            for i in range(T.k):
                vc = pullback(T, code_C, lifts, T.xreps[i])
                imgs.append(v2i(code_C.sig(vc)))
            basis, piv, tags = [], [], []
            kerv = []
            for i, r in enumerate(imgs):
                cur, tag = r, 1 << i
                for bb, pp, tt in zip(basis, piv, tags):
                    if (cur >> pp) & 1:
                        cur ^= bb
                        tag ^= tt
                if cur:
                    basis.append(cur)
                    piv.append((cur & -cur).bit_length() - 1)
                    tags.append(tag)
                else:
                    kerv.append(tag)
            ksigs = []
            for tag in kerv:
                v = np.zeros(T.n, dtype=np.uint8)
                for i in range(T.k):
                    if (tag >> i) & 1:
                        v ^= T.xreps[i]
                ksigs.append(v2i(T.sig(v)))
            kb, _ = rref_ints(ksigs)
            which = None
            for ab, Wb in W.items():
                wb, wp = rref_ints(list(Wb))
                if len(kb) == len(wb) and all(in_span(x, wb, wp)
                                              for x in kb):
                    which = names[ab]
            rec[cname] = dict(cover_k=code_C.k,
                              cover_orders=list(code_C.G.orders),
                              rank_pullback=len(basis), dim_ker=len(kb),
                              ker_equals=which)
            print(lm, cname, rec[cname], flush=True)
            assert which == cname
        out[str(lm)] = rec
        if lm == (12, 12):
            wit = C.a36_witness(T)
            for cname, gens in covers.items():
                code_C, lifts = cover_code(lm, gens)
                vc = pullback(T, code_C, lifts, wit)
                nt = bool(not code_C.is_stab(vc))
                print("  witness pullback to", cname, "weight",
                      int(vc.sum()), "nontrivial:", nt, flush=True)
                out[str(lm)][cname]["witness_pullback_nontrivial"] = nt
                assert nt
    (DATA / "s11_pullback.json").write_text(json.dumps(out, indent=1))


def class_kinds():
    """Sector x class-kind histogram over the classify populations
    (s11_classify.json --keep): pure-x / pure-y / pure-d / mixed.
    Output s11_class_kinds.json."""
    cl = json.loads((DATA / "s11_classify.json").read_text())
    out = {}
    for name, lm in [("gross(12,6)", (12, 6)),
                     ("two-gross(12,12)", (12, 12)), ("bb72(6,6)", (6, 6))]:
        code = C.member_code(*lm)
        W = {}
        for ab, tag in [((0, 1), "x"), ((1, 0), "y"), ((1, 1), "d")]:
            b, _ = C.image_classes(code, ab, K=4)
            W[tag] = rref_ints(list(b))
        hist = {}
        classes = {}
        for o in cl["frames"][name]["objects"]:
            sig = o["class"]
            kind = "mixed"
            for tag, (bb, pp) in W.items():
                if in_span(sig, bb, pp):
                    kind = "pure-" + tag
            key = f"{o['sector']}|{kind}|w{o['w']}"
            hist[key] = hist.get(key, 0) + 1
            classes.setdefault(kind, set()).add(sig)
        rec = dict(hist=dict(sorted(hist.items())),
                   classes={k: len(v) for k, v in classes.items()})
        out[name] = rec
        print(name, rec, flush=True)
    (DATA / "s11_class_kinds.json").write_text(json.dumps(out, indent=1))


def diag_minima():
    """Exact class minima of the 63 diagonal-parity (pure-d) classes of
    gross by census_pass at W = 14, 16 (walk kernel, census-complete
    per class).  Output s11_diag_classes.json."""
    import time
    from bb_lab import cosetbz
    from bb_lab.tower import rep_for, i2v, validate_banked
    from a38_c37xx_freeze import census_pass
    validate_banked(LAB / "data")
    binp = cosetbz.build_kernel()
    gross = C.member_code(12, 6)
    Wd, _ = C.image_classes(gross, (1, 1), K=4)
    pts = {0}
    for b in Wd:
        pts |= {p ^ b for p in pts}
    pts.discard(0)
    classes = sorted(pts)
    assert len(classes) == 63
    out = {}
    for W in (14, 16):
        t0 = time.time()
        mins = {}
        for lo in range(0, 63, 51):
            chunk = classes[lo:lo + 51]
            hits = census_pass(binp, gross,
                               [(f"C{c}", rep_for(gross, c)) for c in chunk],
                               W, f"s11_diag_{W}_{lo}")
            for c in chunk:
                ws = []
                for h in hits[f"C{c}"]:
                    v = i2v(h, gross.n)
                    assert gross.is_cycle(v) and not gross.is_stab(v)
                    assert v2i(gross.sig(v)) == c
                    ws.append(int(v.sum()))
                if ws:
                    mins[c] = min(ws)
        hist = {}
        for c, w in mins.items():
            hist[w] = hist.get(w, 0) + 1
        print(f"W<={W}: {len(mins)}/63 diagonal classes populated; "
              f"min-weight hist {hist}; {time.time() - t0:.1f} s",
              flush=True)
        out[str(W)] = dict(n_populated=len(mins), hist=hist)
    (DATA / "s11_diag_classes.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    if part == "diag":
        diag_minima()
        sys.exit(0)
    if part in ("parity", "both", "all"):
        parity()
    if part in ("pullback", "both", "all"):
        pullback_checks()
    if part in ("kinds", "all"):
        class_kinds()
