"""A22 phase 1: analytic structure of the delta-side (the 94 alpha classes).

P1. m-preimage smallness: every light alpha = Atil*m with m of small site
    support?  (min over the 16 eta0-line preimages)
P2. geometry: active set S subset supp(m*) + TRIANGLE (A-triangle cover).
P3. symmetries of the delta problem: Galois semilinear partner site-map,
    the Phi involution's alpha-action; orbit count of the 94 reps.
P4. |S| histograms (special subsets vs light configs).
P5. sweep orbit reduction: subsets of <= 7 sites up to the site-symmetry
    group -> Lean certificate count.

Usage: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a22_delta_structure.py
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a22_common import *  # noqa: F401,F403
from a22_common import (
    A22_DATA, ATIL, AT_HAT, ETA0, ETA0_ROW, LOG, MU5, SIDX, SITES, THETA,
    TRIANGLE, WP, ZP, apply_theta, canon_alpha, compute_alpha_classes,
    config_cost, active_set, cost_and_h, dft_vec, gfrob, ginv, gmul, gpow,
    idft_vec, site_add, site_sub, site_type,
)

reps = compute_alpha_classes()
print(f"loaded {len(reps)} alpha classes")
assert len(reps) == 94

# ------------------------------------------------------------------ P1
print("\n== P1: m-preimage site-support of the 94 alpha classes ==")
# kernel direction in m-space: mhat = indicator at eta0
KER = idft_vec(np.array([[16 if (p, q) == ETA0 else 0 for q in range(3)]
                         for p in range(5)], dtype=np.int64) % 16)
# fix: indicator value 1 at eta0
KH = np.zeros((5, 3), dtype=np.int64)
KH[ETA0] = 1
KER = idft_vec(KH)
assert np.all(KER != 0), "eta0 character vector should be full-support"

msupp_hist = Counter()
m_star: dict[tuple, np.ndarray] = {}
for rep in reps:
    alpha = np.array(rep, dtype=np.int64)
    ah = dft_vec(alpha)
    assert ah[ETA0] == 0
    mh = np.zeros((5, 3), dtype=np.int64)
    for p in range(5):
        for q in range(3):
            if (p, q) != ETA0:
                mh[p, q] = gmul(int(ah[p, q]), ginv(int(AT_HAT[p, q])))
    m0 = idft_vec(mh)
    best, bs = None, 99
    for c in range(16):
        m = m0.copy()
        if c:
            for k in range(15):
                m[k] ^= gmul(c, int(KER[k]))
        s = int(np.count_nonzero(m))
        if s < bs:
            best, bs = m, s
    # sanity: Atil * best == alpha
    from a22_common import conv_site_gf16
    assert np.array_equal(conv_site_gf16(ATIL, best), alpha)
    msupp_hist[bs] += 1
    m_star[rep] = best
print("   min |supp m| histogram over the 94 classes:", dict(sorted(msupp_hist.items())))

# ------------------------------------------------------------------ P2
print("\n== P2: triangle-cover geometry ==")
cover_ok = 0
for rep in reps:
    alpha = np.array(rep, dtype=np.int64)
    S = set(active_set(alpha))
    m = m_star[rep]
    msites = [SITES[k] for k in range(15) if m[k]]
    covered = set()
    for ms in msites:
        for t in TRIANGLE:
            covered.add(SIDX[site_add(ms, t)])
        # B-side triangle: x * (triangle + msite) pulled back by the pairing
        for t in TRIANGLE:
            covered.add(SIDX[site_add(ms, t)])
    if S <= covered:
        cover_ok += 1
print(f"   S subset (supp m* + TRIANGLE): {cover_ok}/94")
mtypes = Counter()
for rep in reps:
    m = m_star[rep]
    mtypes["".join(sorted(site_type(int(x)) for x in m if x))] += 1
print("   m* value-type profiles:", dict(sorted(mtypes.items(), key=lambda kv: -kv[1])))

# ------------------------------------------------------------------ P3a
print("\n== P3a: Galois semilinear symmetry hunt ==")
# want a site map P and Frobenius power r with
#   THETA[P g, P g'] == Frob^r(THETA[g, g'])  for all g, g'
found_gal = []
for r in (1, 2, 3):
    for c in (1, 2, 3, 4):
        for d in (1, 2):
            for u in range(5):
                for vb in range(3):
                    ok = True
                    Pmap = {}
                    for k, (i, b) in enumerate(SITES):
                        Pmap[k] = SIDX[((c * i + u) % 5, (d * b + vb) % 3)]
                    for kg in range(15):
                        for kp in range(15):
                            if int(THETA[Pmap[kg], Pmap[kp]]) != gfrob(int(THETA[kg, kp]), r):
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        found_gal.append((r, c, d, u, vb))
print(f"   plain equivariance solutions (r,c,d,shift): {found_gal}")

# scalar-twisted version: allow J(alpha)_g = chi(g) * Frob^r(alpha_{P^-1 g})
# equivariance needs THETA[Pg, Pg'] * chi(g') == chi(g) * Frob^r(THETA[g,g'])
# with chi(g) = zeta^{a i(g)} * omega^{e b(g)} * const
found_tw = []
for r in (1, 2, 3):
    for c in (1, 2, 3, 4):
        for d in (1, 2):
            for a_tw in range(5):
                for e_tw in range(3):
                    ok = True
                    for kg in range(15):
                        ig, bg = SITES[kg]
                        Pg = SIDX[((c * ig) % 5, (d * bg) % 3)]
                        chig = gmul(ZP[(a_tw * ig) % 5], WP[(e_tw * bg) % 3])
                        for kp in range(15):
                            ip, bp = SITES[kp]
                            Pp = SIDX[((c * ip) % 5, (d * bp) % 3)]
                            chip = gmul(ZP[(a_tw * ip) % 5], WP[(e_tw * bp) % 3])
                            lhs = gmul(int(THETA[Pg, Pp]), chip)
                            rhs = gmul(chig, gfrob(int(THETA[kg, kp]), r))
                            if lhs != rhs:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        found_tw.append((r, c, d, a_tw, e_tw))
print(f"   twisted equivariance solutions (r,c,d,a,e): {found_tw}")

# whichever symmetries found: check they also fix the eta0 condition and
# permute the 94 reps
def apply_sym(alpha: np.ndarray, r: int, c: int, d: int, a_tw: int, e_tw: int) -> np.ndarray:
    out = np.zeros(15, dtype=np.int64)
    cinv = pow(c, -1, 5)
    dinv = pow(d, -1, 3)
    for kg, (ig, bg) in enumerate(SITES):
        src = SIDX[((cinv * ig) % 5, (dinv * bg) % 3)]
        chig = gmul(ZP[(a_tw * ig) % 5], WP[(e_tw * bg) % 3])
        out[kg] = gmul(chig, gfrob(int(alpha[src]), r))
    return out


rep_set = set(reps)
for sym in found_tw:
    r, c, d, a_tw, e_tw = sym
    img = {canon_alpha(apply_sym(np.array(rp, dtype=np.int64), r, c, d, a_tw, e_tw))
           for rp in reps}
    print(f"   sym {sym}: permutes the 94 reps: {img == rep_set}")

# ------------------------------------------------------------------ P3b
print("\n== P3b: Phi involution action on alpha ==")
# Phi(u|v) = (x^4 y^9 sigma(v) | x y^6 sigma(u)); on sites sigma = id,
# fiber z -> z^-1 (delta conj: Frob^2), sigma twists x^i by zeta^{2i}... derive:
# guess alpha' = ff(beta') with site shift and twist; test candidates vs classes
def phi_alpha(alpha: np.ndarray) -> np.ndarray | None:
    """Hypothesis: alpha'(g) = zeta^{t1 + 2 i(g)} Frob^2(beta'(g + s))  — search
    over shift s and scalar power t1, return the first that maps every rep to
    a rep."""
    return None


beta_of = {rep: apply_theta(np.array(rep, dtype=np.int64)) for rep in reps}
hits = []
for s0 in range(5):
    for s1 in range(3):
        for t1 in range(5):
            for tw in range(5):
                ok = True
                for rep in reps[:8]:
                    alpha = np.array(rep, dtype=np.int64)
                    beta = beta_of[rep]
                    out = np.zeros(15, dtype=np.int64)
                    for kg, (ig, bg) in enumerate(SITES):
                        src = SIDX[((ig + s0) % 5, (bg + s1) % 3)]
                        sc = ZP[(t1 + tw * ig) % 5]
                        out[kg] = gmul(sc, gfrob(int(beta[src]), 2))
                    if canon_alpha(out) not in rep_set:
                        ok = False
                        break
                if ok:
                    hits.append((s0, s1, t1, tw))
print(f"   Phi-candidate (shift, scalarpow, xtwist) forms passing 8 reps: {hits}")
for cand in hits:
    s0, s1, t1, tw = cand
    ok = True
    invol = True
    for rep in reps:
        alpha = np.array(rep, dtype=np.int64)
        beta = beta_of[rep]
        out = np.zeros(15, dtype=np.int64)
        for kg, (ig, bg) in enumerate(SITES):
            src = SIDX[((ig + s0) % 5, (bg + s1) % 3)]
            sc = ZP[(t1 + tw * ig) % 5]
            out[kg] = gmul(sc, gfrob(int(beta[src]), 2))
        if canon_alpha(out) not in rep_set:
            ok = False
            break
    if ok:
        print(f"   Phi form {cand}: permutes ALL 94 reps")

# ------------------------------------------------------------------ P3c
print("\n== P3c: orbit structure of the 94 reps under found symmetries ==")
gens = []
for sym in found_tw:
    r, c, d, a_tw, e_tw = sym
    gens.append(("gal", sym))
for cand in hits:
    gens.append(("phi", cand))


def apply_gen(kind, prm, alpha: np.ndarray) -> np.ndarray:
    if kind == "gal":
        r, c, d, a_tw, e_tw = prm
        return apply_sym(alpha, r, c, d, a_tw, e_tw)
    s0, s1, t1, tw = prm
    beta = apply_theta(alpha)
    out = np.zeros(15, dtype=np.int64)
    for kg, (ig, bg) in enumerate(SITES):
        src = SIDX[((ig + s0) % 5, (bg + s1) % 3)]
        sc = ZP[(t1 + tw * ig) % 5]
        out[kg] = gmul(sc, gfrob(int(beta[src]), 2))
    return out


seen = set()
orbits = []
for rep in reps:
    if rep in seen:
        continue
    orb = {rep}
    frontier = [rep]
    while frontier:
        cur = frontier.pop()
        for kind, prm in gens:
            nxt = canon_alpha(apply_gen(kind, prm, np.array(cur, dtype=np.int64)))
            if nxt in rep_set and nxt not in orb:
                orb.add(nxt)
                frontier.append(nxt)
    orbits.append(orb)
    seen |= orb
print(f"   orbits of the 94 alpha classes under <Galois, Phi>: {len(orbits)}")
print(f"   orbit sizes: {sorted(len(o) for o in orbits)}")

# ------------------------------------------------------------------ P4
print("\n== P4: |S| structure ==")
light_by_S = Counter(len(active_set(np.array(rep, dtype=np.int64))) for rep in reps)
print(f"   |S| histogram of the 94 light alpha classes: {dict(sorted(light_by_S.items()))}")
cost_by = Counter((len(active_set(np.array(r, dtype=np.int64))),
                   config_cost(np.array(r, dtype=np.int64))) for r in reps)
print(f"   (|S|, cost) histogram: {dict(sorted(cost_by.items()))}")

# ------------------------------------------------------------------ P5
print("\n== P5: sweep orbit reduction (site-symmetry on subsets) ==")
# site maps that preserve the sweep problem: shifts (15) x found Galois site
# maps x Phi site action (if site-preserving).  Build the subgroup of S_15.
site_perms = set()
for u in range(5):
    for vb in range(3):
        perm = tuple(SIDX[((i + u) % 5, (b + vb) % 3)] for (i, b) in SITES)
        site_perms.add(perm)
for sym in found_tw:
    r, c, d, a_tw, e_tw = sym
    perm = tuple(SIDX[((c * i) % 5, (d * b) % 3)] for (i, b) in SITES)
    site_perms.add(perm)
# close under composition
changed = True
perms = set(site_perms)
while changed:
    changed = False
    for p1 in list(perms):
        for p2 in list(perms):
            comp = tuple(p1[p2[k]] for k in range(15))
            if comp not in perms:
                perms.add(comp)
                changed = True
print(f"   site-symmetry group order: {len(perms)}")
subsets = []
for size in range(0, 8):
    subsets += list(itertools.combinations(range(15), size))
canon_subs = set()
for S in subsets:
    best = min(tuple(sorted(p[k] for k in S)) for p in perms)
    canon_subs.add(best)
print(f"   subsets |S|<=7: {len(subsets)}; orbits under site symmetry: {len(canon_subs)}")
by_size = Counter(len(s) for s in canon_subs)
print(f"   orbit count by size: {dict(sorted(by_size.items()))}")
