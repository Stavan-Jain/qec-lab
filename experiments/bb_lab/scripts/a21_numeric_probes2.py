#!/usr/bin/env python
"""A21 session 1, probe set 2: structural verification for the analytic attack.

P0  ψ-symmetry: ψ(gx,gy) = (gx, 6gx+4gy) ∈ Aut(Z₅×Z₁₅), B = xy⁶·ψ(A),
    ψ⁻¹(B) = x·A, and the induced chain symmetry
    Φ(u_L,u_R) = (ψ⁻¹u_R, x²y⁶·ψ⁻¹u_L) maps cycles→cycles (split-swap)
    and boundaries→boundaries.
P1  Lemma W: wt(P⋆T) = 3|T| − 2E_P(T) + 4N_P(T) for Sidon |P| = 3.
P2  dA-cluster censuses: max E on 4-sets / 5-sets (translation-reduced,
    on the torus); all 5-sets with (E,N)=(6,0) → images → B-translate test.
P2' 4-sets with (E,N) ∈ {(3,0),(4,0),(5,1)} → images → σ_B-family test.
P3  (2,2) size-4 σ-match: σ_B(δL) vs σ_A(δR) translate test, δL ∈ dB, δR ∈ dA.
P4  (3,3) sub-case probe: weight-3 u_L by (E_B, N_B) class; which classes
    admit ANY weight-3 u_R (census says only E=0 generators).
"""
from __future__ import annotations

import numpy as np
from itertools import combinations
from collections import Counter

L, M = 5, 15
N = L * M

A_SUPP = [(0, 0), (0, 1), (1, 0)]
B_SUPP = [(1, 6), (1, 10), (2, 12)]


def idx(gx, gy):
    return (gx % L) * M + (gy % M)


def unidx(g):
    return divmod(g, M)


def translate(supp_idx, g):
    gx, gy = unidx(g)
    return frozenset(idx(px + gx, py + gy) for (px, py) in
                     (unidx(s) for s in supp_idx))


def conv_matrix(supp):
    C = np.zeros((N, N), dtype=np.uint8)
    for g in range(N):
        gx, gy = divmod(g, M)
        for (px, py) in supp:
            C[g, idx(gx - px, gy - py)] ^= 1
    return C


CA, CB = conv_matrix(A_SUPP), conv_matrix(B_SUPP)
A_IDX = frozenset(idx(*p) for p in A_SUPP)
B_IDX = frozenset(idx(*p) for p in B_SUPP)


def diffset(supp_idx):
    out = set()
    for s in supp_idx:
        sx, sy = unidx(s)
        for t in supp_idx:
            if s == t:
                continue
            tx, ty = unidx(t)
            out.add(idx(sx - tx, sy - ty))
    return frozenset(out)


dA, dB = diffset(A_IDX), diffset(B_IDX)


def psi(g):
    gx, gy = unidx(g)
    return idx(gx, 6 * gx + 4 * gy)


def psi_inv(g):
    gx, gy = unidx(g)
    # solve 6a+4b ≡ gy with a = gx:  4b ≡ gy − 6gx;  4⁻¹ = 4 (mod 15)
    return idx(gx, 4 * (gy - 6 * gx))


def apply_map(vec, fn):
    out = np.zeros(N, dtype=np.uint8)
    for g in np.flatnonzero(vec):
        out[fn(g)] ^= 1
    return out


def conv(P_idx, vec):
    out = np.zeros(N, dtype=np.uint8)
    for g in np.flatnonzero(vec):
        gx, gy = unidx(g)
        for p in P_idx:
            px, py = unidx(p)
            out[idx(gx + px, gy + py)] ^= 1
    return out


def wt(v):
    return int(v.sum())


def EN_count(P_idx, T_list):
    """E = # unordered pairs with difference in dP; N = # mult-3 cells."""
    E = 0
    dP = diffset(P_idx)
    for i in range(len(T_list)):
        for j in range(i + 1, len(T_list)):
            ax, ay = unidx(T_list[i])
            bx, by = unidx(T_list[j])
            if idx(ax - bx, ay - by) in dP:
                E += 1
    Tset = set(T_list)
    Ncnt = 0
    for c in range(N):
        cx, cy = unidx(c)
        if all(idx(cx - px, cy - py) in Tset
               for (px, py) in (unidx(p) for p in P_idx)):
            Ncnt += 1
    return E, Ncnt


rng = np.random.default_rng(21)

# ── P0: the ψ symmetry ────────────────────────────────────────────────
print("── P0: ψ-symmetry ──")
psiA = frozenset(psi(g) for g in A_IDX)
xy6 = idx(1, 6)
print("ψ(A)·xy⁶ == B:", translate(psiA, xy6) == B_IDX)
psi_invB = frozenset(psi_inv(g) for g in B_IDX)
print("ψ⁻¹(B) == x·A:", psi_invB == translate(A_IDX, idx(1, 0)))
print("ψ bijective:", len({psi(g) for g in range(N)}) == N)
print("ψ(dA) == dB:", frozenset(psi(g) for g in dA) == dB)

# Φ on random cycles and boundaries
d2 = np.vstack([CA, CB])
ok_cyc, ok_bd = True, True
x2y6 = idx(2, 6)
for _ in range(30):
    f = (rng.random(N) < 0.3).astype(np.uint8)
    uL, uR = (CA @ f.astype(np.int64)) % 2, (CB @ f.astype(np.int64)) % 2
    # boundary (uL,uR); apply Φ
    vL = apply_map(uR.astype(np.uint8), psi_inv)
    vR_ = apply_map(uL.astype(np.uint8), psi_inv)
    gx, gy = unidx(x2y6)
    vR = np.zeros(N, dtype=np.uint8)
    for g in np.flatnonzero(vR_):
        ggx, ggy = unidx(g)
        vR[idx(ggx + gx, ggy + gy)] = 1
    # cycle check for (vL, vR):
    c1 = (CB @ vL.astype(np.int64) + CA @ vR.astype(np.int64)) % 2
    if c1.any():
        ok_cyc = False
    # boundary check: is (vL,vR) = ∂₂(x·ψ⁻¹f)?
    f2 = apply_map(f, psi_inv)
    f3 = np.zeros(N, dtype=np.uint8)
    for g in np.flatnonzero(f2):
        ggx, ggy = unidx(g)
        f3[idx(ggx + 1, ggy)] = 1
    bL, bR = (CA @ f3.astype(np.int64)) % 2, (CB @ f3.astype(np.int64)) % 2
    if not (np.array_equal(bL, vL) and np.array_equal(bR, vR)):
        ok_bd = False
print("Φ maps boundaries to cycles:", ok_cyc)
print("Φ(∂₂f) == ∂₂(x·ψ⁻¹f):", ok_bd)

# ── P1: Lemma W ───────────────────────────────────────────────────────
print("\n── P1: Lemma W formula ──")
bad = 0
for trial in range(4000):
    sz = rng.integers(1, 9)
    T = list(rng.choice(N, size=sz, replace=False))
    for P_idx, Cm in ((A_IDX, CA), (B_IDX, CB)):
        v = np.zeros(N, dtype=np.uint8)
        v[T] = 1
        w_true = wt((Cm @ v.astype(np.int64)) % 2)
        E, Ncnt = EN_count(P_idx, T)
        if w_true != 3 * len(T) - 2 * E + 4 * Ncnt:
            bad += 1
print(f"Lemma W violations: {bad}/8000")

# ── P2: dA-cluster censuses ───────────────────────────────────────────
print("\n── P2: cluster censuses (translation-reduced, cell 0 ∈ T) ──")
# max E on 4-sets: enumerate 4-sets containing 0
maxE4, maxE5 = 0, 0
quads_dist = Counter()
for rest in combinations(range(1, N), 3):
    T = [0, *rest]
    E, Ncnt = EN_count(A_IDX, T)
    maxE4 = max(maxE4, E)
    quads_dist[(E, Ncnt)] += 1
print(f"4-sets: max E_A = {maxE4}; (E,N) tally "
      f"{dict(sorted(quads_dist.items(), reverse=True)[:6])}")

# 5-sets: E ≥ 6 needs a connected dense cluster; enumerate 5-sets with 0,
# pruning by requiring ≥ 6 edges — use connected growth instead
nbrs = [[idx(gx + dx, gy + dy) for (dx, dy) in
         [unidx(d) for d in dA]] for g in range(N) for (gx, gy) in [unidx(g)]]
found5 = {}
seen = set()
frontier = [frozenset([0])]
allsets = {frozenset([0])}
for size in range(1, 5):
    nxt = set()
    for T in allsets:
        if len(T) != size:
            continue
        cand = set()
        for t in T:
            cand.update(nbrs[t])
        cand -= T
        for c in cand:
            nxt.add(T | {c})
    allsets |= nxt
conn5 = [T for T in allsets if len(T) == 5]
print(f"connected 5-clusters containing 0 (dA-graph): {len(conn5)}")
c5_dist = Counter()
census_15 = []
for T in conn5:
    Tl = sorted(T)
    E, Ncnt = EN_count(A_IDX, Tl)
    maxE5 = max(maxE5, E)
    c5_dist[(E, Ncnt)] += 1
    if E - 2 * Ncnt == 6:
        v = np.zeros(N, dtype=np.uint8)
        v[Tl] = 1
        img = (CA @ v.astype(np.int64)) % 2
        assert wt(img) == 3
        img_idx = frozenset(np.flatnonzero(img))
        isB = any(translate(B_IDX, g) == img_idx for g in range(N))
        census_15.append((frozenset(Tl), img_idx, isB))
print(f"5-clusters: max E_A = {maxE5}; (E,N) tally "
      f"{dict(sorted(c5_dist.items(), reverse=True)[:8])}")
print(f"(1,5) census: {len(census_15)} clusters with E−2N = 6 "
      f"(each 0-anchored; /5 orientations ≈ shapes)")
n_hitB = sum(1 for _, _, isB in census_15 if isB)
print(f"  images that are B-translates: {n_hitB}  (target 0)")
# dedupe up to translation for the shape table
shapes = set()
for T, img, _ in census_15:
    reps = []
    for g in range(N):
        gx, gy = unidx(g)
        Tt = frozenset(idx(tx - gx, ty - gy) for (tx, ty) in
                       (unidx(t) for t in T))
        reps.append(tuple(sorted(Tt)))
    shapes.add(min(reps))
print(f"  distinct shapes up to translation: {len(shapes)}")
for s in sorted(shapes):
    coords = [unidx(g) for g in s]
    E, Ncnt = EN_count(A_IDX, list(s))
    v = np.zeros(N, dtype=np.uint8)
    v[list(s)] = 1
    img = sorted(unidx(g) for g in np.flatnonzero((CA @ v.astype(np.int64)) % 2))
    print(f"    shape {coords} (E={E},N={Ncnt}) → image {img}")

# ── P2': (2,4) census: 4-sets with E−2N ∈ {3,4} ──────────────────────
print("\n── P2': (2,4) census ──")
sigmaB_family = set()   # all σ_B(δ)+t as frozensets, δ ≠ 0
for delta in range(1, N):
    dx, dy = unidx(delta)
    s = set()
    for b in B_IDX:
        bx, by = unidx(b)
        s ^= {idx(bx, by)}
        s ^= {idx(bx + dx, by + dy)}
    fs = frozenset(s)
    for g in range(N):
        sigmaB_family.add(translate(fs, g))
print(f"σ_B family size (all δ, all translates): {len(sigmaB_family)}")
# 4-sets: connected components can be (4) or (3,1); enumerate ALL 4-sets w/ 0
n24_hits = 0
n24_cands = 0
shape24 = Counter()
for rest in combinations(range(1, N), 3):
    T = [0, *rest]
    E, Ncnt = EN_count(A_IDX, T)
    if E - 2 * Ncnt in (3, 4):
        n24_cands += 1
        shape24[(E, Ncnt)] += 1
        v = np.zeros(N, dtype=np.uint8)
        v[T] = 1
        img = (CA @ v.astype(np.int64)) % 2
        img_idx = frozenset(np.flatnonzero(img))
        if img_idx in sigmaB_family:
            n24_hits += 1
print(f"4-sets (0-anchored) with E−2N ∈ {{3,4}}: {n24_cands}, "
      f"by (E,N): {dict(shape24)}")
print(f"  images matching σ_B family: {n24_hits}  (target 0)")

# ── P3: (2,2) size-4 σ-match ─────────────────────────────────────────
print("\n── P3: (2,2) size-4 match ──")
n_match = 0
for dL in dB:
    dLx, dLy = unidx(dL)
    sB = set()
    for b in B_IDX:
        bx, by = unidx(b)
        sB ^= {idx(bx, by)}
        sB ^= {idx(bx + dLx, by + dLy)}
    for dR in dA:
        dRx, dRy = unidx(dR)
        sA = set()
        for a in A_IDX:
            ax, ay = unidx(a)
            sA ^= {idx(ax, ay)}
            sA ^= {idx(ax + dRx, ay + dRy)}
        if any(translate(frozenset(sB), g) == frozenset(sA)
               for g in range(N)):
            n_match += 1
print(f"σ_B(δL) ≅ σ_A(δR) translate matches over δL ∈ dB, δR ∈ dA: "
      f"{n_match}/36  (target 0)")

# ── P4: (3,3) sub-case structure ─────────────────────────────────────
print("\n── P4: (3,3) probe by (E_B, N_B) class of u_L ──")
# classify weight-3 u_L (0-anchored) and test for weight-3 u_R existence
from a21_numeric_session1 import Solver  # reuse GF(2) solver

solA = Solver(CA)
annA = solA.all_kernel_elements()
cls_counter = Counter()
sol_counter = Counter()
for rest in combinations(range(1, N), 2):
    T = [0, *rest]
    E, Ncnt = EN_count(B_IDX, T)
    cls_counter[(E, Ncnt)] += 1
    v = np.zeros(N, dtype=np.uint8)
    v[T] = 1
    c = (CB @ v.astype(np.int64)) % 2
    x0 = solA.solve(c)
    if x0 is None:
        continue
    for z in (annA ^ x0):
        if wt(z) == 3:
            sol_counter[(E, Ncnt)] += 1
print(f"u_L classes (E_B, N_B) tally: {dict(cls_counter)}")
print(f"classes admitting weight-3 u_R: {dict(sol_counter)}")
