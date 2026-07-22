#!/usr/bin/env python
"""A21 probe set 3: the finite tables consumed by the analytic split map.

T0  (3,3)(0,0) sumset-rigidity localization table: u_L ⊆ A⊖B, gates.
T1  σ_B(β) pair-value table (β ∈ dB): which pair differences ∈ dA?
    (kills (3,3)(1,0)×(1,0): A-translate must fit in σ_B(β)+t ⊔ B+t+s.)
T2  {0}∪dB pair-value table: which differences ∈ dA?  (kills (3,1)×(1,0).)
T3  A-cherry image table: #dB-pairs in c = A⋆(A-cherry) (kills (2,0)×(2,0)).
T4  (2,4) size-4 residual: (4,0)-shape images vs σ_B(δ), δ ∈ (1,6)-class.
T5  checks: ψ² = id, ψ(dB) = dA, A⊕B aperiodic, ord(1+x) = 15 on ker ε,
    weight ladder Σ_j |(1+x)^j p| = 40, d(A⊖B) multiplicity profile.
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


def add(g, h):
    gx, gy = unidx(g)
    hx, hy = unidx(h)
    return idx(gx + hx, gy + hy)


def sub(g, h):
    gx, gy = unidx(g)
    hx, hy = unidx(h)
    return idx(gx - hx, gy - hy)


A_IDX = [idx(*p) for p in A_SUPP]
B_IDX = [idx(*p) for p in B_SUPP]


def diffset(S):
    return frozenset(sub(s, t) for s in S for t in S if s != t)


dA, dB = diffset(A_IDX), diffset(B_IDX)


def conv_sets(P, T):
    out = set()
    for p in P:
        for t in T:
            c = add(p, t)
            out ^= {c}
    return frozenset(out)


def pretty(g):
    return str(unidx(g))


def prettys(S):
    return "{" + ", ".join(pretty(g) for g in sorted(S)) + "}"


# ── T5 first (basic identities the others rely on) ────────────────────
print("── T5: base identities ──")


def psi(g):
    gx, gy = unidx(g)
    return idx(gx, 6 * gx + 4 * gy)


print("ψ² = id:", all(psi(psi(g)) == g for g in range(N)))
print("ψ(dB) = dA:", frozenset(psi(g) for g in dB) == dA)
AplusB = conv_sets(A_IDX, B_IDX)
print(f"|A⊕B| = {len(AplusB)}; periods:",
      [pretty(p) for p in range(1, N)
       if frozenset(add(c, p) for c in AplusB) == AplusB] or "none")
AmB = frozenset(sub(a, b) for a in A_IDX for b in B_IDX)
print(f"|A⊖B| = {len(AmB)} (9 ⟺ D2)")
dAmB = Counter(sub(u, v) for u in AmB for v in AmB if u != v)
prof = Counter(dAmB.values())
print(f"d(A⊖B) multiplicity profile: {dict(sorted(prof.items()))}")
hi = sorted(g for g, m in dAmB.items() if m >= 3)
print(f"  mult ≥ 3 values: {len(hi)}; dA∪dB ⊆ them: "
      f"{set(dA | dB) <= set(hi)}; extras: "
      f"{[pretty(g) for g in hi if g not in (dA | dB)]}")

# (1+x) in F2[Z5]: order on ker ε and weight ladder
poly = np.zeros(5, dtype=np.uint8)


def mulx(p, k=1):  # multiply by x^k in F2[Z5]
    return np.roll(p, k)


def times_1px(p):
    return p ^ np.roll(p, 1)


p0 = np.zeros(5, dtype=np.uint8)
p0[0] = p0[1] = 1  # 1 + x
cur = p0.copy()
order = None
for j in range(1, 20):
    cur = times_1px(cur) if False else cur
    break
# order of multiplication-by-(1+x) on ker ε: iterate powers of (1+x)
cur = p0.copy()
for j in range(2, 40):
    cur = cur ^ np.roll(cur, 1)  # cur *= (1+x)
    if np.array_equal(cur, p0):
        order = j - 1
        break
print(f"multiplicative order of (1+x) in F₂[Z₅]* mod... "
      f"(1+x)^k = (1+x) at k = {order} ⟹ ord = {order - 1 if order else '?'}")
# weight ladders: for every even-weight p ≠ 0, Σ_{j=0}^{14} |(1+x)^j p|
tot_ok = True
for mask in range(1, 32):
    p = np.array([(mask >> i) & 1 for i in range(5)], dtype=np.uint8)
    if int(p.sum()) % 2:
        continue
    s, cur = 0, p.copy()
    for j in range(15):
        s += int(cur.sum())
        cur = cur ^ np.roll(cur, 1)
    if not np.array_equal(cur, p):
        tot_ok = False
        print("  (1+x)^15 ≠ 1 on ker ε!!")
    if s != 40:
        tot_ok = False
        print(f"  ladder ≠ 40 for p = {p}: {s}")
print(f"weight ladders all 40 on ker ε ∖ 0: {tot_ok}")

# ── T0: (3,3)(0,0) localization table ────────────────────────────────
print("\n── T0: (0,0) sumset-rigidity table (u_L ⊆ A⊖B, 0 ∈ u_R) ──")
BmA = frozenset(sub(b, a) for b in B_IDX for a in A_IDX)
n_all = n_d1 = n_int = n_final = 0
survivors = []
for u_L in combinations(sorted(AmB), 3):
    n_all += 1
    if diffset(u_L) & dB:
        continue                      # gate 1: direct B-sumset
    n_d1 += 1
    inter = None
    for t in u_L:
        s = frozenset(add(t, v) for v in BmA)
        inter = s if inter is None else (inter & s)
    if len(inter) < 3 or 0 not in inter:
        continue                      # gate 2: room for u_R
    n_int += 1
    for rest in combinations(sorted(inter - {0}), 2):
        u_R = (0, *rest)
        if diffset(u_R) & dA:
            continue                  # gate 3: direct A-sumset
        if conv_sets(B_IDX, u_L) == conv_sets(A_IDX, u_R):
            n_final += 1
            survivors.append((u_L, u_R))
print(f"3-subsets of A⊖B: {n_all}; pass d(u_L)∩dB=∅: {n_d1}; "
      f"pass |∩(t+B⊖A)| ≥ 3: {n_int}; full-equation survivors: {n_final}")
for u_L, u_R in survivors:
    # is it a generator? u_L = A+g, u_R = B+g
    gens = [g for g in range(N)
            if frozenset(add(a, g) for a in A_IDX) == frozenset(u_L)
            and frozenset(add(b, g) for b in B_IDX) == frozenset(u_R)]
    print(f"  u_L {prettys(u_L)}, u_R {prettys(u_R)} "
          f"→ generator g = {[pretty(g) for g in gens]}")

# gate-2 detail for the writeup: |∩| distribution over gate-1 survivors
dist = Counter()
for u_L in combinations(sorted(AmB), 3):
    if diffset(u_L) & dB:
        continue
    inter = None
    for t in u_L:
        s = frozenset(add(t, v) for v in BmA)
        inter = s if inter is None else (inter & s)
    dist[len(inter)] += 1
print(f"|∩(t+B⊖A)| distribution over gate-1 survivors: {dict(sorted(dist.items()))}")

# ── T1: σ_B(β) pair values vs dA ─────────────────────────────────────
print("\n── T1: σ_B(β) pair-difference values ∩ dA (β ∈ dB) ──")
for beta in sorted(dB):
    sig = conv_sets(B_IDX, [0, beta])
    assert len(sig) == 4
    pv = Counter(sub(u, v) for u in sig for v in sig if u != v)
    in_dA = [pretty(g) for g in pv if g in dA]
    in_dB = [pretty(g) for g in pv if g in dB]
    print(f"  β = {pretty(beta)}: pairs ∩ dA = {in_dA or '∅'}, "
          f"pairs ∩ dB = {in_dB}")

# mirror (for (4,2)/(2,2) writeups): σ_A(α) pair values ∩ dB
print("   mirror σ_A(α) pairs ∩ dB (α ∈ dA):")
for alpha in sorted(dA):
    sig = conv_sets(A_IDX, [0, alpha])
    pv = Counter(sub(u, v) for u in sig for v in sig if u != v)
    in_dB = [pretty(g) for g in pv if g in dB]
    print(f"  α = {pretty(alpha)}: pairs ∩ dB = {in_dB or '∅'}")

# ── T2: {0}∪dB pair values vs dA ─────────────────────────────────────
print("\n── T2: S_B = {0}∪dB pair differences ∩ dA ──")
SB = frozenset({0}) | dB
pv = Counter(sub(u, v) for u in SB for v in SB if u != v)
in_dA = [pretty(g) for g in pv if g in dA]
print(f"  pairs of S_B with difference ∈ dA: {in_dA or '∅'}")
SA = frozenset({0}) | dA
pv = Counter(sub(u, v) for u in SA for v in SA if u != v)
in_dB = [pretty(g) for g in pv if g in dB]
print(f"  pairs of S_A with difference ∈ dB: {in_dB or '∅'}")

# ── T3: A-cherry images: #dB-pairs ───────────────────────────────────
print("\n── T3: A-cherries (α₁,α₂ ∈ dA, α₁+α₂ ∉ dA ∪ {0}): "
      "dB-pairs in image ──")
tally3 = Counter()
worst = []
for a1 in sorted(dA):
    for a2 in sorted(dA):
        s = add(a1, a2)
        if s == 0 or s in dA:
            continue
        T = [0, a1, add(a1, a2)]
        if len(set(T)) < 3:
            continue
        c = conv_sets(A_IDX, T)
        assert len(c) == 5, (pretty(a1), pretty(a2), len(c))
        ndB = sum(1 for u in c for v in c if u < v and sub(u, v) in dB
                  or u < v and sub(v, u) in dB)
        # count unordered pairs with either orientation in dB
        ndB = sum(1 for u, v in combinations(sorted(c), 2)
                  if sub(u, v) in dB or sub(v, u) in dB)
        tally3[ndB] += 1
        if ndB >= 2:
            worst.append((pretty(a1), pretty(a2), ndB))
print(f"  #dB-pairs distribution over A-cherries: {dict(sorted(tally3.items()))}")
print(f"  cherries with ≥ 2 dB-pairs (would survive the pair-count gate): "
      f"{worst or 'NONE'}")
# mirror: B-cherries: #dA-pairs in image (for the same kill stated B-side)
tallyB = Counter()
for b1 in sorted(dB):
    for b2 in sorted(dB):
        s = add(b1, b2)
        if s == 0 or s in dB:
            continue
        T = [0, b1, add(b1, b2)]
        if len(set(T)) < 3:
            continue
        c = conv_sets(B_IDX, T)
        if len(c) != 5:
            continue
        ndA = sum(1 for u, v in combinations(sorted(c), 2)
                  if sub(u, v) in dA or sub(v, u) in dA)
        tallyB[ndA] += 1
print(f"  mirror (B-cherries, #dA-pairs): {dict(sorted(tallyB.items()))}")

# ── T4: (2,4) size-4 residual ────────────────────────────────────────
print("\n── T4: (2,4) size-4 residual: (4,0) 4-cluster images vs σ_B ──")
# enumerate connected 4-sets containing 0 with E=4, N=0 (dA-graph)
nbr = {g: [add(g, d) for d in dA] for g in range(N)}
sets4 = set()
grow = {frozenset([0])}
allsets = set(grow)
for size in range(1, 4):
    nxt = set()
    for T in allsets:
        if len(T) != size:
            continue
        cand = set()
        for t in T:
            cand.update(nbr[t])
        for c in cand - T:
            nxt.add(T | {c})
    allsets |= nxt


def EN(P, T):
    dP = diffset(P)
    Tl = sorted(T)
    E = sum(1 for i in range(len(Tl)) for j in range(i + 1, len(Tl))
            if sub(Tl[i], Tl[j]) in dP)
    Ncnt = sum(1 for c in range(N)
               if all(sub(c, p) in T for p in P))
    return E, Ncnt


shapes40 = []
seen = set()
for T in allsets:
    if len(T) != 4:
        continue
    E, Ncnt = EN(A_IDX, T)
    if (E, Ncnt) != (4, 0):
        continue
    canon = min(tuple(sorted(sub(t, t0) for t in T)) for t0 in T)
    if canon in seen:
        continue
    seen.add(canon)
    shapes40.append(canon)
print(f"(4,0) shapes up to translation: {len(shapes40)}")
sigB_dB = {}
for beta in sorted(dB):
    sigB_dB[beta] = conv_sets(B_IDX, [0, beta])
for sh in shapes40:
    c = conv_sets(A_IDX, sh)
    pv = [g for u, v in combinations(sorted(c), 2)
          for g in (sub(u, v),) if g in dB or sub(v, u) in dB]
    match = [pretty(beta) for beta, sig in sigB_dB.items()
             if any(frozenset(add(x, t) for x in c) == sig
                    for t in range(N))]
    print(f"  shape {[pretty(t) for t in sh]} → image "
          f"{[pretty(g) for g in sorted(c)]}; dB-pairs {len(pv)}; "
          f"σ_B matches {match or 'NONE'}")
