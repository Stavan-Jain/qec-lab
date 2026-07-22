#!/usr/bin/env python
"""A21 probe set 4: closing checks for the (0,0) combos + metric facts.

C1  set comparisons for the (3,3)(0,0) ±combos:
      (i)   B⊖A  vs A⊕B   translate?   (u_L = −A+t, u_R = B+s)
      (ii)  −(A⊖B) vs A⊖B translate?   (u_L = −A+t, u_R = −B+s)
      (iii) A⊕B  vs A⊖B   translate?   (u_L = A+t,  u_R = −B+s)
C2  dA-word metric: BFS lengths of all dB elements and 2dA elements
    (consumed by the (2,4)/(2,2) spread bounds).
C3  span-lemma sub-facts: ker(1+x) in F₂[Z₅] = {0, ring};
    B row-support = {6,10,12}: cyclic span 6, min gap 2.
C4  (2,4) size-6 connected re-verification of the spread kill:
    max image spread over connected (3,0)/(5,1) 4-clusters vs 7.
"""
from __future__ import annotations

import numpy as np
from itertools import combinations
from collections import Counter, deque

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
dA = frozenset(sub(s, t) for s in A_IDX for t in A_IDX if s != t)
dB = frozenset(sub(s, t) for s in B_IDX for t in B_IDX if s != t)


def conv_sets(P, T):
    out = set()
    for p in P:
        for t in T:
            out ^= {add(p, t)}
    return frozenset(out)


def translates_of(S):
    return {frozenset(add(s, g) for s in S) for g in range(N)}


AplusB = conv_sets(A_IDX, B_IDX)
AmB = frozenset(sub(a, b) for a in A_IDX for b in B_IDX)
BmA = frozenset(sub(b, a) for b in B_IDX for a in A_IDX)
negAmB = frozenset(sub(0, v) for v in AmB)

print("── C1: set comparisons ──")
print("(i)   B⊖A ∈ translates(A⊕B):", BmA in translates_of(AplusB))
print("(ii)  −(A⊖B) ∈ translates(A⊖B):", negAmB in translates_of(AmB))
print("(iii) A⊕B ∈ translates(A⊖B):", AplusB in translates_of(AmB))

print("\n── C2: dA-word metric lengths ──")
dist = {0: 0}
q = deque([0])
while q:
    g = q.popleft()
    for d in dA:
        h = add(g, d)
        if h not in dist:
            dist[h] = dist[g] + 1
            q.append(h)
for lab, S in (("dB", sorted(dB)), ("2dA", sorted(add(d, d) for d in dA))):
    print(f"  {lab}: " + ", ".join(
        f"{unidx(g)}→{dist[g]}" for g in S))
print(f"  eccentricity of 0: {max(dist.values())}")

print("\n── C3: span-lemma sub-facts ──")
ring_ok = []
for mask in range(32):
    p = [(mask >> i) & 1 for i in range(5)]
    q_ = [(p[i] ^ p[(i - 1) % 5]) for i in range(5)]
    if not any(q_):
        ring_ok.append(mask)
print(f"  ker(1+x) in F₂[Z₅]: masks {ring_ok} (expect [0, 31] = {{0, ring}})")
rows = sorted({gy for (gx, gy) in B_SUPP})
gaps = [(rows[(i + 1) % 3] - rows[i]) % 15 for i in range(3)]
print(f"  B rows {rows}, cyclic gaps {gaps}, span = {15 - max(gaps)}, "
      f"min gap = {min(gaps)}")

print("\n── C4: (2,4) spread re-verification ──")
nbr = {g: [add(g, d) for d in dA] for g in range(N)}
allsets = {frozenset([0])}
for size in range(1, 4):
    nxt = set()
    for T in allsets:
        if len(T) != size:
            continue
        for t in T:
            for c in nbr[t]:
                if c not in T:
                    nxt.add(T | {c})
    allsets |= nxt


def EN(P, T):
    dP = frozenset(sub(s, t) for s in P for t in P if s != t)
    Tl = sorted(T)
    E = sum(1 for i in range(len(Tl)) for j in range(i + 1, len(Tl))
            if sub(Tl[i], Tl[j]) in dP)
    Ncnt = sum(1 for c in range(N) if all(sub(c, p) in T for p in P))
    return E, Ncnt


max_spread = {}
for T in allsets:
    if len(T) != 4:
        continue
    E, Ncnt = EN(A_IDX, T)
    if (E, Ncnt) not in [(3, 0), (5, 1), (4, 0)]:
        continue
    img = conv_sets(A_IDX, sorted(T))
    sp = max(dist[sub(u, v)] for u in img for v in img)
    key = (E, Ncnt)
    max_spread[key] = max(max_spread.get(key, 0), sp)
print(f"  max image spread by connected class: {max_spread}  "
      f"(kill needs < 7 for size-6, and no ±(1,6)-pair for size-4)")
# also verify every full B-translate and every σ_B(β∈dB) has a pair at
# dA-distance exactly the (1,6) length
l16 = dist[idx(1, 6)]
print(f"  |(1,6)| = {l16}; B-internal pair lengths: "
      f"{sorted(dist[sub(u, v)] for u, v in combinations(B_IDX, 2))}")
for beta in sorted(dB):
    sig = conv_sets(B_IDX, [0, beta])
    lens = sorted(dist[sub(u, v)] for u, v in combinations(sorted(sig), 2))
    print(f"  σ_B({unidx(beta)}) pair lengths: {lens}")
