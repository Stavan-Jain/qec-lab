"""Structural tools: Tanner graph, automorphism orbits (nauty), wheel-layer
thickness-2 decomposition + planarity certificate, permutation schedule.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def tanner_edges(H_X: np.ndarray, H_Z: np.ndarray):
    """Vertices: X-checks [0, mX), Z-checks [mX, mX+mZ), qubits [mX+mZ, +n).
    Returns (nv, edges, classes) with classes = (X-check set, Z-check set, qubit set)."""
    mX, n = H_X.shape; mZ = H_Z.shape[0]
    off_q = mX + mZ
    edges = []
    for c, q in zip(*np.nonzero(H_X)):
        edges.append((int(c), off_q + int(q)))
    for c, q in zip(*np.nonzero(H_Z)):
        edges.append((mX + int(c), off_q + int(q)))
    return off_q + n, edges, (set(range(mX)), set(range(mX, off_q)), set(range(off_q, off_q + n)))


def tanner_automorphisms(H_X, H_Z, qubit_blocks=None):
    """Aut group of the Tanner graph with vertex classes X/Z/qubit (optionally
    finer qubit blocks). Returns dict(order=..., orbits on X-checks, ...)."""
    import pynauty
    nv, edges, (Xs, Zs, Qs) = tanner_edges(H_X, H_Z)
    adj = {v: [] for v in range(nv)}
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    coloring = [Xs, Zs, Qs] if qubit_blocks is None else [Xs, Zs] + [set(b) for b in qubit_blocks]
    g = pynauty.Graph(nv, directed=False, adjacency_dict=adj, vertex_coloring=coloring)
    gens, grp1, grp2, orbits, norb = pynauty.autgrp(g)
    orbits = list(orbits)
    x_orbits = len({orbits[v] for v in Xs})
    z_orbits = len({orbits[v] for v in Zs})
    q_orbits = len({orbits[v] for v in Qs})
    order = grp1 * 10 ** grp2
    return {"order": order, "num_generators": len(gens), "orbits_X": x_orbits, "orbits_Z": z_orbits,
            "orbits_qubits": q_orbits, "num_orbits": norb, "orbit_labels": orbits}


def tanner_automorphisms_swap(H_X, H_Z):
    """Aut allowing X↔Z swap (checks colored together) — detects X/Z dualities."""
    import pynauty
    nv, edges, (Xs, Zs, Qs) = tanner_edges(H_X, H_Z)
    adj = {v: [] for v in range(nv)}
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    g = pynauty.Graph(nv, directed=False, adjacency_dict=adj, vertex_coloring=[Xs | Zs, Qs])
    gens, grp1, grp2, orbits, norb = pynauty.autgrp(g)
    return {"order": grp1 * 10 ** grp2, "orbits_checks": len({orbits[v] for v in Xs | Zs})}


@dataclass
class Perms:
    """Permutation description of a two-block code H_X=[A|B], H_Z=[Bᵀ|Aᵀ]
    with A = ΣA_p, B = ΣB_q permutation matrices (m×m each).
    A_p[c] = the L-qubit index hit by X-check c via A_p; likewise B_q[c] (R-qubit)."""
    A: list[np.ndarray]
    B: list[np.ndarray]

    @property
    def m(self):
        return len(self.A[0])

    def inv(self, perm):
        out = np.empty_like(perm); out[perm] = np.arange(len(perm)); return out


def perms_from_matrices(A_mats: list[np.ndarray], B_mats: list[np.ndarray]) -> Perms:
    def to_perm(M):
        M = np.asarray(M)
        assert (M.sum(1) == 1).all() and (M.sum(0) == 1).all(), "not a permutation matrix"
        return np.argmax(M, axis=1)
    return Perms([to_perm(M) for M in A_mats], [to_perm(M) for M in B_mats])


def layer_graphs(P: Perms, a_idx=(1, 2), b_idx=(2,), b_rest=(0, 1), a_rest=(0,)):
    """Build the two BB-style layers as networkx graphs on the 4m Tanner vertices.
    Layer GA = {A_{a_idx}, B_{b_idx}} (2 A-perms + 1 B-perm), GB = the rest.
    Vertex ids: X c -> c ; L l -> m + l ; R r -> 2m + r ; Z z -> 3m + z."""
    import networkx as nx
    m = P.m
    def edges_A(p):  # X-check c — L qubit A_p[c]; Z-check z — R qubit A_p^{-1}[z]
        inv = P.inv(P.A[p])
        return [(c, m + P.A[p][c]) for c in range(m)] + [(3 * m + z, 2 * m + inv[z]) for z in range(m)]
    def edges_B(q):  # X-check c — R qubit B_q[c]; Z-check z — L qubit B_q^{-1}[z]
        inv = P.inv(P.B[q])
        return [(c, 2 * m + P.B[q][c]) for c in range(m)] + [(3 * m + z, m + inv[z]) for z in range(m)]
    GA = nx.Graph(); GB = nx.Graph()
    GA.add_nodes_from(range(4 * m)); GB.add_nodes_from(range(4 * m))
    for p in a_idx: GA.add_edges_from(edges_A(p))
    for q in b_idx: GA.add_edges_from(edges_B(q))
    for p in a_rest: GB.add_edges_from(edges_A(p))
    for q in b_rest: GB.add_edges_from(edges_B(q))
    return GA, GB


def thickness2_certificate(P: Perms):
    """Check all 9 choices of (which A-perm is alone, which B-perm is alone)
    and return the planar ones with component (wheel) statistics."""
    import networkx as nx
    m = P.m
    out = []
    for a_alone in range(3):
        for b_alone in range(3):
            a_idx = tuple(i for i in range(3) if i != a_alone)
            b_rest = tuple(i for i in range(3) if i != b_alone)
            GA, GB = layer_graphs(P, a_idx=a_idx, b_idx=(b_alone,), b_rest=b_rest, a_rest=(a_alone,))
            planA, _ = nx.check_planarity(GA); planB, _ = nx.check_planarity(GB)
            total_edges = GA.number_of_edges() + GB.number_of_edges()
            compsA = [len(c) for c in nx.connected_components(GA)]
            compsB = [len(c) for c in nx.connected_components(GB)]
            out.append({"A_alone": a_alone, "B_alone": b_alone, "planar_A": planA, "planar_B": planB,
                        "edges": total_edges, "degA": sorted(set(d for _, d in GA.degree())),
                        "degB": sorted(set(d for _, d in GB.degree())),
                        "wheelsA": sorted(compsA), "wheelsB": sorted(compsB)})
    return out


def two_block_perms(G, A: dict, a: int, b, code=None) -> "Perms":
    """Permutation decomposition of LP(A, b): H_X = [A | b* I_a].
    A_p = binary expansion of the p-th monomial-permutation part of A (X-check row → L-qubit col),
    B_q = M(conj(b_q)) ⊗ I_a (X-check row → R-qubit col). Verified against the code if given."""
    from .lp import monomial_decomposition, expand, lp_square
    parts = monomial_decomposition(G, A, a)
    A_mats = [expand(G, P, (a, a)) for P in parts]
    bq_parts = [frozenset([g]) for g in sorted(b)]
    B_mats = [expand(G, {(i, i): G.conj(bq) for i in range(a)}, (a, a)) for bq in bq_parts]
    P = perms_from_matrices(A_mats, B_mats)
    if code is None:
        code = lp_square(G, A, a, b)
    m = P.m
    HX = np.zeros_like(code.H_X); HZ = np.zeros_like(code.H_Z)
    for p in P.A:
        HX[np.arange(m), p] ^= 1
        HZ[np.arange(m), m + P.inv(p)] ^= 1     # Aᵀ block: Z-check z ↔ R-qubit A_p^{-1}(z)
    for q in P.B:
        HX[np.arange(m), m + q] ^= 1
        HZ[np.arange(m), P.inv(q)] ^= 1         # Bᵀ block: Z-check z ↔ L-qubit B_q^{-1}(z)
    assert (HX == code.H_X).all(), "A/B permutation parts do not reproduce H_X"
    assert (HZ == code.H_Z).all(), "A/B permutation parts do not reproduce H_Z"
    return P


def _perm_order(p: np.ndarray) -> int:
    n = len(p); seen = np.zeros(n, bool); order = 1
    from math import gcd
    for s in range(n):
        if not seen[s]:
            L = 0; x = s
            while not seen[x]:
                seen[x] = True; x = p[x]; L += 1
            order = order * L // gcd(order, L)
    return order


def toric_layouts(P: Perms, verify: bool = True) -> list[dict]:
    """Generalised Lemma 4 (arXiv:2308.07915).  For i≠j, g≠h let σ = A_j A_i⁻¹ and
    τ = B_g B_h⁻¹ acting on the m L-qubits.  If σ, τ commute and (a,b) ↦ σ^a τ^b ℓ₀ is a
    bijection Z_μ × Z_λ → L (μ = ord σ, λ = ord τ, μλ = m), then placing
      L(σ^a τ^b ℓ₀) at (2a, 2b),  X(A_i⁻¹ ℓ) at (2a+1, 2b),  R(B_g A_i⁻¹ ℓ) at (2a+1, 2b+1),
      Z(B_g ℓ) at (2a, 2b+1)
    makes the Cayley graph of Z_{2μ} × Z_{2λ} a spanning subgraph of the Tanner graph
    (only the commutation of every B_q with every A_p is used).  Returns all such layouts;
    with verify=True each is checked edge-by-edge against the Tanner graph."""
    m = P.m
    out = []
    inv = P.inv
    for i in range(3):
        for j in range(3):
            if i == j: continue
            sigma = P.A[j][inv(P.A[i])]        # σ = A_j ∘ A_i⁻¹  (apply A_i⁻¹ then A_j)
            for g in range(3):
                for h in range(3):
                    if g == h: continue
                    tau = P.B[g][inv(P.B[h])]
                    if not (sigma[tau] == tau[sigma]).all():
                        continue
                    mu = _perm_order(sigma); lam = _perm_order(tau)
                    if mu * lam != m:
                        continue
                    # regularity: orbit map bijective
                    pos = {}
                    ell = 0; ok = True
                    a_pt = 0
                    row_start = 0
                    cur = 0
                    for a in range(mu):
                        x = cur
                        for b in range(lam):
                            if x in pos: ok = False; break
                            pos[x] = (a, b); x = tau[x]
                        if not ok: break
                        cur = sigma[cur]
                    if not ok or len(pos) != m:
                        continue
                    layout = {"i": i, "j": j, "g": g, "h": h, "mu": mu, "lam": lam, "torus": (2 * mu, 2 * lam)}
                    if verify:
                        # build coordinates for all vertices and check the 4 grid edges of each L and R vertex
                        Ai_inv = inv(P.A[i])
                        coordL = {ell: (2 * a, 2 * b) for ell, (a, b) in pos.items()}
                        coordX = {Ai_inv[ell]: (2 * a + 1, 2 * b) for ell, (a, b) in pos.items()}
                        coordR = {P.B[g][Ai_inv[ell]]: (2 * a + 1, 2 * b + 1) for ell, (a, b) in pos.items()}
                        coordZ = {P.B[g][ell]: (2 * a, 2 * b + 1) for ell, (a, b) in pos.items()}
                        assert len(coordX) == m and len(coordR) == m and len(coordZ) == m
                        M1, M2 = 2 * mu, 2 * lam
                        # Tanner adjacency sets
                        XL = set(); XR = set(); ZL = set(); ZR = set()
                        for c in range(m):
                            for p in range(3):
                                XL.add((c, P.A[p][c])); XR.add((c, P.B[p][c]))
                        for z in range(m):
                            for p in range(3):
                                ZR.add((z, inv(P.A[p])[z])); ZL.add((z, inv(P.B[p])[z]))
                        posL = {v: k for k, v in coordL.items()}; posX = {v: k for k, v in coordX.items()}
                        posR = {v: k for k, v in coordR.items()}; posZ = {v: k for k, v in coordZ.items()}
                        good = True
                        for (a, b), ell in posL.items():
                            for da, db in ((1, 0), (-1, 0)):
                                c = posX[((a + da) % M1, b)]
                                good &= (c, ell) in XL
                            for da, db in ((0, 1), (0, -1)):
                                z = posZ[(a, (b + db) % M2)]
                                good &= (z, ell) in ZL
                        for (a, b), r in posR.items():
                            for da in (1, -1):
                                z = posZ[((a + da) % M1, b)]
                                good &= (z, r) in ZR
                            for db in (1, -1):
                                c = posX[(a, (b + db) % M2)]
                                good &= (c, r) in XR
                        layout["verified"] = bool(good)
                        if not good:
                            continue
                    out.append(layout)
    return out
