"""Deck-nontriviality survey across all published Bravyi et al. BB codes.

For each code: k(cover), then k of the quotient code by every order-2 subgroup
<t> of G (x-half, y-half, diagonal where they exist). Quotient codes are built
directly by coset convolution (no (X,Y) re-presentation needed): qubit/check
index = canonical coset rep; H_X = [A|B], H_Z = [B^T|A^T] as group-algebra
convolution matrices over the quotient.

k(quotient) < k(cover)  <=>  deck-nontrivial on that deck (A12).
Literal presentation only — a jump is decisive code-level only after an orbit
sweep, so jumps are marked "(literal)"; k-preservation is per-presentation
evidence the deck is trivial (sufficient: exhibiting one k-preserving descent).
"""
import numpy as np
from itertools import product

CODES = {
    "[[72,12,6]]   (6,6)":   (6, 6,  [(3,0),(0,1),(0,2)], [(0,3),(1,0),(2,0)]),
    "[[90,8,10]]   (15,3)":  (15, 3, [(9,0),(0,1),(0,2)], [(0,0),(2,0),(7,0)]),
    "[[108,8,10]]  (9,6)":   (9, 6,  [(3,0),(0,1),(0,2)], [(0,3),(1,0),(2,0)]),
    "[[144,12,12]] (12,6)":  (12, 6, [(3,0),(0,1),(0,2)], [(0,3),(1,0),(2,0)]),
    "[[288,12,18]] (12,12)": (12, 12,[(3,0),(0,2),(0,7)], [(0,3),(1,0),(2,0)]),
    "[[360,12,24]] (30,6)":  (30, 6, [(9,0),(0,1),(0,2)], [(0,3),(25,0),(26,0)]),
    "[[756,16,34]] (21,18)": (21, 18,[(3,0),(0,10),(0,17)],[(0,5),(3,0),(19,0)]),
}


def rank_f2(M):
    M = M.copy() % 2
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def code_k(elems, add, neg, A, B):
    """k of the BB code over an abstract abelian group given as element list
    with add/neg, and A,B supports (lists of elements)."""
    idx = {g: i for i, g in enumerate(elems)}
    N = len(elems)

    def conv_mat(P):
        M = np.zeros((N, N), dtype=np.uint8)
        for j, g in enumerate(elems):
            for a in P:
                M[idx[add(a, g)], j] ^= 1
        return M

    MA, MB = conv_mat(A), conv_mat(B)
    # transpose polynomial = negated support
    MAT = conv_mat([neg(a) for a in A])
    MBT = conv_mat([neg(b) for b in B])
    HX = np.concatenate([MA, MB], axis=1)
    HZ = np.concatenate([MBT, MAT], axis=1)
    assert not np.any((HX @ HZ.T) % 2), "CSS violated"
    return 2 * N - rank_f2(HX) - rank_f2(HZ)


def quotient(ell, m, t):
    """Elements/add/neg of (Z_ell x Z_m)/<t>, canonical rep = min of coset."""
    def canon(g):
        h = ((g[0] + t[0]) % ell, (g[1] + t[1]) % m)
        return min(g, h)
    elems = sorted({canon((u, v)) for u in range(ell) for v in range(m)})
    add = lambda a, b: canon(((a[0] + b[0]) % ell, (a[1] + b[1]) % m))
    neg = lambda a: canon(((-a[0]) % ell, (-a[1]) % m))
    return elems, add, neg


for name, (ell, m, A, B) in CODES.items():
    full = [(u, v) for u in range(ell) for v in range(m)]
    add0 = lambda a, b: ((a[0]+b[0]) % ell, (a[1]+b[1]) % m)
    neg0 = lambda a: ((-a[0]) % ell, (-a[1]) % m)
    k = code_k(full, add0, neg0, A, B)
    decks = []
    if ell % 2 == 0:
        decks.append(("x", (ell // 2, 0)))
    if m % 2 == 0:
        decks.append(("y", (0, m // 2)))
    if ell % 2 == 0 and m % 2 == 0:
        decks.append(("xy", (ell // 2, m // 2)))
    if not decks:
        print(f"{name}: k={k}; |G| odd — NOT a Z2-cover of anything (no decks)")
        continue
    verdicts = []
    for dname, t in decks:
        elems, add, neg = quotient(ell, m, t)
        canon = lambda g: min(g, ((g[0]+t[0]) % ell, (g[1]+t[1]) % m))
        Aq = [canon(a) for a in A]
        Bq = [canon(b) for b in B]
        if len(set(Aq)) < len(A) or len(set(Bq)) < len(B):
            verdicts.append(f"{dname}: support-collision")
            continue
        kq = code_k(elems, add, neg, Aq, Bq)
        tag = "R-holds" if kq == k else f"JUMP {kq}->{k} (literal)"
        verdicts.append(f"{dname}: k_base={kq} [{tag}]")
    print(f"{name}: k={k};  " + ";  ".join(verdicts))
