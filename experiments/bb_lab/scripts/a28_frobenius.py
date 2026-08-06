"""A28 — Frobenius compression of the census (odd |G|).

Lemma F (note §5).  For odd |G|, squaring sq: F2[G] -> F2[G], f -> f*f is
a WEIGHT-PRESERVING bijection: in char 2, (sum x_i)^2 = sum x_i^2, and
g -> g^2 is a group automorphism (|G| odd), so sq permutes monomials.
It is Frobenius-semilinear for the ring structure: (fg)^2 = f^2 g^2.
Hence sq maps C(A,B) to C(A^2,B^2) = C(sq(A), sq(B)) — and when (A,B) is
FIXED by some power sq^j composed with a group automorphism fixing the
pair, the census inherits that symmetry: classes come in orbits.

For f2a6: check whether some tau in Aut(Z5 x Z15) has (tau(A^2), tau(B^2))
= (A, B) up to translation; then sq*tau acts on the 113 classes — report
its orbit structure (compression factor for any future certified sweep).

Run: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_frobenius.py
"""

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import REGISTRY, load_f2a6_census


def main():
    inst = REGISTRY["f2a6"]
    G = inst.G

    def sq(u):
        out = 0
        for a, b in G.support(u):
            out ^= G.monomial(2 * a, 2 * b)
        return out

    # Lemma F sanity: weight preserved, multiplicative, on random elements
    import random
    rng = random.Random(7)
    for _ in range(50):
        f = rng.getrandbits(G.N)
        g = rng.getrandbits(G.N)
        assert G.weight(sq(f)) == G.weight(f)
        assert sq(G.mul(f, g)) == G.mul(sq(f), sq(g))
    print("Lemma F: sq is weight-preserving and multiplicative ✓")

    # does sq fix (A,B) up to translation + automorphism?
    # automorphisms of Z5 x Z15: (a,b) -> (p*a + 3*r*b?, ...) — for a first
    # pass, try the DIAGONAL monomial maps (a,b) -> (p*a, q*b),
    # gcd(p,5)=1, gcd(q,15)=1, with translations.
    A, B = inst.A, inst.B
    hits = []
    for p in range(1, 5):
        for q in range(1, 15):
            if q % 3 == 0 or q % 5 == 0:
                continue

            def tau(u, p=p, q=q):
                out = 0
                for a, b in G.support(u):
                    out ^= G.monomial(p * a, q * b)
                return out

            TA, TB = tau(sq(A)), tau(sq(B))
            # match up to (independent? no — same translation for both:
            # a boundary map conjugation needs one shift g with
            # TA = x^g A, TB = x^g B... allow independent shifts? The census
            # symmetry needs del-intertwining: tau(sq(del f)) =
            # del'(tau(sq f)) with (A', B') = (TA, TB); classes map to
            # classes of (A', B'); equality of censuses needs (A', B') =
            # translate of (A, B) — SAME translate not required blockwise?
            # u-block: f*A -> translate_g(f*A) works with A' = g*A; blocks
            # can carry different g only if we also translate f — one g per
            # block is fine: (A', B') = (g1*A, g2*B) gives a census bijection
            # only when g1 = g2 (single f).  Actually (g1 A, g2 B) =
            # del'(f) with del' = (·g1 A, ·g2 B): canonical classes of
            # (u,v) -> (g1^-1 u, g2^-1 v)?? — that map does NOT respect
            # joint translation classes unless g1 = g2.  Require g1 == g2.
            for g in G.elements():
                if G.translate(A, *g) == TA and G.translate(B, *g) == TB:
                    hits.append((p, q, g))
    print(f"monomial symmetries with sq*tau(A,B) = translate(A,B): {hits}")
    if not hits:
        print("-> sq does not normalize this (A,B) among diagonal maps; "
              "census-level Frobenius symmetry NOT available here")
        # still: sq gives a census BIJECTION to the squared code — check it
        insq_supp_A = [ (2*a % 5, 2*b % 15) for a,b in inst.A_supp ]
        print(f"   (A^2 supp = {sorted(insq_supp_A)} vs A supp {sorted(inst.A_supp)})")
        return

    p, q, g = hits[0]

    def act(uv):
        u, v = uv
        def one(u):
            out = 0
            for a, b in G.support(u):
                out ^= G.monomial(p * 2 * a, q * 2 * b)
            return out
        gu = G.translate(one(u), -g[0], -g[1])
        gv = G.translate(one(v), -g[0], -g[1])
        return (gu, gv)

    census = load_f2a6_census()
    reps = []
    for row in census:
        u = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 0])
        v = G.from_support([(x, y) for blk, x, y in row["b_support"] if blk == 1])
        reps.append(inst.canonical((u, v)))
    repset = set(reps)
    orbit_sizes = Counter()
    seen = set()
    for r0 in reps:
        if r0 in seen:
            continue
        orb = set()
        cur = r0
        while cur not in orb:
            orb.add(cur)
            uv = inst.unpack(cur)
            cur = inst.canonical(act(uv))
            assert cur in repset, "census not closed under the symmetry!"
        orbit_sizes[len(orb)] += 1
        seen |= orb
    print(f"census closed under sq*tau ✓; orbit sizes {dict(orbit_sizes)}; "
          f"compression 113 -> {sum(orbit_sizes.values())}")


if __name__ == "__main__":
    main()
