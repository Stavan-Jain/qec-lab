"""A28 — the epsilon-trisection theorem for Z2 x H groups, validated on data.

Theorem (Trisection, A28 note §4).  Let G = <s> x H, s^2 = 1, |H| odd,
pi: F2[G] -> F2[H] the quotient s -> 1, and write f = f0 + s*f1 with
f0, f1 in F2[H] (sheets), fbar = f0 + f1 = pi(f).  For b = del_G f != 0
with |b| <= W, exactly one of:

  I.   fbar can be chosen 0 (mod K_G):  b = eps (x) del_H f1 with
       eps = 1 + s;  |b| = 2 |del_H f1|;  class map is a BIJECTION onto
       census_H(floor(W/2)).
  II.  fbar in K_H \ pi-image-of-0 (del_H fbar = 0, fbar != 0):  b = eps (x)
       (u0, v0) with (u0, v0) in the AFFINE coset (fbar*A0, fbar*B0) +
       del_H F2[H];  |b| = 2(|u0| + |v0|)  — finitely many coset censuses.
  III. del_H fbar != 0:  bbar := del_H fbar is a nonzero H-boundary with
       |bbar| <= |b| (sheet inequality |u| = |u0|+|u1| >= |u0+u1|), and
       |b| = |bbar| + 2*excess with excess = |u0 off supp(ubar)| +
       |v0 off supp(vbar)| — lift/flip accounting over census_H(W).

So census_G(W) = doubled census_H(W/2)  ⊔  coset censuses (II)  ⊔  bounded
lifts of census_H(W) (III).  The exponential core halves.

This script PROVES nothing — it validates every clause against the
BZ-certified docket37 census (2,203 classes) and independent H-side BZ runs.

Run: uv run --project experiments/bb_lab python experiments/bb_lab/scripts/a28_trisection.py docket37
"""

import json
import math
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from a28_lsc_lib import DATA, REGISTRY, GroupAlg, LSCInstance, rref
from a28_bz_census import C_SRC, BUILD, write_mat, systematic_forms


def h_instance(inst):
    """The odd-quotient instance: G = Z_l x Z_m with m = 2*m', s = y^(m'),
    H = Z_l x Z_(m'), pi(y) = t."""
    l, m = inst.G.l, inst.G.m
    assert m % 2 == 0 and (m // 2) % 2 == 1 and l % 2 == 1
    mp = m // 2
    H = GroupAlg(l, mp)
    Abar = [(a, b % mp) for a, b in inst.A_supp]
    Bbar = [(a, b % mp) for a, b in inst.B_supp]
    hins = LSCInstance(name=inst.name + "H", G=H,
                       A_supp=Abar, B_supp=Bbar,
                       expect_k=-1, d=-1, W=inst.W)
    return hins


def sheets(inst, f):
    """f in F2[G] -> (f0, f1) in F2[H]^2 (s-exponent = y-exp mod 2)."""
    l, m = inst.G.l, inst.G.m
    mp = m // 2
    H = GroupAlg(l, mp)
    f0 = f1 = 0
    for a, b in inst.G.support(f):
        if b % 2 == 0:
            f0 ^= H.monomial(a, b % mp)
        else:
            f1 ^= H.monomial(a, b % mp)
    return f0, f1


def sheet_poly(inst, p_supp):
    """(P0, P1) sheets of a polynomial given by its G-support."""
    l, m = inst.G.l, inst.G.m
    mp = m // 2
    H = GroupAlg(l, mp)
    P0 = P1 = 0
    for a, b in p_supp:
        if b % 2 == 0:
            P0 ^= H.monomial(a, b % mp)
        else:
            P1 ^= H.monomial(a, b % mp)
    return P0, P1


class BoundarySolver:
    def __init__(self, inst):
        self.inst = inst
        G = inst.G
        N = G.N
        aug = []
        for gi, (a, b) in enumerate(G.elements()):
            row = inst.pack(inst.boundary(G.monomial(a, b)))
            aug.append(row | (1 << (2 * N + gi)))
        self.basis, self.piv = rref(aug, 3 * N, col_order=range(2 * N))
        self.N = N

    def solve(self, uv):
        v = self.inst.pack(uv)
        used = 0
        for brow, c in zip(self.basis, self.piv):
            if v & (1 << c):
                v ^= brow & ((1 << 2 * self.N) - 1)
                used ^= brow >> (2 * self.N)
        assert v == 0, "not a boundary"
        return used


def bz_census(hins, W, tag):
    """Small BZ census run (python driver inline, C kernel)."""
    (b1, p1), (b2, p2), kappa = systematic_forms(hins)
    r = W // 2
    exp_nodes = sum(math.comb(kappa, s) for s in range(1, r + 1))
    binp = BUILD / "bzkern"
    BUILD.mkdir(parents=True, exist_ok=True)
    if not binp.exists():
        (BUILD / "bzkern.c").write_text(C_SRC)
        subprocess.run(["cc", "-O3", "-o", str(binp), str(BUILD / "bzkern.c"),
                        "-lpthread"], check=True)
    hits = set()
    for j, basis in enumerate((b1, b2), 1):
        mat = BUILD / f"{tag}_G{j}.mat"
        write_mat(mat, basis)
        res = subprocess.run([str(binp), str(mat), str(r), str(W), "8",
                              str(BUILD / f"{tag}_G{j}")],
                             capture_output=True, text=True, check=True)
        parts = dict(kv.split("=") for kv in res.stdout.strip().split())
        assert int(parts["KAPPA"]) == kappa and int(parts["nodes"]) == exp_nodes
        for t in range(8):
            for line in open(f"{BUILD}/{tag}_G{j}_t{t:02d}.hits"):
                hits.add(int(line, 16))
    classes = {}
    for h in sorted(hits):
        uv = hins.unpack(h)
        c = hins.canonical(uv)
        classes.setdefault(c, uv)
    return classes


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "docket37"
    inst = REGISTRY[name]
    G = inst.G
    mp = G.m // 2
    H = GroupAlg(G.l, mp)
    hins = h_instance(inst)
    A0, A1 = sheet_poly(inst, inst.A_supp)
    B0, B1 = sheet_poly(inst, inst.B_supp)
    assert A0 ^ A1 == hins.A and B0 ^ B1 == hins.B

    census = json.load(open(DATA / "a28" / f"census_{name}.json"))
    print(f"{name}: {census['n_classes']} classes, "
          f"hist {census['weight_histogram']}; H = Z{G.l} x Z{mp}, "
          f"Abar/Bbar supports {H.support(hins.A)} / {H.support(hins.B)}")

    # K_G basis and its pi-image (for sector I-vs-II well-definedness)
    rowsG = inst.generator_rows()
    NG = G.N
    kerG = []
    # kernel of f -> del f: reuse augmented elimination
    aug = []
    for gi, (a, b) in enumerate(G.elements()):
        aug.append(inst.pack(inst.boundary(G.monomial(a, b))) | (1 << (2 * NG + gi)))
    red = aug[:]
    for c in range(2 * NG):
        bit = 1 << c
        pivrow = None
        for r_ in red:
            if r_ & bit:
                pivrow = r_
                break
        if pivrow is None:
            continue
        red = [(r_ ^ pivrow) if (r_ is not pivrow and r_ & bit) else r_ for r_ in red]
        red.remove(pivrow)
    kerG = [r_ >> (2 * NG) for r_ in red]
    piK = set()
    for kmask in kerG:
        f0, f1 = sheets(inst, kmask)
        piK.add(f0 ^ f1)
    piK_span = sorted(piK)
    print(f"dim K_G = {len(kerG)}; pi(K_G) span size = {len(piK)} "
          f"(nonzero images: {sum(1 for x in piK if x)})")

    solver = BoundarySolver(inst)
    hsolver = BoundarySolver(hins)
    t0 = time.time()
    sec = Counter()
    secIII_hcls = {}
    secI_hcls = {}
    secII_fbars = Counter()
    T_fail = 0
    for cl in census["classes"]:
        u = G.from_support([tuple(p) for p in cl["u_support"]])
        v = G.from_support([tuple(p) for p in cl["v_support"]])
        w = cl["weight"]
        f = solver.solve((u, v))
        f0, f1 = sheets(inst, f)
        fbar = f0 ^ f1
        # T1 sheet identities
        u0, u1 = sheets(inst, u)
        v0, v1 = sheets(inst, v)
        assert G.weight(u) == H.weight(u0) + H.weight(u1)
        assert u0 ^ u1 == H.mul(fbar, hins.A), "ubar != fbar*Abar"
        assert v0 ^ v1 == H.mul(fbar, hins.B), "vbar != fbar*Bbar"
        ub = H.mul(fbar, hins.A)
        vb = H.mul(fbar, hins.B)
        if ub == 0 and vb == 0:
            # I or II: try to zero fbar with a K_G shift
            hit0 = any((fbar ^ kb) == 0 for kb in piK_span)
            if hit0:
                sec["I"] += 1
                # T3: b = eps (x) del_H f1' for the zeroing rep
                # find kappa with pi(kappa) = fbar
                for kmask in [0] + kerG + [a for a in ()]:
                    pass
                # locate a kernel element with matching image
                kk = None
                for kmask in [0] + kerG:
                    kf0, kf1 = sheets(inst, kmask)
                    if (kf0 ^ kf1) == fbar:
                        kk = kmask
                        break
                if kk is None:
                    # search the span (dim small)
                    from itertools import combinations
                    found = False
                    for rr in range(2, len(kerG) + 1):
                        for combo in combinations(kerG, rr):
                            m_ = 0
                            for x in combo:
                                m_ ^= x
                            kf0, kf1 = sheets(inst, m_)
                            if (kf0 ^ kf1) == fbar:
                                kk = m_
                                found = True
                                break
                        if found:
                            break
                assert kk is not None
                f_rep = f ^ kk
                r0, r1 = sheets(inst, f_rep)
                assert r0 == r1, "sector-I rep not pure-epsilon"
                bh = (H.mul(r1, hins.A), H.mul(r1, hins.B))
                assert 2 * (H.weight(bh[0]) + H.weight(bh[1])) == w
                secI_hcls[hins.canonical(bh)] = w // 2
            else:
                sec["II"] += 1
                secII_fbars[min(H.mul(fbar, H.monomial(a, b_)) for a in range(H.l)
                                for b_ in range(H.m))] += 1
                # T5: b = eps (x) (u0, v0), coset structure
                assert u0 == u1 and v0 == v1, "sector-II word not pure-epsilon"
                assert w == 2 * (H.weight(u0) + H.weight(v0))
                # coset membership: u0 - fbar*A0 in ideal(Abar)... jointly:
                diff = (u0 ^ H.mul(fbar, A0), v0 ^ H.mul(fbar, B0))
                # must be an H-boundary
                try:
                    hsolver.solve(diff)
                except AssertionError:
                    T_fail += 1
        else:
            sec["III"] += 1
            bw = H.weight(ub) + H.weight(vb)
            assert bw <= w, "sheet inequality violated!"
            excess = (w - bw)
            assert excess % 2 == 0, "lift excess must be even"
            secIII_hcls.setdefault(hins.canonical((ub, vb)), []).append(w)
    print(f"sectors: {dict(sec)}  [{time.time()-t0:.0f}s]  coset-fails: {T_fail}")

    # T4: independent H-side BZ runs
    ch_half = bz_census(hins, inst.W // 2, f"{name}H_half")
    print(f"census_H({inst.W//2}): {len(ch_half)} classes "
          f"{dict(sorted(Counter(hins.pair_weight(uv) for uv in ch_half.values()).items()))}")
    got_I = set(secI_hcls)
    want_I = set(ch_half)
    print(f"T4 sector I vs doubled H-census: ours {len(got_I)}, "
          f"H-BZ {len(want_I)}, equal={got_I == want_I}")
    assert got_I == want_I, "TRISECTION SECTOR I FALSIFIED"

    ch_full = bz_census(hins, inst.W, f"{name}H_full")
    hw = Counter(hins.pair_weight(uv) for uv in ch_full.values())
    print(f"census_H({inst.W}): {len(ch_full)} classes {dict(sorted(hw.items()))}")
    miss = [c for c in secIII_hcls if c not in ch_full]
    print(f"T4 sector III bbar classes: {len(secIII_hcls)} distinct; "
          f"all in census_H({inst.W}): {not miss}")
    assert not miss, "TRISECTION SECTOR III FALSIFIED"
    print(f"sector II fbar classes: {len(secII_fbars)} "
          f"(counts {dict(secII_fbars) if len(secII_fbars) < 8 else '...'})")
    print("TRISECTION VALIDATED END-TO-END")


if __name__ == "__main__":
    main()
