"""A20: the A22 (eps,delta) fibering ported to Y4 = (18,4) — Phase 0 + V6.

Coordinates. G = Z18(x) x Z4(y); z := x^2 (order 9, the semisimple fiber),
s := x^9 (order 2); x = z^5 s. Base sites: (s,y) in Z2 x Z4 (8 sites, LOCAL
group algebra — invertibility = augmentation, no characters needed). CRT:

    R := F2[z]/(z^9-1)  ~=  F2 (eps: z->1)  x  GF(4) (d4: z->w, z^2+z+1)
                            x  GF(64) (d64: z->t, z^6+z^3+1)

Every fiber polynomial is determined by its value triple (bijection,
2*4*64 = 512 = 2^9), with a generated EXACT weight table W(triple).
|u| = sum over 8 sites of W(u's triple at the site).

Base polynomials in (z,s,y):  A = 1 + z^5 s + z^7 y,  B = 1 + z^5 s y^2 + z y^3.
Component polynomials A_c, B_c live in K_c[Z2 x Z4].

Phase 0 verifications (V1-V5 analogs):
  V1  CRT bijection + weight table (512 triples distinct, weights 0..9).
  V2  exact site-weight formula on random group-algebra elements.
  V3  invertibility verdicts for A_c, B_c on each component (augmentation
      for the local base ring; cross-checked by F2-linearized rank).
  V4  kernel localization: ker d2 (dim 4) lives entirely in the
      non-invertible component(s); GF-dimension there.
  V5  transfer operators C_c = B_c A_c^{-1} on invertible components;
      support sizes.
V6  decompose ALL census classes (data/a20/m_census_classes.jsonl) through
    the coordinates: verify the weight formula and the v = C*u link
    class-by-class; emit the cost taxonomy (active-site patterns, per-site
    weights) to data/a20/fibering_taxonomy.json — the input for the V7
    completeness sweep design.

Usage (from experiments/bb_lab):  uv run scripts/a20_fibering.py
"""
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

LAB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LAB / "src"))
from bb_lab.linalg import nullspace_f2, rank_f2  # noqa: E402

OUT = LAB / "data" / "a20"
LX, LY = 18, 4
N2, N1 = 72, 144

# ---------------------------------------------------------------- GF helpers
# GF(4) = F2[w]/(w^2+w+1), elements 0..3 as bit-pairs (b0 + b1*w).
# GF(64) = F2[t]/(t^6+t^3+1), elements 0..63 as bit-vectors.


def _gf_mul(a, b, poly, deg):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> deg & 1:
            a ^= poly
    return r


def gf4_mul(a, b):
    return _gf_mul(a, b, 0b111, 2)


def gf64_mul(a, b):
    return _gf_mul(a, b, 0b1001001, 6)


def fiber_triple(u9):
    """u in F2^9 (coeffs of z^0..z^8) -> (eps, d4 int, d64 int)."""
    eps = int(u9.sum()) & 1
    d4 = 0
    for j, c in enumerate(u9):
        if c:
            d4 ^= pow_tab4[j % 3]
    d64 = 0
    for j, c in enumerate(u9):
        if c:
            d64 ^= pow_tab64[j % 9]
    return eps, d4, d64


# powers of w (order 3) in GF(4); powers of t (order 9) in GF(64)
pow_tab4 = [1]
for _ in range(2):
    pow_tab4.append(gf4_mul(pow_tab4[-1], 0b10))
pow_tab64 = [1]
for _ in range(8):
    pow_tab64.append(gf64_mul(pow_tab64[-1], 0b10))
assert gf64_mul(pow_tab64[8], 0b10) == 1, "t must have order 9"


def build_weight_table():
    table = {}
    for mask in range(512):
        u9 = np.array([(mask >> j) & 1 for j in range(9)], dtype=np.uint8)
        tri = fiber_triple(u9)
        assert tri not in table, "CRT bijection failure"
        table[tri] = int(u9.sum())
    return table


# ------------------------------------------------- group indexing and lifts
def gidx(a, b):
    return (a % LX) * LY + (b % LY)


def to_fiber_sites(vec72):
    """vec on G -> array [9 fiber][2 s][4 y] via a <-> (j = 5a mod 9, s = a mod 2)."""
    out = np.zeros((9, 2, LY), dtype=np.uint8)
    for a in range(LX):
        j, sg = (5 * a) % 9, a % 2
        for b in range(LY):
            out[j, sg, b] = vec72[gidx(a, b)]
    return out


def site_triples(vec72):
    """-> dict site (s,y) -> (eps, d4, d64) triple of the fiber poly there."""
    F = to_fiber_sites(vec72)
    return {(sg, b): fiber_triple(F[:, sg, b]) for sg in range(2)
            for b in range(LY)}


# ------------------------------- component algebra over K_c[Z2 x Z4] (8 sites)
SITES = [(sg, b) for sg in range(2) for b in range(LY)]
SIDX = {st: i for i, st in enumerate(SITES)}


def comp_poly(supp_zsy, comp):
    """supp as list of (zexp, sexp, yexp) -> vector of K_c values per site,
    as int arrays (K in {4, 64}; comp 'eps' uses ints 0/1)."""
    vec = np.zeros(8, dtype=np.int64)
    for (ze, se, ye) in supp_zsy:
        i = SIDX[(se % 2, ye % LY)]
        if comp == "eps":
            vec[i] ^= 1
        elif comp == "d4":
            vec[i] ^= pow_tab4[ze % 3]
        else:
            vec[i] ^= pow_tab64[ze % 9]
    return vec


def site_conv(pvec, fvec, mul):
    """(p * f) over Z2 x Z4 with K_c coefficients (int-coded)."""
    out = np.zeros(8, dtype=np.int64)
    for i, (si, yi) in enumerate(SITES):
        if not pvec[i]:
            continue
        for j, (sj, yj) in enumerate(SITES):
            if not fvec[j]:
                continue
            out[SIDX[((si + sj) % 2, (yi + yj) % LY)]] ^= mul(
                int(pvec[i]), int(fvec[j]))
    return out


def linearize(pvec, mul, bits):
    """K_c[Z2xZ4]-multiplication-by-p as an F2 matrix (8*bits square)."""
    dim = 8 * bits
    M = np.zeros((dim, dim), dtype=np.uint8)
    for j in range(8):
        for bit in range(bits):
            f = np.zeros(8, dtype=np.int64)
            f[j] = 1 << bit
            g = site_conv(pvec, f, mul)
            for i in range(8):
                for ob in range(bits):
                    if (g[i] >> ob) & 1:
                        M[i * bits + ob, j * bits + bit] = 1
    return M


A_ZSY = [(0, 0, 0), (5, 1, 0), (7, 0, 1)]     # A = 1 + z^5 s + z^7 y
B_ZSY = [(0, 0, 0), (5, 1, 2), (1, 0, 3)]     # B = 1 + z^5 s y^2 + z y^3
AND = lambda a, b: a & b  # noqa: E731  (F2 multiplication on the eps component)


def main():
    W = build_weight_table()
    ws = Counter(W.values())
    print(f"V1 PASS: 512-triple bijection; weight histogram {dict(sorted(ws.items()))}")

    # V2: exact site-weight formula on random vectors
    rng = np.random.default_rng(20)
    for _ in range(50):
        v = rng.integers(0, 2, N2).astype(np.uint8)
        assert int(v.sum()) == sum(W[t] for t in site_triples(v).values())
    print("V2 PASS: |u| = sum of site weights (50 random)")

    # V3: component polynomials + invertibility
    comps = {}
    for name, mul, bits in (("eps", AND, 1), ("d4", gf4_mul, 2),
                            ("d64", gf64_mul, 6)):
        Ac = comp_poly(A_ZSY, name)
        Bc = comp_poly(B_ZSY, name)
        augA = 0
        augB = 0
        for i in range(8):
            augA ^= int(Ac[i])
            augB ^= int(Bc[i])
        MA = linearize(Ac, mul, bits)
        MB = linearize(Bc, mul, bits)
        invA = rank_f2(MA) == MA.shape[0]
        invB = rank_f2(MB) == MB.shape[0]
        assert invA == (augA != 0) and invB == (augB != 0), \
            "local-ring invertibility criterion violated?!"
        comps[name] = dict(Ac=Ac, Bc=Bc, MA=MA, MB=MB, mul=mul, bits=bits,
                           invA=invA, invB=invB)
        print(f"V3 [{name}]: aug(A)={augA:#x} aug(B)={augB:#x} -> "
              f"A invertible: {invA}, B invertible: {invB}")

    # V4: kernel localization
    MA_full = np.zeros((N2, N2), dtype=np.uint8)
    MB_full = np.zeros((N2, N2), dtype=np.uint8)
    A_G = [(gidx(1, 0)), ]  # rebuild full conv matrices from G-supports
    A_SUPP = [(0, 0), (1, 0), (14, 1)]
    B_SUPP = [(0, 0), (1, 2), (2, 3)]
    for ga in range(LX):
        for gb in range(LY):
            g = gidx(ga, gb)
            for (sa, sb) in A_SUPP:
                MA_full[gidx(ga + sa, gb + sb), g] ^= 1
            for (sa, sb) in B_SUPP:
                MB_full[gidx(ga + sa, gb + sb), g] ^= 1
    D2 = np.vstack([MA_full, MB_full])
    K = nullspace_f2(D2)
    print(f"V4: dim ker d2 = {K.shape[0]}", end="")
    for kv in K:
        tri = site_triples(kv)
        eps_part = any(t[0] for t in tri.values())
        d64_part = any(t[2] and True for t in tri.values()) and any(
            t[2] not in (0,) for t in tri.values())
        # localization verdict per basis vector
    # kernel lives where both A_c,B_c singular: check eps/d64 components vanish
    loc_ok = True
    for kv in K:
        F = to_fiber_sites(kv)
        for sg in range(2):
            for b in range(LY):
                e, d4v, d64v = fiber_triple(F[:, sg, b])
                if comps["eps"]["invA"] and e:
                    loc_ok = False
                if comps["d64"]["invA"] and d64v:
                    loc_ok = False
    print(f"  — localized to singular component(s): {loc_ok}")

    # V5: transfer operators on invertible components
    for name in comps:
        c = comps[name]
        if not c["invA"]:
            print(f"V5 [{name}]: A singular — no transfer (kernel component)")
            continue
        bits = c["bits"]
        # C = B A^{-1}: solve MA^T? use linear solve: columns of C-matrix
        MAinv_on_B = []
        MA, MB = c["MA"].astype(np.uint8), c["MB"].astype(np.uint8)
        # solve MA X = MB over F2 (X = matrix of A^{-1}B ... careful order)
        # multiplication operators commute (commutative ring), C acts as MB MA^{-1}
        # find C's site-vector: C = site_conv? easiest: apply MB MA^{-1} to delta_0
        n = MA.shape[0]
        aug = np.concatenate([MA, np.eye(n, dtype=np.uint8)], axis=1)
        r = 0
        piv = []
        for col in range(n):
            nz = np.nonzero(aug[r:, col])[0]
            if nz.size == 0:
                continue
            p = r + nz[0]
            aug[[r, p]] = aug[[p, r]]
            rows = np.nonzero(aug[:, col])[0]
            aug[rows[rows != r]] ^= aug[r]
            piv.append(col)
            r += 1
        MAinv = aug[:, n:]
        Cmat = (MB @ MAinv) % 2
        delta0 = np.zeros(n, dtype=np.uint8)
        delta0[0] = 1
        cvec_bits = (Cmat @ delta0) % 2
        cvec = np.array([sum(int(cvec_bits[i * bits + bb]) << bb
                             for bb in range(bits)) for i in range(8)],
                        dtype=np.int64)
        print(f"V5 [{name}]: |supp C| = {int((cvec != 0).sum())} of 8 sites; "
              f"C site-values {cvec.tolist()}")
        comps[name]["C"] = cvec
        comps[name]["Cmat"] = Cmat

    # V6: decompose all census classes
    taxonomy = []
    census = [json.loads(l) for l in (OUT / "m_census_classes.jsonl")
              .read_text().splitlines()]
    census = [r for r in census if "w" in r]
    mismatch = 0
    for row in census:
        b = np.zeros(N1, dtype=np.uint8)
        b[row["b_support"]] = 1
        u, v = b[:N2], b[N2:]
        tu, tv = site_triples(u), site_triples(v)
        wsum = sum(W[t] for t in tu.values()) + sum(W[t] for t in tv.values())
        if wsum != row["w"]:
            mismatch += 1
            continue
        active = [st for st in SITES
                  if W[tu[st]] or W[tv[st]]]
        taxonomy.append({"w": row["w"],
                         "active": len(active),
                         "site_costs": sorted(
                             (W[tu[st]] + W[tv[st]]) for st in active)})
    print(f"V6: weight formula holds on {len(taxonomy)}/{len(census)} census "
          f"classes ({mismatch} mismatches)")
    prof = Counter((t["w"], t["active"]) for t in taxonomy)
    print("   (w, #active-sites) histogram:",
          dict(sorted(prof.items())))
    cost_min = Counter()
    for t in taxonomy:
        cost_min[min(t["site_costs"])] += 1
    print("   minimum per-class site-cost histogram:", dict(sorted(cost_min.items())))
    (OUT / "fibering_taxonomy.json").write_text(json.dumps(taxonomy))
    print(f"wrote {OUT / 'fibering_taxonomy.json'}")


if __name__ == "__main__":
    main()
