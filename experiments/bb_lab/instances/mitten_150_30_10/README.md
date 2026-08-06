# [[150,30,10]] mitten code — vendored instance data

The smallest of the eight **mitten codes** of Bhardwaj–Ma–Meister–King–
Bluvstein–Preskill–Cain–Xu–Huang, *High-rate qLDPC processors*,
arXiv:2607.28795: non-abelian lifted-product codes LP(A,B) with 1×2 base
matrices over 𝔽₂[G], G = C₅×S₃ = GAP `SmallGroup(30,1)`, base-matrix
data in the paper's Table XIII (element indices into GAP's
`Elements(G)`), rate 1/5, check weight 9.

## Files

Bit-identical copies of the authors' data release —
<https://github.com/a7b/yarn> (MIT), commit
`82fb695a1e403e8a77a2adc4a42898d96d3bf85b` (2026-08-03), path
`processor_codes/mitten/[[150,30,10]]/`:

| file | shape | meaning |
|---|---|---|
| `Hx.npy` | 60×150 | X-type checks (rows weight 9) |
| `Hz.npy` | 60×150 | Z-type checks (rows weight 9) |
| `Lx.npy` | 30×150 | canonical logical-X basis (rows weight 18, ∈ ker Hz) |
| `Lz.npy` | 30×150 | canonical logical-Z basis (rows weight 10, ∈ ker Hx) |

`Lx · Lzᵀ = I₃₀ (mod 2)`; row *i* of each is the same logical qubit.
sha256 sums are pinned in `certificates/mitten_150_30_10_{X,Z}.cert.json`
via the `h_check_sha256` / `logical_space_sha256` fields (computed over
the bit content, not the npy container).

## Why it's here

`scripts/mitten150_tandem.py` independently certifies d = 10 for this
code with the lab's SAT stack (Tandem MaxSAT optimum + CMS native-XOR
UNSAT ladder, both CSS directions) — the paper's five "exact" distances
are Gurobi-IP results with no published checkable artifacts, and this is
the first non-BB, non-abelian code pushed through the lab pipeline. See
`notes/mitten150_tandem_verification.md`.

The matrices are vendored (rather than fetched) so the verification is
reproducible offline and the certificate hashes stay pinned to reviewed
bytes. GAP-free structural validation against the paper's Definition 4 /
Eq. (2) (block pattern, L/R commutation, non-abelian signature, canonical
weights) lives in the script's `validate` mode; re-deriving the matrices
from Table XIII's element indices requires GAP and is deliberately out of
scope here.
