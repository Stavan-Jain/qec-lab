# Weight-six BB family with increasing actual distance efficiency

Date: 2026-09-13. Status: **formal endpoint verified; all required repository checks passed**.

The requested mathematical and Lean endpoint is implemented. For every
positive integer `s`, set `q=2^s`, `N=4^s-1`, and use the BB checks
`A=1+x+x^q`, `B=1+y+y^q` on `Z_N x Z_N`. The exact parameters are

```
[[2*(4^s-1)^2, 2*4^s, d_s]],
```

where `d_s` is the attained minimum Hamming weight among nonzero words of
the explicit recurrence `c_(j+q)=c_j+c_(j+1)`. The quantum minimum is proved
equal to this classical minimum, including an attaining nontrivial Pauli
logical operator. The distance is not a supplied hypothesis.

The formal bounds are `0<d_s<=4^s-1` and `3^s<=2*s*d_s`. On the explicit
subsequence `s_t=16*4^t`, the actual `k*d_s^2/n` is strictly increasing and
unbounded. Every natural threshold `b` is exceeded at selected position
`b+1`. The stabilizer group admits a generating set whose members all
have Pauli weight exactly six.

## Compiled capstones

All declarations below are in `Quantum.Stabilizer.Homological.BB.Fractal`.

- [Code.lean](/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/Code.lean): `stabilizerCodeWithDistance` packages all
  parameters and the exact distance; `chainComplex_pauli_distance` proves
  the least nontrivial logical Pauli weight.
- [Family.lean](/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/Family.lean): `growingCode`,
  `growingEfficiency_eq_kd2_div_n`, `growingEfficiency_strictMono`, and
  `growingEfficiency_unbounded` package the selected family and its metric.
- [LowerBound.lean](/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/LowerBound.lean): `seedSpace_weight_lower`
  connects the proved open Pascal histories to physical cyclic Hamming weight.
- [CheckWeight.lean](/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom/QEC/Stabilizer/Codes/BivariateBicycle/Fractal/CheckWeight.lean):
  `stabilizerCodeWithDistance_has_weight_six_generating_set` proves the physical check-weight claim.

The underlying Lean code index `u` is zero-based, corresponding to
mathematical `s=u+1`. `distanceAt s=seedDistance(s-1)` supplies positive
mathematical indexing for the selected family.

## Scope and evidence

The stronger analytic bounds
`max(2^s,ceil(3^s/(2*s)))<=d_s<=3^s-1` and the sharper subsequence
`s_t=16*2^t` are proved in the
[mathematical note](../../../experiments/bb_lab/notes/BB_FRACTAL_UNBOUNDED_FOM_THEOREM.md).
The lower `2^s` bound, upper witness `3^s-1`, and graph connectedness are
not claimed as part of the current Lean endpoint. The exact formula
`d_s=3^s-1` remains unproved in this development. None of these stronger claims is a
premise of the formal selected-family result.

Complete finite enumeration gives `d_s=2,8,26,80,242` for `s=1,...,5`.
These reproducible calculations are separate from kernel certification and
are unused in the universal proof. The
[independent audit](../../../experiments/bb_lab/notes/BB_FRACTAL_LOWER_BOUND_AUDIT.md)
records the boundary and averaging checks. No novelty claim is made over
Yoshida's fractal-code results or the established hypergraph-product construction.

All six required repository checks and the umbrella check passed. The new
modules have no compiler/linter warnings, and the full axiom audit accepts
only `propext`, `Classical.choice`, and `Quot.sound`. See
[validation.md](validation.md) for commands, timings, and the audit output.
Pipeline status is `pr-ready`.

## Local delivery

The verified Lean development is committed as `1e7f44460ad65e9ec6b1e65a198faf79477f8141` on
`codex/bb-growing-fom`. The companion research branch is
`codex/bb-growing-fom-research`. Both are ready for review locally.
Automatic approval review rejected remote pushes and draft PR creation
because external publication was not explicitly authorized. No push or
PR creation ran; publication awaits the user's approval.
