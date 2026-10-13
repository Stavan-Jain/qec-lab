# Weight-six BB family: completed proof plan and final validation gate

The user requested an explicit family of BB codes with increasing
`k*d^2/n`, distance proofs, and a QECLean `StabilizerCodeWithDistance`
construction, ideally with weight-six checks.

## Research argument

1. PRs 27 and 28 motivate trinomial recurrence seeds. Their primitive-trinomial
   route does not establish infinitely many suitable primitive degrees;
   repeating a fixed square seed keeps its figure of merit constant.
2. Use `p_s(z)=1+z+z^(2^s)` at period `N=4^s-1`, for every positive `s`.
   Frobenius proves `p_s` divides `z^N+1`; the concrete seed has dimension
   `2^s`. No irreducibility or primitive-polynomial assumption is needed.
3. Prove `3^s<=2*s*d_s` by the open Pascal-history bound, a dyadic
   decomposition into at most `s` physical taps, and averaging over all
   cyclic translations. This is the successful unconditional distance proof.
4. Apply the separable BB/hypergraph-product dimension and distance arguments
   to `A=p_s(x), B=p_s(y)`. Retain the actual attained seed minimum `d_s`
   as the exact quantum distance. The natural physical checks have weight six.

These arguments and the exact-distance packaging now compile in Lean.
The stronger analytic upper witness `d_s<=3^s-1` and lower `d_s>=2^s`
are recorded in the mathematical note but are not part of the formal endpoint.
The conjectural equality `d_s=3^s-1` is not required.

## Formalization and growth

1. Work in the isolated QECLean worktree recorded in `state.yaml`; preserve
   unrelated changes and reuse the existing mathlib artifacts.
2. Reuse the homological Pauli bridge and automatic independent CSS
   presentation. Prove both exact code dimension and attained Pauli distance.
3. Package the concrete checks as `stabilizerCodeWithDistance`, then
   `growingCode` at positive mathematical indices `s_t=16*4^t`.
4. Combine only the formal lower bound and `d_s<=N_s` to prove strict increase
   and unboundedness of the actual ratio. Identify that ratio with `k*d^2/n`.
   The sharper analytic subsequence `16*2^t` uses the extra upper witness and
   is separate from the formal selected family.

The capstones, standard-axiom audit, all six required repository checks,
and umbrella validation passed. Pipeline status is `pr-ready`. See
[result.md](result.md), [progress.md](progress.md), and
[validation.md](validation.md) for the exact endpoint and evidence.
