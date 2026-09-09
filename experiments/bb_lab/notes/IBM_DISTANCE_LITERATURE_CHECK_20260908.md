# Focused primary-source distance check — 2026-09-08

Purpose: check whether recent bicycle-code work supplies a missing
exact-distance argument, rather than assuming a new algebraic
description already proves the desired metric theorem. This is a
bounded check of the items below, not a comprehensive novelty review.
No external claim that the original conjecture is solved follows.

## Recent spectral and colon-ideal method

[Sabo, Can, Marquis, arXiv:2608.27565v1](https://arxiv.org/html/2608.27565v1),
especially Theorems 88 and 90 and §VIII.2, was inspected directly.
The general mixed-block colon bound is limited by the sparse checks;
the stabilizer-excluded refinement assumes semisimplicity. Its cover
discussion concerns lower-bound certificates, not exact distance.

For our fixed characteristic-two family, both periods are even, so
the semisimple hypothesis fails. The unrefined mixed term is at most
six for two weight-three polynomials. The subsequent original-family
audit in `A43_S4_block_purity.md` identifies the common single-block
image and proves both exclusions outside it. The further parent
result `IBM_RELATIVE_COLON_OBSTRUCTION.md` shows that even this
class-matched relaxation is exactly six at gross. These are separate
repository results, not statements attributed to the inspected paper.
S11's winding mixedness remains a different notion.

## Covers

[Symons, Rajput, Browne, arXiv:2511.13560v1](https://arxiv.org/html/2511.13560v1),
§4, proves transfer/projection bounds with explicit hypotheses, not a
universal multiplicative distance law. Odd covers have a lifted upper
bound; preserving dimension adds a base-distance lower bound. These
do not coincide in a growing family. The new local parity/Tor proof
in this repository establishes particular even-cover kernels but
does not supply the missing metric equality.

## Exact finite enumeration is not an analytic family theorem

[Lu et al., arXiv:2608.09115v1](https://arxiv.org/html/2608.09115v1),
§IV-B, explicitly excludes stabilizers during exact kernel-weight
enumeration. This is an important distinction from classical additive
code distance, but the enumeration is still the bottleneck. It does
not remove the need for a universal certificate in the current
no-long-computation proof program.

## Frobenius-related bicycle generators

[Rabeti and Mahdavifar, arXiv:2605.14173v1](https://arxiv.org/html/2605.14173v1),
Theorems 3 and 5, provides a logical parametrization and structured
upper bounds for its univariate subclass. The inspected argument
does not give a matching universal lower bound for the original
two-variable pair. Our dyadic witness construction similarly needs
both a smaller attaining representative and an all-class floor.

## Disposition

These sources inform the proof mechanisms and their hypotheses;
none of the inspected statements closes the original IBM target.
The semisimple colon refinement motivated the completed scope audit
and base-member metric obstruction above. None of these results
justifies restarting finite distance searches or treating an upper
bound as exact. No assertion is made that a solution cannot exist
elsewhere in the literature.

## Matching symmetries and repeated-root coding bounds

[Sahay, Williamson, Brown, arXiv:2602.22770v1](https://arxiv.org/html/2602.22770v1),
sections III–IV and appendix D, was checked as a structurally
different route. It constructs matching graphs from code symmetries
and relates truncated symmetries to logical strings. The reported
gross-code decoders still have low-weight failures; section IV.1
does not claim correction of every weight-five error. Appendix D's
string-space calculation determines topological data, not sharp
Hamming minima. Thus neither the decoder results nor the toric
equivalence supplies the missing all-size exact-distance proof.
Retaining joint syndrome information across projections remains
essential; this is consistent with, but does not prove, our local
compatibility obstructions.

[Zeh and Ulmschneider, arXiv:1506.02820v1](https://arxiv.org/html/1506.02820v1),
sections II–III, develops Hasse-derivative bounds for classical
repeated-root cyclic codes. Such a bound requires a concrete
classical code in its physical cyclic coordinates and its full
root conditions. A quotient-homology description alone does not
establish these hypotheses or eliminate stabilizers from the
minimum. No assertion here imports its classical distance into
the original quantum relative distance. The new Cartier descent
probe instead checks actual chain maps and their physical costs.

## Coprime and generalized bicycle constructions

[Wang and Mueller, arXiv:2408.10001v4](https://arxiv.org/html/2408.10001v4),
equations (11)–(14) and Algorithms 1–2, was checked directly. The
cyclic reduction uses coprime group orders and determines dimension
from a polynomial gcd; its construction algorithms still call a
distance-upper-bound routine. Our IBM rectangle orders have a common
factor six, so coprimality of r and r+1 does not itself make this
the same scalar cyclic construction. Neither a prescribed dimension
nor a gcd calculation supplies the needed exact distance.

[Wang and Pryadko, Distance Bounds for Generalized Bicycle Codes](https://www.mdpi.com/2073-8994/14/7/1348)
gives distance bounds and, in particular, a graph-cover treatment of
weight-four codes. That incidence-graph argument does not automatically
apply to the original weight-six generators or the proposed
weight-eight one-check modification. The repository's fixed Fano-minor
obstruction already explains the original pair's distinction from an
ordinary graph-incidence code. No matching all-size metric theorem for
either current pair is imported from these bounds.
