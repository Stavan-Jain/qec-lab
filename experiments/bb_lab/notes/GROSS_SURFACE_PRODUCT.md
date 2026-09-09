# Exact distance for a gross-augmented surface family

2026-09-08. This is an explicit, analytically proved alternative using
the actual gross code as a fixed factor. It is **not the original IBM
torus family, not a scalar two-polynomial BB family, and not a novelty
claim**. The original proof tracks remain active. The result identifies
what is available if the permitted modification includes replacing the
torus enlargement by a surface-code homological product.

## 1. Parameters and the precise modification

**Theorem H.** Let C be the original gross code on `(12,6)`, retaining
all 72 X and 72 Z generators. Its banked exact parameters are
`[[144,12,12]]`, with `d_X=d_Z=12`. The explicit construction below gives,
for every integer `t>=1`, a connected CSS LDPC code with

\[
\boxed{[[144(3t^2-3t+1),\ 12,\ 12t]]},\qquad d_X=d_Z=12t.
\]

Every displayed stabilizer generator has weight at most eight; each
qubit participates in at most six generators of each CSS type. At
`t=1` the construction is exactly the gross code. No search over t is
needed for the distance statement.

The construction preserves the gross seed, 12 encoded qubits, the
target `12t` distance law, and bounded check weights. It changes the
growing geometry and introduces additional qubit/check types and open
boundaries. Asymptotically it uses six times the conjectured IBM b=1
qubit count at equal distance and k. Thus this is a concrete
gross-derived alternative, **not evidence that IBM's original
`[[72r(r+1),12,12r]]` conjecture has been proved**.

## 2. Complete definition of the checks

All vector spaces and tensor products are over F2 with their stated
coordinate bases. Use the lab gross matrices to define the complex

```
C_2=F2^72 --d_2=H_X^T--> C_1=F2^144 --d_1=H_Z--> C_0=F2^72.
```

The 12 linear dependencies among these two sets of generators are
retained. In particular `H_0(C)` and `H_2(C)` need not vanish.

Let `J_t:F2^t -> F2^(t-1)` be the adjacent-difference map

```
(J_t v)_j=v_j+v_(j+1),     0<=j<t-1.
```

Define two complexes in degrees 1 and 0:

```
R_t: F2^t --J_t--> F2^(t-1),
I_t: F2^(t-1) --J_t^T--> F2^t.
```

Empty spaces at `t=1` have dimension zero. The map `J_t` is surjective
and has kernel spanned by the all-one vector `N_t`. Thus R has only
`H_1(R)=F2`, represented by N. The transpose map is injective; I has
only `H_0(I)=F2`, detected by the sum of its t vertex coordinates.

Form the full graded tensor complex

```
E_t=C tensor R_t tensor I_t,
d_E=d_C tensor 1 tensor 1 + 1 tensor J_t tensor 1
                              + 1 tensor 1 tensor J_t^T.
```

Each term acts only on degree-one coordinates in its interval factor;
the C term uses d_1 or d_2 as appropriate. Qubits are coordinate basis
vectors of `(E_t)_2`. X-checks are the columns of `d_E:(E_t)_3->(E_t)_2`;
Z-checks are the rows of `d_E:(E_t)_2->(E_t)_1`. This fixes all checks
for every t, including the boundary checks. Each factor squares to
zero, and mixed terms occur twice, so `d_E^2=0` and the checks commute.

The four qubit summands, listed by factor degrees, are

```
(1,1,0): C_1 tensor R_1 tensor I_0,   dimension 144 t^2,
(2,0,0): C_2 tensor R_0 tensor I_0,   dimension 72 t(t-1),
(0,1,1): C_0 tensor R_1 tensor I_1,   dimension 72 t(t-1),
(1,0,1): C_1 tensor R_0 tensor I_1,   dimension 144 (t-1)^2.
```

Their sum is the stated n. Splitting each finite-dimensional complex
into its homology and contractible two-term summands gives

```
H_2(E_t)=H_1(C) tensor H_1(R_t) tensor H_0(I_t).
```

No other summand survives, because the interval factors have no other
homology. Thus the code encodes exactly 12 qubits despite the gross
check redundancies. This dimension calculation is not used as a
substitute for the metric proof below.

## 3. First interval: an exact one-sided distance multiplication

Put `D=C tensor R_t`. A degree-two chain has coordinates

```
(u_0,...,u_(t-1); v_0,...,v_(t-2)),
u_j in C_1,       v_j in C_2.
```

It is a cycle precisely when

```
d_1 u_j=0 for every j,
u_j+u_(j+1)=d_2 v_j for 0<=j<t-1.
```

Hence every u_j is a base cycle, and all of their homology classes
are equal. If their common class is nonzero, the total physical
weight is at least `t d_X(C)`.

Conversely, if this common class is zero, choose `h_j in C_2` with
`d_2 h_j=u_j` and subtract the boundary of `(h_j)_j in C_2 tensor R_1`.
This eliminates all u_j. The remaining v_j lie in `ker d_2`. Choose
`h'_0=0` and successively `h'_(j+1)=h'_j+v_j`. Then every h'_j lies
in `ker d_2`, and its boundary eliminates the remaining v_j without
restoring any u_j. Thus the cycle is a boundary. This step does not
assume `ker d_2=0`.

A nontrivial degree-two D cycle must therefore have nonzero common
class, proving its homological distance is at least `t d_X(C)`.
If c is a minimum nontrivial base X cycle, `u_j=c` for every j and
`v_j=0` is a nontrivial cycle of weight exactly `t d_X(C)`. We have
proved equality, with an explicit description of the attaining class.

## 4. Second interval: a weight-contracting homology isomorphism

The augmentation `epsilon:I_0->F2` that sums vertex coordinates,
together with zero on I_1, is a chain map. Its inverse on homology
includes F2 at the first vertex. The two maps are chain-homotopy
inverses: connect each other vertex to the first by the unique path
of interval edges. Its boundary is that vertex plus the first one.
Therefore

```
pi=1_D tensor epsilon: E_t -> D
```

is a homology isomorphism in every degree. On a degree-two chain,
pi sums its D_2 components over the I vertices and discards its
`D_1 tensor I_1` component. Hamming weight cannot increase under this
map: each nonzero output coordinate needs at least one nonzero
input coordinate, and these inputs are distinct for distinct outputs.

Every nontrivial E degree-two cycle projects to a nontrivial D
degree-two cycle. Section 3 consequently proves

```
d_X(E_t) >= t d_X(C).
```

Embedding a minimum D cycle at one I vertex attains that weight.
It is nontrivial because its projection is the original D class.
Hence `d_X(E_t)=t d_X(C)=12t`.

## 5. Dual distance and general Pauli operators

Reverse degrees and dualize the full product complex. Since E has
degrees 0 through 4, degree-two homology of its reversed dual is
the degree-two cohomology of E, which measures Z-logical distance
in our convention. This dual complex is

```
C^vee tensor I_t tensor R_t.
```

Here `(C^vee)_j=C_(2-j)^*`; reversing the two interval factors swaps
R and I. Permute the last two tensor factors, a permutation of physical
coordinates preserving Hamming weight. We obtain exactly the form
treated in §§3–4, with base C^vee. Its base homological distance is
`d_Z(C)=12`. Thus `d_Z(E_t)=12t`, including a matching pure Z
representative. No general product-distance conjecture is used.

For a CSS code, any nontrivial Pauli has a nontrivial pure X or pure
Z component. The component's support is contained in the Pauli's
support. Consequently the full Pauli distance equals the smaller
of these two sector distances, proving Theorem H's exact d.

## 6. Bounded checks, participation, and connectedness

The gross H_X and H_Z have row weight six and column weight three.
For X checks, the source summand degrees of E_3 are `(2,1,0)`,
`(2,0,1)`, and `(1,1,1)`. Their weights are bounded respectively by
`6+2`, `6+2`, and `3+2+2`. For Z checks use the degree-one target
summands `(1,0,0)`, `(0,1,0)`, `(0,0,1)`; the bounds are `3+2+2`,
`6+2`, and `6+2`. Removing missing boundary terms only lowers weight.

For the four qubit types in §2, their X-check participation bounds
are respectively `5,4,6,5`; their Z-check participation bounds are
`5,6,4,5`. These follow by counting possible cofaces and faces in
each factor. Thus both generator weight and qubit participation
are uniformly bounded: this is an LDPC family, unlike gross/Shor
concatenation.

The `(1,1,0)` qubits form a t-by-t array of gross-qubit sets. Within
each set the gross X-check incidence graph connects all 144 qubits.
The Z checks of type `(1,0,0)` connect equal gross coordinates at
adjacent R vertices; X checks of type `(1,1,1)` connect them at
adjacent I vertices. Every remaining qubit type shares one of these
or a `(2,1,0)` / `(0,1,0)` check with this array. Hence the full
qubit/check incidence graph is connected for every t. This is not
a growing number of disconnected gross codes with fixed distance.

## 7. Provenance and verification scope

This is a specialization of the known augmented-surface / homological
product construction, not a new general code design. The relevant
primary reference is [Delfosse and Hastings, §§1 and 4](https://arxiv.org/html/2009.14226v2#S4),
which explains distance multiplication for suitable surface products.
Its initial setup removes redundancies from the small factor. Here
the retained gross dependencies are handled directly in §§2–4, using
the interval factors' homology. No decoder theorem with independent-
check hypotheses is imported for this redundant presentation.

[Zeng and Pryadko](https://arxiv.org/abs/2007.12152) treat products with
two-term complexes more generally; their published puncturing lemma
has an [erratum](https://doi.org/10.1103/PhysRevA.104.029902). Our proof
does not invoke that lemma or claim multiplicativity for arbitrary
pairs of quantum complexes.

The fixed gross exact-distance certificate is the only base-distance
input. The construction and all-parameter lower and upper proofs are
explicit above. A short diagnostic checks the generated differentials,
ranks, weights, connectedness, and lifted witnesses at selected small
t; those tests are convention checks, not the universal distance proof.

Acceptance boundary: Theorem H is an exact growing-distance family
derived from gross with bounded checks. It does not supply a distance
theorem for the original trinomials on expanding IBM rectangles. The
active IBM/closely-modified-BB research goal is not silently redefined
as constructing an arbitrary augmented surface code.
