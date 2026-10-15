# Why the Steane Code Has Distance 3

*An explainer for the formal proof in QECLean's
`Codes/Small/Steane7Distance.lean` ([PR #78](https://github.com/Stavan-Jain/QECLean/pull/78)).
Audience: you know what a qubit and a Pauli matrix are, and you have seen the
words "stabilizer code" before. Nothing about the argument itself is assumed;
everything used is built up along the way.*

## What is being claimed

The Steane code stores one logical qubit in seven physical qubits. Its
parameters are written `[[7, 1, 3]]`: seven qubits, one logical qubit, and
**distance 3**. Distance is the number that tells you what the code can do:
a code of distance `d` detects any error touching fewer than `d` qubits and
corrects any error touching fewer than `d/2`. Distance 3 therefore means
"corrects any single-qubit error", which is what makes Steane's code the
standard small example.

The claim `d = 3` is two separate statements, and the proof has to deliver
both:

1. **Lower bound.** No error on one or two qubits can go undetected while
   corrupting the logical qubit.
2. **Upper bound.** Some error on exactly three qubits does go undetected
   while corrupting the logical qubit.

The second half is often forgotten in informal accounts, but without it one
has only shown `d ≥ 3`, a bound rather than a value.

## The vocabulary, from scratch

**Pauli strings.** An operator on seven qubits built as a tensor product of
single-qubit Paulis `I, X, Y, Z`, possibly times a phase `±1, ±i`. We write
them as strings, so `IIIXIXX` means `X` on qubits 3, 5, 6 and identity
elsewhere (qubits are numbered `0` to `6`). The **weight** of a Pauli string
is the number of non-identity letters; `IIIXIXX` has weight 3.

**The commutation rule.** Two single-qubit Paulis commute if they are equal
or one of them is `I`, and anticommute otherwise. Two Pauli strings
therefore either commute or anticommute, and the rule is a parity count:

> `P` and `Q` anticommute exactly when the number of qubit positions at which
> their letters anticommute is **odd**.

This one rule drives the entire proof. It is worth pausing on the two special
cases we will use constantly. If `Q` is `Z` on some set of qubits `S` and `I`
elsewhere, then position `i` contributes to the count exactly when `i ∈ S` and
`P` has an `X` or a `Y` there. If `Q` is `X` on `S`, position `i` contributes
exactly when `i ∈ S` and `P` has a `Z` or a `Y` there.

**Stabilizer group.** A set of commuting Pauli strings (the *checks*, or
*generators*) and everything you can multiply together from them. The code
space is the joint `+1` eigenspace of all of them. Measuring the checks is
how the code detects errors: an error `E` flips the outcome of a check
exactly when `E` anticommutes with it.

**Undetectable errors and logical operators.** An error that commutes with
every check leaves all measurement outcomes unchanged; the code cannot see
it. Such operators form the **centralizer** of the stabilizer group. Some of
them are harmless: the stabilizers themselves act as the identity on the code
space. The dangerous ones are the undetectable operators that are **not**
stabilizers (and not a phase times a stabilizer, which is the same thing
physically). These are the **nontrivial logical operators**: they silently
change the encoded qubit.

**Distance.** The smallest weight of a nontrivial logical operator. This
rephrases the two halves above: the lower bound says every Pauli string of
weight 1 or 2 is either detected or is a stabilizer, and the upper bound
exhibits a weight-3 nontrivial logical.

## The Steane code and the Hamming code

Steane's code is built from a classical code: the `[7, 4, 3]` Hamming code,
whose parity-check matrix `H` has three rows and seven columns.

```
qubit:   0  1  2  3  4  5  6
row 1:   1  1  1  0  1  0  0      r₁ = {0, 1, 2, 4}
row 2:   1  1  0  1  0  1  0      r₂ = {0, 1, 3, 5}
row 3:   1  0  1  1  0  0  1      r₃ = {0, 2, 3, 6}
```

Read the columns top to bottom: they are `111, 110, 101, 011, 100, 010, 001`,
the seven **nonzero** binary vectors of length three, each appearing
**exactly once**. Those two facts, no zero column and no repeated column,
are the whole reason the Hamming code has distance 3, and they will be the
whole reason the Steane code does too.

The Steane code has six checks: for each row `r`, a `Z`-check that is `Z` on
the qubits of `r`, and an `X`-check that is `X` on the qubits of `r`.

```
Z₁ = ZZZIZII    X₁ = XXXIXII
Z₂ = ZZIZIZI    X₂ = XXIXIXI
Z₃ = ZIZZIIZ    X₃ = XIXXIIX
```

The checks all commute with each other. Two `Z`-checks trivially do, as do
two `X`-checks. A `Z`-check on row `r` and an `X`-check on row `s`
anticommute at every qubit in `r ∩ s`, and every pair of rows meets in an
even number of qubits: a row meets itself in 4 qubits, and two different
rows meet in exactly 2 (`r₁ ∩ r₂ = {0, 1}`, `r₁ ∩ r₃ = {0, 2}`,
`r₂ ∩ r₃ = {0, 3}`). Codes whose checks split into a pure-`Z` family and a
pure-`X` family like this are called **CSS codes**.

## The lower bound: nothing of weight 1 or 2 slips through

Take any Pauli string `E` of weight 1 or 2 and suppose it commutes with all
six checks. We show `E` is the identity up to phase, so it is not a nontrivial
logical operator. In fact we show something stronger and cleaner: no Pauli
string of weight 1 or 2 commutes with all six checks at all.

**Step 1: split `E` into an X-part and a Z-part.** Record two subsets of the
qubits:

- `A`, the qubits where `E` has an `X` or a `Y` (the letters with an
  "X-component");
- `B`, the qubits where `E` has a `Z` or a `Y` (the letters with a
  "Z-component").

Each is contained in the support of `E`, so `|A| ≤ 2` and `|B| ≤ 2`. And
`E` is the identity up to phase exactly when `A` and `B` are both empty.

**Step 2: the Z-checks see only `A`.** By the commutation rule, `E`
anticommutes with the `Z`-check on row `r` exactly when `|A ∩ r|` is odd. So
"`E` commutes with all three `Z`-checks" says that `A` meets every row in an
even number of qubits. In matrix language: writing `A` as a `0/1` vector `a`,
this is `H · a = 0` over the two-element field. In other words, **`A` is a
codeword of the Hamming code.**

**Step 3: the Hamming code has no codeword of weight 1 or 2.** This is where
the columns come in.

- If `A = {i}` had one element, `H · a` would be column `i` of `H`. No column
  is zero, so `H · a ≠ 0`.
- If `A = {i, j}` had two elements, `H · a` would be column `i` plus column
  `j`. Two distinct columns are never equal, so their sum is not zero.

So `|A| ≤ 2` and `H · a = 0` force `A = ∅`.

**Step 4: the same for `B` with the X-checks.** The `X`-checks see the
Z-components: `E` anticommutes with the `X`-check on row `r` exactly when
`|B ∩ r|` is odd. The Steane code uses the *same* matrix `H` for its
`X`-checks, so the argument of Steps 2 and 3 repeats verbatim and gives
`B = ∅`.

**Conclusion.** `A = B = ∅`, so `E` has no non-identity letters, contradicting
weight 1 or 2. Every Pauli string of weight 1 or 2 anticommutes with at least
one check, so it is detected, and no such string can be a nontrivial logical
operator.

Two remarks on what just happened. First, the argument never looked at the
specific letters of `E`, only at how many qubits carry an `X`-component and
how many carry a `Z`-component; that is why it covers all 4·7 weight-1 strings
and all 9·21 weight-2 strings without enumerating them. Second, the two
properties of `H` that were used, nonzero columns and distinct columns, are
exactly the statement that the classical Hamming code has distance 3, so the
proof is really "classical distance 3 on both sides gives quantum distance
at least 3".

## The upper bound: a weight-3 logical operator

Now we need an undetectable operator of weight 3 that is not a stabilizer.

**Start from the obvious logical operator.** `X̄ = XXXXXXX`, `X` on all seven
qubits, commutes with every check: with the `X`-checks trivially, and with
each `Z`-check because it anticommutes with it at every qubit of the row,
which is 4 positions, an even number. Likewise `Z̄ = ZZZZZZZ` commutes with
every check. And `X̄` and `Z̄` anticommute with each other, because they
disagree at all 7 positions and 7 is odd.

This already tells us `X̄` is a nontrivial logical operator, by a small but
important argument. Suppose `X̄` were a stabilizer, or a phase times one.
Every stabilizer commutes with `Z̄` (we just checked `Z̄` commutes with all
the checks, hence with everything they generate), and a phase does not affect
commutation. So `X̄` would commute with `Z̄`. It does not. Hence `X̄` is
undetectable and not a stabilizer: a nontrivial logical operator. Its weight,
however, is 7, not 3.

**Lower its weight with a stabilizer.** Multiplying a logical operator by a
stabilizer does not change what it does to the encoded qubit, since the
stabilizer acts as the identity on the code space. Multiply `X̄` by the
`X`-check `X₁ = XXXIXII`:

```
X̄  = X X X X X X X
X₁ = X X X I X I I
──────────────────
     I I I X I X X      =  X on qubits {3, 5, 6}
```

Call this operator `X̄₃`. It has weight 3. It commutes with every check,
being the product of two operators that do. And it still anticommutes with
`Z̄`: it meets `Z̄` at qubits 3, 5, 6, three positions, an odd number. The
same small argument as before then shows `X̄₃` is not a stabilizer, nor a
phase times one. So `X̄₃` is a nontrivial logical operator of weight 3.

Combined with the lower bound, the minimum weight of a nontrivial logical
operator is exactly 3. **The Steane code has distance 3.**

## How the Lean proof follows this outline

The formal proof is organised so that each paragraph above is one
declaration. A few points where the formalisation sharpens the informal
story are worth knowing.

**The two halves are one theorem.** `HasCodeDistance C d` in
`Framework/Core/Logical/CodeDistance.lean` is defined as the conjunction of
the lower bound and the witness, and `hasCodeDistance_of` reduces the lower
bound to "for each `w` with `1 ≤ w < d`, no weight-`w` string is a nontrivial
logical".

**The commutation rule becomes an inner product.** Each single-qubit Pauli
has an X-bit and a Z-bit (`I ↦ 00, X ↦ 10, Z ↦ 01, Y ↦ 11`). The parity count
in the commutation rule is a sum over qubits in the two-element field, the
*symplectic inner product*, and `anticommutes_iff_symplectic_inner_one` is
the formal commutation rule. For a Pauli that is the identity away from one
or two named qubits, the sum collapses to one or two terms; that is
`anticommute_weightOneAt_iff` and `anticommute_weightTwoAt_iff` in
`Framework/Core/CSS/CSSDistance.lean`. Pairing against `Z` reads off the
X-bit and pairing against `X` reads off the Z-bit, which is Step 2 and Step
4 in one line each (`symplecticProductSingle_Z_right`, `_X_right`).

**"`A` is a codeword of weight ≤ 2" becomes a case split on two bits.** The
formal lower bound never builds the set `A`. For a weight-2 string with
letters `P` at `i` and `Q` at `j`, it cases on whether `P` has an X-component
and whether `Q` does:

- both do: pick a row containing exactly one of `i, j` (distinct columns);
  the `Z`-check on it anticommutes;
- exactly one does: pick a row containing that qubit (nonzero columns); the
  `Z`-check on it anticommutes, since the other qubit contributes nothing;
- neither does: both letters are `Z`, and the `X`-check on a row containing
  exactly one of `i, j` anticommutes.

Weight 1 is the same with one letter. This is the content of
`weight_one_anticomm_witness_of_columns` and
`weight_two_anticomm_witness_of_columns`, and their packaging
`hasCodeDistance_three_of_columns` is stated for an arbitrary CSS code: give
it the rows of the `Z`-check matrix and the `X`-check matrix as families of
`Finset`s, proofs that no column is zero (`hz_cover`, `hx_cover`) and that
columns are distinct (`hz_sep`, `hx_sep`), and a weight-3 witness, and it
returns distance 3. Steane passes the same three rows for both sides.

**The concrete facts are decided by the kernel.** The only Steane-specific
inputs are the row family

```lean
def row : Fin 3 → Finset (Fin 7) := ![{0, 1, 2, 4}, {0, 1, 3, 5}, {0, 2, 3, 6}]
```

the two column facts

```lean
lemma row_cover    : ∀ i : Fin 7, ∃ r, i ∈ row r := by decide
lemma row_separate : ∀ i j : Fin 7, i ≠ j → ∃ r, (i ∈ row r ↔ j ∉ row r) := by decide
```

and the witness `logicalXw3 = σ[IIIXIXX]`, with `weight logicalXw3 = 3` and
its anticommutation with `Z̄` also closed by `decide`. `decide` asks Lean's
kernel to evaluate the statement by computation, so the proof term is
checked by the same trusted core as everything else. The earlier, parked
version of this file instead enumerated all 784 weight-2 cases and needed
`native_decide`, which trusts the compiler rather than the kernel and is
banned in this library. The column argument replaces that table with two
small finite checks and reads like the proof above.

**The witness's nontriviality is one lemma.** The small argument "anticommutes
with a centralizer element, so not a stabilizer even up to phase" is
`isNontrivialLogicalOperator_of_anticommute_centralizer` in
`Framework/Core/Logical/LogicalOperators.lean`. It is applied once, to `X̄₃`
and `Z̄`.

The final theorem is a single application:

```lean
theorem code_has_distance_three : HasCodeDistance stabilizerCode 3 :=
  hasCodeDistance_three_of_columns row row generators stabilizerCode stabilizerCode_toSubgroup_eq
    zOn_row_mem xOn_row_mem row_cover row_cover row_separate row_separate
    ⟨logicalXw3, logicalXw3_isNontrivial, logicalXw3_weight⟩
```

and `#print axioms` on it reports only `propext`, `Classical.choice`, and
`Quot.sound`, mathlib's three standard axioms.

## Where the argument stops

The argument is specific to distance 3 because it only has to rule out
weights 1 and 2, where "codeword of weight ≤ 2" reduces to two facts about
columns. For a CSS code of larger distance one would formalise the general
statement that the X-part of an undetected Pauli is a codeword, and then
bound the classical distance directly. And the argument needs the CSS split:
for a code such as the `[[5, 1, 3]]` code, whose checks mix `X` and `Z` on
the same qubit, the Z-checks do not see only the X-part, and the library's
proof for that code is a genuine search over anticommutation witnesses
instead.
