# LP1 theory notes — monomial lifted products LP(A, b⊗I)

Notation.  G a finite abelian group, R = F₂[G] (a commutative Frobenius algebra), N = |G|.
For f ∈ R write M(f) for the N×N binary matrix M(f)[g,h] = [g − h ∈ supp f]; M is a ring
homomorphism with M(f)ᵀ = M(f*), f*(g) = f(−g).  Ann(J) = {v : vJ = 0}.  For an a×a matrix A
over R write A† for the conjugate transpose and M(A) for the block expansion.

**Definition.**  LP(A, b) is the CSS code with H_X = [M(A) | M(b*) ⊗ I_a], H_Z = [M(b) ⊗ I_a | M(A†)]
(qubits: a "L" blocks then a "R" blocks of size N).  CSS: H_X H_Zᵀ = M(A)M(b*)ᵀ + M(b*)M(A†)ᵀ
= M(A)M(b) + M(b)M(A) = 0 since M(b) is a scalar block matrix over the commutative R.  This is the
lifted product LP(A, b) of Panteleev–Kalachev with 1×1 second factor; a = 1 gives BB codes.
Weights: with every row and column of A containing exactly three monomials and |supp b| = 3,
every check has weight 6 and every qubit degree 6 (3 X + 3 Z).

## Lemma 1 (normal form, a = 2)
A ↦ D₁ A D₂ with monomial diagonal D₁, D₂ is a code equivalence (shift X-checks and R-qubits
of copy i by s_i, Z-checks and L-qubits of copy j by t_j; b⊗I is invariant because both ends of
each b-edge move together).  The only connected weight pattern for a = 2 is [[2,1],[1,2]], so
WLOG A = [[p,1],[1,q]]; the residual freedom is (p,q) ↦ (x^s p, x^{−s} q), which fixes p ∋ 1.
Further equivalences: b ↦ x^s b; simultaneous φ ∈ Aut(G) on (p,q,b); (A,b) ↦ (A*, b*) (X↔Z).

## Lemma 2 (k)
For A = [[p,1],[1,q]]:  k = 2·dim_F₂ R/(b*, 1 + pq).
*Proof.*  Left kernel of H_X: z = (z₀, z₁) with z M(A) = 0 and z (M(b*)⊗I) = 0.  Row-vector
multiplication by M(f) is multiplication by f*, so the conditions read p* z₀ + z₁ = 0,
z₀ + q* z₁ = 0, b z₀ = b z₁ = 0.  Hence z₁ = p* z₀ and (1 + p*q*) z₀ = 0, i.e.
ker H_Xᵀ ≅ Ann(b) ∩ Ann(1 + p*q*) = Ann(J) with J = (b, (1+pq)*).  Frobenius duality gives
dim Ann(J) = N − dim J = dim R/J, and conjugation maps J to (b*, 1+pq).  The same computation
for H_Z (with A† in place of A) gives the same dimension, and k = dim ker H_Xᵀ + dim ker H_Zᵀ. ∎
Corollary: k ≤ 2·dim ker b = 2δ_b, with equality iff pq ≡ 1 (mod (b*)).
(Checked numerically on 150 random instances over Z₁₄, Z₂₈, Z₆×Z₆, Z₁₂×Z₃, Z₁₀, Z₉×Z₄.
The unconjugated statement k = 2 dim R/(b, 1+pq) is false for even |G| — e.g. Z₁₄, p=1+x¹⁰,
q=1+x⁹, b=1+x¹¹+x¹²: true k = 0, unconjugated formula 6.)

## Lemma 3 (constituent bound)
If u ∈ ker b (as a column vector, M(b)u = 0) then (u, 0, 0, 0) ∈ ker H_Z; it is a nontrivial
X-logical unless it lies in rowspace(H_X), i.e. unless u = p* z₀ + z₁ … for some z ∈ Ann(b*)².
Generically it is nontrivial, so d ≤ D_b := min{|u| : 0 ≠ u ∈ ker b}, and therefore
kd²/n ≤ 2δ_b · D_b² / (4N) = δ_b D_b² / (2N).  (Also d ≤ 2·min{|v₀| : v₀ ∈ Ann(1+p*q*)}
via the R-only logicals (0,0,v₀,p*v₀).)  D_b is computed exactly by enumerating 2^{δ_b} words.

## Lemma 4 (transitivity ⇒ 2BGA)
Let H_X = [A|B], H_Z = [Bᵀ|Aᵀ] with A = ΣA_p, B = ΣB_q sums of m×m permutation matrices such
that every A_p commutes with every B_q.  Put Γ_A = ⟨A_p⟩, Γ_B = ⟨B_q⟩ ≤ Sym(m).  If Γ_A is
transitive then its centraliser is semiregular, so Γ_B ⊆ C(Γ_A) has |Γ_B| ≤ m; if Γ_B is also
transitive then Γ_B = C(Γ_A) is regular and so is Γ_A, and the code is the 2BGA code of Γ_A
(A acts by left, B by right multiplication).  Consequently every non-2BGA code of this type has
an intransitive Γ_A or Γ_B, and the corresponding constituent code ker A or ker B is a direct
sum over the orbits.  In LP(A, b⊗I) the B-orbits are the a copies (length N = n/(2a)).
*Certificate used in practice:* a 2BGA code's Tanner graph is check-transitive; nauty reports
2 X-check orbits for the codes we report ⇒ not 2BGA for any group ⇒ not BB.

## Proposition 5 (thickness ≤ 2)
Write A = A₁+A₂+A₃, B = B₁+B₂+B₃ with A_p, B_q permutation matrices (König edge colouring of
the (3,3)-regular base multigraph; for a=2, A₁ = diag(p₁,q₁), A₂ = diag(p₂,q₂), A₃ = antidiag(1,1),
B_q = b_q ⊗ I).  Layer G_A carries the edges of A₂, A₃ (X–L), B₃ (X–R), B₃ᵀ (Z–L), A₂ᵀ, A₃ᵀ (Z–R);
layer G_B the rest.  Each is 3-regular.  In G_A the X–L edges form cycles alternating A₃, A₂⁻¹
(outer cycles), the Z–R edges form cycles alternating A₃⁻¹, A₂ (inner cycles), and B₃ matches
outer to inner cycles vertex by vertex *because B₃ commutes with A₂ and A₃*: the 4-cycles
(B₃, A₂, B₃⁻¹, A₂⁻¹) close.  Hence every component is a wheel (two equal cycles joined by
radial edges in cyclic order), which is planar.  This is Lemma 2 of arXiv:2308.07915 verbatim;
the only ingredient is A_pB_q = B_qA_p, which holds for scalar B_q.  (networkx planarity check
passes for all 9 layer choices on every code reported.)

## Proposition 6 (depth-7 syndrome cycle)
Table 5 of arXiv:2308.07915 (CNOT layers A₁(R→Z); A₂(X→L), A₃(R→Z); B₂(X→R), B₁(L→Z);
B₁(X→R), B₂(L→Z); B₃(X→R), B₃(L→Z); A₁(X→L), A₂(R→Z); A₃(X→L); then MeasX/InitZ) measures all
checks and acts trivially on logicals.  The tableau proof uses (A₁+A₃)B + A(B₁+B₂) = A₂B + AB₃
(pure algebra) and, in round 6, B A₂ = A₂ B — true here since B is scalar.  Stim confirms that all
detectors of the noiseless memory circuit are deterministic for every code reported.

## Proposition 7 (toric layout)
Fix i≠j, g≠h and let σ = A_j A_i⁻¹, τ = B_g B_h⁻¹ act on the m L-qubits.  Suppose σ, τ commute
and (a, b) ↦ σ^a τ^b ℓ₀ is a bijection Z_μ × Z_λ → L with μ = ord σ, λ = ord τ.  Place
L(σ^a τ^b ℓ₀) at (2a, 2b), X(A_i⁻¹ ℓ) at (2a+1, 2b), R(B_g A_i⁻¹ ℓ) at (2a+1, 2b+1),
Z(B_g ℓ) at (2a, 2b+1) for ℓ = L(2a,2b).  Then every edge of the Cayley graph of Z_{2μ} × Z_{2λ}
is a Tanner edge: X(2a+1,2b)–L(2a,2b) is an A_i edge, X(2a+1,2b)–L(2a+2,2b) an A_j edge
(A_j A_i⁻¹ ℓ = σℓ), X–R(2a+1,2b+1) a B_g edge, X(2a+1,2b+2)–R(2a+1,2b+1) a B_h edge
(B_h⁻¹B_g A_i⁻¹ ℓ = A_i⁻¹ τ ℓ since B's commute with A_i), Z(2a,2b+1)–L(2a,2b) a B_g edge,
Z–L(2a,2b+2) a B_h edge, Z–R(2a±1,2b+1) A_i / A_j edges (using B_g A_i⁻¹ = A_i⁻¹ B_g).
Only B-with-A commutation is used; A_i need not commute with σ.  For a = 2 with {i,j} = {1,3}:
σ = [[0, q₁⁻¹],[p₁⁻¹, 0]] swaps the copies and σ² = (p₁q₁)⁻¹·I, so the criterion is
G = ⟨p₁q₁⟩ ⊕ ⟨b_g − b_h⟩ (internal direct sum), giving the torus Z_{4·ord(p₁q₁)} × Z_{2·ord(τ)}
in which the two long-range edges per vertex are *copy-dependent* shifts.  Implemented and
verified edge-by-edge in `lp_lab.structure.toric_layouts` (16 layouts recovered for gross).

## Remark (relation to coset codes)
Γ_A = ⟨A_p⟩ ≤ Z_N ≀ S_a acts transitively on [a]×G and the diagonal translations are central, so
LP(A, b⊗I) is an Aydin–Tamo–Barg coset code Q^H_{G'}(a,b) with G' = Z_N ≀ S_a (order a!·N^a),
H = a point stabiliser and b central; their Lemma III.2 gives θ = 2 by the same wheel argument
and their Algorithm V a depth-8 cycle.  Their searches ranged over GAP small groups, not these
wreath products; the toric-layout criterion (Prop. 7) and the polynomial k/cap formulas
(Lemmas 2–3) are specific to the LP presentation.
