# Formalization progress: growing BB figure of merit

Date: 2026-09-13. Worktree:
`/Users/stavanjain/Code/QECLean/.Codex/worktrees/bb-growing-fom`, branch
`codex/bb-growing-fom`.

| Obligation | Current evidence |
|---|---|
| Trinomial divisibility and recurrence algebra | Lean module compiled |
| Concrete recurrence space and dimension `2^s` | Lean module compiled |
| Nonzero open-history weight `>=3^s` | Lean module compiled |
| At most `s` physical taps per truncated history coordinate | Lean module compiled |
| Shift averaging and `3^s<=2*s*d_s` | Lean module compiled; standard axiom audit passed |
| Actual finite minimum, positivity, and attainment | Lean module compiled |
| Separable BB homology dimension and exact Pauli distance | Lean modules and code capstone compiled |
| `StabilizerCodeWithDistance` for the concrete BB checks | Code capstone compiled |
| Physical weight-six stabilizer generating set | Lean module compiled |
| Strict growth on `s_t=16*4^t` and unbounded actual ratio | Family capstone compiled; standard axiom audit passed |
| Required repository validation | All six checks and umbrella check passed; no new diagnostics |

No unproved closed distance formula or primitive-trinomial existence
assumption occurs in the endpoint. The central bound was independently
reviewed for fixed taps, nonzero translated histories, and full cyclic
averaging. The formal distance is the attained minimum of the actual
recurrence space, transferred exactly to nontrivial logical Pauli weight.

The sharper analytic witness `d_s<=3^s-1`, elementary bound `d_s>=2^s`,
and subsequence `s_t=16*2^t` remain outside the formal endpoint. The finite
enumerations are not Lean certificates and supply no universal premise.

See [result.md](result.md) for capstone links and [validation.md](validation.md)
for the completed checks. Pipeline status is `pr-ready`.
