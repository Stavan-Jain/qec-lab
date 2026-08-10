# Build-time budget ledger — mitten_150_30_10

**Constraint (user, 2026-08-10): the certification adds ≤ 5 minutes to a
full QECLean `lake build`. Fallback ceiling 10 minutes — hard cap, never
exceeded.** Measured on the dev machine (Apple silicon, the same class
of box all program reference costs below come from). Wall clock of the
added modules; parallelism across independent files may be counted, but
each allocation below is written as if serial so parallel build is pure
margin.

## Allocations (target tier: 5 min = 300 s)

| Layer | Files | Allocation | Notes |
|---|---|---:|---|
| M1 framework | `LiftedProduct.lean` | 10 s | pure proofs |
| M2 data | `M150/Data.lean` (generated) | 20 s | packed-Nat literals; elaboration only |
| M2 packaging | `M150/{Defs,StabilizerCode}.lean` | 40 s | 2 decoder identities ≈ 5–10 s each (gross analog: ~5 s at 72×144); rest batched |
| M3 witness | `M150/Witness.lean` | 10 s | kernel `decide` |
| M4 floor, ker-H_X side | `M150/FloorZSide*.lean` | 90 s | certificates + ≤ 2 batched native leaves |
| M4 floor, ker-H_Z side | `M150/FloorXSide*.lean` | 90 s | same shape, independent files (parallel margin) |
| M5 assembly | `M150/Distance.lean`, umbrella | 20 s | transport + bundle |
| **Total** | | **280 s** | 20 s slack at the 5-min tier |

The 10-minute fallback, if invoked, is spent **only** on the M4 lines
(→ 240 s/side); every other allocation is firm.

## Reference costs (program history, same machine class)

| Fact | Cost | Lesson |
|---|---:|---|
| A15 sweep layer as filter-form `native_decide` masks (~560 M) | 53 min | banned form |
| Same layer as Gaussian-pivot certificates (`KernelCert` + emitter v3) | 3.9 s | mandatory form |
| Gross decoder identity (72×144, one `native_decide`) | ~5 s | M2 model |
| Z3Z6 bitmask-reworked sweep (per sweep, statement-identical) | ~84 s | ceiling for ONE residual leaf |
| Packed-Nat table walk (RBlock, was 24 min naive) | 0.4 s | table discipline |
| `native_decide` fixed overhead per invocation | ~1–3 s | batch; prefer kernel `decide` under ~1 s |

## Binding design rules

1. Enumeration lives offline (Python emitter, falsify-first); Lean
   checks **certificates** (pivot/rank/left-inverse tables).
2. Few `native_decide` invocations, each a table-driven ∀-goal over
   packed-Nat data; never per-case invocations in a case tree.
3. Kernel `decide` for sub-second facts; `maxRecDepth` bumps local.
4. One heavy leaf per file so lake parallelism turns leaves into margin;
   leaf files import-minimal so they never rebuild on upstream edits.
5. Every session touching Lean re-measures its files and updates the
   ledger below; a PR is not opened while Measured-Total > 300 s
   (> 600 s: work does not merge, full stop — fallback is the
   hypothesis-backed conditional tier for the offending piece).

## Measured ledger

| Date | File | Wall (s) | Method | Notes |
|---|---|---:|---|---|
| — | — | — | — | no Lean landed yet (M0) |

**Measured total: 0 s / 300 s.**

## Measurement protocol

```
cd $QECLEAN_WORKTREE
touch QEC/Stabilizer/Codes/Mitten/M150/*.lean QEC/Stabilizer/Framework/Homological/LiftedProduct.lean
time lake build QEC.Stabilizer.Codes.Mitten
```

Record wall + user CPU; once per phase also measure the cold delta of a
full `lake build` with and without the mitten umbrella imported from
`QEC.lean`.
