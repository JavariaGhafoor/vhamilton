# From Hamilton to VHamilton: implementation changes

## Conceptual diff against original Hamilton

The starting point is Hamilton's original four categories (`Clifford`,
`BareAction`, `Action`, `Proc`) and its existing Clifford operations, diagonal
direct sums, conjugation, and additive parallel composition. **All of those
constructors remain**. The only new language forms are:

```text
Hamilton:   d ::= id | δ_k | -δ_k | δ_k^j | d ⊕ d
VHamilton:  d ::= id | δ_k | -δ_k | δ_k^j
                 | δ_k^{α_j} | -δ_k^{α_j} | d ⊕ d
```

For a valuation `ρ` supplied outside the language, the new bare-action clauses
are `⟦±δ_k^{α_j}⟧_ρ = ±ρ(α_j)/2^k`. The runtime chooses `ρ`, the compiler
computes phase angles, and the VQE layer arranges successive blocks. There is
no new scalar-expression grammar or core sequence constructor.

## Implementation touchpoints relative to Hamilton

| File | Why it changes |
|---|---|
| `src/Syntax.hs` | Add `DeltaParam`/`NegDeltaParam`, parameter discovery, JSON encoding and the external `Valuation` type |
| `src/Semantics.hs` | Extend the existing dyadic-halving cases to parameterized signed actions |
| `src/Compiler.hs` | Read parameter valuations to compute phase angles; retain closed-program entry points |
| `scripts/semantics/denot.py` | Extend the independent matrix denotation for parameterized actions |
| `src/ExamplesVHamilton.hs` | Show one open term denoting a family of Hamiltonians |
| `src/VQE.hs` | Add VQE *outside* the language: block order, parameter grid and objective calls |
| `app/VQEDemo.hs`, `app/ChemVQE.hs` | Express variational generators with the new bare actions |
| `test/` | Exercise valuations, signs, halving and existing Hamilton behaviors |

The fixed target Hamiltonians in these experiments continue to be ordinary
Hamilton programs (`DeltaPow k j`). VHamilton introduces free parameters in
the **ansatz generators**, not in the target Hamiltonian or in a general
scaling operation. See `README.md` and `docs/VHAMILTON.md` for a gradual
Hamilton-first explanation.

---

## Migration from an intermediate scalar-expression prototype

The following notes describe removal of scalar support from an *intermediate
VHamilton prototype*; these scalar forms were **not part of original Hamilton**.
They should not be confused with the minimal Hamilton → VHamilton language diff.

### Minimal bare-action rewrite (2026-10-15)

Updated the previous VHamilton implementation to the revised *minimal* syntax.

- Restored separate Hamilton bare actions: `Id`, `Delta k`, `NegDelta k`, and
  `DeltaPow k j` (integer `j`).
- Added **only** `DeltaParam k alpha` and `NegDeltaParam k alpha` as VHamilton
  primitive extensions.
- Removed `Scalar`, `Const`, `Param`, `Neg`, the generalized `Delta k x`, and
  term-level arbitrary real-scalar instantiation.
- Refactored halving/scheduling, phase-angle computation, and Python JSON
  denotation for the precise six bare-action cases.
- The existing Haskell VQE runner compiles each *open* process with the
  external valuation `rho(alpha)=q/2^r` through `compileTrotterWith`.
- The ideal-unitary checker can also evaluate open terms with
  `compareToIdealWith rho`.
- Re-expressed both the one-qubit calibration and all three DMANH+ chemistry
  VQE generators using `DeltaParam` / `NegDeltaParam`, including negative `ZY`.
- Added regression tests for bare-action semantics, signs, free parameters,
  dyadic schedule segmentation, open compilation, and matrix expectation.
- Updated project language and VQE documentation.

For building and running see `README.md` or `docs/VHAMILTON.md`.

**Validation note:** The available container does not have GHC/Cabal or
Qiskit, so full Haskell build + end-to-end VQE execution must be performed in
a configured project environment. NumPy/SciPy-only denotational regression
checks are included and can be run separately.
