# VQE through minimal VHamilton

VQE ansatz blocks are ordinary open VHamilton `Proc` terms. The object language
contains only Hamilton's original dyadic bare actions plus the two new atoms
`δ_k^{α_j}` and `-δ_k^{α_j}`. In particular, there are **no scalar-expression
constructors** and no process scaling operator.

The optimizer supplies a valuation externally, `ρ(α_j) = q_j/2^r`, where `q_j`
is an integer grid coordinate and `r` is variational precision. The compiler
evaluates phase angles under `ρ`, *without substituting arbitrary real constants
into the object-language syntax*. The ideal denotation likewise uses `ρ`.

Three independent numerical controls remain: target Hamiltonian dyadic
precision `k`, variational grid precision `r`, and Trotter refinement `n`.
Sequential blocks live in `VAnsatz` metadata, not in Hamilton parallel
composition. The same open term can be reused at every parameter precision.

## Experiments

One qubit: `cabal run vqe-demo` searches `q=0..7` with `α=q/16` for
`H=X+Z`, using the open `αY` generator. Expected optimum `q=5`, energy
`-sqrt(2)` (subject to numerical tolerance).

Two qubits: `cabal run chem-vqe -- 10 3 0 exhaustive` executes the
three-parameter DMANH+ ansatz, in the minimal syntax. All existing optimizer
modes (`refine`, `coordinate`, `analytic`, `analytic-direct`, `exhaustive`)
evaluate these same processes, with new valuations for each candidate.

The exact energy oracle in `scripts/vqe/vqe_eval.py` uses the target
Hamiltonian matrix and the circuit's statevector, as before. The outer
experiment drivers use the same Haskell executables, so they now use the
updated syntax automatically.

See `docs/VHAMILTON.md` for the grammar and semantics.
