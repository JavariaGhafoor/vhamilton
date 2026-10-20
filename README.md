# From Hamilton to VHamilton: noiseless VQE

This repository starts with **Hamilton**, a compositional language for constructing
Hamiltonians and compiling their associated time-evolution operators into quantum
circuits. It then makes one **minimal syntactic extension**, called **VHamilton**,
so that a single program can describe a *family* of Hamiltonians indexed by
variational parameters. That family supplies the parameterized time-evolution
blocks used in the noiseless variational quantum eigensolver (VQE) experiments.

**What changes in the language?** Only the set of **bare actions**: VHamilton
adds `δ_k^{α_j}` and `-δ_k^{α_j}`. It does **not** introduce scalar
multiplication, arbitrary arithmetic expressions, or any new process operator.
The original Hamilton grammar and constructors remain available.

## 1. Starting point: Hamilton

Hamilton constructs diagonal bare actions, changes their bases with Clifford
conjugation, and combines Hamiltonians using additive parallel composition.
A compact version of the **original** grammar is:

```text
Clifford     u ::= I | H | S | CX | SWAP | u · u | u ⊗ u
Bare action  d ::= id | δ_k | -δ_k | δ_k^j | d ⊕ d
Action       a ::= d | u a u†
Process      p ::= a | p ∥ p
```

Here `k` is a nonnegative integer and `j` is an integer. The three basic
numerical dyadic values have denotations

```text
⟦id⟧ = 0,    ⟦δ_k⟧ = 1/2^k,    ⟦-δ_k⟧ = -1/2^k,
⟦δ_k^j⟧ = j/2^k.
```

Bare actions denote scalar diagonal entries; `⊕` combines them into a diagonal
matrix. Clifford conjugation changes the basis of a Hamiltonian action:
`⟦u a u†⟧_H = U ⟦a⟧_H U†`. Parallel composition adds Hamiltonians:
`⟦p ∥ q⟧_H = ⟦p⟧_H + ⟦q⟧_H` (for matching dimensions).
The corresponding unitary semantics used by this implementation is
`⟦p⟧_U = exp(-2πi ⟦p⟧_H)`. For a noncommuting sum, the compiler uses a
Trotter schedule to approximate this exponential.

For example, Hamilton can express the **fixed** one-qubit Pauli generators:

```haskell
z = Atom (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1)))
x = Atom (Conj H (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1))))
h = x ∥ z  -- H = X + Z, not a sequential circuit
```

`DeltaPow k j` is just the Haskell name for Hamilton's original `δ_k^j`:
it is **not exponentiation or scalar multiplication syntax**.

## 2. Why extend Hamilton?

For a fixed Hamiltonian, the compiled evolution has fixed phases. In VQE,
we instead need to repeatedly change the durations of chosen evolution blocks,
prepare trial states, and compare their energies. Rebuilding the *program*
for each duration obscures the shared structure of the ansatz.

**VHamilton keeps the Hamiltonian construction fixed and makes the duration
symbolic.** A parameter `α_j` names a time (in the phase convention above):
the optimizer chooses that time later. One open VHamilton term therefore
represents a parameter-indexed family rather than a single closed Hamiltonian.

## 3. Exact syntax difference: Hamilton → VHamilton

Only the **bare-action grammar** is extended:

```diff
- d ::= id | δ_k | -δ_k | δ_k^j | d ⊕ d
+ d ::= id | δ_k | -δ_k | δ_k^j
+     | δ_k^{α_j} | -δ_k^{α_j} | d ⊕ d
```

The Clifford, action, and process grammars from Section 1 are unchanged.
`δ_k^j` (fixed integer `j`) belongs to **Hamilton**, whereas
`δ_k^{α_j}` and `-δ_k^{α_j}` (symbolic parameter name `α_j`) are the **only
new VHamilton primitives**. Neither a `Scalar` type nor a `scale` operation
appears in the object language.

In Haskell (`src/Syntax.hs`), the distinction is explicit:

```haskell
data BareAction
  = Id
  | Delta Natural
  | NegDelta Natural
  | DeltaPow Natural Integer       -- original Hamilton: δ_k^j
  | DeltaParam Natural String     -- new VHamilton:     δ_k^{α_j}
  | NegDeltaParam Natural String  -- new VHamilton:    -δ_k^{α_j}
  | BareAction :⊕: BareAction
```

A **valuation** `ρ` is provided *outside* the language. The only additional
semantic clauses are

```text
⟦δ_k^{α_j}⟧_ρ  =  ρ(α_j)/2^k
⟦-δ_k^{α_j}⟧_ρ = -ρ(α_j)/2^k.
```

Existing direct sums, conjugations, and parallel sums are interpreted exactly
as before, now relative to `ρ` when parameters occur. Dyadic refinement keeps
its original structure: subdividing time by `2^n` increments the index
`k → k+n`; this also applies to the two new parameterized actions.

For instance, the *same* open program expresses the family `α Z`:

```haskell
zVar = Atom (Bare (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha"))
```

With `ρ(alpha)=1/8` it denotes `(1/8)Z`; with `ρ(alpha)=3/8` it denotes
`(3/8)Z`. Neither valuation changes `zVar`'s syntax.

## 4. From parameterized processes to VQE

The **VQE layer is outside the Hamilton/VHamilton grammar**. It supplies an
initial state, a target Hamiltonian, an ordered list of parameterized evolution
blocks, an energy evaluator, and an optimizer. In particular, `∥` still means
*addition of Hamiltonians*; it does **not** mean circuit sequencing.

At VQE grid precision `r`, the optimizer chooses integer coordinates `q_j`
and supplies the valuation `ρ(α_j) = q_j / 2^r`. Compilation then follows

```text
open VHamilton blocks + valuation ρ
       → Hamilton schedules (Trotter refinement n, if needed)
       → quantum circuits, executed block by block
       → trial state |ψ(q)⟩
       → E(q) = ⟨ψ(q)|H_target|ψ(q)⟩
       → classical search over q.
```

The target Hamiltonian remains a **closed Hamilton** process. The ansatz
blocks are **open VHamilton** processes. For the built-in one-qubit experiment,

```text
H_target = X + Z,             |ψ₀⟩ = |0⟩,
G(α) = αY,
|ψ(α)⟩ = exp(-2πi G(α))|0⟩.
```

`Y` is obtained by Clifford-conjugating a parameterized diagonal `Z` action:

```haskell
pauliYVar = Atom (Conj (S · H)
              (Bare (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha")))
```

The demo searches `q = 0..7`, `r = 4` (`α=q/16`). Its ideal minimum is
`-√2`, attained at `q=5` within this search set. The chemistry experiment
uses three open parameterized blocks and a closed dyadically approximated
two-qubit target; their search and compilation machinery uses the **same**
syntax. Details: [`docs/VHAMILTON.md`](docs/VHAMILTON.md) and
[`docs/VQE_README.md`](docs/VQE_README.md).

Three independent precision knobs should not be confused: `k` sets the dyadic
approximation of a *fixed target*; `r` sets the parameter search grid; `n` sets
the Trotter refinement of the compiled evolution.

## 5. Repository layout

```text
app/                 Haskell executables (Hamilton, VHamilton and VQE)
src/Syntax.hs        original bare actions + two VHamilton actions; valuations
src/Semantics.hs     schedules and dyadic halving
src/Compiler.hs      compilation with/without valuations
src/Examples.hs      Hamilton examples
src/ExamplesVHamilton.hs  parameterized examples
src/ExamplesChem.hs  fixed chemistry target and parameterized generators
src/VQE.hs           external VQE ansatz/optimizer layer
test/                Haskell/Python regression tests
scripts/             independent Python semantics, chemistry, VQE, experiments
docs/                language, experiments, chemistry and optimization guides
results/             paper-facing recorded experiment outputs
pih.cabal            Cabal project definition
requirements.txt     Python requirements
```

For a guide to individual files, see [`docs/FILE_GUIDE.md`](docs/FILE_GUIDE.md).
For the step-by-step language change, see
[`docs/VHAMILTON.md`](docs/VHAMILTON.md); for historical implementation
changes, see [`CHANGES_VHAMILTON.md`](CHANGES_VHAMILTON.md).

## 6. Build, tests, and experiments

Run these commands from the project root:

```bash
cabal build all
cabal run vhamilton-demo              # same open term, multiple valuations
cabal run vqe-demo                   # one-qubit H=X+Z benchmark
cabal run qaoa-demo -- 3 4 5         # one-layer QAOA, MaxCut on one edge (writes qaoa_grid_r*.csv)
cabal run trotter-vqe-demo -- 4 0 1 2 # VQE over Z(alpha) || X(beta) at Trotter levels n (writes trotter_grid_r4_n*.csv)
cabal test regression --test-show-details=direct
python3 -m unittest discover -s test -p 'test_*.py' -v

# Two-qubit chemistry VQE (k, r, n, optimizer)
cabal run chem-vqe -- 10 3 0 exhaustive
cabal run chem-vqe -- 14 7 0 direct

# Reproducible recorded experiment sequences
python3 -m scripts.experiments.run_hamilton
python3 -m scripts.experiments.run_hamilton_vqe
python3 -m scripts.experiments.run_all_experiments
```

Output artifacts live in `results/hamilton/` and
`results/hamilton_vqe/`. See
[`docs/EXPERIMENT_RUN_ORDER.md`](docs/EXPERIMENT_RUN_ORDER.md) for the precise
run sequence and destinations.

**Validation status:** Python/NumPy/SciPy denotational tests were run during
the minimal-syntax update. This environment lacked GHC/Cabal and Qiskit, so
Haskell compilation and end-to-end circuit-backed VQE execution require
verification in a configured environment. This documentation revision does
not claim additional successful end-to-end runs.

See `docs/PAPER_FIGURES.md` for the command behind every figure and number in the PLanQC 2027 abstract.
