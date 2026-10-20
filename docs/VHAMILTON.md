# Deriving VHamilton from Hamilton

This document describes the implementation **as a language extension of
Hamilton**, not as a replacement calculus. The defining design constraint is
that the original Hamilton syntax and semantics remain unchanged except for
two new parameterized **bare actions**.

## 1. Hamilton before parameters

Hamilton has four syntactic categories:

```text
u ::= I | H | S | CX | SWAP | u · u | u ⊗ u

d ::= id | δ_k | -δ_k | δ_k^j | d ⊕ d

a ::= d | u a u†

p ::= a | p ∥ p
```

`u` describes Clifford transformations; `d` constructs diagonal Hamiltonians
from one-dimensional dyadic entries; `a` changes their bases; `p` combines
same-dimensional Hamiltonians additively. The relevant typing conditions are:
Clifford composition and conjugation require equal dimensions, Clifford tensor
products multiply dimensions, direct sums add dimensions, and parallel
composition requires equal dimensions. Each `δ_k`, `-δ_k`, and `δ_k^j` has
dimension 1; `H`, `S`, and `I` have dimension 2; `CX` and `SWAP` have dimension 4.

For `k ∈ ℕ` and `j ∈ ℤ`, bare entries denote

```text
⟦id⟧ = 0
⟦δ_k⟧ = 2^{-k}
⟦-δ_k⟧ = -2^{-k}
⟦δ_k^j⟧ = j · 2^{-k}
⟦d₁ ⊕ d₂⟧ = diag(⟦d₁⟧, ⟦d₂⟧)  (block diagonal)
```

Conjugation denotes `U ⟦a⟧_H U†`, and parallel composition denotes
`⟦p ∥ q⟧_H = ⟦p⟧_H + ⟦q⟧_H`. In this codebase, the associated evolution is
`U(p) = exp(-2πi ⟦p⟧_H)`. For noncommuting sums, the existing compiler uses a
Trotter schedule; these notions are already part of Hamilton.

A fixed `Z` Hamiltonian, followed by the basis-changed `X` Hamiltonian, can
be *described* with the original syntax (they are not sequentially executed by
`∥`):

```haskell
z = Atom (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1)))
x = Atom (Conj H (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1))))
h = x ∥ z  -- X + Z
```

## 2. The problem raised by VQE

VQE repeatedly evaluates trial states prepared with variable evolution
durations. Hamilton's fixed integer-dyadic bare actions describe particular
Hamiltonians, but offer no way to give *one term* a free variational time.
VHamilton adds this ability at the leaves, leaving all Hamilton composition
mechanisms intact.

The intuition is **parameters are time**: a symbolic parameter picks the
amount of evolution in a generator. Here time is normalized so a generator
`G` produces `exp(-2πi tG)` when `t` is supplied as its parameter.

## 3. Minimal syntactic change

Old Hamilton:

```text
d ::= id | δ_k | -δ_k | δ_k^j | d ⊕ d
```

New VHamilton:

```text
d ::= id | δ_k | -δ_k | δ_k^j
    | δ_k^{α_j} | -δ_k^{α_j} | d ⊕ d
```

There are exactly **two additional alternatives**. The other three grammar
categories `u`, `a`, and `p` are unchanged. In particular, this language has
*no* general scalar `s`, multiplication `s · p`, arbitrary superscript
expression, negation-expression constructor, or new sequencing operator.

The Haskell constructors in `src/Syntax.hs` are:

| Origin | Paper syntax | Haskell constructor |
|---|---|---|
| Hamilton | `id` | `Id` |
| Hamilton | `δ_k` | `Delta k` |
| Hamilton | `-δ_k` | `NegDelta k` |
| Hamilton | `δ_k^j` | `DeltaPow k j` |
| **VHamilton** | `δ_k^{α_j}` | `DeltaParam k "alpha_j"` |
| **VHamilton** | `-δ_k^{α_j}` | `NegDeltaParam k "alpha_j"` |
| Hamilton | `d₁ ⊕ d₂` | `d1 :⊕: d2` |

`DeltaPow k j` has an integer numerator `j`; the superscript is **not**
Haskell exponentiation. `DeltaParam` instead stores a **name**. This separation
is important: the new parameter action is not a generalized
`Delta k (ScalarExpression)`.

## 4. Semantics under a valuation

A valuation `ρ : ParameterName → ℝ` is supplied *externally* when interpreting
or compiling an open term. The only new scalar-entry clauses are

```text
⟦δ_k^{α_j}⟧_ρ  =  ρ(α_j) / 2^k
⟦-δ_k^{α_j}⟧_ρ = -ρ(α_j) / 2^k.
```

Direct sum, Clifford conjugation, and parallel composition retain their
Hamilton equations, now parameterized by `ρ`. The resulting unitary family is

```text
H_ρ(p) = ⟦p⟧_{H,ρ}
U_ρ(p) = exp(-2πi H_ρ(p)).
```

For example:

```haskell
zVar = Atom (Bare (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha"))
```

is one syntactic process. Under `ρ(alpha)=1/8`, it denotes `(1/8)Z`; under
`ρ(alpha)=3/8`, it denotes `(3/8)Z`. These are two members of the *same*
family, with no syntax rewriting or newly inserted constants.

In the code, `type Valuation = Map String Double`, and missing parameter names
are rejected when their denotations are needed. The compiler entry point for
open terms is `compileTrotterWith rho n p`; the closed Hamilton entry point
`compileTrotter n p` remains available.

## 5. What the compiler must change

The existing pipeline stays intact:

```text
Hamilton/VHamilton process
    → canonical additive schedule
    → divide each scheduled action by 2^n
    → compile Clifford conjugations and diagonal phases
    → repeat the schedule 2^n times
    → circuit.
```

With parameters, halving preserves the name and sign:

```text
δ_k^{α_j} / 2^n  corresponds to  δ_{k+n}^{α_j}
-δ_k^{α_j} / 2^n corresponds to -δ_{k+n}^{α_j}.
```

The compiler evaluates parameterized phases from `ρ(α_j)` using the same
`-2π` phase convention as the original actions. The new branches appear in
`src/Syntax.hs` (syntax, free-variable collection, JSON),
`src/Semantics.hs` (halving), `src/Compiler.hs` (phase angles and parameterized
compilation), and `scripts/semantics/denot.py` (independent denotation).
No changes to the meaning of `⊕`, Clifford conjugation, or `∥` are needed.

## 6. VQE built on top of the language

A VQE experiment is **not** a new VHamilton program constructor. The runner
in `src/VQE.hs` manages an ordered list of `VBlock`s; each block contains an
ordinary `Proc` and a compilation refinement `n`. The runner supplies a
valuation for the current candidate and concatenates the compiled block
circuits in chronological order. Thus:

- `p ∥ q` forms the **sum of Hamiltonian generators**, not two successive
  evolution blocks.
- `VBlock` order describes **successive circuit evolutions** at the VQE layer.
- `ρ(α_j) = q_j/2^r` describes the VQE grid at precision `r`; `r` is not a
  language constructor or an additional dyadic index inside the program.

For the one-qubit calibration:

```haskell
-- Original Hamilton target: H = X + Z
z = Atom (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1)))
x = Atom (Conj H (Bare (DeltaPow 0 1 ⊕ DeltaPow 0 (-1))))
target = x ∥ z

-- New VHamilton ansatz generator: G(alpha) = alpha Y
yVar = Atom (Conj (S · H)
           (Bare (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha")))
```

The VQE layer prepares `|ψ(α)⟩ = exp(-2πi αY)|0⟩` and evaluates
`E(α) = ⟨ψ(α)|(X+Z)|ψ(α)⟩`. Searching `α=q/16` for `q=0,...,7` reaches the
exact ground energy `-√2` at `q=5` (up to floating-point tolerance).

For the two-qubit chemistry ansatz, the implementation constructs the open
generators

```text
G1(α1) = α1 YI
G2(α2) = (α2/2)(IY + ZY)
G3(α3) = (α3/2)(IY - ZY).
```

Here `ZY`'s negative contribution in `G3` uses the **atomic**
`NegDeltaParam`; it is not constructed with general-purpose scalar negation.
The fixed chemistry target continues to be expressible using Hamilton's
`DeltaPow` bare actions. The exact statevector energy evaluator sits outside
the language.

These examples do not imply that arbitrary VQE energies are classically
tractable: the general statevector and exact expectation calculation may be
exponential in the number of qubits. This project uses small noiseless
benchmarks to validate the language and compiler.

## 7. Files and reproduction

Start with `README.md` in the project root. The main implementation files are
`src/Syntax.hs`, `src/Semantics.hs`, `src/Compiler.hs`,
`src/ExamplesVHamilton.hs`, `src/VQE.hs`, `app/VQEDemo.hs`, and
`app/ChemVQE.hs`. Use `docs/VQE_README.md` for experiment parameters and
`docs/EXPERIMENT_RUN_ORDER.md` for the complete experiment order.

```bash
cabal build all
cabal run vhamilton-demo
cabal run vqe-demo
cabal test regression --test-show-details=direct
python3 -m unittest discover -s test -p 'test_*.py' -v
cabal run chem-vqe -- 10 3 0 exhaustive
cabal run chem-vqe -- 14 7 0 direct
```

**Validation note:** Python denotational regressions passed at the minimal
syntax update. GHC/Cabal and Qiskit were unavailable in the validation
environment, so end-to-end Haskell+Qiskit execution must be checked separately.
