# VHamilton in Agda: equivalence of variational programs

`VHamilton.agda` formalizes a fragment of VHamilton with a symbolic
Hamiltonian semantics in the Pauli basis, and proves two variational
programs equivalent for every valuation and every clock level.

Check (Agda 2.6.3, agda-stdlib 1.7.3):

    LC_ALL=C.UTF-8 agda -i . -i /usr/share/agda-stdlib VHamilton.agda

Contents:

* `Pauli`, `PString`: Pauli strings, the basis in which Hamiltonians are read.
* `Cliff`: transformations `I₁ H S CX SWAP`, composition `·`, tensor `⊗`;
  `cj` conjugates a Pauli string, giving a sign and a string (tables checked
  against the matrices with NumPy).
* `Scalar`, `Diag`, `Act`, `Proc`: syntax with fixed clocks `±δ k`, open
  clocks `±δ^ α`, qubit-sized direct sums `⊕`, conjugation `conj u a`, and `∥`.
* `⟦_⟧`: semantics as a list of Pauli terms whose coefficients are affine
  forms in the clock values; `d₀ ⊕ d₁ = 𝟙⊗(d₀+d₁)/2 + Z⊗(d₀−d₁)/2`.
* `_≃_`: equal Pauli coefficients as affine forms (Hamiltonian equivalence); `sound` turns it into
  equality under every valuation `ν` and subscript assignment `k`.
* `ZZ≃CX-IZ-CX`: `ZZ(γ) ≃ CX (I⊗Z(γ)) CX†`, all 16 Pauli strings by four
  `refl` clauses, one per first Pauli; `ZZ-equiv` is the statement for every valuation and level.

Scope: each variable is used at one level; the matrix semantics and the
compiler are not formalized.  Changing one entry of `IZ` makes the proof fail
(`0 != 1 of type ℕ` at the string `𝟙𝟙`), so the check is not vacuous.
