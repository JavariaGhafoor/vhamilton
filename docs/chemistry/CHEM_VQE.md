# Discrete chemistry VQE

This extension moves the VQE prototype from the analytically solved one-qubit
`X + Z` example to the same two-qubit DMANH+ upper block used by the Hamilton
chemistry application.

## What changed

`scripts/chemistry/chem_decompose.py` now retains the identity coefficient as metadata and rounds
it to the same dyadic level `k` as the other Pauli coefficients.  This matters
for VQE because

```
<H + c I> = <H> + c,
```

whereas the identity term may be dropped from time evolution because it only
adds a global phase.

The Haskell chemistry layer now provides two targets:

- `chemProc d`: non-identity Hamiltonian, appropriate for the original dynamics
  / Trotter experiments.
- `chemProcWithIdentity d`: full dyadic Hamiltonian, appropriate for VQE and
  absolute-energy estimation.

## Ansatz

`chem-vqe` uses three open VHamilton process blocks:

```
G1 = YI
G2 = (IY + ZY)/2 = |0><0| ⊗ Y
G3 = (IY - ZY)/2 = |1><1| ⊗ Y
```

with free parameter names `alpha_j` in the two parameterized bare-action forms. At search precision `r`, the VQE layer supplies the valuation

```
rho(alpha_j) = q_j / 2^r.
```

The same VHamilton process terms are therefore reused at every `r`; only the classical valuation grid changes.

Starting from `|00>`, the first block chooses the weight of the first-qubit
branches, and the two conditional blocks independently rotate the second qubit
inside the `|0>` and `|1>` branches.  In the continuous-parameter limit this
family can prepare any real normalized two-qubit state, which is a natural
match for the real-symmetric DMANH+ block.

The terms inside `G2` commute with one another, as do the terms inside `G3`, so
this first chemistry VQE experiment has no Trotter error even at `n=0`.  This is
deliberate: it isolates Hamiltonian discretization `k` and variational
parameter discretization `r` before introducing a noncommuting/Trotterized
ansatz.

## Run

First run the fast independent matrix-semantics reference:

```bash
python -m scripts.vqe.verify_chem_vqe 10 3
```

For `k=10, r=3`, it should report approximately:

```text
best q                  = [2, 0, 1]
full VQE energy scaled  = -0.202148437500
ground H_k scaled       = -0.234070730182
variational gap scaled  = 3.192229268194e-02
VQE energy (mH)         = -10.521725097656
ground H_k (mH)         = -12.183264470605
physical exact (mH)     = -12.162810425974
```

Then run the compiler-backed VHamilton implementation:

```bash
cabal run chem-vqe
```

or explicitly:

```bash
cabal run chem-vqe -- 10 3 0
```

The positional arguments are:

```text
k = Hamiltonian dyadic precision
r = variational-parameter dyadic precision
n = Trotter refinement
```

The default `r=3` searches four values per parameter (`q=0,1,2,3`), or 64
candidate states total.  The compiler-backed version currently crosses the
Haskell/Python boundary repeatedly, so larger grids can take substantially
longer.  Use the Python verifier for quick sweeps and the Haskell executable for
compiler-backed validation until the evaluator is batched.

## Quantities reported

The executable distinguishes:

- **variational gap**: `E_VQE(H_k) - E_0(H_k)`, caused by the finite ansatz/grid;
- **Hamiltonian discretization shift**: `E_0(H_k)/scale - E_0(H_physical)`;
- **total signed VQE error**: `E_VQE/scale - E_0(H_physical)`.

No sampling or hardware noise is present yet; `scripts/vqe/vqe_eval.py` remains an exact
statevector oracle.

## Search modes

The compiler-backed executable now supports five search strategies over the
same discrete ansatz:

```bash
cabal run chem-vqe -- 14 7 0 refine
cabal run chem-vqe -- 14 7 0 coordinate
cabal run chem-vqe -- 14 7 0 analytic
cabal run chem-vqe -- 14 7 0 analytic-direct
cabal run chem-vqe -- 14 7 0 exhaustive
```

`analytic` exploits the known generator spectra to reconstruct each
one-coordinate energy curve from three samples for `G1` or five samples for
`G2/G3`.  It then jumps directly from the fixed `r0=3` seed to the requested
`r` and performs analytic coordinate sweeps there.  `analytic-direct` removes
the coarse exhaustive seed entirely: for `r>=3` it starts deterministically at
`q_j=2^(r-3)`, corresponding to `alpha_j=1/8` and the balanced angle
`beta_j=pi/4`, which prepares `|++>` for this ansatz.  See
`CHEM_VQE_REFINEMENT.md` for the derivation, complexity comparison, and
limitations.


> **Syntax update (2026-10-15):** The former generalized `delta_k^x`/scalar-expression representation has been replaced by the minimal `delta_k^{alpha_j}` and `-delta_k^{alpha_j}` bare actions. Historical statements above about scalars describe the old implementation. See `docs/VHAMILTON.md` (from the project root) for the current syntax.
