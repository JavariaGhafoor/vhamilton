# Pre-VQE precision selection: `scripts/chemistry/chem_decompose_vqe.py`

`scripts/chemistry/chem_decompose_vqe.py` selects the finite dyadic precisions `k` and `r`
**before running VQE**.

It is distinct from `scripts/experiments/chem_vqe_sweep.py`:

- `scripts/chemistry/chem_decompose_vqe.py`: classical preprocessing / prediction; **no VQE grid search**.
- `scripts/experiments/chem_vqe_sweep.py`: exhaustive experiment used afterward to validate the prediction.

## What is predicted

For a requested physical-energy budget `eps_target`, the predictor searches
candidate `(k,r)` values and uses the non-cancelling budget

```text
eps_H(k) + eps_r(k,r) <= eps_target.
```

### Hamiltonian precision `k`

For every candidate `k`, it builds the dyadically rounded Hamiltonian `H_k` and
classically diagonalizes the 4x4 benchmark matrix. It computes

```text
eps_H(k) = |E0(H_k) - E0(H)|.
```

This is ordinary benchmark preprocessing, not VQE.

### Parameter precision `r`

The current ansatz is

```text
G1 = YI
G2 = (IY + ZY)/2 = |0><0| x Y
G3 = (IY - ZY)/2 = |1><1| x Y
```

and its continuous parameters can prepare any real normalized two-qubit state.
The chemistry block is real symmetric, so the exact ground eigenvector of each
`H_k` can be inverted analytically to continuous ansatz parameters `alpha*`.

For each candidate `r`, `alpha*` is rounded once to the nearest available
`q/2^r`. The script evaluates only that **single constructed state**. It never
enumerates the parameter grid and never minimizes over `q`.

Its energy excess

```text
eps_r(k,r) = E(rounded alpha*; H_k) - E0(H_k)
```

is a constructive upper bound on the later discrete-VQE variational gap: a VQE
optimizer can only find an energy less than or equal to this candidate energy.

The recommended pair is the feasible pair with the lowest later exhaustive VQE
candidate count `2^(3(r-1))`, with ties broken by smaller `k` and then `r`.

## Run

For a `0.1 mH` target:

```bash
python -m scripts.chemistry.chem_decompose_vqe --target-error-mh 0.1
```

The prediction is also written to:

```text
chem_vqe_precision.json
```

Then run VQE only at the recommended precisions:

```bash
cabal run chem-vqe -- <k> <r> 0
```

The full sweep remains useful as an *after-the-fact validation* of the
prediction:

```bash
python -m scripts.experiments.chem_vqe_sweep --target-error-mh 0.1
```

## Scope

The near-best `r` predictor deliberately uses the exact 4x4 benchmark ground
state and the exact expressibility of this particular three-parameter real-state
ansatz. This is appropriate for validating the discrete-VQE methodology on the
small DMANH+ block, but it is not a scalable method for large Hamiltonians where
the ground eigenvector is not classically available. A large-system version
would need a priori norm/sensitivity bounds or another source of information
about the ansatz energy landscape.
