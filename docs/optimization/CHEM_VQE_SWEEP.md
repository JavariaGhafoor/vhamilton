# Chemistry VQE precision sweep

`scripts/experiments/chem_vqe_sweep.py` performs the systematic noiseless experiment over the two
finite precision knobs currently present in the chemistry VQE:

- `k`: dyadic precision of the Hamiltonian coefficients,
- `r`: dyadic precision of the three variational parameters.

For each pair `(k,r)` it computes

- the exact ground energy `E0(H_k)` of the dyadically rounded Hamiltonian,
- the discrete VQE optimum `E*_{k,r}` over the full finite parameter grid,
- the Hamiltonian discretization error,
- the variational gap,
- the total physical-energy error,
- the minimizing integer parameter triple `q`.

The current ansatz is the same one used by `chem-vqe`:

```
G1 = YI
G2 = (IY + ZY)/2
G3 = (IY - ZY)/2
alpha_j = q_j / 2^r
q_j = 0, ..., 2^(r-1)-1
```

The number of candidate states is therefore

```
2^(3(r-1)).
```

The script is intentionally a fast matrix-semantics sweep. Candidate states depend
on `r` but not on `k`, so each state grid is generated once and reused for every
Hamiltonian precision. Selected points should still be validated with the
compiler-backed executable, e.g. `cabal run chem-vqe -- 10 3 0`.

## Run

Install the additional plotting dependency once:

```bash
python -m pip install -r requirements.txt
```

Default sweep (`k=3..12`, `r=2..6`):

```bash
python -m scripts.experiments.chem_vqe_sweep
```

Specify a physical-energy tolerance, for example 0.1 mH:

```bash
python -m scripts.experiments.chem_vqe_sweep --target-error-mh 0.1
```

Custom range:

```bash
python -m scripts.experiments.chem_vqe_sweep --k-min 4 --k-max 14 --r-min 2 --r-max 7
```

## Outputs

By default, files are written under `chem_vqe_sweep_results/`:

- `chem_vqe_sweep.csv`: one row for every `(k,r)` pair,
- `summary.json`: lowest-error pair and, when requested, the cheapest pair meeting the tolerance,
- `hamilton_error_vs_k.png`,
- `variational_gap_vs_r.png`,
- `total_error_heatmap.png`.

"Cheapest" is defined explicitly as the smallest number of exhaustive VQE candidate
evaluations; ties are broken by lower `k`, then lower `r`. This is a computational
search-cost definition, not a hardware gate-cost model.

## Error definitions

With physical ground energy `E0(H)` and scale `s`, the script reports

```
Hamilton error    = |E0(H_k)/s - E0(H)|
Variational gap   = E*_{k,r}/s - E0(H_k)/s
Total abs. error  = |E*_{k,r}/s - E0(H)|
```

The Hamiltonian and variational errors are reported separately because their signed
contributions can partially cancel in the total physical-energy error.

### Avoiding accidental error cancellation

The CSV also contains

```
noncancelling_error_budget_mh = Hamilton error + variational gap.
```

This is useful because a coarse `H_k` can shift the exact ground energy in the
opposite direction from the positive variational gap, making the total signed VQE
error look artificially small. The script therefore reports both the pair with the
smallest actual total absolute energy error and the pair with the smallest
non-cancelling component-wise error budget. If `--target-error-mh` is supplied, it
reports cheapest pairs under both criteria.

## Pre-selecting `k` and `r` without VQE

For the intended workflow, run `scripts/chemistry/chem_decompose_vqe.py` first. It predicts a
sufficient near-best `(k,r)` from the requested error budget without enumerating
or optimizing over the VQE grid. `scripts/experiments/chem_vqe_sweep.py` should then be treated as a
validation experiment rather than the precision-selection algorithm.

```bash
python -m scripts.chemistry.chem_decompose_vqe --target-error-mh 0.1
cabal run chem-vqe -- <predicted-k> <predicted-r> 0
```
