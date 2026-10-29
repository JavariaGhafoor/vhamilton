# 13 Error calculus: measured errors against the a priori bounds

Checks the bounds of the paper's error calculus (Sec. 7.4) on the DMANH+ target, and draws its figure.

- `hamiltonian_rounding.csv`: δ_H(k) = E0(H_k) − E0(H̄) for k = 3..14, against the rounding bound M·2^−(k+1)/s, with M = 10 Pauli terms.
- `variational_gap.csv`: δ_r(k,r) = (grid minimum) − E0(H_k) for k ∈ {6, 10, 14} and r = 2..8, computed by exhaustive search over the one-period grid. It is compared against the first-order bound 2πp‖H_k‖2^−r and the second-order bound 2π²p²‖H_k‖4^−r, with p = 3.
- `error_calculus.tex`: the pgfplots figure, copied to the paper's `figures/`.

All values are in mHa. The script asserts that every bound holds and that δ_r is nonincreasing in r. δ_r(10,3) and δ_r(10,5) agree with `06_k_r_error_separation`; δ_H agrees with `07_exhaustive_k_r_sweep`.
