# Dyadic coefficients (VHamilton) versus floating-point angles (Qiskit)

Run: `python3 -m scripts.dyadic_vs_float.compare > results/dyadic_vs_float/run.log`
(NumPy/SciPy for part A, Qiskit 2.6.0 for part B).

## A. Rounding a program to k bits: the cost follows from k

Target: the DMANH+ upper block of `scripts/chemistry/chem_decompose.py`, scaled
as there (largest non-identity coefficient 1/4; 1 mHa = 0.019212 scaled units),
10 Pauli terms including II.  Every coefficient c is rounded to q/2^k.

Bounds known from k before running anything:
* each coefficient moves by at most 2^-(k+1);
* the ground energy moves by at most sum |dc| (Weyl's inequality);
* the unit-time evolution e^{-2πiH} moves by at most 2π · sum_{non-I} |dc| (operator norm).

| k | max dc | 2^-(k+1) | ΔE0 (mHa) | Weyl bound (mHa) | ‖ΔU‖ | bound |
|---|---|---|---|---|---|---|
| 4 | 1.6e-2 | 3.1e-2 | 0.85 | 4.5 | 0.26 | 0.50 |
| 8 | 1.9e-3 | 2.0e-3 | 0.071 | 0.38 | 0.020 | 0.038 |
| 14 | 3.0e-5 | 3.1e-5 | 4.4e-4 | 1.0e-2 | 4.4e-4 | 1.1e-3 |
| 20 | 4.6e-7 | 4.8e-7 | 2.3e-6 | 9.9e-5 | 5.4e-6 | 1.1e-5 |

(all even k from 4 to 20 in `run.log`).  At k = 14, the precision used in the
Sec. 5 experiments, the ground energy moves by 4.4e-4 mHa, far below chemical
accuracy (1.6 mHa); this matches `hamilton_vqe/05_pre_vqe_precision_selection`.

## B. Qiskit: radians as doubles, approximations chosen by the transpiler

* Angles are doubles in radians.  π is irrational, so not even T's angle π/4
  is stored exactly (it is stored as 884279719003555/2^50).
* Transpiling `rz(θ)` to the Clifford+T basis `[h, s, sdg, t, tdg, x, z, cx]`
  (optimization level 1) matches θ against π/4 within a tolerance:
  θ = π/4 + 1e-9 and π/4 + 1e-6 both become a single T, silently dropping
  errors of 5e-10 and 5e-7.  θ = π/4 + 1e-3 is instead synthesised into
  125 T gates with error 2e-13.
* The dyadic angles 3π/8 (coefficient 3/32, k = 5) and 3π/2^13 (3/2^14,
  k = 14) are each synthesised into about 125 T gates (error ≈ 2e-13).

In VHamilton a Pauli rotation e^{-2πi (q/2^k) P} is exactly a Clifford+T
circuit iff q/2^k (in lowest terms) has k ≤ 4 (Z_4 denotes T): this can be
read off the program, instead of being decided by a tolerance at compile time.
