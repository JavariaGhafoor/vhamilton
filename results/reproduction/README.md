# Independent rerun of the results (2026-10-29)

Every results folder was regenerated from a clean checkout of this branch, on a different machine (x86-64 Linux, GHC 9.4.7, Python 3, Qiskit 2.6.0), and compared with the committed outputs. The Haskell executables were built with `ghc` directly; the only change was replacing `Text.show i` with `Text.pack (show i)`, because this machine's `text` package is older than 2.1.2. `chem_vqe_sweep` was run with a no-op `matplotlib` stub, because it only plots.

## Identical, or equal to floating-point roundoff

| Folder | Result |
|---|---|
| 01 | One-qubit demo: same energies and optimum (q = 5, −√2). Only the printed header differs; it comes from an older version of the demo. |
| 03 | Target decomposition: identical JSON. |
| 04 | (10,3) exhaustive at n = 0, 2: identical. |
| 05 | Precision selection, 6 targets: same (k, r, q); values within 4e-15. |
| 06, 10 | Step 5 at (10,5), (14,7), (14,10): same valuations, energies, sweeps and queries (29, 38, 53). Fit residuals differ at the 1e-16 level. |
| 07 | Exhaustive (k,r) sweep, 50 rows: within 2e-14. |
| 08 | Step 1 (exhaustive) and step 4 (fitted search from a coarse seed, 128 queries): identical. |
| 11, 12 | Validation (30 settings) and scaling (21 runs): within 6e-17. |
| `qaoa_one_edge` | All three grids: identical. |
| `trotter_vqe` | Levels 0–3: all 256 energies identical per level. Level 4 was stopped at 158/256 by a time limit; all 80 logged energies match. An independent NumPy model of the Trotter circuit reproduces all five committed grids (levels 0–4) to 5e-13, the CSV rounding. |
| `trotter_pair` | `run.log` and all four CSVs: bit-identical. |

## Differences caused by exact ties

These come from valuations that prepare the same state. When φ1 = 0 or π/2, one branch of the first qubit is empty, so the clock behind it has no effect, and which tied valuation is reported depends on the last bits of the Qiskit energies.

- `04_k10_r3_n1_rerun.log`: same energy (−0.2021484375), but the reported valuation is (2,2,1) instead of (2,0,1). All four values of α2 tie.
- `08_full_dyadic_beam_rerun.log` (step 2): same result, (28,30,19), but 325 queries instead of 315. A tie at r = 3 is broken differently, so the beam keeps different valuations.
- `08_coordinate_beam_rerun.log` (step 3): same result, but 198 queries instead of 189, from the same tie.
