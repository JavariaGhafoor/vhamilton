# VQE over a noncommuting block, compiled at Trotter levels n

Ansatz: the single open block `Z(α) ∥ X(β)` on |0⟩ (`app/TrotterVQEDemo.hs`);
target H = X + Y + Z, ground energy −√3 ≈ −1.732.  Grid r = 4 (q ∈ {0,…,15}
per clock, 256 valuations), exhaustive search of the compiled energy E_n at
each level n = 0…4.

Reproduce:

    cabal run trotter-vqe-demo -- 4 0 1 2 3 4                 # summary.log, trotter_grid_r4_n*.csv
    python3 scripts/figures/trotter_vs_exact.py results/trotter_vqe/trotter_grid_r4_n*.csv
                                                             # trotter_vs_exact.txt

## The program's own optimum

The program denotes e^{-2πi(αZ+βX)}: one rotation about an axis between +X and
+Z (both clocks forward).  Its best energy is −√2 ≈ −1.414, on the grid and in
a continuous scan alike (`trotter_vs_exact.txt`).

## What the search reports per level (`summary.log`, `trotter_vs_exact.txt`)

| n | compiled best E_n | at q | program's E there |
|---|---|---|---|
| 0 | −1.414 | (2,11)* | −1.414 |
| 1 | **−1.714** | (2,11) | **−1.016** |
| 2–4 | −1.414 | (0,11)* | −1.414 |

\* ties: several valuations reach −√2; `trotter_vs_exact.txt` reports (0,3).

At n = 1 the compiled circuit √A √B √A √B (A, B the full Z and X turns) is a
different motion from the program's single turn and reaches below anything
the program denotes.

## Re-running the n = 1 optimum at other levels

The compiled energy of ν = (2,11) at each level, read from the Haskell grids
(`grep '^2,11,' trotter_grid_r4_n*.csv`):

| n | 0 | 1 | 2 | 3 | 4 | program (exact) |
|---|---|---|---|---|---|---|
| E_n(2,11) | −1.414 | −1.714 | −1.357 | −1.148 | −1.072 | −1.016 |

From n = 1 on the energies close in on the program's own value at that
valuation, so re-running a returned valuation at higher levels exposes −1.714
as a property of the level-1 circuit rather than of the program.
