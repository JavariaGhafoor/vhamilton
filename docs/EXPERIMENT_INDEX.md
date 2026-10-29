# Experiments behind each recorded figure and number

Every recorded number and data-driven figure is produced by the commands below (run from this
repository's root; the energy oracle needs the Python packages in
`requirements.txt`).

| Item | Command | Output |
|---|---|---|
| Fixed-clock identities (Z_k ladder; ZZ_3 = CZ then S⊗S up to phase) | `python3 scripts/figures/verify_entangler_identity.py` | `results/figures/entangler_identity_check.txt` |
| QAOA run, precision ladder r = 3, 4, 5 | `cabal run qaoa-demo -- 3 4 5` | `results/qaoa_one_edge/qaoa_grid_r{3,4,5}.csv`, `summary.log` |
| QAOA ladder TikZ source | `python3 scripts/figures/qaoa_ladder_tikz.py --size 1.75 --gap 0.4 results/qaoa_one_edge/qaoa_grid_r{3,4,5}.csv` | `results/figures/qaoa_ladder_small.tex` (full-size version: `qaoa_ladder.tex`) |
| Closed form E = sin(4πγ) sin(8πβ) | `python3 scripts/figures/verify_qaoa_closed_form.py results/qaoa_one_edge/qaoa_grid_r*.csv` | `results/qaoa_one_edge/closed_form_check.txt` |
| Earlier single-panel heatmap (superseded) | `python3 scripts/figures/qaoa_heatmap_tikz.py results/qaoa_one_edge/qaoa_grid_r5.csv` | `results/figures/qaoa_heatmap.tex` |
| Trotterized VQE table | `cabal run trotter-vqe-demo -- 4 0 1 2 3 4` | `results/trotter_vqe/trotter_grid_r4_n*.csv`, `summary.log` |
| Compiled vs exact landscapes for the Trotter run | `python3 scripts/figures/trotter_vs_exact.py results/trotter_vqe/trotter_grid_r4_n*.csv` | `results/trotter_vqe/trotter_vs_exact.txt` |
| Chemistry search staircase: exhaustive / refine / coordinate / analytic-direct | see `results/hamilton_vqe/08_search_method_comparison/COMMAND_*.txt` | `search_method_comparison.csv` and logs |
| Chemistry search staircase: analytic with coarse seed (128 queries) | `cabal run chem-vqe -- 14 7 0 analytic` | `results/hamilton_vqe/08_search_method_comparison/analytic.log` |
| Chemistry precision scaling and validation numbers | (existing experiments) | `results/hamilton_vqe/10_*`, `11_*`, `12_*` |

## New executables

* `qaoa-demo` (`app/QAOADemo.hs`): one-layer QAOA for MaxCut on a single edge.
  Closed preparation block (YI+IY)/8, open cost block γ ZZ, open mixer
  β(XI+IX); exhaustive grid over one period at each precision r given on the
  command line.
* `trotter-vqe-demo` (`app/TrotterVQEDemo.hs`): one open noncommuting block
  Z(α) ∥ X(β) on |0⟩ with target X+Y+Z, compiled at each Trotter level n given
  on the command line, so the optimizer only ever sees compiled energies.
  At n = 1 the compiled optimum (−1.714) lies below anything the program's own
  trajectories reach (−√2), while that valuation's exact energy is −1.016;
  the mean gap between compiled and exact landscapes falls 1.04, 0.64, 0.53,
  0.24, 0.10 for n = 0..4.

## Build note

The library uses `Data.Text.show`, which needs `text >= 2.1.2` (GHC ≥ 9.6's
boot package). With an older `text`, replace `Text.show i` by
`pack (show i)` in `src/Syntax.hs`; the results above were produced that way
with GHC 9.4.7, with no other change to the code.
