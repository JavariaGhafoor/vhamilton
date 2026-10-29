# Experiment run order: Hamilton, then Hamilton-VQE

Run all commands from the project root with the Python environment activated.
The canonical runners persist every recorded numerical result under `results/`.

## A. Hamilton experiments

Run all Hamilton experiments:

```bash
python -m scripts.experiments.run_hamilton
```

Individual order:

1. **ZX first-order Trotter convergence**, with its convergence table.
   ```bash
   python -m scripts.experiments.run_hamilton --only 01_zx_trotter_convergence
   ```
   Results: `results/hamilton/01_zx_trotter_convergence/`
   - `table_zx_conv.tex`
   - `circuit_zx_n0.tex`
   - `circuit_zx_n1.tex`
   - logs and exact commands

2. **ZX schedule-order equivalence**: two admissible orderings at `n=1`.
   ```bash
   python -m scripts.experiments.run_hamilton --only 02_zx_schedule_orderings
   ```
   Results: `results/hamilton/02_zx_schedule_orderings/table_zx_ord.tex`

3. **Asymmetric ZXY schedule-order ranking**: six orderings at `n=1`.
   ```bash
   python -m scripts.experiments.run_hamilton --only 03_zxy_schedule_orderings
   ```
   Results: `results/hamilton/03_zxy_schedule_orderings/table_zxy_ord.tex`

4. **Exact commuting ZZ interaction**
   ```bash
   python -m scripts.experiments.run_hamilton --only 04_zz_commuting_exact
   ```
   Results: `results/hamilton/04_zz_commuting_exact/`
   - `frob_zzparity.tex`
   - `circuit_zzparity_n0.tex`

5. **DMANH+ dyadic representation error vs `k`**
   ```bash
   python -m scripts.experiments.run_hamilton --only 05_chemistry_dyadic_error
   ```
   Results: `results/hamilton/05_chemistry_dyadic_error/`
   - `table_dyadic.tex`
   - decompositions for `k=3,4,6,8,10,12` and reference `k=20`

6. **DMANH+ Trotter convergence vs `n` at `k=12`**
   ```bash
   python -m scripts.experiments.run_hamilton --only 06_chemistry_trotter_convergence
   ```
   Results: `results/hamilton/06_chemistry_trotter_convergence/`
   - `table_chem_trotter.tex`
   - saved `k=12` decomposition

Every value is therefore regenerated into the
same numbered result folders rather than existing only in terminal output.

---

## B. Hamilton-VQE experiments, in reproducibility order

0. **VHamilton minimal language examples** — open process instantiation, Trotter refinement, and a commuting exact case.
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 00_vhamilton_minimal_examples
   ```
   Results: `results/hamilton_vqe/00_vhamilton_minimal_examples/vhamilton_demo.log`

Run the full Hamilton-VQE sequence:

```bash
python -m scripts.experiments.run_hamilton_vqe
```

1. **One-qubit `X+Z` calibration**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 01_one_qubit_calibration
   ```
   Expected checkpoint: `q*=5`, `E=-sqrt(2)`.
   Results: `results/hamilton_vqe/01_one_qubit_calibration/vqe_demo.log`

2. **Compiler/regression validation**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 02_compiler_regression
   ```
   Results: `results/hamilton_vqe/02_compiler_regression/regression.log`

3. **DMANH+ chemistry target, Pauli decomposition, identity metadata**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 03_chemistry_target_and_identity
   ```
   Results: `results/hamilton_vqe/03_chemistry_target_and_identity/`
   - `physical_pauli_decomposition.csv`
   - `chemistry_target_metadata.json`
   - `dyadic_decomposition_k14.json`

4. **Current ansatz has no active Trotter error**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 04_current_ansatz_trotter_invariance
   ```
   Runs `(k,r)=(10,3)` for `n=0,1,2` with exhaustive search and stores all logs plus
   `trotter_invariance_summary.csv`.

5. **Pre-VQE precision selection**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 05_pre_vqe_precision_selection
   ```
   Targets: `1, 0.1, 0.05, 0.01, 0.005, 0.001 mH`.
   Results: `results/hamilton_vqe/05_pre_vqe_precision_selection/`
   - one full predictor JSON per target
   - one log per target
   - `precision_selection_table.csv` containing the selected precisions

6. **Separate Hamiltonian precision `k` from parameter precision `r`**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 06_k_r_error_separation
   ```
   Runs `(10,3)` exhaustive and `(10,5)` analytic-direct.
   Results: `results/hamilton_vqe/06_k_r_error_separation/`
   - raw compiler-backed logs
   - `k_r_separation_summary.csv`

7. **Full exhaustive `(k,r)` benchmark sweep**: data for the error figures.
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 07_exhaustive_k_r_sweep
   ```
   Range: `k=3..12`, `r=2..6`, conservative target `0.1 mH`.
   Results: `results/hamilton_vqe/07_exhaustive_k_r_sweep/`
   - `chem_vqe_sweep.csv`
   - `summary.json`
   - `hamilton_error_vs_k.png`
   - `variational_gap_vs_r.png`
   - `total_error_heatmap.png`
   - `sweep.log`

8. **Search-method comparison at `(k,r)=(14,7)`**
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 08_search_method_comparison
   ```
   Runs independent exhaustive NumPy reference plus `refine`, `coordinate`, and
   `direct` compiler-backed searches.
   Results: `results/hamilton_vqe/08_search_method_comparison/`
   - raw logs for every method
   - `search_method_comparison.csv`
   - `search_query_comparison.png`

9. **Harmonic-fit conditioning** — conditioning analysis preceding the experimental results.
   ```bash
   python -m scripts.experiments.run_hamilton_vqe --only 09_harmonic_fit_conditioning
   ```
   Results: `results/hamilton_vqe/09_harmonic_fit_conditioning/harmonic_fit_conditioning.json`

10. **Representative analytic-direct runs**
    ```bash
    python -m scripts.experiments.run_hamilton_vqe --only 10_analytic_direct_representative
    ```
    Runs `(10,5)`, `(14,7)`, `(14,10)`.
    Results: `results/hamilton_vqe/10_analytic_direct_representative/`
    - three compiler-backed logs
    - `analytic_direct_representative.csv`

11. **Systematic analytic-direct vs exhaustive validation**
    ```bash
    python -m scripts.experiments.run_hamilton_vqe --only 11_analytic_direct_global_validation
    ```
    Range: `k=3..12`, `r=4,5,6` (30 cases).
    Results: `results/hamilton_vqe/11_analytic_direct_global_validation/`
    - `direct_validation_matrix.csv`
    - `validation_summary_table.csv`
    - `direct_validation_heatmap.png`
    - `validation.log`

12. **Observed direct-search query/sweep scaling**
    ```bash
    python -m scripts.experiments.run_hamilton_vqe --only 12_analytic_direct_precision_scaling
    ```
    Range: `k in {8,10,14}`, `r=4..10`.
    Results: `results/hamilton_vqe/12_analytic_direct_precision_scaling/`
    - `direct_scaling_results.csv`
    - `direct_queries_vs_r.png`
    - `direct_sweeps_vs_r.png`
    - `scaling.log`

### Full experiment sequence

To run Hamilton first and then Hamilton-VQE in one shot:

```bash
python -m scripts.experiments.run_all_experiments
```

This is deliberately the reproducibility order, not the shortest runtime order.
