# Results

Outputs are split by experiment suite and numbered in run order:

```text
results/
├── hamilton/
│   ├── 01_zx_trotter_convergence/
│   ├── 02_zx_schedule_orderings/
│   ├── 03_zxy_schedule_orderings/
│   ├── 04_zz_commuting_exact/
│   ├── 05_chemistry_dyadic_error/
│   └── 06_chemistry_trotter_convergence/
└── hamilton_vqe/
    ├── 01_one_qubit_calibration/
    ├── 02_compiler_regression/
    ├── 03_chemistry_target_and_identity/
    ├── 04_current_ansatz_trotter_invariance/
    ├── 05_pre_vqe_precision_selection/
    ├── 06_k_r_error_separation/
    ├── 07_exhaustive_k_r_sweep/
    ├── 08_search_method_comparison/
    ├── 09_harmonic_fit_conditioning/
    ├── 10_analytic_direct_representative/
    ├── 11_analytic_direct_global_validation/
    └── 12_analytic_direct_precision_scaling/
```

Use the runners in `scripts/experiments/`. Each runner stores the exact command,
raw log, and structured/plot output needed by the corresponding result.
See `docs/EXPERIMENT_RUN_ORDER.md`.
