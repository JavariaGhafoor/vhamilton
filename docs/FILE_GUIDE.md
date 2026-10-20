# File guide

## Haskell core

- `src/Syntax.hs` — Hamilton syntax plus two parameterized bare actions, external valuations, dimensional typing, and JSON encoding.
- `src/Semantics.hs` — semantic definitions.
- `src/Compiler.hs` — Hamilton-to-circuit compiler.
- `src/Qasm.hs` — QASM generation and diagonal-gate lowering interface.
- `src/Denot.hs` — bridge to independent Python semantic checks.
- `src/Examples.hs` — generic examples.
- `src/ExamplesChem.hs` — DMANH+ chemistry Hamiltonian bridge and parameterized Pauli constructors.
- `src/ExamplesVHamilton.hs` — minimal open-process VHamilton examples.
- `src/VQE.hs` — VQE sequencing/valuation layer over open VHamilton processes.
- `src/Quantikz.hs` — circuit rendering helpers.

## Haskell executables

- `app/VHamiltonDemo.hs` — minimal VHamilton language examples.
- `app/ChemVQE.hs` — chemistry VQE executable; includes exhaustive, refine,
  coordinate, analytic, and analytic-direct search modes.
- `app/VQEDemo.hs` — small VQE demonstration.
- `app/GenTables.hs` — table generation.
- `app/GenCircuits.hs` — circuit generation.

## Python semantics

- `scripts/semantics/denot.py` — independent Hamiltonian/unitary semantics and QASM interpretation.
- `scripts/semantics/compare_ideals.py` — compare ideal unitaries.
- `scripts/semantics/pauli_semantics_check.py` — regression helper for Pauli semantics.

## Python compiler utilities

- `scripts/compiler/diag_to_qasm3.py` — lower diagonal unitary gates to elementary QASM-compatible gates.
- `scripts/compiler/runq.py` — load an OpenQASM 3 file and print its unitary.

## Python chemistry

- `scripts/chemistry/chem_decompose.py` — DMANH+ Pauli decomposition and dyadic rounding.
- `scripts/chemistry/chem_decompose_vqe.py` — pre-VQE k/r precision selector for this benchmark.
- `scripts/chemistry/verify_chem.py` — independent chemistry semantics check.

## Python VQE

- `scripts/vqe/vqe_eval.py` — exact statevector energy oracle used by the Haskell VQE code.
- `scripts/vqe/verify_chem_vqe.py` — independent exhaustive matrix-semantic VQE reference.
- `scripts/vqe/verify_chem_vqe_analytic.py` — independent analytic/analytic-direct VQE reference.

## Experiments

- `scripts/experiments/chem_vqe_sweep.py` — fast matrix-semantic k/r sweep.
- `scripts/experiments/validate_analytic_direct.py` — global-optimum validation matrix and query-scaling experiments.
- `results/hamilton/` — Hamilton-paper outputs, numbered in paper order.
- `results/hamilton_vqe/` — Hamilton-VQE outputs, numbered in document order.

## Documentation

- `docs/VHAMILTON.md` — VHamilton syntax, semantics, examples, and VQE mapping.
- `docs/VQE_README.md` — VQE implementation overview.
- `docs/chemistry/CHEM_VQE.md` — chemistry VQE details.
- `docs/chemistry/CHEM_DECOMPOSE_VQE.md` — pre-VQE precision selection.
- `docs/optimization/CHEM_VQE_REFINEMENT.md` — dyadic refinement, coordinate, and analytic search.
- `docs/optimization/CHEM_VQE_SWEEP.md` — k/r sweep experiment.
- `docs/testing/REGRESSION_TESTS.md` — compiler/semantics regression suite.

## Document experiment runners

- `scripts/experiments/run_hamilton.py` — reproduces the numerical/figure artifacts from the Hamilton paper in paper order and writes them to `results/hamilton/`.
- `scripts/experiments/run_hamilton_vqe.py` — reproduces the experiments/results from the Hamilton-VQE writeup in document order and writes them to `results/hamilton_vqe/`.
- `scripts/experiments/run_all_experiments.py` — runs the Hamilton sequence first and then the Hamilton-VQE sequence.
- `scripts/experiments/_experiment_utils.py` — shared logging, command-recording, CSV and compiler-output parsing helpers.
- `docs/EXPERIMENT_RUN_ORDER.md` — canonical ordered checklist with commands and result destinations.
- `results/README.md` — result-folder map.
