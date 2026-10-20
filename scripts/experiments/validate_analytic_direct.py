#!/usr/bin/env python3
"""Independent matrix-semantic validation for the noiseless Hamilton VQE study.

This script reproduces the two extra experiment tables used in the comprehensive
writeup:

  1. direct_validation_matrix.csv
     Compare analytic-direct against the exact exhaustive discrete optimum for
     k = 3..12 and r = 4..6.

  2. direct_scaling_results.csv
     Measure analytic-direct query count and sweep count for
     k in {8,10,14} and r = 4..10, without enumerating the huge target grid.

The implementation is deliberately independent of the Haskell/QASM compiler:
it works directly with the 4x4 DMANH+ matrix and the closed-form three-parameter
ansatz state.  This makes it a useful semantic cross-check of the compiler-backed
runs.

Run:
    python -m scripts.experiments.validate_analytic_direct

or choose one experiment:
    python -m scripts.experiments.validate_analytic_direct --experiment validation
    python -m scripts.experiments.validate_analytic_direct --experiment scaling

Choose an output directory with:
    python -m scripts.experiments.validate_analytic_direct --out-dir results
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


# -----------------------------------------------------------------------------
# DMANH+ upper 4x4 block (milliHartree), as used by chem_decompose.py.
# The tiny published asymmetry is symmetrized before Pauli decomposition.
# -----------------------------------------------------------------------------

UPPER_BLOCK_MH = np.array(
    [
        [37.68, -7.537, 0.7182, -0.0024],
        [-7.5367, 7.182, -7.752, 0.8407],
        [0.7182, -7.752, 0.2725, -6.909],
        [-0.0024, 0.8407, -6.909, -7.460],
    ],
    dtype=float,
)
UPPER_BLOCK_MH = 0.5 * (UPPER_BLOCK_MH + UPPER_BLOCK_MH.T)

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI_1Q = {"I": I, "X": X, "Y": Y, "Z": Z}


def pauli_matrix(label: str) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for ch in label:
        out = np.kron(out, PAULI_1Q[ch])
    return out


def pauli_decompose_two_qubit(H: np.ndarray, tol: float = 1e-12) -> Dict[str, float]:
    coeffs: Dict[str, float] = {}
    for a in "IXYZ":
        for b in "IXYZ":
            label = a + b
            P = pauli_matrix(label)
            c = np.trace(P.conj().T @ H) / 4.0
            if abs(c.real) > tol:
                coeffs[label] = float(c.real)
    return coeffs


PHYSICAL_COEFFS = pauli_decompose_two_qubit(UPPER_BLOCK_MH)

# Same scaling convention as the chemistry code: map the largest NON-IDENTITY
# coefficient to 0.25 in the dimensionless Hamiltonian.
MAX_NONIDENTITY = max(abs(c) for p, c in PHYSICAL_COEFFS.items() if p != "II")
SCALE = 0.25 / MAX_NONIDENTITY


def dyadic_hamiltonian(k: int) -> np.ndarray:
    """Return full H_k, INCLUDING the identity coefficient.

    c_P^(k) = round(scale * c_P * 2^k) / 2^k.
    """
    Hk = np.zeros((4, 4), dtype=complex)
    denom = 2**k
    for label, coeff_mh in PHYSICAL_COEFFS.items():
        num = round(SCALE * coeff_mh * denom)
        Hk += (num / denom) * pauli_matrix(label)
    return Hk


# -----------------------------------------------------------------------------
# Closed-form ansatz
#
# G1 = YI
# G2 = (IY + ZY)/2 = |0><0| x Y
# G3 = (IY - ZY)/2 = |1><1| x Y
#
# beta_j = 2*pi*q_j/2^r
# |psi> = [ cos b1 cos b2,
#           cos b1 sin b2,
#           sin b1 cos b3,
#           sin b1 sin b3 ]^T
# -----------------------------------------------------------------------------


def ansatz_state(q: Sequence[int], r: int) -> np.ndarray:
    beta = 2.0 * np.pi * np.asarray(q, dtype=float) / (2**r)
    b1, b2, b3 = beta
    return np.array(
        [
            np.cos(b1) * np.cos(b2),
            np.cos(b1) * np.sin(b2),
            np.sin(b1) * np.cos(b3),
            np.sin(b1) * np.sin(b3),
        ],
        dtype=complex,
    )


def exact_energy(Hk: np.ndarray, q: Sequence[int], r: int) -> float:
    psi = ansatz_state(q, r)
    return float(np.real(np.vdot(psi, Hk @ psi)))


# -----------------------------------------------------------------------------
# Exact exhaustive discrete reference, used only at feasible r <= 6.
# -----------------------------------------------------------------------------


def exhaustive_reference(Hk: np.ndarray, r: int) -> Tuple[Tuple[int, int, int], float]:
    """Return the exact minimum over q_j in {0,...,2^(r-1)-1}."""
    m = 2 ** (r - 1)
    q1, q2, q3 = np.meshgrid(
        np.arange(m), np.arange(m), np.arange(m), indexing="ij"
    )

    b1 = 2.0 * np.pi * q1 / (2**r)
    b2 = 2.0 * np.pi * q2 / (2**r)
    b3 = 2.0 * np.pi * q3 / (2**r)

    psi = np.stack(
        [
            np.cos(b1) * np.cos(b2),
            np.cos(b1) * np.sin(b2),
            np.sin(b1) * np.cos(b3),
            np.sin(b1) * np.sin(b3),
        ],
        axis=-1,
    )

    # H_k is real for this benchmark. einsum evaluates every candidate at once.
    energies = np.einsum(
        "...i,ij,...j->...", psi, Hk.real, psi, optimize=True
    )
    idx = np.unravel_index(int(np.argmin(energies)), energies.shape)
    qstar = tuple(int(x) for x in idx)
    return qstar, float(energies[idx])


# -----------------------------------------------------------------------------
# Analytic-direct coordinate optimizer.
#
# The energy along one coordinate is an exact finite trigonometric polynomial.
# G1 has spectral gaps {0, +-2}:
#   E_1(beta) = c + a2 cos(2 beta) + b2 sin(2 beta).
#
# G2,G3 have spectral gaps {0, +-1, +-2}:
#   E_j(beta) = c + a1 cos(beta) + b1 sin(beta)
#                 + a2 cos(2 beta) + b2 sin(2 beta).
#
# We reconstruct these curves from 3 or 5 exact energy queries, then evaluate
# the fitted model classically on every allowed q along that coordinate.  The
# latter does NOT count as a VQE energy query.  Only calls to exact_energy do.
# -----------------------------------------------------------------------------

BETA_SAMPLES_G1 = np.array([0.0, np.pi / 4.0, np.pi / 2.0])
V_G1 = np.column_stack(
    [
        np.ones(3),
        np.cos(2.0 * BETA_SAMPLES_G1),
        np.sin(2.0 * BETA_SAMPLES_G1),
    ]
)
V_G1_INV = np.linalg.inv(V_G1)

BETA_SAMPLES_G23 = np.array(
    [0.0, np.pi / 4.0, np.pi / 2.0, 3.0 * np.pi / 4.0, 7.0 * np.pi / 8.0]
)
V_G23 = np.column_stack(
    [
        np.ones(5),
        np.cos(BETA_SAMPLES_G23),
        np.sin(BETA_SAMPLES_G23),
        np.cos(2.0 * BETA_SAMPLES_G23),
        np.sin(2.0 * BETA_SAMPLES_G23),
    ]
)
V_G23_INV = np.linalg.inv(V_G23)


@dataclass
class DirectResult:
    q: Tuple[int, int, int]
    energy: float
    queries: int
    sweeps: int


def analytic_direct(Hk: np.ndarray, r: int, max_sweeps: int = 12) -> DirectResult:
    if r < 4:
        raise ValueError("analytic-direct validation uses r >= 4")

    max_q = 2 ** (r - 1)
    beta_grid = 2.0 * np.pi * np.arange(max_q, dtype=float) / (2**r)

    # Deterministic maximin balanced seed:
    # beta_j = pi/4 -> alpha_j = 1/8 -> q_j = 2^(r-3).
    q = [2 ** (r - 3)] * 3

    cache: Dict[Tuple[int, int, int], float] = {}
    query_count = 0

    def query(qvec: Sequence[int]) -> float:
        nonlocal query_count
        key = tuple(int(x) for x in qvec)
        if key not in cache:
            cache[key] = exact_energy(Hk, key, r)
            query_count += 1
        return cache[key]

    # Count the balanced starting-energy evaluation.
    current_energy = query(q)

    # These integer q values correspond exactly to the beta sample points above.
    q_samples_g1 = [0, 2 ** (r - 3), 2 ** (r - 2)]
    q_samples_g23 = [
        0,
        2 ** (r - 3),
        2 ** (r - 2),
        3 * 2 ** (r - 3),
        7 * 2 ** (r - 4),
    ]

    sweeps_done = 0
    for sweep in range(1, max_sweeps + 1):
        q_at_start = tuple(q)

        # ----- q1 / G1: 3-coefficient harmonic model -----
        e_samples: List[float] = []
        for sample_q in q_samples_g1:
            test = q.copy()
            test[0] = sample_q
            e_samples.append(query(test))

        c, a2, b2 = V_G1_INV @ np.asarray(e_samples)
        model = c + a2 * np.cos(2.0 * beta_grid) + b2 * np.sin(2.0 * beta_grid)
        min_model = float(np.min(model))
        q[0] = int(np.flatnonzero(model <= min_model + 1e-12)[0])
        current_energy = query(q)

        # ----- q2 and q3 / G2,G3: 5-coefficient harmonic model -----
        for j in (1, 2):
            e_samples = []
            for sample_q in q_samples_g23:
                test = q.copy()
                test[j] = sample_q
                e_samples.append(query(test))

            c, a1, b1, a2, b2 = V_G23_INV @ np.asarray(e_samples)
            model = (
                c
                + a1 * np.cos(beta_grid)
                + b1 * np.sin(beta_grid)
                + a2 * np.cos(2.0 * beta_grid)
                + b2 * np.sin(2.0 * beta_grid)
            )
            min_model = float(np.min(model))
            q[j] = int(np.flatnonzero(model <= min_model + 1e-12)[0])
            current_energy = query(q)

        sweeps_done = sweep
        if tuple(q) == q_at_start:
            break

    return DirectResult(tuple(q), current_energy, query_count, sweeps_done)


# -----------------------------------------------------------------------------
# Experiments used in the writeup
# -----------------------------------------------------------------------------

VALIDATION_FIELDS = [
    "k",
    "r",
    "direct_q1",
    "direct_q2",
    "direct_q3",
    "ref_q1",
    "ref_q2",
    "ref_q3",
    "q_match",
    "direct_energy_scaled",
    "ref_energy_scaled",
    "energy_diff_scaled",
    "energy_diff_abs_scaled",
    "sweeps",
    "queries",
    "grid_size",
    "reduction",
]

SCALING_FIELDS = [
    "k",
    "r",
    "q1",
    "q2",
    "q3",
    "energy_scaled",
    "queries",
    "sweeps",
    "grid_size",
]


def run_validation_matrix(path: Path) -> None:
    rows = []
    print("\nValidation matrix: analytic-direct vs exhaustive reference")
    print("k = 3..12, r = 4..6\n")

    by_r = {4: [0, 0], 5: [0, 0], 6: [0, 0]}  # [matches,total]

    for k in range(3, 13):
        Hk = dyadic_hamiltonian(k)
        for r in range(4, 7):
            direct = analytic_direct(Hk, r)
            ref_q, ref_e = exhaustive_reference(Hk, r)

            diff = direct.energy - ref_e
            match = direct.q == ref_q
            grid_size = (2 ** (r - 1)) ** 3

            by_r[r][0] += int(match)
            by_r[r][1] += 1

            rows.append(
                {
                    "k": k,
                    "r": r,
                    "direct_q1": direct.q[0],
                    "direct_q2": direct.q[1],
                    "direct_q3": direct.q[2],
                    "ref_q1": ref_q[0],
                    "ref_q2": ref_q[1],
                    "ref_q3": ref_q[2],
                    "q_match": match,
                    "direct_energy_scaled": direct.energy,
                    "ref_energy_scaled": ref_e,
                    "energy_diff_scaled": diff,
                    "energy_diff_abs_scaled": abs(diff),
                    "sweeps": direct.sweeps,
                    "queries": direct.queries,
                    "grid_size": grid_size,
                    "reduction": grid_size / direct.queries,
                }
            )

            print(
                f"k={k:2d} r={r}: direct={direct.q} ref={ref_q} "
                f"match={str(match):5s} Q={direct.queries:2d} L={direct.sweeps}"
            )

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=VALIDATION_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    total_matches = sum(v[0] for v in by_r.values())
    total_cases = sum(v[1] for v in by_r.values())
    print("\nExact global recovery by r:")
    for r in (4, 5, 6):
        print(f"  r={r}: {by_r[r][0]}/{by_r[r][1]}")
    print(f"  overall: {total_matches}/{total_cases}")
    print(f"\nWrote {path}")


def run_scaling(path: Path) -> None:
    rows = []
    print("\nAnalytic-direct scaling experiment")
    print("k in {8,10,14}, r = 4..10\n")

    for k in (8, 10, 14):
        Hk = dyadic_hamiltonian(k)
        for r in range(4, 11):
            direct = analytic_direct(Hk, r)
            grid_size = (2 ** (r - 1)) ** 3

            rows.append(
                {
                    "k": k,
                    "r": r,
                    "q1": direct.q[0],
                    "q2": direct.q[1],
                    "q3": direct.q[2],
                    "energy_scaled": direct.energy,
                    "queries": direct.queries,
                    "sweeps": direct.sweeps,
                    "grid_size": grid_size,
                }
            )

            print(
                f"k={k:2d} r={r:2d}: q={direct.q} "
                f"E={direct.energy:+.12f} Q={direct.queries:2d} "
                f"L={direct.sweeps} grid={grid_size}"
            )

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SCALING_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--experiment",
        choices=("all", "validation", "scaling"),
        default="all",
        help="which experiment to run",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("../../Downloads"),
        help="directory for output CSV files",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"scale = {SCALE:.12f}")
    print(f"cond(V_G1)  = {np.linalg.cond(V_G1):.12f}")
    print(f"cond(V_G23) = {np.linalg.cond(V_G23):.12f}")

    if args.experiment in ("all", "validation"):
        run_validation_matrix(args.out_dir / "direct_validation_matrix.csv")

    if args.experiment in ("all", "scaling"):
        run_scaling(args.out_dir / "direct_scaling_results.csv")


if __name__ == "__main__":
    main()
