#!/usr/bin/env python3
"""Predict dyadic precisions (k, r) for the DMANH+ discrete VQE before VQE.

This script performs only classical preprocessing.  It does *not* enumerate the
VQE parameter grid and does *not* run a variational optimizer.

For a requested physical-energy budget eps_target it searches candidate
Hamiltonian precisions k and parameter precisions r and recommends the cheapest
pair whose non-cancelling predicted error budget fits inside eps_target.

The two error contributions are estimated as follows.

1. Hamiltonian dyadic error

   Construct H_k by rounding every scaled Pauli coefficient to a dyadic with
   denominator 2^k, diagonalize the small 4x4 benchmark matrix, and compute

       eps_H(k) = |E0(H_k)/s - E0(H)|.

   This is exact classical benchmark preprocessing, not VQE.

2. Variational parameter discretization

   The current three-generator ansatz is

       G1 = YI
       G2 = (IY + ZY)/2 = |0><0| x Y
       G3 = (IY - ZY)/2 = |1><1| x Y,

       U(alpha) = exp(-2pi i alpha3 G3)
                  exp(-2pi i alpha2 G2)
                  exp(-2pi i alpha1 G1).

   For real amplitudes this continuous ansatz can prepare any normalized
   two-qubit real state (up to global phase).  Since the benchmark H_k is real
   symmetric, its ground eigenvector can therefore be converted analytically to
   continuous parameters alpha*.  For each r, alpha* is rounded independently
   to the nearest available dyadic q/2^r.  The energy of this *single
   constructed state* gives a sufficient upper bound on the best discrete-VQE
   gap:

       eps_r(k,r) = E(round_r(alpha*); H_k)/s - E0(H_k)/s >= 0.

   No grid search is performed.  A later VQE can only do at least as well as
   this constructed candidate, so eps_r is a conservative constructive bound
   for this benchmark ansatz.

The recommendation minimizes the number of candidate states that a later full
exhaustive VQE would need, 2^(3(r-1)), subject to

    eps_H(k) + eps_r(k,r) <= eps_target.

Ties are broken by lower k and then lower r.

Important scope note
--------------------
The near-best r prediction exploits exact diagonalization of this 4x4 benchmark
and the fact that this particular ansatz is exactly expressive for real
2-qubit states.  It is therefore a benchmark precision-selection method, not a
scalable replacement for VQE on large Hamiltonians.  The key property for the
current research workflow is that it chooses k and r *before* the VQE search.

Example
-------
    python -m scripts.chemistry.chem_decompose_vqe --target-error-mh 0.1

Then validate the predicted point with the compiler-backed implementation:

    cabal run chem-vqe -- <k> <r> 0
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.chemistry.chem_decompose import HTILDE, PAULI, pauli_decompose


def pmat(label: str) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for c in label:
        out = np.kron(out, PAULI[c])
    return out


PAULI_MATS = {a + b: pmat(a + b) for a in "IXYZ" for b in "IXYZ"}


def physical_block(block: str) -> np.ndarray:
    h = HTILDE[:4, :4] if block == "upper" else HTILDE[4:, 4:]
    return (h + h.T) / 2.0


def exact_decomposition(block: str, target: float):
    h = physical_block(block)
    coeffs = pauli_decompose(h)
    identity = float(coeffs.get("II", 0.0))
    nonidentity = {lab: float(c) for lab, c in coeffs.items() if lab != "II"}
    maxc = max(abs(v) for v in nonidentity.values())
    scale = target / maxc
    evals = np.linalg.eigvalsh(h)
    physical_ground = float(evals[0].real)
    return h, identity, nonidentity, scale, physical_ground


def dyadic_hamiltonian(
    k: int,
    identity: float,
    nonidentity: dict[str, float],
    scale: float,
):
    denom = 2**k
    identity_num = round(scale * identity * denom)
    h = (identity_num / denom) * PAULI_MATS["II"]
    rounded = {"II": identity_num / denom}

    for label, coeff in nonidentity.items():
        num = round(scale * coeff * denom)
        value = num / denom
        if num != 0:
            h = h + value * PAULI_MATS[label]
            rounded[label] = value
        else:
            rounded[label] = 0.0

    return h, rounded


def make_real_ground_vector(v: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    """Remove arbitrary global phase and return a real normalized eigenvector."""
    v = np.asarray(v, dtype=complex)
    pivot = int(np.argmax(np.abs(v)))
    if abs(v[pivot]) > tol:
        v = v * np.exp(-1j * np.angle(v[pivot]))

    if np.linalg.norm(v.imag) > tol:
        raise ValueError(
            "Ground eigenvector is not real up to global phase; the current "
            "three-rotation constructive predictor does not apply."
        )

    out = v.real
    out = out / np.linalg.norm(out)
    return out


def state_to_continuous_alphas(v: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """Invert the current real-state ansatz analytically.

    With beta_j = 2*pi*alpha_j, the prepared amplitudes are

        [cos b1 cos b2,
         cos b1 sin b2,
         sin b1 cos b3,
         sin b1 sin b3].

    beta_j are chosen in [0, pi], hence alpha_j in [0, 1/2].
    """
    v = np.asarray(v, dtype=float)
    v = v / np.linalg.norm(v)

    # We are free to choose global sign.  Make the |11> amplitude nonnegative
    # whenever possible because sin(beta1), sin(beta3) are nonnegative on
    # beta in [0, pi].
    if v[3] < -tol:
        v = -v

    a, b, c, d = v
    n0 = float(np.hypot(a, b))
    n1 = float(np.hypot(c, d))

    if n0 < tol:
        beta1 = np.pi / 2.0
        beta2 = 0.0
    else:
        # Choose the signed first-branch amplitude so that sin(beta2) >= 0.
        if abs(b) > tol:
            c1 = float(np.copysign(n0, b))
        elif abs(a) > tol:
            c1 = float(np.copysign(n0, a))
        else:
            c1 = n0

        beta1 = float(np.arctan2(n1, c1))
        beta2 = float(np.arctan2(b / c1, a / c1))
        if beta2 < 0:
            beta2 += 2.0 * np.pi

    if n1 < tol:
        beta3 = 0.0
    else:
        beta3 = float(np.arctan2(d / n1, c / n1))
        if beta3 < 0:
            beta3 += 2.0 * np.pi

    alphas = np.array([beta1, beta2, beta3], dtype=float) / (2.0 * np.pi)

    if np.any(alphas < -1e-10) or np.any(alphas > 0.5 + 1e-10):
        raise ValueError(f"Constructed alpha outside [0, 1/2]: {alphas}")

    return np.clip(alphas, 0.0, 0.5)


def state_from_alphas(alphas: np.ndarray) -> np.ndarray:
    """Construct the ansatz state analytically, without Qiskit or VQE."""
    b1, b2, b3 = 2.0 * np.pi * np.asarray(alphas, dtype=float)
    return np.array(
        [
            np.cos(b1) * np.cos(b2),
            np.cos(b1) * np.sin(b2),
            np.sin(b1) * np.cos(b3),
            np.sin(b1) * np.sin(b3),
        ],
        dtype=complex,
    )


def round_alphas_to_grid(alphas: np.ndarray, r: int):
    if r < 1:
        raise ValueError("r must be at least 1")

    denom = 2**r
    grid_size = 2 ** (r - 1)  # q = 0, ..., 2^(r-1)-1 => alpha in [0, 1/2)

    q = np.rint(np.asarray(alphas) * denom).astype(int)
    q = np.clip(q, 0, grid_size - 1)
    rounded = q.astype(float) / denom
    return q, rounded, grid_size


def energy(psi: np.ndarray, h: np.ndarray) -> float:
    return float(np.real(np.vdot(psi, h @ psi)))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Predict k and r from an energy budget without running VQE."
    )
    ap.add_argument("--block", choices=["upper", "lower"], default="upper")
    ap.add_argument(
        "--target",
        type=float,
        default=0.25,
        help="scale max non-identity Pauli coefficient to this fraction (< 1/2)",
    )
    ap.add_argument(
        "--target-error-mh",
        type=float,
        default=0.1,
        help="maximum non-cancelling physical-energy error budget in milliHartree",
    )
    ap.add_argument("--k-min", type=int, default=3)
    ap.add_argument("--k-max", type=int, default=16)
    ap.add_argument("--r-min", type=int, default=2)
    ap.add_argument("--r-max", type=int, default=10)
    ap.add_argument(
        "--output-json",
        type=Path,
        default=Path("chem_vqe_precision.json"),
        help="where to save the recommendation and all predictor rows",
    )
    args = ap.parse_args()

    if not (0.0 < args.target < 0.5):
        ap.error("--target must lie strictly between 0 and 0.5")
    if args.target_error_mh <= 0:
        ap.error("--target-error-mh must be positive")
    if args.k_min < 0 or args.k_max < args.k_min:
        ap.error("invalid k range")
    if args.r_min < 1 or args.r_max < args.r_min:
        ap.error("invalid r range")

    _, identity, nonidentity, scale, physical_ground = exact_decomposition(
        args.block, args.target
    )

    rows = []

    for k in range(args.k_min, args.k_max + 1):
        hk, rounded = dyadic_hamiltonian(k, identity, nonidentity, scale)
        evals, evecs = np.linalg.eigh(hk)
        e0_scaled = float(evals[0].real)
        e0_mh = e0_scaled / scale
        ham_shift_mh = e0_mh - physical_ground
        ham_error_mh = abs(ham_shift_mh)

        ground = make_real_ground_vector(evecs[:, 0])
        alpha_star = state_to_continuous_alphas(ground)
        reconstructed = state_from_alphas(alpha_star)
        fidelity = float(abs(np.vdot(ground, reconstructed)) ** 2)
        if fidelity < 1.0 - 1e-9:
            raise RuntimeError(
                f"Continuous ansatz reconstruction failed at k={k}: fidelity={fidelity}"
            )

        for r in range(args.r_min, args.r_max + 1):
            q, alpha_rounded, grid_size = round_alphas_to_grid(alpha_star, r)
            psi = state_from_alphas(alpha_rounded)
            candidate_scaled = energy(psi, hk)
            candidate_mh = candidate_scaled / scale

            # Since E0 is the minimum eigenvalue this should be nonnegative;
            # clamp tiny numerical negative values to zero.
            parameter_gap_mh = max(0.0, candidate_mh - e0_mh)
            budget_mh = ham_error_mh + parameter_gap_mh
            signed_candidate_error_mh = candidate_mh - physical_ground

            rows.append(
                {
                    "k": k,
                    "r": r,
                    "grid_size": grid_size,
                    "later_exhaustive_candidate_count": grid_size**3,
                    "q": [int(x) for x in q],
                    "alpha_star": [float(x) for x in alpha_star],
                    "alpha_rounded": [float(x) for x in alpha_rounded],
                    "continuous_ansatz_fidelity": fidelity,
                    "ground_hk_scaled": e0_scaled,
                    "ground_hk_mh": e0_mh,
                    "hamilton_shift_mh": ham_shift_mh,
                    "hamilton_error_mh": ham_error_mh,
                    "constructed_candidate_scaled": candidate_scaled,
                    "constructed_candidate_mh": candidate_mh,
                    "parameter_gap_upper_bound_mh": parameter_gap_mh,
                    "predicted_noncancelling_budget_mh": budget_mh,
                    "constructed_signed_error_vs_physical_mh": signed_candidate_error_mh,
                    "rounded_pauli_coefficients": rounded,
                }
            )

    feasible = [
        row
        for row in rows
        if row["predicted_noncancelling_budget_mh"] <= args.target_error_mh
    ]

    recommendation = None
    if feasible:
        recommendation = min(
            feasible,
            key=lambda row: (
                row["later_exhaustive_candidate_count"],
                row["k"],
                row["r"],
            ),
        )

    # Also expose the numerically lowest predicted budget in case the requested
    # target is tighter than the scanned ranges can reach.
    lowest_budget = min(
        rows,
        key=lambda row: (
            row["predicted_noncancelling_budget_mh"],
            row["later_exhaustive_candidate_count"],
            row["k"],
            row["r"],
        ),
    )

    result = {
        "method": "classical constructive precision prediction; no VQE grid search",
        "scope": (
            "DMANH+ 4x4 benchmark; r prediction uses exact ground-state "
            "diagonalization and analytic inversion of the current real-state ansatz"
        ),
        "block": args.block,
        "scale": scale,
        "target_error_mh": args.target_error_mh,
        "physical_ground_mh": physical_ground,
        "k_range": [args.k_min, args.k_max],
        "r_range": [args.r_min, args.r_max],
        "recommendation": recommendation,
        "lowest_predicted_budget_in_scan": lowest_budget,
        "rows": rows,
    }

    args.output_json.write_text(json.dumps(result, indent=2) + "\n")

    print("DMANH+ pre-VQE precision prediction")
    print("  VQE search performed          = no")
    print("  target non-cancelling budget = %.12g mH" % args.target_error_mh)
    print("  physical exact ground        = %.12f mH" % physical_ground)
    print("  k scan                       = %d..%d" % (args.k_min, args.k_max))
    print("  r scan                       = %d..%d" % (args.r_min, args.r_max))

    if recommendation is None:
        print("\nNo (k,r) in the requested ranges is predicted to meet the target.")
        print("Best predicted budget in the scanned ranges:")
        rec = lowest_budget
    else:
        print("\nRecommended precisions before running VQE")
        rec = recommendation

    print(f"  k                            = {rec['k']}")
    print(f"  r                            = {rec['r']}")
    print(f"  predicted q near optimum     = {rec['q']}")
    print(
        "  later exhaustive candidates = {:,}".format(
            rec["later_exhaustive_candidate_count"]
        )
    )
    print(f"  Hamilton dyadic error        = {rec['hamilton_error_mh']:.12e} mH")
    print(
        "  parameter rounding bound     = "
        f"{rec['parameter_gap_upper_bound_mh']:.12e} mH"
    )
    print(
        "  non-cancelling budget        = "
        f"{rec['predicted_noncancelling_budget_mh']:.12e} mH"
    )
    print(
        "  constructed candidate error  = "
        f"{rec['constructed_signed_error_vs_physical_mh']:+.12e} mH"
    )
    print("  continuous alpha*            = [" + ", ".join(f"{x:.9f}" for x in rec["alpha_star"]) + "]")
    print("  rounded alpha                = [" + ", ".join(f"{x:.9f}" for x in rec["alpha_rounded"]) + "]")
    print(f"  details written to           = {args.output_json}")

    if recommendation is not None:
        print("\nRun the actual compiler-backed VQE only after this selection:")
        print(f"  cabal run chem-vqe -- {rec['k']} {rec['r']} 0")
        print("\nOptional validation experiment (not needed for selection):")
        print("  python -m scripts.experiments.chem_vqe_sweep --target-error-mh %.12g" % args.target_error_mh)


if __name__ == "__main__":
    main()
