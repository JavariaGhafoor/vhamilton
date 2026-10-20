#!/usr/bin/env python3
"""Systematic noiseless discrete-VQE sweep for the DMANH+ upper block.

This is the fast matrix-semantics experiment corresponding to app/ChemVQE.hs.
It sweeps Hamilton dyadic precision k and variational dyadic precision r,
saves every result to CSV, identifies the lowest-error pair, optionally finds
the cheapest pair that meets a requested physical-energy tolerance, and writes
convergence plots.

The ansatz is the same as ChemVQE.hs:
    G1 = YI
    G2 = (IY + ZY)/2 = |0><0| x Y
    G3 = (IY - ZY)/2 = |1><1| x Y
    alpha_j = q_j / 2^r,
    q_j in {0, ..., 2^(r-1)-1}.

Examples:
    python -m scripts.experiments.chem_vqe_sweep
    python -m scripts.experiments.chem_vqe_sweep --k-min 3 --k-max 12 --r-min 2 --r-max 6
    python -m scripts.experiments.chem_vqe_sweep --target-error-mh 0.1
    python -m scripts.experiments.chem_vqe_sweep --output-dir sweep_results

The exhaustive-search cost grows as 2^(3(r-1)); r=6 has 32,768 candidates.
The script precomputes each variational state once per r and evaluates all k
values against those states vectorially, so it is much faster than launching
Qiskit once per candidate.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm

from scripts.chemistry.chem_decompose import HTILDE, PAULI, pauli_decompose


def pmat(label: str) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for c in label:
        out = np.kron(out, PAULI[c])
    return out


PAULI_MATS = {"".join([a, b]): pmat(a + b) for a in "IXYZ" for b in "IXYZ"}

G1 = PAULI_MATS["YI"]
G2 = 0.5 * (PAULI_MATS["IY"] + PAULI_MATS["ZY"])
G3 = 0.5 * (PAULI_MATS["IY"] - PAULI_MATS["ZY"])
GENERATORS = (G1, G2, G3)
PSI0 = np.array([1, 0, 0, 0], dtype=complex)


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
    physical_ground = float(np.min(np.linalg.eigvalsh(h)).real)
    return h, identity, nonidentity, scale, physical_ground


def dyadic_hamiltonian(k: int, identity: float, nonidentity: dict[str, float], scale: float):
    denom = 2**k
    identity_num = round(scale * identity * denom)
    h = (identity_num / denom) * PAULI_MATS["II"]
    rounded = {"II": identity_num / denom}
    for label, coeff in nonidentity.items():
        num = round(scale * coeff * denom)
        if num != 0:
            value = num / denom
            h = h + value * PAULI_MATS[label]
            rounded[label] = value
    return h, rounded


def states_for_r(r: int):
    if r < 1:
        raise ValueError("r must be at least 1")

    grid_size = 2 ** (r - 1)
    denom = 2**r

    # Only 3 * grid_size small 4x4 exponentials are needed.  Candidate states
    # are then built from these cached unitaries.
    unitaries = []
    for g in GENERATORS:
        unitaries.append([
            expm(-2j * np.pi * (q / denom) * g)
            for q in range(grid_size)
        ])

    candidate_count = grid_size**3
    states = np.empty((candidate_count, 4), dtype=complex)
    qs = np.empty((candidate_count, 3), dtype=np.int32)

    idx = 0
    for q1 in range(grid_size):
        psi1 = unitaries[0][q1] @ PSI0
        for q2 in range(grid_size):
            psi2 = unitaries[1][q2] @ psi1
            for q3 in range(grid_size):
                states[idx] = unitaries[2][q3] @ psi2
                qs[idx] = (q1, q2, q3)
                idx += 1

    return states, qs, grid_size


def evaluate_states(states: np.ndarray, h: np.ndarray):
    # E_b = <psi_b|H|psi_b> for every candidate b.
    energies = np.einsum("bi,ij,bj->b", states.conj(), h, states).real
    best_idx = int(np.argmin(energies))
    return float(energies[best_idx]), best_idx


def write_csv(rows, path: Path):
    fields = [
        "k", "r", "grid_size", "candidate_count", "best_q1", "best_q2", "best_q3",
        "scale", "vqe_scaled", "ground_hk_scaled", "physical_ground_mh",
        "vqe_mh", "ground_hk_mh", "hamilton_shift_mh", "hamilton_error_mh",
        "variational_gap_scaled", "variational_gap_mh",
        "noncancelling_error_budget_mh",
        "total_signed_error_mh", "total_abs_error_mh",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def plot_hamilton_error(rows, k_values, output_dir: Path):
    by_k = {}
    for row in rows:
        by_k.setdefault(row["k"], row)
    ys = [by_k[k]["hamilton_error_mh"] for k in k_values]

    fig, ax = plt.subplots()
    ax.plot(k_values, ys, marker="o")
    ax.set_xlabel("Hamiltonian dyadic precision k")
    ax.set_ylabel("|E0(H_k) - E0(H)| (mH)")
    ax.set_title("Hamiltonian discretization convergence")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_dir / "hamilton_error_vs_k.png", dpi=180)
    plt.close(fig)


def plot_variational_gap(rows, k_values, r_values, output_dir: Path):
    fig, ax = plt.subplots()
    for k in k_values:
        subset = sorted((row for row in rows if row["k"] == k), key=lambda x: x["r"])
        ax.plot(
            [row["r"] for row in subset],
            [row["variational_gap_mh"] for row in subset],
            marker="o",
            label=f"k={k}",
        )
    ax.set_xlabel("Variational dyadic precision r")
    ax.set_ylabel("E*_k,r - E0(H_k) (mH)")
    ax.set_title("Discrete variational convergence")
    ax.grid(True, alpha=0.3)
    if len(k_values) <= 12:
        ax.legend(fontsize="small", ncols=2)
    fig.tight_layout()
    fig.savefig(output_dir / "variational_gap_vs_r.png", dpi=180)
    plt.close(fig)


def plot_total_error(rows, k_values, r_values, output_dir: Path):
    matrix = np.full((len(r_values), len(k_values)), np.nan)
    lookup = {(row["k"], row["r"]): row for row in rows}
    for i, r in enumerate(r_values):
        for j, k in enumerate(k_values):
            matrix[i, j] = lookup[(k, r)]["total_abs_error_mh"]

    fig, ax = plt.subplots()
    im = ax.imshow(matrix, aspect="auto", origin="lower")
    ax.set_xticks(range(len(k_values)), labels=k_values)
    ax.set_yticks(range(len(r_values)), labels=r_values)
    ax.set_xlabel("Hamiltonian dyadic precision k")
    ax.set_ylabel("Variational dyadic precision r")
    ax.set_title("Absolute physical VQE error (mH)")
    fig.colorbar(im, ax=ax, label="|E*_k,r - E0(H)| (mH)")
    fig.tight_layout()
    fig.savefig(output_dir / "total_error_heatmap.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k-min", type=int, default=3)
    ap.add_argument("--k-max", type=int, default=12)
    ap.add_argument("--r-min", type=int, default=2)
    ap.add_argument("--r-max", type=int, default=6)
    ap.add_argument("--block", choices=["upper", "lower"], default="upper")
    ap.add_argument("--target", type=float, default=0.25,
                    help="same scaling target used by scripts.chemistry.chem_decompose")
    ap.add_argument("--target-error-mh", type=float, default=None,
                    help="optional absolute physical-energy tolerance used to choose a cheapest pair")
    ap.add_argument("--output-dir", default="chem_vqe_sweep_results")
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    if args.k_min < 0 or args.k_max < args.k_min:
        ap.error("require 0 <= k-min <= k-max")
    if args.r_min < 1 or args.r_max < args.r_min:
        ap.error("require 1 <= r-min <= r-max")
    if not (0.0 < args.target < 0.5):
        ap.error("--target must lie strictly between 0 and 0.5")
    if args.target_error_mh is not None and args.target_error_mh <= 0:
        ap.error("--target-error-mh must be positive")

    k_values = list(range(args.k_min, args.k_max + 1))
    r_values = list(range(args.r_min, args.r_max + 1))
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    _, identity, nonidentity, scale, physical_ground = exact_decomposition(args.block, args.target)

    # H_k and exact E0(H_k) depend only on k.
    hk_data = {}
    for k in k_values:
        h, rounded = dyadic_hamiltonian(k, identity, nonidentity, scale)
        e0_scaled = float(np.min(np.linalg.eigvalsh(h)).real)
        e0_mh = e0_scaled / scale
        hk_data[k] = {
            "h": h,
            "ground_scaled": e0_scaled,
            "ground_mh": e0_mh,
            "hamilton_shift_mh": e0_mh - physical_ground,
            "hamilton_error_mh": abs(e0_mh - physical_ground),
            "rounded": rounded,
        }

    rows = []
    print(f"DMANH+ {args.block} block discrete-VQE sweep")
    print(f"k = {args.k_min}..{args.k_max}, r = {args.r_min}..{args.r_max}")
    print(f"scale = {scale:.12f}, exact physical ground = {physical_ground:.12f} mH")
    print()

    # Candidate states are independent of k, so enumerate once for each r.
    for r in r_values:
        states, qs, grid_size = states_for_r(r)
        candidate_count = len(states)
        print(f"r={r}: {candidate_count:,} candidate states")

        for k in k_values:
            hkd = hk_data[k]
            best_scaled, idx = evaluate_states(states, hkd["h"])
            best_q = qs[idx].tolist()
            vqe_mh = best_scaled / scale
            var_gap_scaled = best_scaled - hkd["ground_scaled"]
            var_gap_mh = vqe_mh - hkd["ground_mh"]
            total_signed = vqe_mh - physical_ground

            row = {
                "k": k,
                "r": r,
                "grid_size": grid_size,
                "candidate_count": candidate_count,
                "best_q1": int(best_q[0]),
                "best_q2": int(best_q[1]),
                "best_q3": int(best_q[2]),
                "scale": scale,
                "vqe_scaled": best_scaled,
                "ground_hk_scaled": hkd["ground_scaled"],
                "physical_ground_mh": physical_ground,
                "vqe_mh": vqe_mh,
                "ground_hk_mh": hkd["ground_mh"],
                "hamilton_shift_mh": hkd["hamilton_shift_mh"],
                "hamilton_error_mh": hkd["hamilton_error_mh"],
                "variational_gap_scaled": var_gap_scaled,
                "variational_gap_mh": var_gap_mh,
                "noncancelling_error_budget_mh": hkd["hamilton_error_mh"] + max(var_gap_mh, 0.0),
                "total_signed_error_mh": total_signed,
                "total_abs_error_mh": abs(total_signed),
            }
            rows.append(row)

    csv_path = output_dir / "chem_vqe_sweep.csv"
    write_csv(rows, csv_path)

    best = min(rows, key=lambda row: (row["total_abs_error_mh"], row["candidate_count"], row["k"]))
    best_budget = min(rows, key=lambda row: (row["noncancelling_error_budget_mh"], row["candidate_count"], row["k"]))

    summary = {
        "block": args.block,
        "k_values": k_values,
        "r_values": r_values,
        "scale": scale,
        "physical_ground_mh": physical_ground,
        "lowest_total_error": best,
        "lowest_noncancelling_error_budget": best_budget,
        "target_error_mh": args.target_error_mh,
        "cheapest_meeting_total_target": None,
        "cheapest_meeting_budget_target": None,
        "cheapest_definition": "minimum exhaustive candidate evaluations; ties broken by lower k, then lower r",
    }

    print("\nLowest absolute physical-energy error in the sweep")
    print(f"  k,r                  = ({best['k']}, {best['r']})")
    print(f"  best q               = [{best['best_q1']},{best['best_q2']},{best['best_q3']}]")
    print(f"  candidates           = {best['candidate_count']:,}")
    print(f"  VQE energy           = {best['vqe_mh']:.12f} mH")
    print(f"  Hamilton error       = {best['hamilton_error_mh']:.12e} mH")
    print(f"  variational gap      = {best['variational_gap_mh']:.12e} mH")
    print(f"  total absolute error = {best['total_abs_error_mh']:.12e} mH")
    print("  note                  = total error can benefit from cancellation")

    print("\nLowest non-cancelling error budget in the sweep")
    print(f"  k,r                  = ({best_budget['k']}, {best_budget['r']})")
    print(f"  best q               = [{best_budget['best_q1']},{best_budget['best_q2']},{best_budget['best_q3']}]")
    print(f"  Hamilton error       = {best_budget['hamilton_error_mh']:.12e} mH")
    print(f"  variational gap      = {best_budget['variational_gap_mh']:.12e} mH")
    print(f"  error budget (sum)   = {best_budget['noncancelling_error_budget_mh']:.12e} mH")

    if args.target_error_mh is not None:
        feasible = [row for row in rows if row["total_abs_error_mh"] <= args.target_error_mh]
        if feasible:
            cheapest = min(feasible, key=lambda row: (row["candidate_count"], row["k"], row["r"]))
            summary["cheapest_meeting_total_target"] = cheapest
            print(f"\nCheapest pair meeting |total error| <= {args.target_error_mh:g} mH")
            print("  cost definition      = minimum exhaustive candidate evaluations; ties -> lower k")
            print(f"  k,r                  = ({cheapest['k']}, {cheapest['r']})")
            print(f"  best q               = [{cheapest['best_q1']},{cheapest['best_q2']},{cheapest['best_q3']}]")
            print(f"  candidates           = {cheapest['candidate_count']:,}")
            print(f"  total absolute error = {cheapest['total_abs_error_mh']:.12e} mH")
        else:
            print(f"\nNo pair in this sweep reaches |total error| <= {args.target_error_mh:g} mH")

        feasible_budget = [row for row in rows if row["noncancelling_error_budget_mh"] <= args.target_error_mh]
        if feasible_budget:
            cheapest_budget = min(feasible_budget, key=lambda row: (row["candidate_count"], row["k"], row["r"]))
            summary["cheapest_meeting_budget_target"] = cheapest_budget
            print(f"\nCheapest pair meeting non-cancelling error budget <= {args.target_error_mh:g} mH")
            print(f"  k,r                  = ({cheapest_budget['k']}, {cheapest_budget['r']})")
            print(f"  candidates           = {cheapest_budget['candidate_count']:,}")
            print(f"  Hamilton error       = {cheapest_budget['hamilton_error_mh']:.12e} mH")
            print(f"  variational gap      = {cheapest_budget['variational_gap_mh']:.12e} mH")
            print(f"  error budget (sum)   = {cheapest_budget['noncancelling_error_budget_mh']:.12e} mH")
        else:
            print(f"\nNo pair reaches non-cancelling error budget <= {args.target_error_mh:g} mH")

    with (output_dir / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)

    if not args.no_plots:
        plot_hamilton_error(rows, k_values, output_dir)
        plot_variational_gap(rows, k_values, r_values, output_dir)
        plot_total_error(rows, k_values, r_values, output_dir)

    print(f"\nSaved CSV:     {csv_path}")
    print(f"Saved summary: {output_dir / 'summary.json'}")
    if not args.no_plots:
        print(f"Saved plots:   {output_dir / 'hamilton_error_vs_k.png'}")
        print(f"               {output_dir / 'variational_gap_vs_r.png'}")
        print(f"               {output_dir / 'total_error_heatmap.png'}")


if __name__ == "__main__":
    main()
