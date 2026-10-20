#!/usr/bin/env python3
"""Run the Hamilton-VQE experiment suite in its intended reproducibility order.

All experiment outputs needed for reproducibility are persisted below
results/hamilton_vqe/<experiment>/ rather than existing only in terminal output.

Run everything:
    python -m scripts.experiments.run_hamilton_vqe

Run one stage:
    python -m scripts.experiments.run_hamilton_vqe --only 08_search_method_comparison
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scripts.chemistry.chem_decompose import HTILDE, pauli_decompose
from scripts.experiments import validate_analytic_direct as vad
from scripts.experiments._experiment_utils import (
    ensure_dir,
    parse_chem_vqe,
    project_root,
    run_logged,
    write_rows_csv,
)

ROOT = project_root()
DEFAULT_RESULTS = ROOT / "results" / "hamilton_vqe"

EXPERIMENTS = [
    "00_vhamilton_minimal_examples",
    "01_one_qubit_calibration",
    "02_compiler_regression",
    "03_chemistry_target_and_identity",
    "04_current_ansatz_trotter_invariance",
    "05_pre_vqe_precision_selection",
    "06_k_r_error_separation",
    "07_exhaustive_k_r_sweep",
    "08_search_method_comparison",
    "09_harmonic_fit_conditioning",
    "10_analytic_direct_representative",
    "11_analytic_direct_global_validation",
    "12_analytic_direct_precision_scaling",
]


def chem_run(k: int, r: int, n: int, mode: str, out: Path, stem: str) -> dict:
    text = run_logged(
        ["cabal", "run", "chem-vqe", "--", str(k), str(r), str(n), mode],
        out,
        log_name=f"{stem}.log",
        command_name=f"COMMAND_{stem}.txt",
    )
    row = {"k": k, "r": r, "n": n, "mode": mode, **parse_chem_vqe(text)}
    return row


def write_physical_pauli_table(out: Path) -> None:
    h = (HTILDE[:4, :4] + HTILDE[:4, :4].T) / 2.0
    coeffs = pauli_decompose(h)
    rows = [{"pauli": lab, "coefficient_mH": coeffs[lab]} for lab in sorted(coeffs)]
    write_rows_csv(out / "physical_pauli_decomposition.csv", rows)
    meta = {
        "physical_ground_mH": float(np.linalg.eigvalsh(h)[0]),
        "identity_coefficient_mH": float(coeffs.get("II", 0.0)),
    }
    (out / "chemistry_target_metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    run_logged(
        [sys.executable, "-m", "scripts.chemistry.chem_decompose", "--block", "upper", "--level", "14"],
        out,
        log_name="dyadic_decomposition_k14.json",
        command_name="COMMAND_dyadic_decomposition_k14.txt",
    )


def plot_search_comparison(rows: list[dict], out: Path) -> None:
    names = [r["method"] for r in rows]
    q = [r["queries"] for r in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.bar(names, q)
    ax.set_yscale("log")
    ax.set_ylabel("VQE energy queries (log scale)")
    ax.set_title("Search comparison at k=14, r=7")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(out / "search_query_comparison.png", dpi=180)
    plt.close(fig)


def validation_plots(csv_path: Path, out: Path) -> None:
    rows = list(csv.DictReader(csv_path.open()))
    ks = sorted({int(r["k"]) for r in rows})
    rs = sorted({int(r["r"]) for r in rows})
    mat = np.zeros((len(rs), len(ks)))
    scale = float(vad.SCALE)
    for row in rows:
        i = rs.index(int(row["r"]))
        j = ks.index(int(row["k"]))
        mat[i, j] = float(row["energy_diff_abs_scaled"]) / scale
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    im = ax.imshow(mat, aspect="auto")
    ax.set_xticks(range(len(ks)), ks)
    ax.set_yticks(range(len(rs)), rs)
    ax.set_xlabel("Hamiltonian precision k")
    ax.set_ylabel("parameter precision r")
    ax.set_title("Analytic-direct excess above exhaustive optimum (mH)")
    for i in range(len(rs)):
        for j in range(len(ks)):
            val = mat[i, j]
            label = "0" if val < 1e-12 else f"{val:.3g}"
            ax.text(j, i, label, ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="energy excess (mH)")
    fig.tight_layout()
    fig.savefig(out / "direct_validation_heatmap.png", dpi=180)
    plt.close(fig)

    summary = []
    for r in rs:
        rr = [x for x in rows if int(x["r"]) == r]
        exact = sum(str(x["q_match"]).lower() == "true" for x in rr)
        summary.append({
            "r": r,
            "exact_recoveries": f"{exact}/{len(rr)}",
            "mean_queries": float(np.mean([int(x["queries"]) for x in rr])),
            "mean_sweeps": float(np.mean([int(x["sweeps"]) for x in rr])),
            "max_excess_mH": max(float(x["energy_diff_abs_scaled"]) / scale for x in rr),
        })
    write_rows_csv(out / "validation_summary_table.csv", summary)


def scaling_plots(csv_path: Path, out: Path) -> None:
    rows = list(csv.DictReader(csv_path.open()))
    ks = sorted({int(r["k"]) for r in rows})
    for field, ylabel, filename in [
        ("queries", "objective-query count Q", "direct_queries_vs_r.png"),
        ("sweeps", "coordinate sweeps L", "direct_sweeps_vs_r.png"),
    ]:
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        for k in ks:
            rr = sorted((x for x in rows if int(x["k"]) == k), key=lambda x: int(x["r"]))
            ax.plot([int(x["r"]) for x in rr], [int(x[field]) for x in rr], marker="o", label=f"k={k}")
        ax.set_xlabel("parameter precision r")
        ax.set_ylabel(ylabel)
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / filename, dpi=180)
        plt.close(fig)


def run_one(name: str, base: Path) -> None:
    out = ensure_dir(base / name)

    if name == "00_vhamilton_minimal_examples":
        run_logged(["cabal", "run", "vhamilton-demo"], out, log_name="vhamilton_demo.log")

    elif name == "01_one_qubit_calibration":
        run_logged(["cabal", "run", "vqe-demo"], out, log_name="vqe_demo.log")

    elif name == "02_compiler_regression":
        run_logged(["cabal", "test", "regression", "--test-show-details=direct"], out,
                   log_name="regression.log")

    elif name == "03_chemistry_target_and_identity":
        write_physical_pauli_table(out)

    elif name == "04_current_ansatz_trotter_invariance":
        rows = [chem_run(10, 3, n, "exhaustive", out, f"k10_r3_n{n}") for n in (0, 1, 2)]
        write_rows_csv(out / "trotter_invariance_summary.csv", rows)

    elif name == "05_pre_vqe_precision_selection":
        targets = [1.0, 0.1, 0.05, 0.01, 0.005, 0.001]
        rows = []
        for target in targets:
            tag = str(target).replace(".", "p")
            json_path = out / f"precision_target_{tag}_mH.json"
            run_logged(
                [sys.executable, "-m", "scripts.chemistry.chem_decompose_vqe",
                 "--target-error-mh", str(target), "--output-json", str(json_path)],
                out,
                log_name=f"precision_target_{tag}_mH.log",
                command_name=f"COMMAND_precision_target_{tag}_mH.txt",
            )
            data = json.loads(json_path.read_text())
            rec = data["recommendation"] or data["lowest_predicted_budget_in_scan"]
            rows.append({
                "target_mH": target,
                "k": rec["k"], "r": rec["r"],
                "hamilton_error_mH": rec["hamilton_error_mh"],
                "parameter_rounding_bound_mH": rec["parameter_gap_upper_bound_mh"],
                "budget_mH": rec["predicted_noncancelling_budget_mh"],
                "reference_grid_size": rec["later_exhaustive_candidate_count"],
                "predicted_q": rec["q"],
            })
        write_rows_csv(out / "precision_selection_table.csv", rows)

    elif name == "06_k_r_error_separation":
        rows = [
            chem_run(10, 3, 0, "exhaustive", out, "k10_r3_exhaustive"),
            chem_run(10, 5, 0, "direct", out, "k10_r5_direct"),
        ]
        write_rows_csv(out / "k_r_separation_summary.csv", rows)

    elif name == "07_exhaustive_k_r_sweep":
        run_logged(
            [sys.executable, "-m", "scripts.experiments.chem_vqe_sweep",
             "--k-min", "3", "--k-max", "12", "--r-min", "2", "--r-max", "6",
             "--target-error-mh", "0.1", "--output-dir", str(out)],
            out,
            log_name="sweep.log",
        )

    elif name == "08_search_method_comparison":
        ref_text = run_logged([sys.executable, "-m", "scripts.vqe.verify_chem_vqe", "14", "7"], out,
                              log_name="exhaustive_reference.log", command_name="COMMAND_exhaustive_reference.txt")
        m = re.search(r"best q\s*=\s*(\[[^\n]+\])", ref_text)
        ref_q = json.loads(m.group(1)) if m else None
        rows = [{"method": "exhaustive_reference", "queries": 262144, "best_q": ref_q,
                 "reference_grid": 262144, "reduction": 1.0}]
        for mode, label in [("refine", "full_dyadic_beam"), ("coordinate", "coordinate_beam"), ("direct", "analytic_direct")]:
            row = chem_run(14, 7, 0, mode, out, label)
            rows.append({"method": label, "queries": row["queries"], "best_q": row["best_q"],
                         "reference_grid": row["reference_grid"], "reduction": row["reduction"]})
        write_rows_csv(out / "search_method_comparison.csv", rows)
        plot_search_comparison(rows, out)

    elif name == "09_harmonic_fit_conditioning":
        data = {
            "cond_V_G1_2norm": float(np.linalg.cond(vad.V_G1)),
            "cond_V_G23_2norm": float(np.linalg.cond(vad.V_G23)),
            "G1_sample_count": int(vad.V_G1.shape[0]),
            "G23_sample_count": int(vad.V_G23.shape[0]),
        }
        (out / "harmonic_fit_conditioning.json").write_text(json.dumps(data, indent=2) + "\n")
        (out / "COMMAND.txt").write_text("computed directly from validate_analytic_direct.V_G1 and V_G23\n")
        print(json.dumps(data, indent=2))

    elif name == "10_analytic_direct_representative":
        rows = [
            chem_run(10, 5, 0, "direct", out, "k10_r5_direct"),
            chem_run(14, 7, 0, "direct", out, "k14_r7_direct"),
            chem_run(14, 10, 0, "direct", out, "k14_r10_direct"),
        ]
        write_rows_csv(out / "analytic_direct_representative.csv", rows)

    elif name == "11_analytic_direct_global_validation":
        run_logged(
            [sys.executable, "-m", "scripts.experiments.validate_analytic_direct",
             "--experiment", "validation", "--out-dir", str(out)],
            out,
            log_name="validation.log",
        )
        validation_plots(out / "direct_validation_matrix.csv", out)

    elif name == "12_analytic_direct_precision_scaling":
        run_logged(
            [sys.executable, "-m", "scripts.experiments.validate_analytic_direct",
             "--experiment", "scaling", "--out-dir", str(out)],
            out,
            log_name="scaling.log",
        )
        scaling_plots(out / "direct_scaling_results.csv", out)

    else:
        raise ValueError(name)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--only", choices=EXPERIMENTS, default=None)
    args = ap.parse_args()
    ensure_dir(args.results_dir)
    chosen = [args.only] if args.only else EXPERIMENTS
    for name in chosen:
        print(f"\n=== Hamilton-VQE experiment: {name} ===")
        run_one(name, args.results_dir)
    print(f"\nHamilton-VQE experiment results saved under: {args.results_dir}")


if __name__ == "__main__":
    main()
