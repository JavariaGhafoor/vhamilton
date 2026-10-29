"""Error calculus for the DMANH+ chemistry VQE: measured errors against the
a priori bounds of the paper (Sec. 7.4).

  delta_H(k)   = E0(H_k) - E0(H_bar)           vs  Weyl bound  M 2^-(k+1)
  delta_r(k,r) = min_grid E_{H_k} - E0(H_k)     vs  first order 2 pi p ||H_k|| 2^-r
                                               and second order 2 pi^2 p^2 ||H_k|| 4^-r

delta_r is computed by exhaustive search over the one-period grid
q_j in {0, ..., 2^(r-1)-1}.  All energies are reported in mHa.

Writes a CSV per panel and a pgfplots figure to OUT_DIR.

Usage:
  python -m scripts.experiments.error_calculus_bounds \
    --out-dir results/hamilton_vqe/13_error_calculus_bounds
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from scripts.experiments.validate_analytic_direct import (
    PHYSICAL_COEFFS, SCALE, UPPER_BLOCK_MH, dyadic_hamiltonian)

P = 3                       # clocks in the ansatz
M = len(PHYSICAL_COEFFS)    # Pauli terms in the target, identity included


def grid_minimum(hk: np.ndarray, r: int) -> float:
    phi = 2.0 * np.pi * np.arange(2 ** (r - 1)) / 2 ** r
    b1, b2, b3 = np.meshgrid(phi, phi, phi, indexing="ij")
    psi = np.stack([np.cos(b1) * np.cos(b2), np.cos(b1) * np.sin(b2),
                    np.sin(b1) * np.cos(b3), np.sin(b1) * np.sin(b3)], axis=-1)
    return float(np.einsum("...i,ij,...j->...", psi, hk.real, psi).min())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--k-values", type=int, nargs="+", default=[6, 10, 14])
    ap.add_argument("--r-max", type=int, default=8)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    e0_phys = float(np.linalg.eigvalsh(UPPER_BLOCK_MH)[0])

    ham_rows = []
    for k in range(3, 15):
        e0 = float(np.linalg.eigvalsh(dyadic_hamiltonian(k))[0]) / SCALE
        ham_rows.append({"k": k, "delta_H_mh": e0 - e0_phys,
                         "weyl_bound_mh": M * 2.0 ** -(k + 1) / SCALE})
        assert abs(e0 - e0_phys) <= ham_rows[-1]["weyl_bound_mh"]

    var_rows = []
    for k in args.k_values:
        hk = dyadic_hamiltonian(k)
        evals = np.linalg.eigvalsh(hk)
        norm = float(max(abs(evals)))
        prev = None
        for r in range(2, args.r_max + 1):
            gap = (grid_minimum(hk, r) - float(evals[0])) / SCALE
            first = 2 * np.pi * P * norm * 2.0 ** -r / SCALE
            second = 2 * np.pi ** 2 * P ** 2 * norm * 4.0 ** -r / SCALE
            assert gap <= first and gap <= second
            assert prev is None or gap <= prev + 1e-12   # nested grids
            prev = gap
            var_rows.append({"k": k, "r": r, "delta_r_mh": gap,
                             "first_order_bound_mh": first,
                             "second_order_bound_mh": second})

    for name, rows in (("hamiltonian_rounding.csv", ham_rows),
                       ("variational_gap.csv", var_rows)):
        with open(args.out_dir / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    def coords(pairs):
        return " ".join(f"({x},{y:.6e})" for x, y in pairs)

    marks = {6: "square*", 10: "triangle*", 14: "*"}
    shades = {6: "black!35", 10: "black!60", 14: "violet"}
    kmax = max(args.k_values)
    lines = [
        r"\begin{tikzpicture}",
        r"\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.3cm},",
        r"  width=0.38\linewidth, height=5cm, ymode=log, grid=major, grid style={gray!20},",
        r"  tick label style={font=\scriptsize}, label style={font=\footnotesize},",
        r"  legend style={font=\scriptsize, draw=none, fill=white, fill opacity=0.8, text opacity=1}, legend cell align=left]",
        r"\nextgroupplot[xlabel={Hamiltonian precision $k$}, ylabel={mHa}, xtick={4,6,...,14},",
        r"  title={\footnotesize (a) rounding the target}, ymin=1e-6, ymax=1e2, legend pos=south west]",
        r"\addplot[thick, violet, mark=*, mark size=1.6pt] coordinates {"
        + coords((x["k"], abs(x["delta_H_mh"])) for x in ham_rows) + "};",
        r"\addlegendentry{$|\delta_H(k)|$}",
        r"\addplot[thick, dashed, black!60] coordinates {"
        + coords((x["k"], x["weyl_bound_mh"]) for x in ham_rows) + "};",
        r"\addlegendentry{bound $M2^{-(k+1)}/s$}",
        r"\nextgroupplot[xlabel={clock precision $r$}, xtick={2,...," + str(args.r_max) + "},",
        r"  title={\footnotesize (b) the grid of clock values}, ymin=1e-4, ymax=1e3,",
        r"  legend style={at={(1.04,0.5)}, anchor=west}]",
    ]
    for k in args.k_values:
        lines.append(rf"\addplot[thick, {shades[k]}, mark={marks[k]}, mark size=1.6pt] coordinates {{"
                     + coords((x["r"], x["delta_r_mh"]) for x in var_rows if x["k"] == k) + "};")
        lines.append(rf"\addlegendentry{{$\delta_r({k},r)$}}")
    for key, style, label in (("first_order_bound_mh", "dotted", "first-order bound"),
                              ("second_order_bound_mh", "dashed", "second-order bound")):
        lines.append(rf"\addplot[thick, {style}, black!60] coordinates {{"
                     + coords((x["r"], x[key]) for x in var_rows if x["k"] == kmax) + "};")
        lines.append(rf"\addlegendentry{{{label}}}")
    lines += [r"\end{groupplot}", r"\end{tikzpicture}", ""]
    (args.out_dir / "error_calculus.tex").write_text("\n".join(lines))

    for x in ham_rows:
        print(f"k={x['k']:2d} delta_H={x['delta_H_mh']:+.6f} bound={x['weyl_bound_mh']:.4f}")
    for x in var_rows:
        print(f"k={x['k']:2d} r={x['r']} delta_r={x['delta_r_mh']:.6g} "
              f"first={x['first_order_bound_mh']:.4g} second={x['second_order_bound_mh']:.4g}")


if __name__ == "__main__":
    main()
