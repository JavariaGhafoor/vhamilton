#!/usr/bin/env python3
"""Run the numerical experiments appearing in the Hamilton paper, in paper order.

Every result used in the paper is written below results/hamilton/<experiment>/.
Run from the project root:

    python -m scripts.experiments.run_hamilton
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from scripts.experiments._experiment_utils import ensure_dir, project_root, run_logged

ROOT = project_root()
DEFAULT_RESULTS = ROOT / "results" / "hamilton"

EXPERIMENTS = [
    "01_zx_trotter_convergence",
    "02_zx_schedule_orderings",
    "03_zxy_schedule_orderings",
    "04_zz_commuting_exact",
    "05_chemistry_dyadic_error",
    "06_chemistry_trotter_convergence",
]


def save_decomp(k: int, out: Path) -> None:
    cmd = [sys.executable, "-m", "scripts.chemistry.chem_decompose", "--block", "upper", "--level", str(k)]
    text = run_logged(cmd, out, log_name=f"chem_decompose_k{k}.log", command_name=f"COMMAND_decomp_k{k}.txt")
    # stdout is JSON; keep a directly consumable copy too.
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return
    (out / f"chem_decomposition_k{k}.json").write_text(json.dumps(data, indent=2) + "\n")


def run_one(name: str, base: Path) -> None:
    out = ensure_dir(base / name)
    if name == "01_zx_trotter_convergence":
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "zx-conv"], out,
                   log_name="zx_trotter.log", command_name="COMMAND_table.txt")
        run_logged(["cabal", "run", "gen-circuits", "--", str(out), "zx"], out,
                   log_name="zx_circuits.log", command_name="COMMAND_circuits.txt")
    elif name == "02_zx_schedule_orderings":
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "zx-orderings"], out)
    elif name == "03_zxy_schedule_orderings":
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "zxy-orderings"], out)
    elif name == "04_zz_commuting_exact":
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "zz-parity"], out,
                   log_name="zz_exact.log", command_name="COMMAND_table.txt")
        run_logged(["cabal", "run", "gen-circuits", "--", str(out), "zz"], out,
                   log_name="zz_circuit.log", command_name="COMMAND_circuit.txt")
    elif name == "05_chemistry_dyadic_error":
        # These are exactly the precisions used in the Hamilton-paper table; k=20 is its reference.
        for k in (3, 4, 6, 8, 10, 12, 20):
            save_decomp(k, out)
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "dyadic-error"], out)
    elif name == "06_chemistry_trotter_convergence":
        save_decomp(12, out)
        run_logged(["cabal", "run", "gen-tables", "--", str(out), "chem-trotter"], out)
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
        print(f"\n=== Hamilton experiment: {name} ===")
        run_one(name, args.results_dir)
    print(f"\nHamilton-paper results saved under: {args.results_dir}")


if __name__ == "__main__":
    main()
