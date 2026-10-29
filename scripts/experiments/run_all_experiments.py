#!/usr/bin/env python3
"""Run the Hamilton experiments first, then the Hamilton-VQE experiments."""
from scripts.experiments.run_hamilton import main as run_hamilton
from scripts.experiments.run_hamilton_vqe import main as run_vqe

if __name__ == "__main__":
    # These mains parse argv, so this convenience script is intentionally for the no-argument full run.
    import sys
    if len(sys.argv) != 1:
        raise SystemExit("run_all_experiments takes no arguments; use the two individual runners for --only")
    run_hamilton()
    run_vqe()
