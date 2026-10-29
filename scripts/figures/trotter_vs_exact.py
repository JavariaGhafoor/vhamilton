#!/usr/bin/env python3
"""Compare trotter-vqe-demo output with the program's exact trajectories.

Target H = X + Y + Z, ansatz e^{-2 pi i (alpha Z + beta X)} |0>, valuations q / 2^r.
For each compiled grid (one Trotter level n) report the compiled optimum, the
exact energy of that same valuation, the gap between compiled and exact
landscapes, and how many compiled energies fall below the exact family's
minimum (which only compilation error can produce).

Usage: python3 trotter_vs_exact.py results/trotter_vqe/trotter_grid_r4_n*.csv
"""
import csv
import re
import sys

import numpy as np
from scipy.linalg import expm

X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]])
Z = np.diag([1, -1]).astype(complex)
H = X + Y + Z
psi0 = np.array([1, 0], complex)


def exact(a, b):
    s = expm(-2j * np.pi * (a * Z + b * X)) @ psi0
    return float(np.real(s.conj() @ H @ s))


grid = np.linspace(0, 1, 401)
bound = min(exact(a, b) for a in grid for b in grid)
print(f"exact family minimum (continuous scan): {bound:.6f}; ground energy: {np.linalg.eigvalsh(H)[0]:.6f}")
for path in sys.argv[1:]:
    m = re.search(r"_r(\d+)_n(\d+)\.csv$", path)
    r, n = int(m.group(1)), int(m.group(2))
    rows = list(csv.DictReader(open(path)))
    N = 2 ** r
    best = min(rows, key=lambda x: float(x["energy"]))
    qa, qb = int(best["q_alpha"]), int(best["q_beta"])
    gaps = [abs(float(x["energy"]) - exact(int(x["q_alpha"]) / N, int(x["q_beta"]) / N)) for x in rows]
    below = sum(float(x["energy"]) < bound - 1e-9 for x in rows)
    print(f"n={n}: compiled best {float(best['energy']):+.6f} at q=({qa},{qb}); exact there {exact(qa / N, qb / N):+.6f}; "
          f"mean |compiled-exact| {np.mean(gaps):.3f}, max {max(gaps):.3f}; below exact minimum: {below}")
