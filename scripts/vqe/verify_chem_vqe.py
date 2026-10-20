#!/usr/bin/env python3
"""Fast matrix-semantics reference for app/ChemVQE.hs.

This script does not compile QASM.  It independently evaluates the same dyadic
Hamiltonian and the same three-parameter discrete ansatz directly as matrices,
so it is useful for checking the Haskell/compiler run.

Usage:
    python -m scripts.vqe.verify_chem_vqe [k] [r]

Defaults: k=10, r=3.
"""

import itertools
import json
import subprocess
import sys

import numpy as np
from scipy.linalg import expm

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
P = {"I": I, "X": X, "Y": Y, "Z": Z}


def pmat(label):
    out = np.array([[1]], dtype=complex)
    for c in label:
        out = np.kron(out, P[c])
    return out


def fetch(k):
    p = subprocess.run(
        [sys.executable, "-m", "scripts.chemistry.chem_decompose", "--block", "upper", "--level", str(k)],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(p.stdout)


def build_h(d):
    h = (d["identity_num"] / 2 ** d["identity_k"]) * pmat("II")
    for t in d["terms"]:
        h += (t["num"] / 2 ** t["k"]) * pmat(t["pauli"])
    return h


def main():
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    r = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    if r < 1:
        raise SystemExit("r must be at least 1")

    d = fetch(k)
    h = build_h(d)
    e0 = float(np.min(np.linalg.eigvalsh(h)).real)

    g1 = pmat("YI")
    g2 = 0.5 * (pmat("IY") + pmat("ZY"))
    g3 = 0.5 * (pmat("IY") - pmat("ZY"))
    generators = [g1, g2, g3]

    grid = range(2 ** (r - 1))
    psi0 = np.array([1, 0, 0, 0], dtype=complex)

    best = None
    for qs in itertools.product(grid, repeat=3):
        psi = psi0.copy()
        for q, g in zip(qs, generators):
            alpha = q / 2 ** r
            psi = expm(-2j * np.pi * alpha * g) @ psi
        e = float(np.real(np.vdot(psi, h @ psi)))
        if best is None or e < best[0]:
            best = (e, qs)

    e, qs = best
    scale = d["scale"]
    print(f"k={k}, r={r}, candidates={len(grid)**3}")
    print(f"best q                  = {list(qs)}")
    print(f"full VQE energy scaled  = {e:.12f}")
    print(f"ground H_k scaled       = {e0:.12f}")
    print(f"variational gap scaled  = {e-e0:.12e}")
    print(f"VQE energy (mH)         = {e/scale:.12f}")
    print(f"ground H_k (mH)         = {e0/scale:.12f}")
    print(f"physical exact (mH)     = {d['physical_ground_mh']:.12f}")


if __name__ == "__main__":
    main()
