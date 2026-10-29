#!/usr/bin/env python3
"""
chem_decompose.py — Pauli-decompose the DMANH+ Hamiltonian block.

Reads the published DMANH+ block, Pauli-decomposes it, scales the Hamiltonian,
and rounds every Pauli coefficient to a dyadic rational at level k.

For Hamiltonian-simulation experiments the identity term can be omitted from the
process because it contributes only a global phase.  For VQE, however, the
identity coefficient shifts every energy and therefore must be retained as
metadata.  The JSON emitted here contains both the non-identity Hamilton process
and the dyadic identity offset needed to reconstruct absolute energies.

Output on stdout:
  {
    "scale": float,
    "level": int,
    "identity_coeff_mh": float,
    "scaled_identity": float,
    "identity_num": int,
    "identity_k": int,
    "physical_ground_mh": float,
    "scaled_physical_ground": float,
    "terms": [{"pauli": "ZZ", "num": int, "k": int}, ...]
  }

The non-identity term (num / 2^k) * P is represented directly by Hamilton.
The full dyadic Hamiltonian used for absolute-energy reporting is

    H_k(full) = (identity_num / 2^k) I + sum_P (num_P / 2^k) P.

The original physical energy scale is recovered by dividing scaled energies by
"scale".
"""

import sys
import json
import argparse
from itertools import product
import numpy as np

I = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I, "X": X, "Y": Y, "Z": Z}

# Block-diagonalized Htilde (mH). Upper and lower 4x4 diagonal blocks.
HTILDE = np.array([
    [37.68,  -7.537,  0.7182, -0.0024, 0, 0, 0, 0],
    [-7.5367, 7.182, -7.752,  0.8407,  0, 0, 0, 0],
    [0.7182, -7.752,  0.2725, -6.909,  0, 0, 0, 0],
    [-0.0024, 0.8407, -6.909, -7.460,  0, 0, 0, 0],
    [0, 0, 0, 0,  7.4960, -8.0473, 0.2976, 0.5454],
    [0, 0, 0, 0, -8.0473, -0.2705, -7.204, 0.4200],
    [0, 0, 0, 0,  0.2976, -7.204,  6.884, -7.420],
    [0, 0, 0, 0,  0.5454,  0.4200, -7.420, 37.64],
], dtype=float)


def pauli_matrix(label):
    m = np.array([[1]], dtype=complex)
    for ch in label:
        m = np.kron(m, PAULI[ch])
    return m


def pauli_decompose(h, tol=1e-6):
    n = int(np.log2(h.shape[0]))
    out = {}
    for labels in product("IXYZ", repeat=n):
        lab = "".join(labels)
        c = np.trace(pauli_matrix(lab).conj().T @ h) / (2 ** n)
        if abs(c.real) > tol:
            out[lab] = float(c.real)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block", choices=["upper", "lower"], default="upper")
    ap.add_argument("--level", type=int, default=10,
                    help="dyadic precision level k")
    ap.add_argument("--target", type=float, default=0.25,
                    help="map max non-identity |coeff| to this fraction (< 1/2)")
    args = ap.parse_args()

    if args.level < 0:
        ap.error("--level must be non-negative")
    if not (0.0 < args.target < 0.5):
        ap.error("--target must lie strictly between 0 and 0.5")

    blk = HTILDE[:4, :4] if args.block == "upper" else HTILDE[4:, 4:]
    blk = (blk + blk.T) / 2.0  # symmetrize tiny published asymmetry

    coeffs = pauli_decompose(blk)
    identity_coeff = float(coeffs.get("II", 0.0))
    nonidentity = {lab: c for lab, c in coeffs.items() if lab != "II"}

    maxc = max(abs(v) for v in nonidentity.values())
    scale = args.target / maxc

    # Round the identity term too.  It is not put into `terms`, because doing so
    # would add only a global phase to compiled dynamics, but VQE adds this
    # dyadic value back to every reported absolute energy.
    identity_scaled = scale * identity_coeff
    identity_num = round(identity_scaled * 2 ** args.level)

    terms = []
    for label, coef in sorted(nonidentity.items(), key=lambda kv: -abs(kv[1])):
        num = round(scale * coef * 2 ** args.level)
        if num != 0:
            terms.append({"pauli": label, "num": int(num), "k": args.level})

    physical_ground = float(np.min(np.linalg.eigvalsh(blk)).real)

    json.dump({
        "scale": scale,
        "level": args.level,
        "identity_coeff_mh": identity_coeff,
        "scaled_identity": identity_scaled,
        "identity_num": int(identity_num),
        "identity_k": args.level,
        "physical_ground_mh": physical_ground,
        "scaled_physical_ground": scale * physical_ground,
        "terms": terms,
    }, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
