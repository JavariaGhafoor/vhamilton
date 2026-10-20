#!/usr/bin/env python3
"""Independent Pauli-label semantic regression helper.

Reads JSON on stdin:
  {"proc": <Hamilton Proc>, "pauli": "XI", "coeff": 0.125}

Computes the Hamiltonian denotation of the process using scripts.semantics.denot and compares
it with coeff * (P_first kron P_second) in the standard |00>,|01>,|10>,|11>
basis.  This explicitly checks label/wire ordering such as XI vs IX and XY vs
YX; it is intentionally independent of the QASM compiler.
"""

import json
import sys

import numpy as np

from scripts.semantics.denot import hamiltonian

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I, "X": X, "Y": Y, "Z": Z}


def pauli_matrix(label: str) -> np.ndarray:
    if not label:
        raise ValueError("empty Pauli label")
    out = PAULI[label[0]]
    for ch in label[1:]:
        out = np.kron(out, PAULI[ch])
    return out


def main() -> None:
    payload = json.loads(sys.stdin.read())
    label = payload["pauli"]
    coeff = float(payload["coeff"])

    actual = hamiltonian(payload["proc"])
    expected = coeff * pauli_matrix(label)
    err = float(np.linalg.norm(actual - expected, ord="fro"))

    sys.stdout.write(json.dumps({"frob": err}))


if __name__ == "__main__":
    main()
