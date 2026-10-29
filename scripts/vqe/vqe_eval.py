#!/usr/bin/env python3
"""Exact VQE energy oracle for the Hamilton codebase.

Input on stdin:
  {"hamiltonian": <Hamilton Proc JSON>, "qasm": "<OpenQASM 3>"}

Output:
  {"energy": float, "ground_energy": float}

This is intentionally an exact statevector oracle.  It is the noiseless
reference against which later shot noise, gate noise, and Trotter effects can
be studied separately.
"""

import json
import sys
import numpy as np
from qiskit.qasm3 import loads
from qiskit.quantum_info import Statevector

# Reuse the Hamiltonian denotation already used by the Hamilton implementation.
from scripts.semantics.denot import hamiltonian


def main() -> None:
    payload = json.loads(sys.stdin.read())
    H = hamiltonian(payload["hamiltonian"])

    qc = loads(payload["qasm"])
    psi = Statevector.from_instruction(qc).data

    energy = float(np.real(np.vdot(psi, H @ psi)))
    ground = float(np.min(np.linalg.eigvalsh(H)).real)

    sys.stdout.write(json.dumps({
        "energy": energy,
        "ground_energy": ground,
    }))


if __name__ == "__main__":
    main()
