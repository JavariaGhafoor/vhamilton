#!/usr/bin/env python3
"""
runq.py — load an OpenQASM 3 file and print its unitary matrix.

Usage: python3 -m scripts.compiler.runq <file.qasm>
"""

import sys
import numpy as np
import qiskit.qasm3 as qasm3
from qiskit.quantum_info import Operator

EPS = 5e-3
COL = 12  # column width, matches Haskell formatting

def fmt(z: complex) -> str:
    re = 0.0 if abs(z.real) < EPS else float(z.real)
    im = 0.0 if abs(z.imag) < EPS else float(z.imag)
    if re == 0.0 and im == 0.0:
        s = "0"
    elif im == 0.0:
        s = f"{re:.2f}"
    elif re == 0.0:
        s = f"{im:.2f}i"
    else:
        sign = "+" if im > 0 else ""
        s = f"{re:.2f}{sign}{im:.2f}i"
    return f"{s:>{COL}}"

def print_matrix(U: np.ndarray, header: str) -> None:
    print(header)
    for row in U:
        print("[ " + " ".join(fmt(z) for z in row) + " ]")
    print()

def main(argv):
    if len(argv) != 2:
        print("usage: python3 -m scripts.compiler.runq <file.qasm>")
        return 2
    with open(argv[1], encoding="utf-8") as f:
        qc = qasm3.loads(f.read())
    if qc.num_qubits == 0:
        print("(no qubits)")
        return 0
    qc_no_meas = qc.remove_final_measurements(inplace=False)
    U = Operator(qc_no_meas).data
    print_matrix(U, "=== U_circ ===")
    print("=== circuit ===")
    print(qc.draw(output="text"))
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
