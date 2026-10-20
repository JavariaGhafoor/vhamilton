#!/usr/bin/env python3
"""
diag_to_qasm3.py — lower diagonal unitary gates to elementary gates (§Implementation).

Part of the compilation pipeline (fig:circuits): decomposes QDiagonalAt gates
emitted by Qasm.hs into {h, s, x, cx, swap, p, cp} via Qiskit DiagonalGate.

Reads JSON from stdin:
  {"off": int, "n": int, "phis": [float]}

  off:  qubit offset (diagonal acts on qubits [off .. off+n-1])
  n:    number of qubits (0 means global phase)
  phis: phase angles — the diagonal entries are exp(i*phi)

Writes JSON to stdout:
  {"gates": [{"name": str, "params": [float], "qubits": [int]}]}
"""

import sys
import json
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import DiagonalGate

SUPPORTED = {"h", "s", "x", "cx", "swap", "p", "cp"}

def main():
    payload = json.loads(sys.stdin.read())
    off  = int(payload["off"])
    n    = int(payload["n"])
    phis = payload["phis"]

    if n == 0:
        if len(phis) != 1:
            raise SystemExit("0-qubit diagonal must have exactly one phase")
        print(json.dumps({
            "gates": [{"name": "gphase",
                       "params": [float(phis[0])],
                       "qubits": []}]
        }))
        return

    expected = 1 << n
    if len(phis) != expected:
        raise SystemExit(f"expected {expected} phases, got {len(phis)}")

    qc = QuantumCircuit(off + n)
    qc.append(DiagonalGate(np.exp(1j * np.array(phis, dtype=float))),
              list(range(off, off + n)))

    tqc = transpile(qc, basis_gates=list(SUPPORTED), optimization_level=0)

    gates = []
    if abs(tqc.global_phase) > 1e-12:
        gates.append({"name": "gphase",
                      "params": [float(tqc.global_phase)],
                      "qubits": []})

    for inst, qargs, _ in tqc.data:
        if inst.name not in SUPPORTED:
            raise SystemExit(f"unsupported gate after transpile: {inst.name}")
        gates.append({
            "name":   inst.name,
            "params": [float(p) for p in inst.params],
            "qubits": [tqc.find_bit(q).index for q in qargs],
        })

    print(json.dumps({"gates": gates}))

if __name__ == "__main__":
    main()
