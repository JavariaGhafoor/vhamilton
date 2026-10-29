#!/usr/bin/env python3
"""
denot.py — Hamiltonian and unitary semantics of Hamilton processes.

Reads JSON from stdin:
  {"proc": <Proc>, "qasm": "<OpenQASM 3 string>"}

Writes JSON to stdout:
  {"frob_aligned": float, "ideal": [[[re, im]]], "dim": int}

Implements the semantics:
  hamiltonian0  -- ⟦d⟧_H  for bare actions (diagonal matrices):
    ⟦id⟧_H       =  0
    ⟦δ_k^j⟧_H    = j/2^k
    ⟦±δ_k^alpha⟧_H,rho = ±rho(alpha)/2^k
    ⟦d ⊕ d'⟧_H   =  ⟦d⟧_H ⊕ ⟦d'⟧_H   (block diagonal)
  hamiltonian1  -- ⟦a⟧_H  for actions: ⟦u a u†⟧_H  =  U ⟦a⟧_H U†
  hamiltonian   -- ⟦p⟧_H  for processes: ⟦p ∥ p'⟧_H  =  ⟦p⟧_H + ⟦p'⟧_H
  clifford      -- ⟦u⟧_U  for Clifford transformations (unitary matrices)

  Unitary semantics: ⟦p⟧_U = exp(-2πi ⟦p⟧_H)  (computed via scipy.linalg.expm)
"""

import sys
import json
import numpy as np
from math import pi
from scipy.linalg import expm, block_diag

# ---------------------------------------------------------------------------
# Clifford matrices
# ---------------------------------------------------------------------------

CLIFFORD = {
    "I":    lambda _: np.eye(2, dtype=complex),
    "H":    lambda _: np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2),
    "S":    lambda _: np.array([[1,0],[0,1j]], dtype=complex),
    "CX":   lambda _: np.array([[1,0,0,0],[0,1,0,0],
                                 [0,0,0,1],[0,0,1,0]], dtype=complex),
    "SWAP": lambda _: np.array([[1,0,0,0],[0,0,1,0],
                                 [0,1,0,0],[0,0,0,1]], dtype=complex),
}

def clifford(ast):
    tag = ast["tag"]
    if tag in CLIFFORD:
        return CLIFFORD[tag](None)
    c = ast["contents"]
    if tag == "Tensor":
        # u ⊗ v: u on lower qubits → kron(v, u) in little-endian
        return np.kron(clifford(c[1]), clifford(c[0]))
    if tag == "Compose":
        # u · v  builds the matrix product clifford(u) @ clifford(v).
        # (Must match the compiler's composition order; the reversed order
        #  flips the sign of Y-basis conjugations, e.g. S·H, producing -Y
        #  instead of Y and a spurious Trotter-convergence floor.)
        return clifford(c[0]) @ clifford(c[1])
    raise ValueError(f"unknown Clifford tag: {tag!r}")

# ---------------------------------------------------------------------------
# Hamiltonian of bare actions, actions, and processes
# ---------------------------------------------------------------------------

def parameter_value(rho, name):
    if name not in rho:
        raise ValueError(f"unbound VHamilton parameter: {name!r}")
    return float(rho[name])


def hamiltonian0(ast, rho=None):
    """Hamiltonian denotation of a bare action under an external valuation."""
    if rho is None:
        rho = {}
    tag = ast["tag"]
    if tag == "Id":
        return np.array([[0.0]], dtype=complex)
    c = ast["contents"]
    if tag == "Delta":
        value = 1.0 / 2**int(c)
    elif tag == "NegDelta":
        value = -1.0 / 2**int(c)
    elif tag == "DeltaPow":
        k, j = int(c[0]), int(c[1])
        value = j / 2**k
    elif tag == "DeltaParam":
        k, name = int(c[0]), c[1]
        value = parameter_value(rho, name) / 2**k
    elif tag == "NegDeltaParam":
        k, name = int(c[0]), c[1]
        value = -parameter_value(rho, name) / 2**k
    elif tag == "Oplus":
        return block_diag(hamiltonian0(c[0], rho),
                          hamiltonian0(c[1], rho)).astype(complex)
    else:
        raise ValueError(f"unknown bare action tag: {tag!r}")
    return np.array([[value]], dtype=complex)


def hamiltonian1(ast, rho=None):
    """Hamiltonian denotation of an action."""
    tag = ast["tag"]
    c = ast["contents"]
    if tag == "Base":
        return hamiltonian0(c, rho)
    if tag == "Conj":
        U = clifford(c[0])
        H = hamiltonian1(c[1], rho)
        return U @ H @ U.conj().T
    raise ValueError(f"unknown action tag: {tag!r}")


def hamiltonian(ast, rho=None):
    """Hamiltonian denotation of a process under an external valuation."""
    tag = ast["tag"]
    c = ast["contents"]
    if tag == "Base":
        return hamiltonian1(c, rho)
    if tag == "Par":
        return hamiltonian(c[0], rho) + hamiltonian(c[1], rho)
    raise ValueError(f"unknown process tag: {tag!r}")

# ---------------------------------------------------------------------------
# Frobenius distance with global phase alignment
# ---------------------------------------------------------------------------

def frob(A):
    return float(np.linalg.norm(A, ord="fro"))

def align_phase(Uc, Ui):
    """Return e^{-iφ} Uc minimising ‖e^{-iφ} Uc - Ui‖_F."""
    inner = np.trace(Ui.conj().T @ Uc)
    return Uc if abs(inner) < 1e-14 else Uc / (inner / abs(inner))

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    from qiskit.qasm3 import loads
    from qiskit.quantum_info import Operator

    payload = json.loads(sys.stdin.read())
    H  = hamiltonian(payload["proc"], payload.get("valuation", {}))
    Ui = expm(-2j * pi * H)
    Uc = Operator(loads(payload["qasm"])).data
    Uc_aligned = align_phase(Uc, Ui)
    sys.stdout.write(json.dumps({
        "frob_aligned": frob(Uc_aligned - Ui),
        "ideal": [[[float(z.real), float(z.imag)] for z in row]
                  for row in Ui],
        "dim": int(Ui.shape[0]),
    }))

if __name__ == "__main__":
    main()
