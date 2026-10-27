"""Dyadic coefficients (VHamilton) versus floating-point angles (Qiskit).

Part A: the chemistry target of Sec. 5 (DMANH+ upper block, scaled as in
chem_decompose.py), with every Pauli coefficient rounded to k bits.  For each k
we report the measured coefficient, ground-energy and unit-time unitary errors
next to the bounds that follow from k alone:

    |c - q/2^k| <= 2^-(k+1)               per coefficient
    |E0(H) - E0(H_k)| <= sum |dc|         (Weyl)
    ||e^{-2 pi i H} - e^{-2 pi i H_k}|| <= 2 pi sum_{non-I} |dc|

Part B: how Qiskit (2.x) treats the same kind of angles: radians are doubles,
the transpiler snaps angles to Clifford+T within a tolerance, and other angles
are approximated by Solovay-Kitaev with an error the program does not state.

Run:  python3 -m scripts.dyadic_vs_float.compare > results/dyadic_vs_float/run.log
"""
from math import pi

import numpy as np
from scipy.linalg import expm

from scripts.chemistry.chem_decompose import HTILDE, pauli_decompose, pauli_matrix


def chemistry(target=0.25):
    blk = HTILDE[:4, :4]
    blk = (blk + blk.T) / 2.0
    coeffs = pauli_decompose(blk)
    nonid = {l: c for l, c in coeffs.items() if l != "II"}
    scale = target / max(abs(c) for c in nonid.values())
    return {l: scale * c for l, c in coeffs.items()}, scale


def build(coeffs):
    return sum(c * pauli_matrix(l) for l, c in coeffs.items())


def part_a():
    coeffs, scale = chemistry()
    H = build(coeffs)
    E0 = np.linalg.eigvalsh(H)[0]
    U = expm(-2j * pi * build({l: c for l, c in coeffs.items() if l != "II"}))
    L = len(coeffs)
    print("Part A: chemistry target, %d Pauli terms (incl. II), 1 mHa = %.6f scaled" % (L, scale))
    print("%3s %10s %10s %12s %12s %12s %12s %8s" % (
        "k", "max|dc|", "2^-(k+1)", "dE0 (mHa)", "Weyl (mHa)", "||dU||", "2pi*sum", "bits"))
    for k in range(4, 21, 2):
        ck = {l: round(c * 2 ** k) / 2 ** k for l, c in coeffs.items()}
        dc = {l: abs(c - ck[l]) for l, c in coeffs.items()}
        Ek = np.linalg.eigvalsh(build(ck))[0]
        Uk = expm(-2j * pi * build({l: c for l, c in ck.items() if l != "II"}))
        dU = np.linalg.norm(U - Uk, 2)
        bnd_U = 2 * pi * sum(v for l, v in dc.items() if l != "II")
        bits = sum(int(abs(round(c * 2 ** k))).bit_length() + 1 for c in coeffs.values())
        print("%3d %10.2e %10.2e %12.2e %12.2e %12.2e %12.2e %8d" % (
            k, max(dc.values()), 2.0 ** -(k + 1), abs(E0 - Ek) / scale,
            sum(dc.values()) / scale, dU, bnd_U, bits))


def part_b():
    from fractions import Fraction
    from qiskit import QuantumCircuit, transpile
    from qiskit.quantum_info import Operator

    print("\nPart B: Qiskit angles")
    f = Fraction(pi / 4)
    print("pi/4 is stored as the double %d/%d (a dyadic number of radians; pi/4 is irrational)" % (
        f.numerator, f.denominator))
    basis = ["h", "s", "sdg", "t", "tdg", "x", "z", "cx"]

    def dist(a, b):  # operator-norm distance, minimised over global phase
        A, B = Operator(a).data, Operator(b).data
        ph = np.vdot(B.flatten(), A.flatten())
        ph = ph / abs(ph) if abs(ph) > 0 else 1
        return np.linalg.norm(A - ph * B, 2)

    print("%-24s %-28s %12s" % ("rz angle", "Clifford+T result", "error"))
    for name, th in [("pi/4  (T, k=4)", pi / 4),
                     ("pi/4 + 1e-9", pi / 4 + 1e-9),
                     ("pi/4 + 1e-6", pi / 4 + 1e-6),
                     ("pi/4 + 1e-3", pi / 4 + 1e-3),
                     ("3pi/8 (q=3, k=5)", 3 * pi / 8),
                     ("3pi/2^13 (q=3, k=14)", 3 * pi / 2 ** 13)]:
        qc = QuantumCircuit(1)
        qc.rz(th, 0)
        t = transpile(qc, basis_gates=basis, optimization_level=1)
        ops = dict(t.count_ops())
        print("%-24s %-28s %12.2e" % (name, ops, dist(qc, t)))


if __name__ == "__main__":
    part_a()
    part_b()
