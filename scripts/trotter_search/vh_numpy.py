"""Exact NumPy replica of VHamilton's compiled energies for small ansatze.

Mirrors Compiler.compileTrotterWith: a block p = a_1 ∥ ... ∥ a_m compiled at
level n is the canonical schedule [a_1/2^n, ..., a_m/2^n] repeated 2^n times,
each atom run for unit time, so the block unitary is

    U_n = (e^{-2πi H_m/2^n} ··· e^{-2πi H_1/2^n})^{2^n}.

Blocks run in sequence on |0...0>.  Every term here is c·α·P for a Pauli
string P (c a fixed dyadic factor, α an open clock), so each factor has the
closed form cos(θ) I - i sin(θ) P with θ = 2π c α / 2^n.  'schedule="sym"'
alternates the canonical and reversed schedule between repetitions (a
compiler choice that keeps the same atoms and the same count).

Validated against the Haskell pipeline (Compiler + Qiskit statevector) by
validate_against_haskell.py.
"""
import json
import subprocess
import numpy as np

PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.diag([1.0, -1.0]).astype(complex),
}


def pauli(label):
    m = np.array([[1.0 + 0j]])
    for ch in label:          # first character = high (slow) qubit, as in ExamplesChem
        m = np.kron(m, PAULI[ch])
    return m


def chem_target(k, block="upper"):
    out = subprocess.run(
        ["python3", "-m", "scripts.chemistry.chem_decompose", "--block", block, "--level", str(k)],
        check=True, capture_output=True, text=True).stdout
    d = json.loads(out)
    h = sum(t["num"] / 2 ** t["k"] * pauli(t["pauli"]) for t in d["terms"])
    h = h + d["identity_num"] / 2 ** d["identity_k"] * np.eye(4)
    return h, d


class Ansatz:
    """blocks: list of blocks; a block is a list of terms (label, clock, coeff)
    meaning coeff·α_clock·P_label, in schedule (∥) order."""

    def __init__(self, blocks, nclocks, nqubits, trotterized=None):
        self.blocks = blocks
        self.nclocks = nclocks
        self.dim = 2 ** nqubits
        self.mats = [[(pauli(l), j, c) for (l, j, c) in b] for b in blocks]
        # Which blocks are compiled at the requested level n; the others
        # (blocks whose atoms commute) are compiled exactly at level 0.
        self.trotterized = [True] * len(blocks) if trotterized is None else trotterized

    def atoms(self):
        return sum(len(b) for b in self.blocks)

    def levels(self, n):
        return [n if t else 0 for t in self.trotterized]

    def atoms_per_query(self, n):
        return sum(len(b) * 2 ** lv for b, lv in zip(self.blocks, self.levels(n)))

    def states(self, alphas, n, schedule="canon"):
        """alphas: (B, nclocks) clock times.  Returns (B, dim) compiled states."""
        alphas = np.atleast_2d(np.asarray(alphas, dtype=float))
        bsz = alphas.shape[0]
        psi = np.zeros((bsz, self.dim), dtype=complex)
        psi[:, 0] = 1.0
        eye = np.eye(self.dim, dtype=complex)
        for block, lv in zip(self.mats, self.levels(n)):
            reps = 2 ** lv
            facs = []
            for (p, j, c) in block:
                th = 2 * np.pi * c * alphas[:, j] / reps
                facs.append(np.cos(th)[:, None, None] * eye - 1j * np.sin(th)[:, None, None] * p)
            for rep in range(reps):
                order = facs if (schedule == "canon" or rep % 2 == 0) else facs[::-1]
                for f in order:
                    psi = np.einsum("bij,bj->bi", f, psi)
        return psi

    def exact_states(self, alphas):
        """The program's own trajectories: each block runs e^{-2πi Σ_j H_j}."""
        alphas = np.atleast_2d(np.asarray(alphas, dtype=float))
        bsz = alphas.shape[0]
        psi = np.zeros((bsz, self.dim), dtype=complex)
        psi[:, 0] = 1.0
        for block in self.mats:
            g = np.zeros((bsz, self.dim, self.dim), dtype=complex)
            for (p, j, c) in block:
                g += (c * alphas[:, j])[:, None, None] * p
            w, v = np.linalg.eigh(g)
            u = np.einsum("bij,bj,bkj->bik", v, np.exp(-2j * np.pi * w), v.conj())
            psi = np.einsum("bij,bj->bi", u, psi)
        return psi


def energies(h, psi):
    return np.real(np.einsum("bi,ij,bj->b", psi.conj(), h, psi))
