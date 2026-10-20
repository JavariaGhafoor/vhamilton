"""No Qiskit/Haskell needed: check the *new* VHamilton JSON denotation.

Run from the project root:
    python3 -m unittest discover -s test -p 'test_*.py' -v
"""

import unittest
import numpy as np
from scipy.linalg import expm
from scripts.semantics.denot import hamiltonian, hamiltonian0


def bare(tag, contents=None):
    return {"tag": tag, **({"contents": contents} if contents is not None else {})}


def diagonal(actions):
    """Four-dimensional diagonal with the exact Haskell Oplus nesting."""
    return bare("Oplus", [bare("Oplus", actions[:2]), bare("Oplus", actions[2:])])


def pauli_diagonal(label, k, alpha, sign=1):
    result = []
    for a, b in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        parity = sum(bit for bit, char in zip((a, b), label) if char != "I") % 2
        plus = sign * (1 if parity == 0 else -1) > 0
        result.append(bare("DeltaParam" if plus else "NegDeltaParam", [k, alpha]))
    return diagonal(result)


class MinimalSyntaxTests(unittest.TestCase):
    def test_hamilton_original_bare_actions(self):
        for tag, c, expected in [
            ("Id", None, 0), ("Delta", 3, 0.125),
            ("NegDelta", 3, -0.125),
            ("DeltaPow", [3, -7], -7/8),
        ]:
            with self.subTest(tag=tag):
                np.testing.assert_allclose(hamiltonian0(bare(tag, c)), [[expected]])

    def test_only_two_parameterized_atoms_and_sign(self):
        rho = {"alpha": 3/16}
        p = bare("DeltaParam", [2, "alpha"])
        m = bare("NegDeltaParam", [2, "alpha"])
        np.testing.assert_allclose(hamiltonian0(p, rho), [[3/64]])
        np.testing.assert_allclose(hamiltonian0(m, rho), [[-3/64]])
        with self.assertRaisesRegex(ValueError, 'unbound'):
            hamiltonian0(p)
        # Generalized scalar expressions are explicitly not accepted.
        with self.assertRaisesRegex(ValueError, 'unknown bare action'):
            hamiltonian0(bare('DeltaX', [2, {"tag": "Const", "contents": 0.5}]))

    def test_parameterized_y_generator(self):
        # S·H rotates Z into +Y; coefficient selected by the valuation.
        rho = {"alpha": 5/16}
        z_core = bare("Oplus", [bare("DeltaParam", [0, "alpha"]),
                                  bare("NegDeltaParam", [0, "alpha"])])
        cliff = {"tag": "Compose", "contents": [{"tag": "S"}, {"tag": "H"}]}
        term = {"tag": "Base", "contents": {"tag": "Conj", "contents": [
            cliff, {"tag": "Base", "contents": z_core}]}}
        y = np.array([[0, -1j], [1j, 0]])
        np.testing.assert_allclose(hamiltonian(term, rho), (5/16)*y, atol=1e-14)

    def test_dmanh_three_generators(self):
        # First isolate IY and ZY diagonal cores. S·H on the low bit rotates Z -> Y.
        rho = {"alpha2": 1/4, "alpha3": 3/16}
        iy = pauli_diagonal('IY', 1, 'alpha2')
        zy = pauli_diagonal('ZY', 1, 'alpha2')
        minus_zy = pauli_diagonal('ZY', 1, 'alpha3', sign=-1)
        # The sign pattern corresponds to conditional Y generators after
        # the existing Clifford basis conjugations.
        np.testing.assert_allclose(np.diag(hamiltonian0(iy, rho)),
                                   [1/8, -1/8, 1/8, -1/8])
        np.testing.assert_allclose(np.diag(hamiltonian0(zy, rho)),
                                   [1/8, -1/8, -1/8, 1/8])
        np.testing.assert_allclose(np.diag(hamiltonian0(minus_zy, rho)),
                                   [-3/32, 3/32, 3/32, -3/32])

    def test_one_qubit_noiseless_vqe_grid(self):
        """Check the unchanged calibration optimum using *new* action nodes."""
        zcore = bare("Oplus", [bare("Delta", 0), bare("NegDelta", 0)])
        z = {"tag": "Base", "contents": {"tag": "Base", "contents": zcore}}
        x = {"tag": "Base", "contents": {"tag": "Conj", "contents": [
            {"tag": "H"}, {"tag": "Base", "contents": zcore}]}}
        target = {"tag": "Par", "contents": [x, z]}
        H = hamiltonian(target)
        ycore = bare("Oplus", [bare("DeltaParam", [0, "alpha"]),
                                  bare("NegDeltaParam", [0, "alpha"])])
        y = {"tag": "Base", "contents": {"tag": "Conj", "contents": [
            {"tag": "Compose", "contents": [{"tag": "S"}, {"tag": "H"}]},
            {"tag": "Base", "contents": ycore}]}}
        psi0 = np.array([1.0, 0.0], dtype=complex)
        energies = []
        for q in range(8):
            G = hamiltonian(y, {"alpha": q/16})
            psi = expm(-2j*np.pi*G) @ psi0
            energies.append(float(np.real(np.vdot(psi, H @ psi))))
        self.assertEqual(int(np.argmin(energies)), 5)
        self.assertAlmostEqual(energies[5], -np.sqrt(2), places=12)
        self.assertAlmostEqual(float(np.min(np.linalg.eigvalsh(H))), -np.sqrt(2))

    def test_halving_preserves_parameters_and_sign(self):
        # A refinement n adds n to k; repeating 2^n segments restores the angle.
        rho = {"a": 3/8}
        for tag in ("DeltaParam", "NegDeltaParam"):
            for n in range(5):
                original = hamiltonian0(bare(tag, [1, "a"]), rho)[0, 0]
                segmented = hamiltonian0(bare(tag, [1+n, "a"]), rho)[0, 0]
                np.testing.assert_allclose(2**n*segmented, original)


if __name__ == '__main__':
    unittest.main()
