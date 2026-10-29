#!/usr/bin/env python3
"""Check two families of fixed-clock identities.

* Z_k := delta_k (+) (-delta_k) has unit-time evolution e^{-2 pi i Z / 2^k},
  equal up to global phase to I, -I, Z, S, T for k = 0..4.
* ZZ_3 has unit-time evolution e^{-i pi ZZ/4} = e^{-i pi/4} (S (x) S) CZ.
"""
import numpy as np
from scipy.linalg import expm

Z = np.diag([1, -1]).astype(complex)
S = np.diag([1, 1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])
I = np.eye(2)


def equal_up_to_phase(a, b):
    k = np.argmax(np.abs(b))
    ph = a.flat[k] / b.flat[k]
    return abs(abs(ph) - 1) < 1e-12 and np.allclose(a, ph * b)


for k, g, name in [(0, I, "I"), (1, -I, "-I"), (2, Z, "Z"), (3, S, "S"), (4, T, "T")]:
    U = expm(-2j * np.pi * Z / 2 ** k)
    print(f"Z_{k}: unit-time evolution = {name} up to phase: {equal_up_to_phase(U, g)}"
          + (f"; exactly: {np.allclose(U, g)}" if k < 2 else ""))

ZZ = np.kron(Z, Z)
CZ = np.diag([1, 1, 1, -1])
U = expm(-1j * np.pi / 4 * ZZ)
print("ZZ_3: e^{-i pi ZZ/4} == e^{-i pi/4} (S x S) CZ:",
      np.allclose(U, np.exp(-1j * np.pi / 4) * np.kron(S, S) @ CZ))
