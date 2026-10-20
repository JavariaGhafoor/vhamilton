# Regression tests

Run these from the project root with the Python virtual environment activated:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
cabal test regression --test-show-details=direct
```

The suite freezes the corrected compiler semantics and the first discrete VQE example.

It checks:

- one-qubit `X/8`, `Y/8`, and `Z/8` processes against
  `exp(-2*pi*i [[p]]_H)` up to global phase;
- two-qubit `XI/8`, `IX/8`, `XY/8`, `YX/8`, and `ZZ/8` compiled circuits against
  the same Hamiltonian semantics;
- independently, that the labels `XI`, `IX`, `XY`, `YX`, and `ZZ` denote
  `P_first kron P_second` in the standard `|00>,|01>,|10>,|11>` basis;
- the VQE baseline `H = X + Z`, `U(q)=exp(-2*pi*i*(q/16)Y)`, `q=0..7`, which must return
  `q*=5` and `E*=-sqrt(2)`.

The Pauli compiler tests use coefficient `1/8`, rather than coefficient `1`, because
`exp(-2*pi*i P)=I` for a Pauli operator `P` with eigenvalues `+/-1`; coefficient `1`
would therefore be too weak to detect basis/conjugation errors.

The compiler fixes being protected are:

1. `u · v` emits `v` and then `u`, so chronological circuit execution has matrix `UV`.
2. `Conj u a` emits `u†`, then `a`, then `u`, so its matrix is `U A U†`.
