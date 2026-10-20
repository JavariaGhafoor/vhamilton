#!/usr/bin/env python3
"""compare_ideals.py -- compare two processes' ideal unitaries (§Application).

Protocol (stdin):  {"proc1": <Proc>, "proc2": <Proc>}
Output (stdout):   {"frob": float}

Used to measure dyadic approximation error (§Application, dyadic error table):
  ‖⟦p_k⟧_U - ⟦p_ref⟧_U‖_F  for varying dyadic precision levels k.

Computes ‖exp(-2πi H1) - exp(-2πi H2)‖_F with no phase alignment,
since the distance between two ideal unitaries is a genuine physical quantity.
The unitary semantics ⟦·⟧_U is from fig:sem (implemented in scripts.semantics.denot).
"""

import sys
import json
from scipy.linalg import expm
from math import pi

from scripts.semantics.denot import hamiltonian, frob

def main():
    data = json.loads(sys.stdin.read())
    H1 = hamiltonian(data["proc1"])
    H2 = hamiltonian(data["proc2"])
    U1 = expm(-2j * pi * H1)
    U2 = expm(-2j * pi * H2)
    sys.stdout.write(json.dumps({"frob": frob(U1 - U2)}))

if __name__ == "__main__":
    main()
