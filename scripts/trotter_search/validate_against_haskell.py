"""Compare vh_numpy energies with the Haskell pipeline (trotter-chem-energies:
VHamilton compiler -> OpenQASM -> Qiskit statevector) on random valuations.

Usage: python3 -m scripts.trotter_search.validate_against_haskell BIN [POINTS] [controlled|pair]
"""
import subprocess
import sys
import numpy as np

from scripts.trotter_search.vh_numpy import chem_target, energies
from scripts.trotter_search.experiment import parallel_ansatz, pair_ansatz

K, R = 14, 7
binary = sys.argv[1]
npts = int(sys.argv[2]) if len(sys.argv) > 2 else 8
name = sys.argv[3] if len(sys.argv) > 3 else "controlled"
h, _ = chem_target(K)
ans = parallel_ansatz() if name == "controlled" else pair_ansatz()
rng = np.random.default_rng(7)
worst = 0.0
for n in range(5):
    for sched in ("canon", "sym"):
        qs = rng.integers(0, 2 ** (R - 1), size=(npts, 3))
        qs[0] = (28, 22, 16)
        inp = "\n".join(" ".join(map(str, q)) for q in qs) + "\n"
        out = subprocess.run([binary, name, str(K), str(R), str(n), sched], input=inp,
                             capture_output=True, text=True, check=True).stdout.split("\n")
        hs = np.array([float(l.split()[-1]) for l in out if l.strip()])
        ns = energies(h, ans.states(qs / 2 ** R, n, sched))
        d = np.abs(hs - ns).max()
        worst = max(worst, d)
        print(f"n={n} {sched:5s} points={npts}  max |haskell - numpy| = {d:.2e}")
print(f"worst difference {worst:.2e}")
sys.exit(0 if worst < 1e-9 else 1)
