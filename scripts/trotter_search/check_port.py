"""Check the Python search ports against the Haskell ChemVQE runs
(results/hamilton_vqe/08_search_method_comparison): k = 14, r = 7, n = 0."""
import itertools
import sys
import numpy as np

from scripts.trotter_search.vh_numpy import Ansatz, chem_target, energies
from scripts.trotter_search import search as S

h, _ = chem_target(14)
seq = Ansatz([[("YI", 0, 1.0)], [("IY", 1, 0.5), ("ZY", 1, 0.5)], [("IY", 2, 0.5), ("ZY", 2, -0.5)]], 3, 2)


def efn(q, r):
    return float(energies(h, seq.states(np.array(q) / 2 ** r, 0))[0])


legacy1 = lambda m: [0, m // 4, m // 2]
legacy2 = lambda m: [0, m // 4, m // 2, 3 * m // 4, 7 * m // 8]
specs = [S.FitSpec([2.0], legacy1), S.FitSpec([1.0, 2.0], legacy2), S.FitSpec([1.0, 2.0], legacy2)]
expect = {"refine": (315, (28, 30, 19)), "coordinate": (189, (28, 30, 19)),
          "analytic": (128, (28, 30, 19)), "analytic_direct": (38, (28, 30, 19))}
ok = True
for name, run in [("refine", lambda o: S.refine(o, 3, 7)), ("coordinate", lambda o: S.coordinate(o, 3, 7)),
                  ("analytic", lambda o: S.analytic(o, 3, 7, specs)),
                  ("analytic_direct", lambda o: S.analytic_direct(o, 3, 7, specs))]:
    o = S.Oracle(efn)
    e, q = run(o)
    good = (o.queries, tuple(q)) == expect[name]
    ok &= good
    print(f"{name:16s} queries={o.queries:4d} q={tuple(q)} E={e:.12f}  haskell={expect[name]}  {'OK' if good else 'MISMATCH'}")
sys.exit(0 if ok else 1)
