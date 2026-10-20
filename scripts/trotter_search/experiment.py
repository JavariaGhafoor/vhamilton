"""Search methods of Sec. 5 on a Trotterized ansatz, and search in the Trotter
dimension.

Usage: python3 -m scripts.trotter_search.experiment [controlled|pair]

Example "pair" (results/trotter_pair): one generator per clock,

    [ α1 YI ∥ α2 ZY ]  at level n     (YI and ZY anticommute)
    [ α3 IY ]          at level 0     (IY commutes with both: exact, run once)

which spans the same three Paulis as the default example below.

Default example "controlled" (results/trotter_search):

Ansatz: the three generators of the Sec. 5 ansatz written as ONE parallel
VHamilton block (see app/TrotterChemEnergies.hs),

    G(α) = α1 YI ∥ α2 (IY + ZY)/2 ∥ α3 (IY − ZY)/2,

on the DMANH+ target at k = 14, searched at r = 7 (2^18 valuations).  The
program denotes e^{-2πi G(α)}; level n compiles it as 2^n repetitions of its
schedule.  At n = 0 that is exactly the sequential Sec. 5 circuit.

Writes results/trotter_search/{landscapes,methods,certify,continuation}.csv
and prints a log.
"""
import csv
import itertools
import math
import os
import sys
import numpy as np

from scripts.trotter_search.vh_numpy import Ansatz, chem_target, energies
from scripts.trotter_search import search as S

K, R, P = 14, 7, 3
NMAX = 5


def parallel_ansatz():
    return Ansatz([[("YI", 0, 1.0), ("IY", 1, 0.5), ("ZY", 1, 0.5),
                    ("IY", 2, 0.5), ("ZY", 2, -0.5)]], 3, 2)


def pair_ansatz():
    return Ansatz([[("YI", 0, 1.0), ("ZY", 1, 1.0)], [("IY", 2, 1.0)]], 3, 2,
                  trotterized=[True, False])


LEG1 = lambda m: [0, m // 4, m // 2]
LEG2 = lambda m: [0, m // 4, m // 2, 3 * m // 4, 7 * m // 8]
EXAMPLES = {
    # name: (ansatz, per-clock generator eigenvalues, per-clock "repeated 2^n times", output dir)
    "controlled": (parallel_ansatz, [[-1, 1], [-1, 0, 1], [-1, 0, 1]], [True, True, True],
                   "results/trotter_search"),
    "pair": (pair_ansatz, [[-1, 1], [-1, 1], [-1, 1]], [True, True, False],
             "results/trotter_pair"),
}
NAME = sys.argv[1] if __name__ == "__main__" and len(sys.argv) > 1 else "controlled"
ANSATZ, EIGS, REPEATED, OUT = EXAMPLES[NAME]


def naive_spec(e):
    return S.FitSpec([2.0], LEG1) if len(e) == 2 else S.FitSpec([1.0, 2.0], LEG2)


def specs_naive():
    """The n = 0 budget of ChemVQE, used blindly at every level."""
    return [naive_spec(e) for e in EIGS]


def specs_aware(n):
    """Budget read off the program and the level: a clock whose block is
    compiled at level n has 2^n times as many frequencies; a clock in an
    exact (level-0) block keeps its n = 0 budget."""
    if n == 0:
        return specs_naive()
    return [S.FitSpec(S.clock_frequencies(e, n)) if rep else naive_spec(e)
            for e, rep in zip(EIGS, REPEATED)]


def main():
    os.makedirs(OUT, exist_ok=True)
    h, d = chem_target(K)
    scale = d["scale"]                      # scaled energy = scale * mHa
    mh = lambda x: x / scale
    ans = ANSATZ()
    m = 2 ** (R - 1)
    Q = np.array(list(itertools.product(range(m), repeat=P)))
    A = Q / 2 ** R
    exact = energies(h, ans.exact_states(A))
    emin = exact.min()
    qx = tuple(int(x) for x in Q[exact.argmin()])
    ground = np.linalg.eigvalsh(h)[0]
    idx = lambda q: (q[0] * m + q[1]) * m + q[2]
    print(f"target k={K}: ground {ground:.9f};  program optimum on the r={R} grid "
          f"{emin:.9f} at {qx}  (gap to ground {mh(emin - ground):.4f} mHa)")

    tables = {}
    for n in range(NMAX + 2):
        for sch in ("canon", "sym"):
            tables[(n, sch)] = energies(h, ans.states(A, n, sch))

    def oracle(n, sch):
        t = tables[(n, sch)]

        def f(q, r):
            lift = 2 ** (R - r)
            return float(t[idx(tuple(x * lift for x in q))])
        return S.Oracle(f, cost_per_query=ans.atoms_per_query(n))

    true = lambda q: float(exact[idx(q)])

    # ---------------------------------------------------------------- landscapes
    # window of +-8 grid steps per clock around the program's optimum
    near = np.all(np.abs(Q - np.array(qx)) <= 8, axis=1)
    rows = []
    for n in range(NMAX + 1):
        for sch in ("canon", "sym"):
            t = tables[(n, sch)]
            i = int(t.argmin())
            q = tuple(int(x) for x in Q[i])
            rows.append(dict(n=n, schedule=sch, atoms_per_query=ans.atoms_per_query(n),
                             compiled_min=t[i], q=q, true_energy=exact[i],
                             regret_mHa=mh(exact[i] - emin),
                             below_program_mHa=mh(min(0.0, t[i] - emin)),
                             mean_err_mHa=mh(np.abs(t - exact).mean()),
                             max_err_mHa=mh(np.abs(t - exact).max()),
                             near_err_mHa=mh(np.abs(t - exact)[near].mean())))
    write(f"{OUT}/landscapes.csv", rows)
    print("\nLandscapes (exhaustive over 262,144 valuations)")
    for r_ in rows:
        print(f"  n={r_['n']} {r_['schedule']:5s} compiled min {r_['compiled_min']:.6f} at {r_['q']}"
              f"  true {r_['true_energy']:.6f}  regret {r_['regret_mHa']:.3f} mHa"
              f"  below program {r_['below_program_mHa']:.3f} mHa"
              f"  mean|E_n-E| {r_['mean_err_mHa']:.3f} mHa  near optimum {r_['near_err_mHa']:.4f} mHa")

    # ---------------------------------------------------------------- methods per level
    rows = []
    print("\nSearch methods at each level (queries / atom-applications / q / true regret)")
    for n in range(NMAX + 1):
        for sch in ("canon", "sym"):
            runs = [("exhaustive", None), ("refine", lambda o: S.refine(o, P, R)),
                    ("coordinate", lambda o: S.coordinate(o, P, R)),
                    ("analytic_naive", lambda o: S.analytic(o, P, R, specs_naive(), log)),
                    ("direct_naive", lambda o: S.analytic_direct(o, P, R, specs_naive(), log=log)),
                    ("analytic_aware", lambda o: S.analytic(o, P, R, specs_aware(n), log)),
                    ("direct_aware", lambda o: S.analytic_direct(o, P, R, specs_aware(n), log=log))]
            for name, run in runs:
                log = []
                o = oracle(n, sch)
                if name == "exhaustive":
                    t = tables[(n, sch)]
                    i = int(t.argmin())
                    e, q, nq = float(t[i]), tuple(int(x) for x in Q[i]), len(t)
                else:
                    e, q = run(o)
                    nq = o.queries
                res = [abs(x[2]) for x in log if x[0] == "residual"]
                conds = [x[2] for x in log if x[0] == "cond"]
                infeasible = sorted({x[1] for x in log if x[0] == "infeasible"})
                rows.append(dict(n=n, schedule=sch, method=name, queries=nq,
                                 atom_applications=nq * ans.atoms_per_query(n),
                                 q=q, compiled_energy=e, true_energy=true(q),
                                 regret_mHa=mh(true(q) - emin),
                                 max_fit_residual=max(res) if res else "",
                                 max_fit_cond=max(conds) if conds else "",
                                 infeasible_clocks=" ".join(map(str, infeasible))))
                rr = rows[-1]
                extra = f"  fit residual {rr['max_fit_residual']:.1e}" if res else ""
                extra += f"  cond {rr['max_fit_cond']:.1e}" if conds else ""
                extra += f"  UNFITTABLE clocks {infeasible}" if infeasible else ""
                print(f"  n={n} {sch:5s} {name:15s} {nq:7d} q {str(q):13s} {nq * ans.atoms_per_query(n):9d} atoms"
                      f"  regret {rr['regret_mHa']:.3f} mHa{extra}")
    write(f"{OUT}/methods.csv", rows)

    # ---------------------------------------------------------------- certify with one query at n+1
    rows = []
    print("\nCertificate: re-evaluate the level-n optimum once at level n+1")
    for n in range(NMAX + 1):
        for sch in ("canon", "sym"):
            t, t1 = tables[(n, sch)], tables[(n + 1, sch)]
            i = int(t.argmin())
            q = tuple(int(x) for x in Q[i])
            jump = mh(t1[i] - t[i])
            rows.append(dict(n=n, schedule=sch, q=q, E_n=t[i], E_n1=t1[i], jump_mHa=jump,
                             true_regret_mHa=mh(exact[i] - emin)))
            print(f"  n={n} {sch:5s} q {q}  E_n {t[i]:.6f}  E_(n+1) {t1[i]:.6f}  jump {jump:+.3f} mHa"
                  f"  true regret {mh(exact[i] - emin):.3f} mHa")
    write(f"{OUT}/certify.csv", rows)

    # ---------------------------------------------------------------- climb the Trotter levels
    rows = []
    print("\nContinuation in n: analytic-direct at n=0, then warm-started local search at n=1,2,...")
    for sch in ("canon", "sym"):
        for local in ("axial", "aware"):
            o0 = oracle(0, sch)
            e, q = S.analytic_direct(o0, P, R, specs_naive())
            total_q, total_cost = o0.queries, o0.cost
            rows.append(dict(schedule=sch, local=local, n=0, queries=o0.queries, cum_queries=total_q,
                             cum_atoms=total_cost, q=q, compiled_energy=e, true_energy=true(q),
                             regret_mHa=mh(true(q) - emin), step_mHa=""))
            prev_q = q
            for n in range(1, NMAX + 1):
                o = oracle(n, sch)
                e0 = o(prev_q, R)
                if local == "axial":
                    e, q = S.coordinate_level(o, R, [(e0, prev_q)])[0]
                else:
                    e, q = S.analytic_direct(o, P, R, specs_aware(n), start=prev_q,
                                             cache={prev_q: e0})
                total_q += o.queries
                total_cost += o.cost
                step = mh(e0 - tables[(n - 1, sch)][idx(prev_q)])
                rows.append(dict(schedule=sch, local=local, n=n, queries=o.queries, cum_queries=total_q,
                                 cum_atoms=total_cost, q=q, compiled_energy=e, true_energy=true(q),
                                 regret_mHa=mh(true(q) - emin), step_mHa=step))
                prev_q = q
            for r_ in rows[-(NMAX + 1):]:
                print(f"  {sch:5s} {local:5s} n={r_['n']} +{r_['queries']:3d} queries (total {r_['cum_queries']:4d},"
                      f" {r_['cum_atoms']:6d} atoms)  q {str(r_['q']):13s} regret {r_['regret_mHa']:.3f} mHa"
                      + (f"  E_n(prev) - E_(n-1)(prev) = {r_['step_mHa']:+.3f} mHa" if r_['n'] else ""))
    write(f"{OUT}/continuation.csv", rows)


def write(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (" ".join(map(str, v)) if isinstance(v, tuple) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
