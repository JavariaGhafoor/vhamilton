"""Exhaustive minimum of the chemistry VQE energy over the full half-range grid
q_j in {0, ..., 2^(r-1)-1}, j = 1, 2, 3, evaluated slice by slice (one q1 at a
time, all (q2, q3) at once).  Reports the three lowest valuations and the margin
between the best and the runner-up.  Two further checks compare double
precision with higher precision: the three lowest valuations are re-evaluated in
NumPy long double (80-bit x87 extended precision on x86-64, 64-bit significand)
and with mpmath at 50 significant digits, and a random sample of valuations
(seed 0) is evaluated in both double and long double.  These are numerical
checks, not a formal error bound.

Usage:
  python -m scripts.experiments.exhaustive_grid_certificate --k 14 --r 10
"""
import argparse

import mpmath
import numpy as np

from scripts.experiments.validate_analytic_direct import SCALE, dyadic_hamiltonian


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=14)
    ap.add_argument("--r", type=int, default=10)
    ap.add_argument("--samples", type=int, default=200_000)
    args = ap.parse_args()
    k, r = args.k, args.r

    hk = dyadic_hamiltonian(k).real
    m = 2 ** (r - 1)
    phi = 2.0 * np.pi * np.arange(m) / 2 ** r
    c, s = np.cos(phi), np.sin(phi)

    candidates = []
    for q1 in range(m):
        psi = np.stack([np.broadcast_to(c[q1] * c[:, None], (m, m)),
                        np.broadcast_to(c[q1] * s[:, None], (m, m)),
                        np.broadcast_to(s[q1] * c[None, :], (m, m)),
                        np.broadcast_to(s[q1] * s[None, :], (m, m))], axis=-1)
        e = np.einsum("abi,ij,abj->ab", psi, hk, psi).ravel()
        for i in np.argpartition(e, 3)[:3]:
            q2, q3 = np.unravel_index(i, (m, m))
            candidates.append((float(e[i]), (q1, int(q2), int(q3))))
    candidates.sort()

    e0 = float(np.linalg.eigvalsh(hk)[0])
    print(f"k={k} r={r} grid={m ** 3:,} valuations")
    for e, q in candidates[:3]:
        print(f"  q={q}  E={e:.14f} (scaled)  E-E0(H_k)={(e - e0) / SCALE:.9f} mHa")
    print(f"  margin best -> runner-up = {(candidates[1][0] - candidates[0][0]) / SCALE:.3e} mHa"
          f" = {candidates[1][0] - candidates[0][0]:.3e} scaled")

    # Higher-precision re-evaluation of the three lowest valuations.  H_k has
    # dyadic entries, exact in double, so only the trigonometry and the
    # arithmetic change precision.
    mpmath.mp.dps = 50
    hmp = mpmath.matrix([[mpmath.mpf(float(x)) for x in row] for row in hk])

    def energy_ld(q):
        p = 2 * np.arccos(np.longdouble(-1)) * np.asarray(q, dtype=np.longdouble) / np.longdouble(2 ** r)
        c, s = np.cos(p), np.sin(p)
        psi = np.array([c[0] * c[1], c[0] * s[1], s[0] * c[2], s[0] * s[2]])
        return psi @ hk.astype(np.longdouble) @ psi

    def energy_mp(q):
        p = [2 * mpmath.pi * x / 2 ** r for x in q]
        c, s = [mpmath.cos(x) for x in p], [mpmath.sin(x) for x in p]
        psi = mpmath.matrix([c[0] * c[1], c[0] * s[1], s[0] * c[2], s[0] * s[2]])
        return (psi.T * hmp * psi)[0]

    print("  higher precision (scaled):")
    mp_vals = []
    for e, q in candidates[:3]:
        emp = energy_mp(q)
        mp_vals.append(emp)
        print(f"    q={q}  double={e:.17g}  long double={float(energy_ld(q)):.17g}  "
              f"mpmath50={mpmath.nstr(emp, 20)}  |double-mpmath|={float(abs(e - emp)):.1e}")
    print(f"  margin at 50 digits = {mpmath.nstr(mp_vals[1] - mp_vals[0], 12)} scaled")

    # Random sample: double against long double.
    rng = np.random.default_rng(0)
    qs = rng.integers(0, m, size=(args.samples, 3))

    def energies(qs, dt, two_pi):
        p = two_pi * qs.astype(dt) / dt(2 ** r)
        c, s = np.cos(p), np.sin(p)
        psi = np.stack([c[:, 0] * c[:, 1], c[:, 0] * s[:, 1], s[:, 0] * c[:, 2], s[:, 0] * s[:, 2]], -1)
        return np.einsum("bi,ij,bj->b", psi, hk.astype(dt), psi)

    d = np.abs(energies(qs, np.float64, 2 * np.pi).astype(np.longdouble)
               - energies(qs, np.longdouble, 2 * np.arccos(np.longdouble(-1))))
    print(f"  sample of {args.samples:,} valuations (seed 0): max |double - long double| = {float(d.max()):.1e} scaled")


if __name__ == "__main__":
    main()
