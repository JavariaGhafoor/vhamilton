"""For every setting of the global validation where step 5 (analytic-direct)
misses the exhaustive optimum, check by direct enumeration that its stopping
point is a coordinatewise minimum: for each clock j and every grid value q,
E(nu*) <= E(nu* with clock j set to q).  No fitted curve is involved.

Usage:
  python -m scripts.experiments.check_coordinate_minima \
    --matrix results/hamilton_vqe/11_analytic_direct_global_validation/direct_validation_matrix.csv
"""
import argparse
import csv

from scripts.experiments.validate_analytic_direct import SCALE, dyadic_hamiltonian, exact_energy


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    args = ap.parse_args()

    misses = [x for x in csv.DictReader(open(args.matrix)) if x["q_match"] != "True"]
    ok_all = True
    for x in misses:
        k, r = int(x["k"]), int(x["r"])
        q = [int(x["direct_q1"]), int(x["direct_q2"]), int(x["direct_q3"])]
        ref = (int(x["ref_q1"]), int(x["ref_q2"]), int(x["ref_q3"]))
        hk = dyadic_hamiltonian(k)
        e = exact_energy(hk, q, r)
        m = 2 ** (r - 1)
        # smallest improvement available by changing one clock (negative = improvement)
        best_change = min(exact_energy(hk, q[:j] + [v] + q[j + 1:], r) - e
                          for j in range(3) for v in range(m) if v != q[j])
        ok = best_change >= -1e-12
        ok_all &= ok
        print(f"k={k:2d} r={r} stop={tuple(q)} optimum={ref} "
              f"excess={(e - exact_energy(hk, ref, r)) / SCALE:.6f} mHa "
              f"best single-clock change={best_change / SCALE:+.6f} mHa  "
              f"{'coordinatewise minimum' if ok else 'NOT a coordinatewise minimum'}")
    print(f"{len(misses)} misses; all coordinatewise minima: {ok_all}")


if __name__ == "__main__":
    main()
