#!/usr/bin/env python3
"""Independent matrix-semantics check for chem-vqe's analytic search mode.

This mirrors the search logic in app/ChemVQE.hs but evaluates the two-qubit
ansatz directly as 4x4 NumPy/SciPy matrices.  It is intended as a fast check of
both the final q tuple and the number of VQE objective evaluations.

Usage:
    python -m scripts.vqe.verify_chem_vqe_analytic [k] [r] [analytic|direct]

Examples:
    python -m scripts.vqe.verify_chem_vqe_analytic 14 7 analytic
    python -m scripts.vqe.verify_chem_vqe_analytic 14 7 direct
"""

from __future__ import annotations

import itertools
import json
import math
import subprocess
import sys

import numpy as np
from scipy.linalg import expm

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
P = {"I": I, "X": X, "Y": Y, "Z": Z}

BEAM_WIDTH = 3
COARSE_R = 3
MAX_SWEEPS = 12
SCAN_BINS = 256
GOLDEN_ITERS = 70


def pmat(label: str) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for c in label:
        out = np.kron(out, P[c])
    return out


G1 = pmat("YI")
G2 = 0.5 * (pmat("IY") + pmat("ZY"))
G3 = 0.5 * (pmat("IY") - pmat("ZY"))
GENERATORS = (G1, G2, G3)
PSI0 = np.array([1, 0, 0, 0], dtype=complex)


def fetch(k: int) -> dict:
    proc = subprocess.run(
        [sys.executable, "-m", "scripts.chemistry.chem_decompose", "--block", "upper", "--level", str(k)],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(proc.stdout)


def build_h(d: dict) -> np.ndarray:
    h = (d["identity_num"] / 2 ** d["identity_k"]) * pmat("II")
    for term in d["terms"]:
        h += (term["num"] / 2 ** term["k"]) * pmat(term["pauli"])
    return h


def exact_energy(h: np.ndarray, qs: tuple[int, int, int], r: int) -> float:
    psi = PSI0.copy()
    for q, g in zip(qs, GENERATORS):
        alpha = q / 2**r
        psi = expm(-2j * np.pi * alpha * g) @ psi
    return float(np.real(np.vdot(psi, h @ psi)))


def phase(q: int, m: int) -> float:
    return math.pi * q / m


def trig_value(coeff: np.ndarray, phi: float) -> float:
    c, a1, b1, a2, b2 = coeff
    return float(
        c
        + a1 * math.cos(phi)
        + b1 * math.sin(phi)
        + a2 * math.cos(2 * phi)
        + b2 * math.sin(2 * phi)
    )


def sample_qs(j: int, m: int) -> list[int]:
    if j == 0:
        return [0, m // 4, m // 2]
    return [0, m // 4, m // 2, 3 * m // 4, 7 * m // 8]


def fit_model(j: int, m: int, qs: list[int], energies: list[float]) -> np.ndarray:
    if j == 0:
        e0, e_quarter, e_half = energies
        c = 0.5 * (e0 + e_half)
        a2 = 0.5 * (e0 - e_half)
        b2 = e_quarter - c
        return np.array([c, 0.0, 0.0, a2, b2])

    design = []
    for q in qs:
        ph = phase(q, m)
        design.append([1.0, math.cos(ph), math.sin(ph), math.cos(2 * ph), math.sin(2 * ph)])
    return np.linalg.solve(np.asarray(design), np.asarray(energies))


def golden_min(coeff: np.ndarray, a: float, b: float) -> float:
    rho = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - rho * (b - a)
    d = a + rho * (b - a)
    fc = trig_value(coeff, c)
    fd = trig_value(coeff, d)
    for _ in range(GOLDEN_ITERS):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - rho * (b - a)
            fc = trig_value(coeff, c)
        else:
            a, c, fc = c, d, fd
            d = a + rho * (b - a)
            fd = trig_value(coeff, d)
    return 0.5 * (a + b)


def minimum_phases(coeff: np.ndarray) -> list[float]:
    xs = np.linspace(0.0, math.pi, SCAN_BINS + 1)
    ys = np.array([trig_value(coeff, x) for x in xs])
    out = [0.0, math.pi]
    for i in range(1, SCAN_BINS):
        if ys[i] <= ys[i - 1] and ys[i] <= ys[i + 1]:
            out.append(golden_min(coeff, float(xs[i - 1]), float(xs[i + 1])))
    return out


def best_dyadic_q(coeff: np.ndarray, m: int) -> int:
    candidates = {0, m - 1}
    for ph in minimum_phases(coeff):
        q_cont = m * ph / math.pi
        q_floor = math.floor(q_cont)
        q_ceil = math.ceil(q_cont)
        q_round = round(q_cont)
        for q in (q_floor, q_ceil, q_round - 1, q_round, q_round + 1):
            candidates.add(max(0, min(m - 1, int(q))))
    return min(candidates, key=lambda q: trig_value(coeff, phase(q, m)))


def analytic_search(h: np.ndarray, r: int):
    if r <= COARSE_R:
        m = 2 ** (r - 1)
        candidates = list(itertools.product(range(m), repeat=3))
        values = [(exact_energy(h, q, r), q) for q in candidates]
        return min(values), len(values), 0

    coarse_m = 2 ** (COARSE_R - 1)
    coarse = []
    for q in itertools.product(range(coarse_m), repeat=3):
        coarse.append((exact_energy(h, q, COARSE_R), q))
    coarse.sort()
    beam = coarse[:BEAM_WIDTH]
    eval_count = len(coarse)

    lift = 2 ** (r - COARSE_R)
    beam = [(e, tuple(lift * x for x in q)) for e, q in beam]
    cache = {q: e for e, q in beam}
    m = 2 ** (r - 1)

    def get(q: tuple[int, int, int]) -> float:
        nonlocal eval_count
        if q not in cache:
            cache[q] = exact_energy(h, q, r)
            eval_count += 1
        return cache[q]

    for sweep in range(1, MAX_SWEEPS + 1):
        old_qs = [q for _, q in beam]
        updated = []
        for _, seed in beam:
            q = list(seed)
            e = get(tuple(q))
            for j in range(3):
                sq = sample_qs(j, m)
                vals = []
                for x in sq:
                    qt = q.copy()
                    qt[j] = x
                    vals.append(get(tuple(qt)))
                coeff = fit_model(j, m, sq, vals)
                q[j] = best_dyadic_q(coeff, m)
                e = get(tuple(q))
            updated.append((e, tuple(q)))

        merged = {q: e for e, q in beam + updated}
        beam = sorted((e, q) for q, e in merged.items())[:BEAM_WIDTH]
        if [q for _, q in beam] == old_qs:
            return beam[0], eval_count, sweep

    return beam[0], eval_count, MAX_SWEEPS



def balanced_seed(r: int) -> tuple[int, int, int]:
    if r < 3:
        raise ValueError("exact balanced alpha=1/8 seed requires r >= 3")
    q = 2 ** (r - 3)
    return (q, q, q)


def analytic_direct_search(h: np.ndarray, r: int):
    if r < 3:
        m = 2 ** (r - 1)
        candidates = list(itertools.product(range(m), repeat=3))
        values = [(exact_energy(h, q, r), q) for q in candidates]
        return min(values), len(values), 0

    q0 = balanced_seed(r)
    cache: dict[tuple[int, int, int], float] = {}
    eval_count = 0
    m = 2 ** (r - 1)

    def get(q: tuple[int, int, int]) -> float:
        nonlocal eval_count
        if q not in cache:
            cache[q] = exact_energy(h, q, r)
            eval_count += 1
        return cache[q]

    q = q0
    e = get(q)
    for sweep in range(1, MAX_SWEEPS + 1):
        old_q = q
        q_work = list(q)
        e_work = e
        for j in range(3):
            sq = sample_qs(j, m)
            vals = []
            for x in sq:
                qt = q_work.copy()
                qt[j] = x
                vals.append(get(tuple(qt)))
            coeff = fit_model(j, m, sq, vals)
            q_work[j] = best_dyadic_q(coeff, m)
            e_work = get(tuple(q_work))

        q_new = tuple(q_work)
        if e_work <= e + 1e-10:
            q, e = q_new, e_work
        else:
            return (e, q), eval_count, sweep
        if q == old_q:
            return (e, q), eval_count, sweep

    return (e, q), eval_count, MAX_SWEEPS

def main() -> None:
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    r = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    mode = sys.argv[3].lower() if len(sys.argv) > 3 else "analytic"
    if r < 1:
        raise SystemExit("r must be at least 1")
    if mode not in {"analytic", "direct", "analytic-direct", "balanced"}:
        raise SystemExit("mode must be analytic or direct")

    d = fetch(k)
    h = build_h(d)
    if mode == "analytic":
        (best_e, best_q), evals, sweeps = analytic_search(h, r)
        mode_label = "analytic coarse-seeded"
    else:
        (best_e, best_q), evals, sweeps = analytic_direct_search(h, r)
        mode_label = "analytic-direct balanced-start"
    ground = float(np.min(np.linalg.eigvalsh(h)).real)
    exhaustive = (2 ** (r - 1)) ** 3

    print(f"k={k}, r={r}")
    print(f"mode                    = {mode_label}")
    if mode != "analytic" and r >= 3:
        print(f"balanced start q        = {list(balanced_seed(r))}")
    print(f"best q                  = {list(best_q)}")
    print(f"analytic VQE scaled     = {best_e:.12f}")
    print(f"ground H_k scaled       = {ground:.12f}")
    print(f"variational gap scaled  = {best_e-ground:.12e}")
    print(f"analytic sweeps         = {sweeps}")
    print(f"energy evaluations      = {evals}")
    print(f"exhaustive target grid  = {exhaustive}")
    print(f"evaluation reduction    = {exhaustive/evals:.1f}x")
    print(f"VQE energy (mH)         = {best_e/d['scale']:.12f}")
    print(f"physical exact (mH)     = {d['physical_ground_mh']:.12f}")


if __name__ == "__main__":
    main()
