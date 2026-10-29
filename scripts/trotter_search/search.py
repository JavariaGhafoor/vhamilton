"""Python ports of the search modes of app/ChemVQE.hs, generalised to a
Trotter-aware frequency budget.

Every mode talks to the landscape only through an Oracle, which counts each
distinct energy query exactly as the Haskell code does (refine mode re-queries
inherited centres; the other modes cache).  With n = 0 and the ChemVQE ansatz
these ports reproduce the Haskell query counts 315 / 189 / 128 / 38 and the
optimum q = (28, 30, 19); see check_port.py.
"""
import itertools
import math
import numpy as np

COARSE_R = 3
BEAM = 3
COORD_MAX_ROUNDS = 20
ANALYTIC_MAX_SWEEPS = 12
SCAN_BINS = 256
GOLDEN_IT = 70


class Oracle:
    """Energy queries at one fixed compilation (ansatz, n, schedule).
    energy_fn(q_tuple, r) -> float.  Counts every call (no caching here)."""

    def __init__(self, energy_fn, cost_per_query=1):
        self.f = energy_fn
        self.queries = 0
        self.cost_per_query = cost_per_query

    def __call__(self, q, r):
        self.queries += 1
        return self.f(tuple(q), r)

    @property
    def cost(self):
        return self.queries * self.cost_per_query


def grid(r):
    return range(2 ** (r - 1))


def best_beam(results):
    # results: list of (E, q); stable sort by energy like Data.List.sortBy
    return sorted(results, key=lambda t: t[0])[:BEAM]


def _coarse(oracle, p):
    res = [(oracle(q, COARSE_R), q) for q in itertools.product(grid(COARSE_R), repeat=p)]
    return res


# ---------------------------------------------------------------- refine (zoom)
def refine(oracle, p, r):
    res = _coarse(oracle, p)
    beam = best_beam(res)
    best = min(res, key=lambda t: t[0])
    for rr in range(COARSE_R + 1, r + 1):
        up = 2 ** (rr - 1) - 1
        cands = []
        for _, q in beam:
            axes = [[x for x in range(2 * qi - 1, 2 * qi + 2) if 0 <= x <= up] for qi in q]
            for c in itertools.product(*axes):
                if c not in cands:
                    cands.append(c)
        res = [(oracle(c, rr), c) for c in cands]
        beam = best_beam(res)
        best = min(res, key=lambda t: t[0])
    return best


# ---------------------------------------------------------------- coordinate
def axial(q, r):
    up = 2 ** (r - 1) - 1
    out = []
    for j, qj in enumerate(q):
        for x in (qj - 1, qj + 1):
            if 0 <= x <= up:
                out.append(q[:j] + (x,) + q[j + 1:])
    return out


def coordinate_level(oracle, r, inherited):
    beam = list(inherited)
    cache = {q: e for e, q in inherited}
    order = [q for _, q in inherited]          # cache order (newest first in Haskell)
    for _ in range(COORD_MAX_ROUNDS):
        proposed = []
        for _, q in beam:
            for c in axial(q, r):
                if c not in proposed:
                    proposed.append(c)
        unseen = [c for c in proposed if c not in cache]
        if not unseen:
            break
        new = [(oracle(c, r), c) for c in unseen]
        for e, c in new:
            cache[c] = e
        order = [c for _, c in new] + order
        nxt = best_beam([(cache[c], c) for c in order])
        stable = [q for _, q in nxt] == [q for _, q in beam]
        beam = nxt
        if stable:
            break
    return best_beam([(cache[c], c) for c in order])


def coordinate(oracle, p, r):
    res = _coarse(oracle, p)
    beam = best_beam(res)
    for rr in range(COARSE_R + 1, r + 1):
        inherited = [(e, tuple(2 * x for x in q)) for e, q in beam]
        beam = coordinate_level(oracle, rr, inherited)
    return beam[0]


# ---------------------------------------------------------------- analytic fits
class TrigModel:
    """E(phi) = c + sum_f a_f cos(f phi) + b_f sin(f phi), phi = 2π α."""

    def __init__(self, freqs, coeffs):
        self.freqs = freqs
        self.coeffs = coeffs

    def __call__(self, phi):
        phi = np.asarray(phi, dtype=float)
        v = self.coeffs[0] * np.ones_like(phi)
        for i, f in enumerate(self.freqs):
            v = v + self.coeffs[1 + 2 * i] * np.cos(f * phi) + self.coeffs[2 + 2 * i] * np.sin(f * phi)
        return v


def design(freqs, phis):
    cols = [np.ones_like(phis)]
    for f in freqs:
        cols += [np.cos(f * phis), np.sin(f * phis)]
    return np.stack(cols, axis=1)


def sample_qs(freqs, m, legacy):
    if legacy is not None:
        return legacy
    s = 1 + 2 * len(freqs)
    if s > m:
        return None                     # more unknowns than grid points on this axis
    return sorted({(i * m) // s for i in range(s)})


def golden(model, a, b):
    rho = (math.sqrt(5) - 1) / 2
    c, d = b - rho * (b - a), a + rho * (b - a)
    fc, fd = model(c), model(d)
    for _ in range(GOLDEN_IT):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - rho * (b - a)
            fc = model(c)
        else:
            a, c, fc = c, d, fd
            d = a + rho * (b - a)
            fd = model(d)
    return (a + b) / 2


def best_q(model, m):
    phis = math.pi * np.arange(SCAN_BINS + 1) / SCAN_BINS
    v = model(phis)
    mins = [0.0, math.pi]
    for i in range(1, SCAN_BINS):
        if v[i] <= v[i - 1] and v[i] <= v[i + 1]:
            mins.append(golden(model, phis[i - 1], phis[i + 1]))
    cands = [0, m - 1]
    for ph in mins:
        qc = m * ph / math.pi
        for x in (math.floor(qc), math.ceil(qc), round(qc) - 1, round(qc), round(qc) + 1):
            x = max(0, min(m - 1, x))
            if x not in cands:
                cands.append(x)
    return min(cands, key=lambda x: float(model(math.pi * x / m)))


class FitSpec:
    """Per-clock frequency sets (in units of phi) and optional legacy sample
    positions (the Haskell n = 0 sample design)."""

    def __init__(self, freqs, legacy=None):
        self.freqs = freqs
        self.legacy = legacy


def analytic_pass(oracle, r, specs, q, e_cur, cache, log):
    m = 2 ** (r - 1)
    for j, spec in enumerate(specs):
        legacy = spec.legacy(m) if spec.legacy else None
        sq = sample_qs(spec.freqs, m, legacy)
        if sq is None:
            log.append(("infeasible", j))
            continue
        tuples = [q[:j] + (x,) + q[j + 1:] for x in sq]
        es = []
        for t in tuples:
            if t not in cache:
                cache[t] = oracle(t, r)
            es.append(cache[t])
        A = design(spec.freqs, math.pi * np.array(sq, dtype=float) / m)
        coef = np.linalg.lstsq(A, np.array(es), rcond=None)[0]
        log.append(("cond", j, float(np.linalg.cond(A))))
        model = TrigModel(spec.freqs, coef)
        qb = best_q(model, m)
        t = q[:j] + (qb,) + q[j + 1:]
        if t not in cache:
            cache[t] = oracle(t, r)
        log.append(("residual", j, cache[t] - float(model(math.pi * qb / m))))
        q, e_cur = t, cache[t]
    return q, e_cur


def analytic(oracle, p, r, specs, log=None):
    log = [] if log is None else log
    res = _coarse(oracle, p)
    beam = best_beam(res)
    lift = 2 ** (r - COARSE_R)
    beam = [(e, tuple(lift * x for x in q)) for e, q in beam]
    cache = {q: e for e, q in beam}
    for _ in range(ANALYTIC_MAX_SWEEPS):
        updated = [analytic_pass(oracle, r, specs, q, e, cache, log) for e, q in beam]
        merged, seen = [], set()
        for q, e in updated + [(q, e) for e, q in beam]:
            if q not in seen:
                seen.add(q)
                merged.append((e, q))
        nxt = best_beam(merged)
        stable = [q for _, q in nxt] == [q for _, q in beam]
        beam = nxt
        if stable:
            break
    return beam[0]


def analytic_direct(oracle, p, r, specs, start=None, log=None, cache=None):
    log = [] if log is None else log
    q = tuple(start) if start is not None else (2 ** (r - 3),) * p
    cache = {} if cache is None else cache
    if q not in cache:
        cache[q] = oracle(q, r)
    e = cache[q]
    for _ in range(ANALYTIC_MAX_SWEEPS):
        q2, e2 = analytic_pass(oracle, r, specs, q, e, cache, log)
        if not e2 <= e + 1e-10:
            break
        if q2 == q:
            break
        q, e = q2, e2
    return (e, q)


# ---------------------------------------------------------------- frequency budget
def clock_frequencies(eigs, n):
    """Frequencies (units of phi = 2π α) of the energy along one clock whose
    per-repetition generator has eigenvalues `eigs` (integers after scaling),
    when the block is compiled at Trotter level n: the clock's factor appears
    2^n times with eigenphases λ φ / 2^n, so the state carries phases s φ / 2^n
    with s a sum of 2^n eigenvalues, and the energy carries their differences."""
    sums = {0}
    for _ in range(2 ** n):
        sums = {s + l for s in sums for l in eigs}
    diffs = sorted({abs(a - b) for a in sums for b in sums} - {0})
    return [d / 2 ** n for d in diffs]
