# Search on a Trotterized ansatz, and search in the Trotter dimension

Code: `scripts/trotter_search/` (driver `experiment.py`), Haskell check
`app/TrotterChemEnergies.hs`.  Raw output: `run.log`, `landscapes.csv`,
`methods.csv`, `certify.csv`, `continuation.csv`, `validate.log`.

Reproduce:

```
python3 -m scripts.trotter_search.check_port          # ports vs Sec. 5 Haskell runs
python3 -m scripts.trotter_search.experiment > results/trotter_search/run.log
python3 -m scripts.trotter_search.validate_against_haskell <trotter-chem-energies binary>
```

## Setup

Same target as Sec. 5 (DMANH+ upper block, k = 14, r = 7, 2^18 valuations),
but the three generators are written as **one parallel block**

    G(α) = α1 YI ∥ α2 (IY + ZY)/2 ∥ α3 (IY − ZY)/2 .

The program denotes e^{-2πi G(α)}.  YI anticommutes with ZY, so the compiled
circuit depends on the Trotter level n: 2^n repetitions of the 5-atom
schedule, 5·2^n atoms per energy query.  At n = 0 the schedule runs the three
generators one after another, which is *exactly* the Sec. 5 sequential
ansatz.  `canon` is Hamilton's compiler; `sym` alternates the schedule and its
reverse between repetitions (same atoms, same count, symmetric product).

"Regret" is the program's own energy E(ν) at the returned valuation minus the
program's optimum on the grid, E(28, 22, 16) = −0.233566 (0.006 mHa above the
ground energy of H_k).  1 mHa = 0.019212 in the scaled units of the logs.

Energies come from a NumPy replica of the compiled circuits
(`vh_numpy.py`); `validate.log` compares it with the Haskell pipeline
(compiler → OpenQASM → Qiskit).  The search modes are Python ports of
`app/ChemVQE.hs`; `check_port.py` shows analytic-direct reproduces the
Haskell run exactly (38 queries) and the other modes reach the same optimum
(28, 30, 19) with 307 / 185 / 133 queries against Haskell's 315 / 189 / 128:
their beams break ties between energies equal to ~1e-12 differently.  All
numbers below are from the ports, so they compare like with like.

## 1. The level changes the answer

| n | sched | compiled optimum ν_n | true regret (mHa) | mean ∣E_n − E∣ (mHa) |
|---|---|---|---|---|
| 0 | either | (28,30,19) | **2.115** | 4.94 |
| 1 | canon | (28,24,18) | 0.125 | 1.37 |
| 2 | canon | (28,23,17) | 0.018 | 0.58 |
| 3 | canon | (28,23,17) | 0.018 | 0.28 |
| 4 | canon | (28,22,17) | 0.006 | 0.14 |
| 5 | canon | (28,22,16) | 0 | 0.07 |
| 1 | sym | (28,21,17) | 0.070 | 5.84 |
| 2 | sym | (28,22,17) | 0.006 | 0.79 |
| 3 | sym | (28,22,16) | **0** | 0.17 |
| 4 | sym | (28,22,16) | 0 | 0.04 |
| 5 | sym | (28,22,16) | 0 | 0.01 |

* The Sec. 5 answer (28,30,19) is the optimum of the n = 0 circuit; as a
  valuation of the parallel program it is 2.1 mHa off, more than chemical
  accuracy (1.6 mHa).
* At canon n = 1, 2 the compiled optimum lies 0.005 mHa *below* anything the
  program denotes (the Sec. 4 effect, small here).
* With the symmetric schedule the landscape error falls ≈4× per level instead
  of ≈2×, and the program's optimum is reached at n = 3 instead of n = 5.

## 2. Sec. 5 methods at a fixed level

Refine (zoom) and coordinate search recover the exhaustive optimum of E_n at
every level and schedule, with an almost level-independent query count
(≈316 and 163–185).  Their cost in atoms doubles per level.

The analytic fits are where Trotterization bites:

* **naive** (the n = 0 frequency budget of ChemVQE used at every level): the
  model is now wrong, fit residuals ≈ 1e-3 (the Haskell code would print its
  warning), and the result is off by up to 0.02 mHa (e.g. sym n = 3–5 returns
  (29,22,17) while the level optimum is (28,22,16)).
* **aware** (budget read off the program and n): a clock whose factor has
  eigenvalues Λ appears 2^n times, so its energy curve carries frequencies
  (s − s')/2^n with s, s' sums of 2^n elements of Λ: 2^n frequencies for YI,
  2^{n+1} for the controlled clocks, i.e. 2^{n+1}+1 and 2^{n+2}+1 samples.
  The fits are then exact (residual 1e-16) and recover the level optimum, but
  analytic-direct costs 38 → 77 → 113 → ≈300 queries for n = 0…3, the
  design matrices become badly conditioned (cond 2e6 at n = 1, 1e18 at n ≥ 2),
  and at n ≥ 4 the controlled clocks need more samples (65) than the axis has
  grid points (64), so they cannot be fitted at all and the method fails
  (regret 1.3–11 mHa).

So precision r stays free, but the frequency budget, and with it the analytic
method's cost, grows as 2^n.

## 3. Optimizing in the Trotter dimension

**Certificate (1 query).**  Re-evaluating the level-n optimum at level n+1
flags a premature level: the jump is +0.87 (canon) / +3.6 (sym) mHa at
n = 0, +0.04 at n = 1, ≤ 0.007 from n = 2, and 0.000 for sym from n = 3.
It is of the order of the true regret (from 1.7× above it, sym n = 0, to
about 3× below it, canon n = 1–2) and shrinks with it, so
"increase n until the jump is below τ" is a usable stopping rule.

**Continuation in n.**  Valuations do not depend on n (the valuation picks
the trajectory, n only compiles it), so the optimum at level n−1 is a valid
start at level n, exactly as q ↦ 2q nests precisions.  Analytic-direct at
n = 0, then warm-started axial search at each next level:

| schedule | reaches regret 0 at | total queries | total atoms |
|---|---|---|---|
| canon | n = 5 | 193 | 6,410 |
| sym | n = 3 | 171 | 2,260 |

against a cold start at the level that is needed:

| method at canon n = 5 | queries | atoms |
|---|---|---|
| exhaustive | 262,144 | 41,943,040 |
| refine | 316 | 50,560 |
| coordinate | 177 | 28,320 |
| analytic-direct (naive) | 37 | 5,920 (regret 0.020) |

Most queries are spent at cheap levels; the expensive levels only polish.
Continuation with the symmetric schedule reaches the program's optimum with
12.5× fewer atom applications than coordinate search at canon n = 5, and
18,500× fewer than exhaustive search there.  Warm starts with the aware fit
instead of axial moves cost more (275 / 331 queries) and, for canon, stop at
(28,23,17) because the fit is unavailable at n ≥ 4.
