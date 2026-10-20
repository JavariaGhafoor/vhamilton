# One generator per clock: a noncommuting pair and a free clock

Run: `python3 -m scripts.trotter_search.experiment pair > results/trotter_pair/run.log`.
Haskell check: `trotter-chem-energies pair …` via
`python3 -m scripts.trotter_search.validate_against_haskell BIN 6 pair` (`validate.log`).
Same target, grid and search ports as `results/trotter_search` (k = 14, r = 7,
2^18 valuations); see that README for the conventions (regret, atoms, ports).

## Ansatz

    [ α1 YI ∥ α2 ZY ]   compiled at level n   (YI, ZY anticommute: [YI, ZY] = 2i XY)
    [ α3 IY ]           compiled at level 0   (IY commutes with both: exact, run once)

Because IY commutes with the pair, e^{-2πi(α1 YI + α2 ZY + α3 IY)} =
e^{-2πi α3 IY} e^{-2πi(α1 YI + α2 ZY)} exactly, so the program is the single
parallel block over all three generators; only the pair is Trotterized.
Atoms per energy query: 2·2^n + 1 (the controlled example: 5·2^n).

It spans the same Paulis as the controlled example of `results/trotter_search`
(α2 G2 + α3 G3 = ((α2+α3)/2) IY + ((α2−α3)/2) ZY), so the program reaches the
same optimum, −0.233566 (0.006 mHa above the ground energy of H_k), here at
ν = (28, 3, 19).  There IY was hidden inside both controlled clocks, so both
clashed with YI; written with one generator per clock, the free clock is
visible in the syntax.  (YI ∥ XY instead of ZY only reaches 0.33 mHa.)

## Results

**Level changes the answer, more sharply.**  The n = 0 circuit's optimum
(36, 37, 57) is 33.1 mHa off as a valuation of the program (the large clock
times make the n = 0 product far from the simultaneous evolution).  Regret of
the level optimum:

| n | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| canon | 33.1 | 0.125 | 0.018 | 0.018 | **0** | 0 |
| sym | 33.1 | 0.070 | **0** | 0 | 0 | 0 |

At canon n = 1, 2 the compiled optimum is again 0.005 mHa below anything the
program denotes.  Averaged over the whole grid the symmetric schedule is not
uniformly better (it loses where clocks are large); near the optimum (small
α1, α2) its error is 2× (n = 2) to 8× (n = 4) smaller, and it reaches the
program's optimum at n = 2 instead of n = 4.

**Sec. 5 methods per level.**  Refine and coordinate recover each level's
optimum (≈ 290 and ≈ 165 queries for n ≥ 1).  At n = 0 they disagree with
exhaustive search: they stop at (28, 5, 25), 2.1 mHa off, while the n = 0
optimum is the 33 mHa trap.

Analytic fits: the free clock keeps its n = 0 budget (3 samples) at every
level; only α1, α2 need 2^{n+1} + 1 samples.

* aware budget: exact fits (residual 1e-15), recovers each level's optimum;
  analytic-direct costs 50, 52, 86–88, 135–152, 243 queries for n = 0…4.
  Feasible up to n = 4 (33 samples on a 64-point axis); fails at n = 5 (65).
  The controlled example already failed at n = 4.
* naive budget (n = 0 counts at every level): residuals 3e-2 to 4e-1, and the
  result misses the level optimum at every n ≥ 1 except sym n = 1
  (true regret up to 0.57 mHa).

**Search in the Trotter dimension.**  One query at n + 1: jumps of
+26 / +21 mHa at n = 0 (canon / sym), +0.045 at n = 1, ≤ 0.007 from n = 2.

Continuation (analytic-direct at n = 0, then warm-started axial search per level):

| | reaches regret 0 at | queries | atoms |
|---|---|---|---|
| canon climb | n = 4 | 181 | 1,709 |
| sym climb | n = 2 | 175 | **867** |
| coordinate, cold, canon n = 4 | – | 163 | 5,379 |
| coordinate, cold, canon n = 5 | – | 163 | 10,595 |
| exhaustive, canon n = 4 | – | 262,144 | 8,650,752 |

Compared with the controlled example (2,260 atoms for the sym climb), writing
one generator per clock and compiling the free clock once costs 2.6× fewer
atoms for the same program optimum.
