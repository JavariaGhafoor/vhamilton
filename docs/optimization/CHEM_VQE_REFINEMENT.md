# Dyadic VQE search modes

`chem-vqe` now provides five classical search modes over the same discrete
variational family

```text
alpha_j = q_j / 2^r,
q_j in {0,...,2^(r-1)-1}.
```

The modes are:

- `refine`: full local dyadic beam refinement;
- `coordinate`: coordinate-wise dyadic beam refinement;
- `analytic`: analytic/trigonometric coordinate minimization with a coarse beam
  seed and a direct jump to the requested dyadic precision;
- `analytic-direct`: the same analytic coordinate minimization started from a
  deterministic balanced interior point, with no coarse exhaustive search;
- `exhaustive`: full target-grid enumeration for validation.

All modes return an integer tuple `q` on the requested dyadic grid.  The
non-exhaustive modes differ only in how they choose which VQE energy values to
query.

## 1. Exhaustive reference

With `p` parameters and `2^(r-1)` values per parameter, exhaustive search uses

```text
N_exhaustive = (2^(r-1))^p = O(2^(p r)).
```

For the current `p=3` ansatz this is `2^(3(r-1))` candidates.

Run:

```bash
cabal run chem-vqe -- 10 5 0 exhaustive
```

This remains the global finite-grid reference when the grid is small enough.

## 2. Full-neighborhood dyadic refinement (`refine`)

At precision `r`, a parameter is `alpha=q/2^r`.  Moving to precision `r+1`,
`2q` represents exactly the same physical parameter:

```text
q / 2^r = (2q) / 2^(r+1).
```

The search therefore begins with a small exhaustive grid at `r0=min(3,r)`,
keeps a beam of the best three points, doubles them at the next precision, and
checks the full Cartesian `+/-1` neighborhood

```text
q_j' in {2q_j-1, 2q_j, 2q_j+1}.
```

For `p` parameters this gives at most `3^p` points per beam center and roughly

```text
O(B * 3^p * r)
```

energy evaluations after including all refinement levels.

Run:

```bash
cabal run chem-vqe -- 10 5 0 refine
```

## 3. Coordinate dyadic refinement (`coordinate`)

Coordinate mode keeps the same nested-dyadic idea, but tests only one-coordinate
moves

```text
q -> q +/- e_j.
```

There are at most `2p` new axial points per beam center per local round rather
than `3^p`.  Inherited doubled centers reuse their already-known energies.
The resulting energy-query complexity is approximately

```text
O(B * p * L * r),
```

where `L` is the number of axial rounds required at each precision level.

Run:

```bash
cabal run chem-vqe -- 10 5 0 coordinate
# shorthand
cabal run chem-vqe -- 10 5 0 coord
```

For the `k=14,r=7` compiler-backed experiment, both `refine` and `coordinate`
returned `q=[28,30,19]`.  The full-neighborhood method used 315 energy
evaluations, whereas coordinate refinement used 189, compared with 262,144
points in the exhaustive target grid.

## 4. Analytic dyadic coordinate minimization (`analytic`)

The current ansatz has additional algebraic structure that can remove the
multiplicative target-precision factor from the expensive energy-query count.
Define

```text
phi = 2*pi*alpha = pi*q / 2^(r-1).
```

For one coordinate at a time, with the other parameters fixed, the VQE energy
is a low-degree trigonometric polynomial determined by the generator spectrum.

For

```text
G1 = YI,
```

whose eigenvalues are `{-1,+1}`, only the frequency difference 2 occurs, so

```text
E_1(phi) = c + a2 cos(2 phi) + b2 sin(2 phi).
```

Three VQE energy values determine this curve exactly.

For

```text
G2 = (IY+ZY)/2,
G3 = (IY-ZY)/2,
```

whose eigenvalues are `{-1,0,+1}`, frequency differences 1 and 2 can occur:

```text
E_j(phi) = c
         + a1 cos(phi) + b1 sin(phi)
         + a2 cos(2 phi) + b2 sin(2 phi).
```

Five VQE energy values determine this curve.  `analytic` fits that curve,
finds its continuous local minima using only classical arithmetic, projects
those minima back to the allowed target dyadic grid, and queries the predicted
point.  Coordinates are updated sequentially and the beam is swept until it no
longer changes.

Unlike `refine` and `coordinate`, it does **not** traverse

```text
r0 -> r0+1 -> ... -> r_target.
```

It lifts each coarse seed directly from `r0=3` to the target precision using

```text
q_target = 2^(r_target-r0) q_coarse,
```

which represents exactly the same physical parameter, and then performs the
analytic coordinate sweeps only at the target `r`.

Run:

```bash
cabal run chem-vqe -- 14 7 0 analytic
# alias
cabal run chem-vqe -- 14 7 0 rotosolve
```

Its expensive VQE-energy-query complexity is approximately

```text
O(2^(p*r0) + B*p*L),
```

for fixed harmonic degree.  Since `r0=3` is fixed, this has no multiplicative
`r_target` factor in the energy-query count.  There is still a small fixed
amount of classical trigonometric fitting/minimization per coordinate.

The harmonic model is exact for the three generators used in the present
benchmark.  The code reports a model residual at each analytic jump; a large
residual is a warning that the assumed low-degree harmonic structure no longer
matches the compiled generator.  This becomes important if the ansatz is later
changed to noncommuting Trotterized blocks.


## 5. Deterministic balanced analytic search (`analytic-direct`)

The remaining exponential-looking term in `analytic`, `2^(p*r0)`, comes only
from the coarse exhaustive seed.  It is not required by the one-coordinate
trigonometric fit.  For the present three-parameter state-preparation ansatz we
can construct a nondegenerate interior seed directly.

Write

```text
beta_j = 2*pi*alpha_j,
alpha_j = q_j / 2^r.
```

Choose the balanced branch angle

```text
beta_j = pi/4.
```

Then `alpha_j=1/8`, so for every `r>=3` the exact dyadic integer is

```text
q_j = 2^(r-3).
```

For the current ansatz this gives

```text
|psi_0> = [1,1,1,1]^T / 2 = |++>,
```

so the first-qubit `|0>` and `|1>` branches are both populated and neither
conditional generator begins on an empty branch.  Constructing the tuple
requires only one integer value per variational parameter: `O(p)` ordinary
arithmetic and **zero search evaluations**.  The executable evaluates the
initial seed once to record its objective value, then performs the same exact
one-coordinate trigonometric fits directly at target precision.

Run:

```bash
cabal run chem-vqe -- 14 7 0 analytic-direct
# aliases
cabal run chem-vqe -- 14 7 0 direct
cabal run chem-vqe -- 14 7 0 balanced
```

For `k=14,r=7`, the independent matrix verifier gives

```text
balanced start q        = [16,16,16]
best q                  = [28,30,19]
energy evaluations      = 38
exhaustive target grid  = 262144
evaluation reduction    = 6898.5x
```

The expensive energy-query complexity is approximately

```text
O(p * L)
```

for fixed harmonic degree, where `L` is the number of analytic coordinate
sweeps.  Thus neither an exponential coarse seed nor a multiplicative target
precision factor remains in the VQE-query count.

This does **not** prove polynomial-time global VQE optimization.  The balanced
seed is constructed in polynomial time and avoids the obvious empty-branch
degeneracy of this ansatz, but analytic coordinate descent can still converge
to a nonglobal basin on a general nonconvex VQE landscape.

## Complexity summary

```text
exhaustive:  O(2^(p r))
refine:      O(B * 3^p * r)
coordinate:  O(B * p * L * r)
analytic:        O(2^(p*r0) + B * p * L)   [energy queries]
analytic-direct: O(p * L)                    [energy queries]
```

The four non-exhaustive modes are heuristic global optimizers.  `analytic` and
`analytic-direct` exactly minimize each fitted one-coordinate curve for this
ansatz, but sequential coordinate minimization (with or without a finite beam)
is not a general proof of the global multi-parameter optimum.  Use `exhaustive`
at feasible `r` and the independent NumPy verifier to validate representative
points.
