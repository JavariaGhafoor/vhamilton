# 14 Exhaustive certificate at (k, r) = (14, 10)

Evaluates the chemistry VQE energy at all 2^27 = 134,217,728 valuations of the half-range grid at r = 10, against H_14. This certifies that step 5 of the search (`10_analytic_direct_representative`, `12_analytic_direct_precision_scaling`) returns the grid optimum, not just a valuation within the constructive bound.

`certificate.log` lists the three lowest valuations:

- The optimum, (223, 238, 155), is 0.000119 mHa above E0(H_14).
- The runner-up, (223, 239, 155), is 5.0e-5 mHa higher than the optimum, far above floating-point error.

Two further checks are numerical, not a formal error bound:

- The three lowest valuations are re-evaluated in NumPy `longdouble` (80-bit x87 extended precision on x86-64) and with mpmath at 50 significant digits. The double-precision energies are within 6e-17 of the 50-digit values, and the order and margin are unchanged.
- 200,000 random valuations (seed 0) are evaluated in both double and `longdouble`. The largest difference is 7.3e-16 (scaled).

The run takes about a minute; the exhaustive pass alone takes about 10 s in NumPy.
