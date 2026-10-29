#!/usr/bin/env python3
"""Check qaoa-demo output against the closed form E = sin(4*pi*gamma) * sin(8*pi*beta).

Usage: python3 verify_qaoa_closed_form.py results/qaoa_one_edge/qaoa_grid_r*.csv
The grid file name gives the precision r; valuations are q / 2^r.
"""
import csv
import math
import re
import sys

for path in sys.argv[1:]:
    r = int(re.search(r"_r(\d+)\.csv$", path).group(1))
    rows = list(csv.DictReader(open(path)))
    err = max(abs(float(x["energy"])
                  - math.sin(4 * math.pi * int(x["q_gamma"]) / 2 ** r)
                  * math.sin(8 * math.pi * int(x["q_beta"]) / 2 ** r)) for x in rows)
    print(f"{path}: {len(rows)} valuations, max |E - closed form| = {err:.3e}, "
          f"max |E| = {max(abs(float(x['energy'])) for x in rows):.6f}")
