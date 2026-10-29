#!/usr/bin/env python3
"""Render the QAOA grids at several precisions as side-by-side TikZ heatmaps.

Usage: python3 qaoa_ladder_tikz.py [--size CM] [--gap CM] qaoa_grid_r3.csv qaoa_grid_r4.csv qaoa_grid_r5.csv > qaoa_ladder.tex
Each panel covers one period, gamma, beta in [0, 1/2), at the same physical size,
so a finer grid shows as smaller cells. The search's reported optimum is starred.
"""
import csv
import re
import sys

SIZE = 2.4    # cm per panel side (--size)
GAP = 0.75    # cm between panels (--gap)
args = sys.argv[1:]
while args and args[0].startswith("--"):
    opt, val, args = args[0], float(args[1]), args[2:]
    if opt == "--size":
        SIZE = val
    elif opt == "--gap":
        GAP = val
    else:
        sys.exit("unknown option " + opt)


def colour(e):
    return (r"blue!%d!white" % round(-e * 85)) if e < 0 else (r"red!%d!white" % round(e * 85))


out = [r"\begin{tikzpicture}[font=\scriptsize]"]
for i, path in enumerate(args):
    r = int(re.search(r"_r(\d+)\.csv$", path).group(1))
    rows = list(csv.DictReader(open(path)))
    n = 1 + max(int(x["q_gamma"]) for x in rows)
    c = SIZE / n
    x0 = i * (SIZE + GAP)
    best = min(rows, key=lambda x: float(x["energy"]))
    for x in rows:
        g, b, e = int(x["q_gamma"]), int(x["q_beta"]), float(x["energy"])
        out.append(r"\fill[%s] (%.4f,%.4f) rectangle ++(%.4f,%.4f);" % (colour(e), x0 + g * c, b * c, c, c))
    if n <= 4:  # coarse grid: show the cells, so an all-zero panel reads as data
        out.append(r"\draw[gray!35, step=%.4f, xshift=%.4fcm] (0,0) grid (%.4f,%.4f);" % (c, x0, SIZE, SIZE))
    out.append(r"\draw[gray] (%.4f,0) rectangle ++(%.4f,%.4f);" % (x0, SIZE, SIZE))
    if float(best["energy"]) < -0.5:
        bg, bb = int(best["q_gamma"]), int(best["q_beta"])
        out.append(r"\node[star, star points=5, fill=yellow, draw=black, inner sep=1.3pt] at (%.4f,%.4f) {};"
                   % (x0 + (bg + 0.5) * c, (bb + 0.5) * c))
    out.append(r"\node[above] at (%.4f,%.4f) {$r=%d$: $%d\times%d$};" % (x0 + SIZE / 2, SIZE, r, n, n))
    out.append(r"\node[below] at (%.4f,0) {$\gamma$};" % (x0 + SIZE / 2))
    if i == 0:
        out.append(r"\node[left] at (%.4f,%.4f) {$\beta$};" % (x0, SIZE / 2))
# colour bar
xb = len(args) * (SIZE + GAP) - GAP + 0.3
for k in range(16):
    e = -1 + 2 * (k + 0.5) / 16
    out.append(r"\fill[%s] (%.4f,%.4f) rectangle ++(0.22,%.4f);" % (colour(e), xb, k * SIZE / 16, SIZE / 16))
out.append(r"\draw[gray] (%.4f,0) rectangle ++(0.22,%.4f);" % (xb, SIZE))
out.append(r"\node[right] at (%.4f,0.08) {$-1$}; \node[right] at (%.4f,%.4f) {$0$}; \node[right] at (%.4f,%.4f) {$+1$};"
           % (xb + 0.22, xb + 0.22, SIZE / 2, xb + 0.22, SIZE - 0.08))
out.append(r"\node[above] at (%.4f,%.4f) {$E$};" % (xb + 0.11, SIZE))
out.append(r"\end{tikzpicture}")
print("\n".join(out))
