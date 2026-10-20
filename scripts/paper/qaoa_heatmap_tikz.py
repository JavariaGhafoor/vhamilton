#!/usr/bin/env python3
"""Render qaoa_grid_r5.csv (from `qaoa-demo 5`) as a TikZ heatmap for the paper.

Usage: python3 qaoa_heatmap_tikz.py qaoa_grid_r5.csv > qaoa_heatmap.tex
Overlays the r=3 grid (circles) and the r=4 optimum q=(2,3) (star).
"""
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1])))
n = 1 + max(int(r["q_gamma"]) for r in rows)          # 16 at r=5
cell = 0.16                                            # cm
out = [r"\begin{tikzpicture}[x=%gcm,y=%gcm,font=\footnotesize]" % (cell, cell)]
for r in rows:
    g, b, e = int(r["q_gamma"]), int(r["q_beta"]), float(r["energy"])
    col = (r"blue!%d!white" % round(-e * 85)) if e < 0 else (r"red!%d!white" % round(e * 85))
    out.append(r"\fill[%s] (%d,%d) rectangle ++(1,1);" % (col, g, b))
out.append(r"\draw[gray] (0,0) rectangle (%d,%d);" % (n, n))
# r=3 grid: q3 in {0..3}  ->  r=5 index 4*q3
for g in range(4):
    for b in range(4):
        out.append(r"\draw[thick] (%g,%g) circle (3.2pt);" % (4 * g + 0.5, 4 * b + 0.5))
# r=4 optimum (2,3) -> r=5 index (4,6)
out.append(r"\node[star, star points=5, fill=yellow, draw=black, inner sep=1.6pt] at (4.5,6.5) {};")
out.append(r"\node[below] at (%g,0) {$\gamma=\nu(\gamma)/32\in[0,\tfrac12)$};" % (n / 2))
out.append(r"\node[rotate=90, above] at (0,%g) {$\beta\in[0,\tfrac12)$};" % (n / 2))
# colour bar
out.append(r"\begin{scope}[shift={(%d,0)}]" % (n + 2))
for i in range(16):
    e = -1 + 2 * (i + 0.5) / 16
    col = (r"blue!%d!white" % round(-e * 85)) if e < 0 else (r"red!%d!white" % round(e * 85))
    out.append(r"\fill[%s] (0,%d) rectangle ++(1,1);" % (col, i))
out.append(r"\draw[gray] (0,0) rectangle (1,16);")
out.append(r"\node[right] at (1,0.5) {$-1$}; \node[right] at (1,8) {$0$}; \node[right] at (1,15.5) {$+1$};")
out.append(r"\node[above] at (0.5,16) {$E$};")
out.append(r"\end{scope}")
out.append(r"\end{tikzpicture}")
print("\n".join(out))
