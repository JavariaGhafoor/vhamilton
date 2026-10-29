"""pgfplots figure of analytic-direct cost against clock precision r, from
results/hamilton_vqe/12_analytic_direct_precision_scaling/direct_scaling_results.csv:
(a) energy queries Q against the size 2^{3(r-1)} of the one-period grid;
(b) coordinate sweeps L, with the per-sweep query bound Q <= 1 + 16 L checked.

Usage:
  python -m scripts.experiments.direct_scaling_figure \
    --csv results/hamilton_vqe/12_analytic_direct_precision_scaling/direct_scaling_results.csv \
    --out results/hamilton_vqe/12_analytic_direct_precision_scaling/direct_scaling.tex
"""
import argparse
import csv
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    for x in rows:
        assert int(x["queries"]) <= 1 + 16 * int(x["sweeps"])
        assert int(x["grid_size"]) == 2 ** (3 * (int(x["r"]) - 1))
    ks = sorted({int(x["k"]) for x in rows})
    rs = sorted({int(x["r"]) for x in rows})
    style = {8: ("black!35", "square*"), 10: ("black!60", "triangle*"), 14: ("violet", "*")}

    def series(k, key):
        return " ".join(f"({x['r']},{x[key]})" for x in rows if int(x["k"]) == k)

    out = [
        r"\begin{tikzpicture}",
        r"\begin{groupplot}[group style={group size=2 by 1, horizontal sep=1.5cm},",
        r"  width=0.42\linewidth, height=4.6cm, grid=major, grid style={gray!20},",
        r"  xlabel={clock precision $r$}, xtick={" + ",".join(map(str, rs)) + "},",
        r"  tick label style={font=\scriptsize}, label style={font=\footnotesize},",
        r"  legend style={font=\scriptsize, draw=none, fill=white}, legend cell align=left]",
        r"\nextgroupplot[ymode=log, ymin=5, ymax=1e9, title={\footnotesize (a) energy queries $Q$},",
        r"  legend pos=north west]",
        r"\addplot[thick, dashed, black!60] coordinates {"
        + " ".join(f"({r},{2 ** (3 * (r - 1))})" for r in rs) + "};",
        r"\addlegendentry{grid size $2^{3(r-1)}$}",
    ]
    for k in ks:
        col, mk = style[k]
        out += [rf"\addplot[thick, {col}, mark={mk}, mark size=1.6pt] coordinates {{{series(k, 'queries')}}};",
                rf"\addlegendentry{{$Q$, $k={k}$}}"]
    out += [r"\nextgroupplot[ymin=0, ymax=6, ytick={0,...,6}, title={\footnotesize (b) coordinate sweeps $s$},",
            r"  legend pos=south east]"]
    for k in ks:
        col, mk = style[k]
        out += [rf"\addplot[thick, {col}, mark={mk}, mark size=1.6pt] coordinates {{{series(k, 'sweeps')}}};",
                rf"\addlegendentry{{$k={k}$}}"]
    out += [r"\end{groupplot}", r"\end{tikzpicture}", ""]
    args.out.write_text("\n".join(out))


if __name__ == "__main__":
    main()
