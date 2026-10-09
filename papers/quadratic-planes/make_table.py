"""Write graphs-table.tex (the table of section 2) from data/quadratic_planes/q*.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "..", "data", "quadratic_planes")
rows = []
for d in sorted(int(f[1:-5]) for f in os.listdir(D) if f.startswith("q") and f.endswith(".json")):
    g = json.load(open(os.path.join(D, f"q{d}.json"))); P = g["points"]; E = g["edges"]; n = len(P)
    deg = [0] * n; dirs = set()
    for a, b in E:
        deg[a] += 1; deg[b] += 1
        v = tuple(P[b][i] - P[a][i] for i in range(4)); dirs.add(max(v, tuple(-x for x in v)))
    rows.append((d, g["D"], n, len(E), min(deg), max(deg), 2 * len(dirs)))
num = lambda x: f"{x // 1000}\\,{x % 1000:03d}" if x >= 1000 else str(x)
with open(os.path.join(HERE, "graphs-table.tex"), "w") as f:
    f.write("% written by make_table.py from data/quadratic_planes/\n")
    f.write("\\begin{tabular}{rrrrcr}\n\\hline\n$d$ & $D$ & vertices & edges & degrees & directions used\\\\\n\\hline\n")
    for d, Dd, n, m, lo, hi, k in rows:
        f.write(f"{d} & {Dd} & {num(n)} & {num(m)} & {lo}--{hi} & {k}\\\\\n")
    f.write("\\hline\n\\end{tabular}\n")
print("wrote graphs-table.tex:", [(r[0], r[2]) for r in rows])
