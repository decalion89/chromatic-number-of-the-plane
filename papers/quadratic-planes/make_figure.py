"""Write q11-graph.tex (Figure 1: the graph over Q(sqrt 11) in TikZ) from data/quadratic_planes/q11.json. The points are
placed at their real coordinates, rounded to five decimals for the drawing; each segment has length exactly 1. The
marker shape of a vertex is its colour in the stored 4-colouring."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
g = json.load(open(os.path.join(HERE, "..", "..", "data", "quadratic_planes", "q11.json")))
d, D = g["d"], g["D"]
r = math.sqrt(d)
P = [((a + b * r) / D, (c + e * r) / D) for a, b, c, e in g["points"]]
E, col = g["edges"], [int(c) for c in g["four_colouring"]]
assert all(abs(math.dist(P[a], P[b]) - 1) < 1e-9 for a, b in E)
mark = {0: ("circle, minimum size=3.6pt", "black"), 1: ("rectangle, minimum size=3.3pt", "white"),
        2: ("regular polygon, regular polygon sides=3, minimum size=5pt", "black"),
        3: ("diamond, minimum size=4.6pt", "white")}
out = ["% written by make_figure.py from data/quadratic_planes/q11.json",
       r"\begin{tikzpicture}[x=3.2cm, y=3.2cm]",
       r"\draw[gray!55, densely dashed, line width=0.35pt] (0,0) circle (1);"]
out += [r"\draw[gray!65, line width=0.22pt] (%.5f,%.5f) -- (%.5f,%.5f);" % (P[a] + P[b]) for a, b in E]
out += [r"\node[draw=black, line width=0.35pt, fill=%s, inner sep=0pt, %s] at (%.5f,%.5f) {};"
        % (mark[col[v]][1], mark[col[v]][0], x, y) for v, (x, y) in enumerate(P)]
out.append(r"\end{tikzpicture}")
with open(os.path.join(HERE, "q11-graph.tex"), "w") as f:
    f.write("\n".join(out) + "\n")
print(f"wrote q11-graph.tex: {len(P)} vertices, {len(E)} edges")
