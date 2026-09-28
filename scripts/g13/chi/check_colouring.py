"""Independent check of a colouring of G_13 (does not import g13.py).

Input: a file with 169 integers (colours of the points (x, y), index 13 x + y), or a kissat model ('v ...' lines)
of a formula whose first 845 variables are x(v, c) = 1 + 5 v + c.  Every pair of points is tested:
z ~ w iff (x - x')^2 - 2 (y - y')^2 = 1 mod 13.  Prints the number of colours, the number of edges checked
(must be 1183) and the monochromatic edges (must be none).
usage: check_colouring.py FILE
"""
import sys

text = open(sys.argv[1]).read()
if any(line.startswith("v ") for line in text.splitlines()):
    lits = [int(t) for line in text.splitlines() if line.startswith("v ") for t in line.split()[1:]]
    true = {l for l in lits if l > 0}
    col = []
    for v in range(169):
        cs = [c for c in range(5) if 1 + 5 * v + c in true]
        assert len(cs) == 1, (v, cs)
        col.append(cs[0])
else:
    col = [int(t) for t in text.split()]
assert len(col) == 169
edges = bad = 0
for v in range(169):
    for w in range(v + 1, 169):
        dx, dy = (v // 13 - w // 13) % 13, (v % 13 - w % 13) % 13
        if (dx * dx - 2 * dy * dy) % 13 == 1:
            edges += 1
            if col[v] == col[w]:
                bad += 1
                print("monochromatic edge", (v // 13, v % 13), (w // 13, w % 13))
print(f"{len(set(col))} colours, {edges} edges checked, {bad} monochromatic")
sys.exit(0 if edges == 1183 and bad == 0 else 1)
