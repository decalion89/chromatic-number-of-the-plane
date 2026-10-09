"""certify_z24.py W.json OUT.cnf: as certify4.py (the formula 'proper 4-colouring, origin coloured 0, no listed cycle
tight'), for points of Q(zeta_24) (Q(sqrt2, sqrt3)^2); edges: all unit pairs (float prefilter, exact check)."""
import sys, json, itertools
import numpy as np
from z24 import is_unit, complex_of
W = json.load(open(sys.argv[1])); out = sys.argv[2]
P = [tuple(p) for p in W["points"]]; n = len(P); D = W["D"]; idx = {p: i for i, p in enumerate(P)}
Z = complex_of(P, D)
edges = set()
for i in range(n):
    d = np.abs(Z[i + 1:] - Z[i])
    for k in np.nonzero(np.abs(d - 1.0) < 1e-9)[0]:
        j = i + 1 + int(k)
        if is_unit([b - a for a, b in zip(P[i], P[j])], D):
            edges.add((i, j))
adj = [set() for _ in range(n)]
for i, j in edges:
    adj[i].add(j); adj[j].add(i)
print(f"{n} points, {len(edges)} unit pairs", flush=True)
v = lambda i, c: 4 * i + c + 1
cls = []
for i in range(n):
    cls.append([v(i, c) for c in range(4)])
    for c1, c2 in itertools.combinations(range(4), 2):
        cls.append([-v(i, c1), -v(i, c2)])
for i, j in sorted(edges):
    for c in range(4):
        cls.append([-v(i, c), -v(j, c)])
cls.append([v(idx[(0,) * 8], 0)])
ncyc = 0
for cyc in W["cycles"]:
    m = len(cyc)
    for a, b in zip(cyc, cyc[1:] + cyc[:1]):
        assert b in adj[a], "a listed cycle uses a non-edge"
    assert len(set(cyc)) == m
    if m % 4:
        continue
    for c in range(4):
        cls.append([-v(cyc[k], (c + k) % 4) for k in range(m)])
    ncyc += 1
with open(out, "w") as fh:
    fh.write(f"p cnf {4 * n} {len(cls)}\n")
    for c in cls:
        fh.write(" ".join(map(str, c)) + " 0\n")
print(f"wrote {out}: {4 * n} variables, {len(cls)} clauses ({ncyc} cycles of length 0 mod 4)", flush=True)
