"""Second, independent check of the finite witness over Q(sqrt11) (written from the file format only, sharing no code
with check_witness.py).  Rebuilds H = A + A from data/quadratic_planes/q11.json, recomputes every unit-distance pair
from the finite list of unit vectors with denominator 30, compares with the witness, checks the (7,2)-colouring and
the cycles, and writes a CNF with its own variable numbering:
  x(v,k) = k*n + v + 1  (vertex v has colour k)            -- different layout from the share
  y(a,b) = 7n + index   ("arc a->b is NOT tight")
  clauses: each vertex some colour; no edge with colour difference 0 or +-1 mod 7;
           y(a,b) -> not (x(a,k) and x(b,k+2)) for every k;  each listed cycle: OR of y over its arcs.
A model gives a (7,2)-colouring in which every listed cycle has a non-tight arc; conversely such a colouring gives a
model.  So UNSAT <=> every (7,2)-colouring has a tight listed cycle.

usage: python3 verify_independent.py [OUT.cnf]
(then, for instance, kissat OUT.cnf OUT.drat; drat-trim OUT.cnf OUT.drat -L OUT.lrat; cake_lpr OUT.cnf OUT.lrat)"""
import json, gzip, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
WIT = os.path.join(HERE, "witness_q11sum.json.gz")
OUT = sys.argv[1] if len(sys.argv) > 1 else None
q = json.load(open(os.path.join(REPO, 'data', 'quadratic_planes', 'q11.json')))
W = json.load(gzip.open(WIT, 'rt'))
D = 30
assert q['D'] == D and q['d'] == 11 and W['denominator'] == D
A = [tuple(p) for p in q['points']]
H = sorted({tuple(x + y for x, y in zip(p, r)) for p in A for r in A})
print('A + A:', len(H), 'points')
P = [tuple(p) for p in W['points']]
assert len(set(P)) == len(P)
assert set(P) == set(H), 'point set differs from A + A'
print('witness points = A + A: yes')
# unit vectors with denominator 30: a^2 + 11 b^2 + c^2 + 11 e^2 = 900, a b + c e = 0
U = []
for b in range(-9, 10):
    for e in range(-9, 10):
        r = D * D - 11 * (b * b + e * e)
        if r < 0:
            continue
        for a in range(-30, 31):
            c2 = r - a * a
            if c2 < 0:
                continue
            c = int(round(c2 ** 0.5))
            for cc in {c, -c}:
                if cc * cc == c2 and a * b + cc * e == 0:
                    U.append((a, b, cc, e))
U = sorted(set(U))
print('unit vectors with denominator 30:', len(U))
idx = {p: i for i, p in enumerate(P)}
E = set()
for i, p in enumerate(P):
    for u in U:
        j = idx.get(tuple(x + y for x, y in zip(p, u)))
        if j is not None and j != i:
            E.add((min(i, j), max(i, j)))
WE = {(min(i, j), max(i, j)) for i, j in W['edges']}
assert len(WE) == len(W['edges'])
print('unit pairs:', len(E), '; witness edges:', len(WE), '; equal:', E == WE)
assert E == WE
col = W['colouring']
assert len(col) == len(P) and all(isinstance(c, int) and 0 <= c < 7 for c in col)
assert all((col[j] - col[i]) % 7 in (2, 3, 4, 5) for i, j in E)
print('(7,2)-colouring: yes')
cyc = W['cycles']
for C in cyc:
    assert len(C) == len(set(C)) >= 3
    for k in range(len(C)):
        a, b = C[k], C[(k + 1) % len(C)]
        assert (min(a, b), max(a, b)) in E
print('cycles:', len(cyc), 'simple cycles along edges; lengths', sorted({len(C) for C in cyc}))
# tight cycles of the given colouring: each listed cycle's tightness under col (for information)
tight = sum(all((col[C[(k + 1) % len(C)]] - col[C[k]]) % 7 == 2 for k in range(len(C))) for C in cyc)
print('listed cycles tight under the given colouring:', tight)
n = len(P)
x = lambda v, k: k * n + v + 1
arcs = {}
for C in cyc:
    for k in range(len(C)):
        a, b = C[k], C[(k + 1) % len(C)]
        if (a, b) not in arcs:
            arcs[(a, b)] = 7 * n + len(arcs) + 1
cl = []
for v in range(n):
    cl.append([x(v, k) for k in range(7)])
for i, j in sorted(E):
    for k in range(7):
        for dlt in (0, 1, -1):
            cl.append([-x(i, k), -x(j, (k + dlt) % 7)])
for (a, b), y in arcs.items():
    for k in range(7):
        cl.append([-y, -x(a, k), -x(b, (k + 2) % 7)])
for C in cyc:
    cl.append([arcs[(C[k], C[(k + 1) % len(C)])] for k in range(len(C))])
nv = 7 * n + len(arcs)
print('CNF:', nv, 'variables,', len(cl), 'clauses')
if OUT:
    with open(OUT, 'w') as f:
        f.write(f'p cnf {nv} {len(cl)}\n')
        for c in cl:
            f.write(' '.join(map(str, c)) + ' 0\n')
    print('written to', OUT)
