"""Second, independent check of the finite witnesses over Q(sqrt d) (written from the file format only, sharing no
code with check_witness.py; d = 11 unless the witness file says).  Compares the points of the witness with the vertex
set A of data/quadratic_planes/q<d>.json (witness_q11sum: H = A + A, rebuilt here; the grown and minimised witnesses
keep part of A), recomputes every unit-distance pair from the finite list of unit vectors with the denominator D of the witness
(all integer solutions of a^2 + d b^2 + c^2 + d e^2 = D^2, a b + c e = 0), compares with the witness, checks the
(7,2)-colouring and the cycles, and writes a CNF with its own variable numbering:
  x(v,k) = k*n + v + 1  (vertex v has colour k)            -- different layout from the share
  y(a,b) = 7n + index   ("arc a->b is NOT tight")
  clauses: each vertex some colour; no edge with colour difference 0 or +-1 mod 7;
           y(a,b) -> not (x(a,k) and x(b,k+2)) for every k;  each listed cycle: OR of y over its arcs;
           x(v0,0) if the witness fixes a vertex v0 ("fixed_vertex": rotating the colours by -c(v0) keeps a
           (7,2)-colouring and all its colour differences, so this changes nothing).
A model gives a (7,2)-colouring in which every listed cycle has a non-tight arc; conversely such a colouring gives a
model.  So UNSAT <=> every (7,2)-colouring has a tight listed cycle.

usage: python3 verify_independent.py [OUT.cnf [WITNESS.json.gz]]     (default witness: witness_q11sum.json.gz)
(then, for instance, kissat OUT.cnf OUT.drat; drat-trim OUT.cnf OUT.drat -L OUT.lrat; cake_lpr OUT.cnf OUT.lrat)"""
import json, gzip, math, os, sys


def _req(ok, *msg):
    """An explicit check (not assert, so that python -O cannot skip it)."""
    if not ok:
        print('REJECTED:', *msg, file=sys.stderr)
        sys.exit(1)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", ".."))
OUT = sys.argv[1] if len(sys.argv) > 1 else None
WIT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "witness_q11sum.json.gz")
W = json.load(gzip.open(WIT, 'rt'))
d = W.get('d', 11); D = W['denominator']
q = json.load(open(os.path.join(REPO, 'data', 'quadratic_planes', f'q{d}.json')))
_req(q['D'] == D and q['d'] == d, 'the witness and q<d>.json differ in d or D')
print(f'field Q(sqrt{d}), denominator {D}')
A = [tuple(p) for p in q['points']]
H = sorted({tuple(x + y for x, y in zip(p, r)) for p in A for r in A})
print('A + A:', len(H), 'points')
P = [tuple(p) for p in W['points']]
_req(all(len(p) == 4 and all(isinstance(t, int) for t in p) for p in P), 'bad points')
_req(len(set(P)) == len(P), 'repeated point')
print('witness:', os.path.basename(WIT), '-', len(P), 'points')
sumset = set(P) == set(H)
print('witness points = A + A:', 'yes' if sumset else 'no')
print(f'points of A in the witness: {len(set(A) & set(P))} of {len(set(A))}')
if os.path.basename(WIT) == 'witness_q11sum.json.gz':
    _req(sumset, 'point set differs from A + A')
# unit vectors with denominator 30: a^2 + 11 b^2 + c^2 + 11 e^2 = 900, a b + c e = 0
U = []
B = math.isqrt(D * D // d)
for b in range(-B, B + 1):
    for e in range(-B, B + 1):
        r = D * D - d * (b * b + e * e)
        if r < 0:
            continue
        for a in range(-D, D + 1):
            c2 = r - a * a
            if c2 < 0:
                continue
            c = math.isqrt(c2)
            for cc in {c, -c}:
                if cc * cc == c2 and a * b + cc * e == 0:
                    U.append((a, b, cc, e))
U = sorted(set(U))
print(f'unit vectors with denominator {D}:', len(U))
idx = {p: i for i, p in enumerate(P)}
E = set()
for i, p in enumerate(P):
    for u in U:
        j = idx.get(tuple(x + y for x, y in zip(p, u)))
        if j is not None and j != i:
            E.add((min(i, j), max(i, j)))
WE = {(min(i, j), max(i, j)) for i, j in W['edges']}
_req(len(WE) == len(W['edges']), 'repeated edge')
print('unit pairs:', len(E), '; witness edges:', len(WE), '; equal:', E == WE)
_req(E == WE, 'the witness edges are not the unit pairs')
col = W['colouring']
_req(len(col) == len(P) and all(isinstance(c, int) and 0 <= c < 7 for c in col), 'bad colouring')
_req(all((col[j] - col[i]) % 7 in (2, 3, 4, 5) for i, j in E), 'not a (7,2)-colouring')
print('(7,2)-colouring: yes')
cyc = W['cycles']
for C in cyc:
    _req(len(C) == len(set(C)) >= 3, 'a cycle is not simple', C)
    for k in range(len(C)):
        a, b = C[k], C[(k + 1) % len(C)]
        _req((min(a, b), max(a, b)) in E, 'cycle uses a non-edge', C)
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
if 'fixed_vertex' in W:
    v0 = W['fixed_vertex']
    _req(isinstance(v0, int) and 0 <= v0 < n, 'bad fixed_vertex')
    cl.append([x(v0, 0)])
    print('fixed vertex:', v0, '(colour 0)')
nv = 7 * n + len(arcs)
print('CNF:', nv, 'variables,', len(cl), 'clauses')
if OUT:
    with open(OUT, 'w') as f:
        f.write(f'p cnf {nv} {len(cl)}\n')
        for c in cl:
            f.write(' '.join(map(str, c)) + ' 0\n')
    print('written to', OUT)
