"""Exact check of a circular colouring c(x) = floor(5 * frac(phi(x))) on a module and on a grown graph.

phi is given by its values on the LLL basis of circgate.py (rationals, e.g. 6-decimal floats).
Checks, in exact rational arithmetic:
  1. frac(phi(u)) in [1/5, 4/5] for every unit vector u (so the colouring is proper on the whole
     module graph along those directions), with the exact least slack;
  2. on a checkpoint's points: the colouring has no monochromatic unit step.
usage: python3 circverify.py <checkpoint.json> <phi_1> ... <phi_r>"""
import sys, json, os, time
from fractions import Fraction as Fr
from math import gcd, floor
import numpy as np
import sympy
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
t0 = time.time()
path = sys.argv[1]; phi = [Fr(x) for x in sys.argv[2:]]
d = json.load(open(path))
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(p, q) for p, q in xy[0]]), F.element([Fr(p, q) for p, q in xy[1]]))
U = [mk(xy) for xy in d["units"]]
cc = [Fr(0)] * F.dim; cc[1] = Fr(1); S3 = F.element(cc)
turn = lambda p: Point(p.x * Fr(1, 2) - p.y * S3 * Fr(1, 2), p.x * S3 * Fr(1, 2) + p.y * Fr(1, 2))
seen = {(round(u.fx, 7), round(u.fy, 7)) for u in U}; n0 = len(U)
for u in list(U):
    w = u
    for _ in range(5):
        w = turn(w); k_ = (round(w.fx, 7), round(w.fy, 7))
        if k_ not in seen: seen.add(k_); U.append(w)
Ev = [tuple(u.x.c) + tuple(u.y.c) for u in U]
V = [mk(xy) for xy in d["points"]]; A = d["A"]
raw = [tuple((p - V[A]).x.c) + tuple((p - V[A]).y.c) for p in V]
den = 1
for v in Ev + raw:
    for x in v: den = den * Fr(x).denominator // gcd(den, Fr(x).denominator)
# NOTE: same basis construction as circgate.py (echelon of the edge directions, LLL in the trace form)
E = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(x) * den) for x in v) for v in Ev[:n0])})
B = echelon(E); r = len(B)
src = open(os.path.join(HERE, "allunits.py")).read()
exec(src[src.index("def lll_exact"):src.index("Gr, Ur = lll_exact(G)")])
wt = list(F._prod) * 2
G0 = [[Fr(sum(w * p * q for w, p, q in zip(wt, B[i], B[j]))) for j in range(r)] for i in range(r)]
_, T = lll_exact(G0)
Ti = sympy.Matrix(T).inv()
def lcoords(v):
    c = coords(B, tuple(int(Fr(x) * den) for x in v))
    assert all(Fr(t).denominator == 1 for t in c)
    return [int(x) for x in (sympy.Matrix([[int(t) for t in c]]) * Ti)]
assert len(phi) == r, f"need {r} values of phi"
CU = [lcoords(v) for v in Ev]
fr = lambda q: q - floor(q)
vals = [fr(sum(Fr(c) * p for c, p in zip(cu, phi))) for cu in CU]
slack = min(min(v - Fr(1, 5), Fr(4, 5) - v) for v in vals)
print(f"{os.path.basename(path)}: {len(U)} units, rank {r}; exact least slack of frac(phi(u)) in [1/5, 4/5]: {slack} = {float(slack):.6f}")
ok_units = slack >= 0
col = [floor(5 * fr(sum(Fr(c) * p for c, p in zip(lcoords(v), phi)))) for v in raw]
idx = {tuple(lcoords(v)): i for i, v in enumerate(raw)} if False else None
# edges of the checkpoint graph: pairs differing by a unit (hash the exact coordinates)
pc = [tuple(lcoords(v)) for v in raw]
pos = {c: i for i, c in enumerate(pc)}
bad = edges = 0
for i, c in enumerate(pc):
    for cu in CU:
        j = pos.get(tuple(a + b for a, b in zip(c, cu)))
        if j is not None and j > i:
            edges += 1; bad += col[i] == col[j]
print(f"  units all proper: {ok_units};  checkpoint graph: {len(pc)} points, {edges} unit edges, monochromatic: {bad}   [{time.time()-t0:.0f}s]")
print(f"  colour class sizes: {[col.count(k) for k in range(5)]}")
