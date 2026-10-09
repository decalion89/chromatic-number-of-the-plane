"""exactcheck.py -- independent exact check of a unit-distance graph over Q(sqrt R) and of its 3-colouring CNF.
usage: python3 exactcheck.py R STATE.json PTS.npy CNF"""
import sys, json, numpy as np
from math import sqrt
R = int(sys.argv[1]); st = json.load(open(sys.argv[2])); P = np.load(sys.argv[3]).tolist(); cnf = sys.argv[4]
D, U = st["D"], st["U"]
def mul(x, y):
    a, b, c, d = x; e, f, g, h = y
    return (a*e + R*b*f - c*g - R*d*h, a*f + b*e - c*h - d*g, a*g + c*e + R*b*h + R*d*f, a*h + d*e + b*g + c*f)
def conj(x): return (x[0], x[1], -x[2], -x[3])
for u in U:
    assert mul(tuple(u), conj(tuple(u))) == (D * D, 0, 0, 0), u
pts = [tuple(p) for p in P]; assert len(set(pts)) == len(pts)
idx = {p: i for i, p in enumerate(pts)}
E = set()
for i, p in enumerate(pts):
    for u in U:
        j = idx.get(tuple(a + b for a, b in zip(p, u)))
        if j is not None and i < j: E.add((i, j))
Ec = set(); nv = None; units = []
for line in open(cnf):
    if line.startswith("p"): nv = int(line.split()[2]); continue
    t = [int(x) for x in line.split()[:-1]]
    if len(t) == 2 and t[0] < 0 and t[1] < 0:
        a, b = (-t[0] - 1) // 3, (-t[1] - 1) // 3
        assert (-t[0] - 1) % 3 == (-t[1] - 1) % 3
        Ec.add((min(a, b), max(a, b)))
    elif len(t) == 1: units.append(t[0])
alo = sum(1 for line in open(cnf) if not line.startswith("p") and len(line.split()) == 4 and all(int(x) > 0 for x in line.split()[:3]))
s = sqrt(R)
z = lambda p: complex((p[0] + p[1] * s) / D, (p[2] + p[3] * s) / D)
md = max(abs(abs(z(pts[j]) - z(pts[i])) - 1) for i, j in E)
ok = Ec == E and nv == 3 * len(pts) and alo == len(pts)
print(f"Q(sqrt{R}): {len(U)} directions with u*conj(u) = D^2 (D = {D}); {len(pts)} distinct points; {len(E)} unit edges; "
      f"CNF: {nv} vars, {alo} at-least-one clauses, edge clauses match: {Ec == E}; units {units}; "
      f"max float error {md:.1e}; {'OK' if ok else 'MISMATCH'}")
