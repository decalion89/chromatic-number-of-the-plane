"""Structure of a colouring of the blocked module: split the points into sheets (cosets of the
803 module M_1, which has index 5^8 in the blocked module) and measure, sheet by sheet, the best
coset fit against M_1's 1728 admissible psi.  'Coset per sheet' would be a colouring family of the
blocked module that no quotient gate could see."""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, json, itertools
from fractions import Fraction as Fr
from math import gcd
from collections import defaultdict
import numpy as np
import sympy
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
exec(open(os.path.join(HN_DIR, "scripts", "gate.py")).read().split("def gate(g, label):")[0])
d = json.load(open(sys.argv[1])); col = d["colouring"]
F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(a, b) for a, b in xy[0]]), F.element([Fr(a, b) for a, b in xy[1]]))
V = [mk(xy) for xy in d["points"]]
d1 = json.load(open(HN_DIR + "/data/five_247_c.json"))
P1 = [mk(xy) for xy in d1["points"]]
g1 = build_graph(P1)
U1 = [tuple((P1[b] - P1[a]).x.c) + tuple((P1[b] - P1[a]).y.c) for a, b in g1.edges()]
raw = [tuple((p - V[0]).x.c) + tuple((p - V[0]).y.c) for p in V]
# rational basis of M_1 (echelon over integers after a common scaling)
den = 1
for v in U1:
    for q in v: den = den * Fr(q).denominator // gcd(den, Fr(q).denominator)
I1 = sorted({max(w, tuple(-x for x in w)) for w in (tuple(int(Fr(q) * den) for q in v) for v in U1)})
B1 = echelon(I1); r = len(B1)
Bm = sympy.Matrix(B1).T                      # 16 x 8, columns = basis of M_1 (scaled by den)
pinv = (Bm.T * Bm).inv() * Bm.T              # exact left inverse on the span
def m1coords(v):
    x = pinv * sympy.Matrix([Fr(q) * den for q in v])
    return [Fr(int(sympy.fraction(t)[0]), int(sympy.fraction(t)[1])) for t in x]
co = [m1coords(v) for v in raw]
sheet = defaultdict(list)
for i, c in enumerate(co):
    sheet[tuple(x - (x.numerator // x.denominator) for x in c)].append(i)
sizes = sorted((len(v) for v in sheet.values()), reverse=True)
print(f"{sys.argv[1]}: n={len(V)}; sheets (cosets of M_1) met: {len(sheet)}; largest sheets {sizes[:8]}", flush=True)
# admissible psi of M_1
C1 = [coords(B1, e) for e in I1]
adm = [psi for psi in itertools.product(range(5), repeat=r) if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in C1)]
perms = list(itertools.permutations(range(5)))
print(f"  M_1: rank {r}, {len(adm)} admissible psi", flush=True)
for key_, idx in sorted(sheet.items(), key=lambda kv: -len(kv[1]))[:6]:
    base = co[idx[0]]
    ints = np.array([[int(a - b) for a, b in zip(co[i], base)] for i in idx], dtype=object)
    I5 = np.array([[int(x) % 5 for x in row] for row in ints], dtype=np.int64)
    cv = np.array([col[i] for i in idx], dtype=np.int64)
    best = 0
    for psi in adm:
        pv = (I5 @ np.array(psi, dtype=np.int64)) % 5
        T = np.zeros((5, 5), dtype=np.int64); np.add.at(T, (cv, pv), 1)
        best = max(best, max(sum(T[c][p[c]] for c in range(5)) for p in perms))
    print(f"  sheet of {len(idx)} points: best coset fit {best}/{len(idx)} = {best/len(idx):.3f}", flush=True)
