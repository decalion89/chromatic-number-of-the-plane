"""All the unit vectors of the lambda-closure M' = M1 + lambda M1 that come for free:
U1 (units of M1, closed under the 60-degree turn), lambda U1, and the Exoo-Ismailescu units
w_e = 5(1 - lambda) e (|w_e| = 5 |1 - lambda| = 1, and w_e = 5e - 5 lambda e lies in M').
Writes a module file (units + one dummy point) for circgate.py / circextend.py.

usage: python3 lamclosure_units.py <M1.json> <extra_units.json or -> <out.json>"""
import sys, json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
exec(open(os.path.join(HERE, "gate.py")).read().split("def gate(g, label):")[0])
d = json.load(open(sys.argv[1])); F = Field(tuple(d["field_generators"]))
mk = lambda xy: Point(F.element([Fr(p, q) for p, q in xy[0]]), F.element([Fr(p, q) for p, q in xy[1]]))
P = [mk(xy) for xy in d["points"]]; g = build_graph(P); Ud = {}
key = lambda p: (round(p.fx, 7), round(p.fy, 7))
for i, j in g.edges():
    for w in (P[j] - P[i], P[i] - P[j]): Ud.setdefault(key(w), w)
def el(rad):
    c = [Fr(0)] * F.dim
    for m in range(F.dim):
        if F._prod[m] in rad: c[m] = Fr(rad[F._prod[m]])
    return F.element(c)
def rot(cos_c, sin_rad):
    cth = F.rational(cos_c); sth = el(sin_rad)
    return lambda p: Point(p.x * cth - p.y * sth, p.x * sth + p.y * cth)
om = rot(Fr(1, 2), {3: Fr(1, 2)}); lam = rot(Fr(49, 50), {11: Fr(3, 50)})
def close60(D):
    for u in list(D.values()):
        w = u
        for _ in range(5):
            w = om(w); D.setdefault(key(w), w)
    return D
U1 = close60(dict(Ud))
U2 = {key(lam(u)): lam(u) for u in U1.values()}
U3 = {}
for u in U1.values():
    w = (u - lam(u)) * 5 if hasattr(u, "__mul__") else None
    w = Point((u.x - lam(u).x) * 5, (u.y - lam(u).y) * 5)
    assert abs(w.fx ** 2 + w.fy ** 2 - 1) < 1e-9
    U3[key(w)] = w
allU = dict(U1); allU.update(U2); allU.update(U3)
if sys.argv[2] != "-":
    ex = json.load(open(sys.argv[2]))
    for xy in ex["units"]: u = mk(xy); allU.setdefault(key(u), u)
allU = close60(allU)
print(f"U1 {len(U1)}, lambda U1 {len(U2)}, E-I units w_e {len(U3)} (new beyond U1, lambda U1: {len(set(U3) - set(U1) - set(U2))}); total {len(allU)}")
enc = lambda u: [[[x.numerator, x.denominator] for x in u.x.c], [[x.numerator, x.denominator] for x in u.y.c]]
json.dump({"field_generators": list(F.gens), "units": [enc(u) for u in allU.values()],
           "points": [enc(Point(F.rational(0), F.rational(0)))], "A": 0, "B": 0,
           "note": "lambda-closure of " + os.path.basename(sys.argv[1]) + ": U1, lambda U1, 5(1-lambda) U1"}, open(sys.argv[3], "w"))
