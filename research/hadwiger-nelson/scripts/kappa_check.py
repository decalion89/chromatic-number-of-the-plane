"""kappa = omega * rho7 = (-11 + 5 sqrt(-3)) / 14 is congruent to 1 mod 5 in Z[omega, 1/7].
Check on the module: for every admissible psi and every unit u with kappa^k u also a unit
of the module (k = +-1, and also rho7, omega alone), compare psi(kappa u) with psi(u)."""
import sys, json, itertools
from fractions import Fraction as Fr
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "gate.py")).read().split("def gate(g, label):")[0])
from hn.geometry import Rotation
ROOT = "/home/user/darwin-50/research/hadwiger-nelson"
d = json.load(open(f"{ROOT}/data/{sys.argv[1]}"))
F = Field(tuple(d["field_generators"]))
P = [Point(F.element([Fr(a, b) for a, b in x]), F.element([Fr(a, b) for a, b in y])) for x, y in d["points"]]
g = build_graph(P)
units = {}
for a, b in g.edges():
    for u in (P[b] - P[a], P[a] - P[b]): units[u] = True
units = list(units)
E = edge_vectors(g); B = echelon(E); r = len(B)
den = None
# coordinates of every unit in the echelon basis (edge_vectors scales by a common denominator)
from math import gcd
out = [tuple(u.x.c) + tuple(u.y.c) for u in units]
D = 1
for v in out:
    for q in v: D = D * Fr(q).denominator // gcd(D, Fr(q).denominator)
vec = lambda p: tuple(int(Fr(q) * D) for q in tuple(p.x.c) + tuple(p.y.c))
def inmod(p):
    v = vec(p) if all((Fr(q) * D).denominator == 1 for q in tuple(p.x.c) + tuple(p.y.c)) else None
    if v is None: return None
    try: return coords(B, v)
    except Exception: return None
Cu = [inmod(u) for u in units]
assert all(c is not None for c in Cu)
# admissible psi mod 5
adm = []
for psi in itertools.product(range(5), repeat=r):
    if all(sum(a * b for a, b in zip(psi, c)) % 5 for c in Cu): adm.append(psi)
print(f"{sys.argv[1]}: {len(units)} units, rank {r}, {len(adm)} admissible psi")
s3 = F.element([Fr(0)] * 0 + [Fr(0), Fr(1)] + [Fr(0)] * (len(F.element([Fr(1)]).c) - 2)) if False else None
# sqrt3 as a field element: basis order follows Field; find it by squaring candidates
one = F.one()
basis = []
nb = len(one.c)
for i in range(nb):
    e = F.element([Fr(1) if j == i else Fr(0) for j in range(nb)])
    basis.append(e)
sq3 = next(e for e in basis if e * e == F.rational(3))
rots = {"omega": Rotation(F.rational(Fr(1, 2)), sq3 * F.rational(Fr(1, 2))),
        "rho7": Rotation(F.rational(Fr(1, 7)), sq3 * F.rational(Fr(4, 7))),
        "kappa": Rotation(F.rational(Fr(-11, 14)), sq3 * F.rational(Fr(5, 14)))}
uset = set(units)
for name, R in rots.items():
    for R2, nm in ((R, name), (R.inverse(), name + "^-1")):
        pairs = [(i, units.index(R2(u))) for i, u in enumerate(units) if R2(u) in uset]
        agree = sum(1 for psi in adm for i, j in pairs
                    if sum(a * b for a, b in zip(psi, Cu[i])) % 5 == sum(a * b for a, b in zip(psi, Cu[j])) % 5)
        # scalar test: is psi(R u) = k psi(u) for one k per psi?
        scal = 0
        for psi in adm:
            ks = {(sum(a * b for a, b in zip(psi, Cu[j])) * pow(sum(a * b for a, b in zip(psi, Cu[i])), -1, 5)) % 5 for i, j in pairs}
            if len(ks) == 1: scal += 1
        print(f"  {nm:9s}: {len(pairs)} units u with R u a unit; psi(Ru) = psi(u) in {agree}/{len(adm)*len(pairs)}; "
              f"psi o R a scalar multiple of psi for {scal}/{len(adm)} psi")
