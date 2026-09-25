"""Unit sets over L = Q(sqrt3, sqrt5) for searches and for split_gate_q35.py.

A unit set is written as the zeta-orbit (zeta = e^{i pi/6}), closed under complex conjugation, of all
products g_1^{e_1} ... g_k^{e_k} with |e_j| <= t_j of named generators:
  tau  = (2 + i sqrt5)/3        not a unit above 3            (chain35.json)
  phi2 = (1 + i sqrt15)/4       not a unit above 2, valuation +-2 there; phi2 + conj(phi2) = 1/2
  dg   = (7 + i sqrt15)/8       de Grey's distance-2 spindle rotation, valuation +-4 above 2
  o1   = ((sqrt5 - sqrt3) - i (sqrt3 + sqrt5))/4   odd valuation (+-1) above 2
  g5   = (3 + 4i)/5,  rho7 = (1 + 4 i sqrt3)/7,  t7 = (2 + 3 i sqrt5)/7
  c1..c4                        c1 + c2 + c3 + c4 = 1/sqrt3: a chain of four unit rhombi (13 vertices)
The json has "points" = {0} u U and "units" = U, the format of grow_lean.py and module_gate.py.

usage: units_q35.py tau:1,phi2:1,o1:1 out.json
"""
import os, sys, json
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from hn.field import Field
from hn.geometry import Point

F = Field((3, 5))
r = F.rational
s3, s5 = F.sqrt(3), F.sqrt(5)
s15 = s3 * s5


def mul(a, b):
    return Point(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x)


def conj(a):
    return Point(a.x, -a.y)


ONE = Point(r(1), r(0))
ZETA = Point(s3 * r(Fr(1, 2)), r(Fr(1, 2)))
GEN = {
    "tau": Point(r(Fr(2, 3)), s5 * r(Fr(1, 3))),
    "phi2": Point(r(Fr(1, 4)), s15 * r(Fr(1, 4))),
    "dg": Point(r(Fr(7, 8)), s15 * r(Fr(1, 8))),
    "o1": Point(s5 * r(Fr(1, 4)) - s3 * r(Fr(1, 4)), s3 * r(Fr(-1, 4)) - s5 * r(Fr(1, 4))),
    "g5": Point(r(Fr(3, 5)), r(Fr(4, 5))),
    "rho7": Point(r(Fr(1, 7)), s3 * r(Fr(4, 7))),
    "t7": Point(r(Fr(2, 7)), s5 * r(Fr(3, 7))),
    "c1": Point(s3 * r(Fr(-1, 6)) - s15 * r(Fr(1, 6)), s3 * r(Fr(-1, 6)) + s15 * r(Fr(1, 6))),
    "c2": Point(s3 * r(Fr(-1, 6)) + s15 * r(Fr(1, 6)), s3 * r(Fr(1, 6)) + s15 * r(Fr(1, 6))),
    "c3": Point(s3 * r(Fr(1, 3)) - s5 * r(Fr(1, 6)), r(Fr(-1, 3)) - s15 * r(Fr(1, 6))),
    "c4": Point(s3 * r(Fr(1, 3)) + s5 * r(Fr(1, 6)), r(Fr(1, 3)) - s15 * r(Fr(1, 6))),
}


def units(spec):
    """spec: list of (name, t)"""
    base = [ONE]
    for name, t in spec:
        g, gi = GEN[name], conj(GEN[name])
        new = []
        for b in base:
            p, q = b, b
            new.append(b)
            for _ in range(t):
                p, q = mul(p, g), mul(q, gi)
                new += [p, q]
        base = new
    out = set()
    for b in base:
        p = b
        for _ in range(12):
            out.add(p)
            out.add(conj(p))
            p = mul(p, ZETA)
    return sorted(out, key=lambda p: (float(p.x), float(p.y)))


def check():
    for g in GEN.values():
        assert g.x * g.x + g.y * g.y == F.one()
    total = Point(sum((GEN[c].x for c in ("c1", "c2", "c3", "c4")), F.zero()),
                  sum((GEN[c].y for c in ("c1", "c2", "c3", "c4")), F.zero()))
    assert total == Point(s3 * r(Fr(1, 3)), F.zero())          # 1/sqrt3


if __name__ == "__main__":
    check()
    spec = [(s.split(":")[0], int(s.split(":")[1])) for s in sys.argv[1].split(",")]
    U = units(spec)
    enc = lambda e: [[c.numerator, c.denominator] for c in e.c]
    O = Point(F.zero(), F.zero())
    json.dump({"field_generators": [3, 5],
               "points": [[enc(O.x), enc(O.y)]] + [[enc(u.x), enc(u.y)] for u in U],
               "units": [[enc(u.x), enc(u.y)] for u in U],
               "note": f"units_q35.py {sys.argv[1]}"}, open(sys.argv[2], "w"))
    print(len(U), "units ->", sys.argv[2])
