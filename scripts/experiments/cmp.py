import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, cmath, math, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from hn.field import Field
from hn.geometry import Point, Rotation

FLD = Field((3, 11))
ONE, ZERO = FLD.rational(1), FLD.zero()
ROT60 = Rotation(FLD.rational(Fraction(1, 2)), FLD.sqrt(3) * FLD.rational(Fraction(1, 2)))
CH = Rotation(FLD.sqrt(33) * FLD.rational(Fraction(1, 6)),
              FLD.sqrt(3) * FLD.rational(Fraction(1, 6)))

seeds, q = [], Point(ONE, ZERO)
for _ in range(3):
    seeds.append(q)
    q = CH(q)
hexes = []
for s in seeds:
    h, p = [], s
    for _ in range(6):
        h.append(p)
        p = ROT60(p)
    hexes.append(h)

aux = {}
for a in range(3):
    for b in range(a, 3):
        for i, u in enumerate(hexes[a]):
            for j, v in enumerate(hexes[b]):
                if (a, i) >= (b, j):
                    continue
                s = Point(u.x + v.x, u.y + v.y)
                if s.norm2() in (ONE, ZERO):
                    continue
                aux.setdefault((s.x, s.y), set()).add((a, i))
                aux[(s.x, s.y)].add((b, j))
pts = list(aux)
E = sum(1 for i in range(len(pts)) for j in range(i + 1, len(pts))
        if Point(*pts[i]).is_unit_apart(Point(*pts[j])))
print(f"EXACT: {len(pts)} auxiliaries, {E} edges")

W = cmath.exp(1j * math.pi / 3)
HEXF = [W ** k for k in range(6)]
th = math.acos(5 / 6)
fh = [[cmath.exp(1j * a) * h for h in HEXF] for a in (0.0, th / 2, th)]
fa = {}
for a in range(3):
    for b in range(a, 3):
        for i, u in enumerate(fh[a]):
            for j, v in enumerate(fh[b]):
                if (a, i) >= (b, j):
                    continue
                s = u + v
                if abs(abs(s) - 1) < 1e-9 or abs(s) < 1e-9:
                    continue
                fa.setdefault((round(s.real, 6), round(s.imag, 6)), set()).add((a, i))
fp = [complex(*k) for k in fa]
for tol in (1e-9, 1e-7, 1e-6, 1e-5):
    e = sum(1 for i in range(len(fp)) for j in range(i + 1, len(fp))
            if abs(abs(fp[i] - fp[j]) - 1) < tol)
    print(f"FLOAT tol {tol:g}: {len(fp)} auxiliaries, {e} edges")
