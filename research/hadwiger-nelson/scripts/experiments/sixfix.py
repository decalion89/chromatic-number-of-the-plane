"""Five-colour forcing in de Grey's G, with the closable set stated correctly.

The first pass tested squarefree(1 - 4D) against products of {-3,-7,-11,-15},
which is only half the square classes his field has.  His coordinates are real
and lie in F = Q(sqrt3, sqrt5, sqrt7, sqrt11); a rotation by theta is the pair
(cos, sin) with cos = 1 - 1/(2D) rational, so the only question is whether

    sin = sqrt(4D - 1) / (2D)

lands in F, i.e. whether the squarefree part of 4D - 1 divides 1155 = 3.5.7.11.
That admits D = 2, 9, 14, 34, 44, ... which the first pass skipped, and it
needs D >= 1/4, since 2d sin(theta/2) = 1 is unachievable below.

Checked against the chain: D = 1, 3, 4, 16 give 4D - 1 = 3, 11, 15, 63 = 9.7,
so the radicands 3, 11, 15, 7 -- de Grey's field, in the order he uses them.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import isqrt
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD, Rotation
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD
RAD = {1}
for p in (3, 5, 7, 11):
    RAD |= {r * p for r in RAD}


def squarefree(m):
    s, d = (-1 if m < 0 else 1), abs(m)
    f = 2
    while f * f <= d:
        e = 0
        while d % f == 0:
            d //= f
            e += 1
        if e % 2:
            s *= f
        f += 1
    return s * d


def closable(D):
    D = Fr(D)
    if D <= Fr(1, 4):
        return False
    r = 4 * D - 1
    return squarefree(r.numerator * r.denominator) in RAD


def rot_closing(D):
    """cos = 1 - 1/(2D), sin = sqrt(4D-1)/(2D), both in F."""
    D = Fr(D)
    c = 1 - Fr(1, 2) / D
    r = 4 * D - 1
    sq = squarefree(r.numerator * r.denominator)
    if sq not in RAD:
        return None
    k2 = (1 - c * c) / sq
    num, den = k2.numerator, k2.denominator
    t = isqrt(num * den)
    if t * t != num * den:
        return None
    rot = Rotation(F.rational(c), F.sqrt(sq) * F.rational(Fr(t, den)))
    assert rot.cos * rot.cos + rot.sin * rot.sin == F.rational(1)
    return rot


print(f"chain check: "
      f"{[(D, 4 * D - 1, squarefree(4 * D - 1)) for D in (1, 3, 4, 16)]}",
      flush=True)
print(f"integer D up to 100 that de Grey's field closes: "
      f"{[D for D in range(1, 101) if closable(D)]}", flush=True)

G = build_G(F, as_graph=False)
n = len(G)
zs = [(float(p.x), float(p.y)) for p in G]
cell = {}
for i, (a, b) in enumerate(zs):
    cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
E = []
for i, (a, b) in enumerate(zs):
    cx, cy = int(a // 1), int(b // 1)
    for da in (-1, 0, 1):
        for db in (-1, 0, 1):
            for j in cell.get((cx + da, cy + db), ()):
                if j <= i:
                    continue
                if abs((a - zs[j][0]) ** 2 + (b - zs[j][1]) ** 2 - 1) > 1e-7:
                    continue
                d = G[i] - G[j]
                if d.x * d.x + d.y * d.y == F.rational(1):
                    E.append((i, j))
print(f"G: {n} vertices, {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

KC = 5
cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
for a, b in E:
    for c in range(KC):
        cls.append([-(1 + a * KC + c), -(1 + b * KC + c)])
sv = Solver(name="cd19", bootstrap_with=cls)
print(f"  5-colourable? {sv.solve()}  [{time.time()-t0:.0f}s]", flush=True)

cand, occ = [], Counter()
for i in range(n):
    ai, bi = zs[i]
    for j in range(i + 1, n):
        v = (ai - zs[j][0]) ** 2 + (bi - zs[j][1]) ** 2
        if v > 36.0:
            continue
        D = Fr(round(v * 1584), 1584)
        if abs(float(D) - v) > 1e-7 or D == 1 or not closable(D):
            continue
        cand.append((i, j, D))
        occ[D] += 1
print(f"  {len(cand)} pairs at a closable distance, "
      f"{[(str(D), c) for D, c in occ.most_common(10)]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

hits, hard = [], []
for k, (i, j, D) in enumerate(cand):
    sv.conf_budget(40000)
    r = sv.solve_limited(assumptions=[1 + i * KC, -(1 + j * KC)])
    if r is False:
        hits.append((i, j, D))
        print(f"  *** FORCED AT FIVE COLOURS: {i},{j} D = {D}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
    elif r is None:
        hard.append((i, j, D))
    if k and k % 5000 == 0:
        print(f"  ... {k}/{len(cand)}, {len(hard)} hard, {len(hits)} forced  "
              f"[{time.time()-t0:.0f}s]", flush=True)
for i, j, D in hard:
    if not sv.solve(assumptions=[1 + i * KC, -(1 + j * KC)]):
        hits.append((i, j, D))
        print(f"  *** FORCED (full run): {i},{j} D = {D}", flush=True)
print(f"  {len(hits)} forced pairs at five colours in G "
      f"({len(hard)} needed a full run)  [{time.time()-t0:.0f}s]", flush=True)
