"""Grow a rotationally symmetric set greedily, instead of guessing a seed.

S is 39 points with eighteen unit edges and SEVENTEEN isolated vertices, and
25 of its 39 points sit at irrational radius from the origin.  So the seed has
almost no structure of its own -- everything in Sa comes from how the twelve
rotated copies interlock.  That is the design criterion, and it can be built
towards rather than guessed at.

Grow: keep a D6-symmetric point set; the candidates are the points at distance
exactly one from TWO points already present, which are where unit-circle pairs
meet; add the whole orbit of whichever candidate creates the most new edges;
repeat.  Each step is the locally densest rotationally symmetric move.

The intersection of the unit circles about A and B at squared distance D is
(A+B)/2 +- sqrt(4-D)/2 * n, with n the unit normal to B-A.  Both sqrt(D) and
sqrt(4-D) must lie in the field for that point to be constructible there --
the same closability condition the spindle needs, arriving from the other
side.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from collections import Counter
from hn.field import Field
from hn.geometry import Point, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

t0 = time.time()
# sqrt(D) AND sqrt(4-D) must both lie in the field for an
# intersection to be constructible, and de Grey's four radicals
# allow almost no pairs -- the growth stalled after three steps
# with six candidates.  Adding sqrt(2) and sqrt(13) widens it.
K = Field((2, 3, 5, 7, 11, 13))
ZERO = Point(K.zero(), K.zero())
W = _rot60(K)


SQUARE_CLASSES = []
for _bits in range(1 << 6):
    _v = 1
    for _i, _p in enumerate((2, 3, 5, 7, 11, 13)):
        if _bits >> _i & 1:
            _v *= _p
    SQUARE_CLASSES.append(_v)
SQUARE_CLASSES.sort()


def sqrt_in_field(e):
    """A square root of the field element e, if the field has one."""
    if all(x == 0 for x in e.c[1:]):
        v = e.c[0]
        if v < 0:
            return None
        num, den = v.numerator, v.denominator
        for r in SQUARE_CLASSES:
            # v = q^2 * r  =>  sqrt(v) = q sqrt(r)
            q2 = Fr(num, den) / r
            if q2 <= 0:
                continue
            n2, d2 = q2.numerator, q2.denominator
            rn, rd = int(round(n2 ** 0.5)), int(round(d2 ** 0.5))
            if rn * rn == n2 and rd * rd == d2:
                s = K.rational(Fr(rn, rd))
                return s if r == 1 else K.sqrt(r) * s
    return None


def orbit(p):
    out, q = [], p
    for _ in range(6):
        out.append(q)
        out.append(Point(q.x, -q.y))
        q = W(q)
    return out


def candidates(P):
    """Points at distance one from two members of P, constructible here."""
    out = set()
    n = len(P)
    for i in range(n):
        for j in range(i + 1, n):
            A, B = P[i], P[j]
            D = A.dist2(B)
            if D == 0:
                continue
            f = float(D)
            if f > 4 or f < 0.05:
                continue
            sD = sqrt_in_field(D)
            s4 = sqrt_in_field(K.rational(4) - D)
            if sD is None or s4 is None:
                continue
            inv = K.rational(1) / sD
            half = K.rational(Fr(1, 2))
            mx = (A.x + B.x) * half
            my = (A.y + B.y) * half
            nx = -(B.y - A.y) * inv * s4 * half
            ny = (B.x - A.x) * inv * s4 * half
            out.add(Point(mx + nx, my + ny))
            out.add(Point(mx - nx, my - ny))
    return [p for p in out if p not in set(P) and float(p.norm2()) < 9]


def edges(P):
    g = build_graph(P)
    return set((min(a, b), max(a, b)) for a, b in g.edges())


def chrom(P, hi=5):
    E = sorted(edges(P))
    n = len(P)
    for k in range(2, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        sv = Solver(name="cd15", bootstrap_with=cls)
        ok = sv.solve()
        sv.delete()
        if ok:
            return k, len(E)
    return f">{hi}", len(E)


P = orbit(Point(K.one(), K.zero())) + [ZERO]
P = list(dict.fromkeys(P))
print(f"start: {len(P)} points (the unit orbit and the origin)"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for step in range(1, 26):
    cand = candidates(P)
    if not cand:
        print("   no constructible candidates left", flush=True)
        break
    cur = len(edges(P))
    best, bestgain = None, -1
    for c in cand[:250]:
        trial = list(dict.fromkeys(P + orbit(c)))
        gain = len(edges(trial)) - cur
        if gain > bestgain:
            best, bestgain = c, gain
    P = list(dict.fromkeys(P + orbit(best)))
    k, m = chrom(P)
    print(f"  step {step}: +{len(orbit(best))} -> {len(P)} pts, {m} edges, "
          f"{m/len(P):.2f} per vertex, chi = {k}  (gain {bestgain}, "
          f"{len(cand)} candidates)  [{time.time()-t0:.0f}s]", flush=True)
    if k == ">5":
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
