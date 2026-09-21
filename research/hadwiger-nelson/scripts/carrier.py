"""Which carrier is the right one?  The question the cap tests never asked.

de Grey's H is 4-chromatic and capped at four colours.  Carrier and question
carry the same number.  That is not decoration: with k colours on a graph that
only needs k-1, one colour is free everywhere, and a ring can almost always be
shown all k by permuting inside the slack.  A cap needs the colouring to be
tight, and tightness is exactly chromaticity.

So every cap test at five colours has to be run on a carrier that actually
needs five.  Sa and Y need four.  Their universes -- Sa or Y plus every
unit-circle intersection -- might need more, since adding points can only
raise the chromatic number, and if they do the five-colour tests on them were
in the right regime after all.  If they do not, those tests were asking a
3-chromatic graph for a 4-cap, and their zeros say nothing.

One SAT call settles it per carrier.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_G, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

SEED = sys.argv[1] if len(sys.argv) > 1 else "Sa"
RADIUS = float(sys.argv[2]) if len(sys.argv) > 2 else 4.5
KS = [int(x) for x in (sys.argv[3] if len(sys.argv) > 3 else "4").split(",")]
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))


def rsqrt(e):
    if any(x for x in e.c[1:]):
        return None
    v = e.c[0]
    if v <= 0:
        return None
    for rr in CLASSES:
        q = v / rr
        nn, dd = q.numerator, q.denominator
        rn, rd = int(round(nn ** .5)), int(round(dd ** .5))
        if rn * rn == nn and rd * rd == dd:
            s = K.rational(Fr(rn, rd))
            return s if rr == 1 else K.sqrt(rr) * s
    return None


builders = {"Sa": build_Sa, "Y": build_Y,
            "G": lambda f: build_G(f, as_graph=False)}
P = list(builders[SEED](K))
have, fresh = set(P), []
n0 = len(P)
for i in range(n0):
    A = P[i]
    for j in range(i + 1, n0):
        B = P[j]
        D = A.dist2(B)
        if not .05 < float(D) < 3.99:
            continue
        sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
        if sD is None or s4 is None:
            continue
        inv = K.rational(1) / sD
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * inv * s4 * half
        ny = (B.x - A.x) * inv * s4 * half
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q in have or float(q.x) ** 2 + float(q.y) ** 2 > RADIUS ** 2:
                continue
            have.add(q)
            fresh.append(q)
P = P + fresh
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{SEED} universe (r<={RADIUS}): {n} points, {len(E)} edges, "
      f"seed {n0}  [{time.time()-t0:.0f}s]", flush=True)
for k in KS:
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    for col in range(k):          # break the colour symmetry on vertex 0
        if col:
            cls.append([-(1 + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    print(f"  k={k}: {'colourable' if ok else 'NOT COLOURABLE -- chi > '+str(k)}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
