"""The centre of a ring need not be a vertex.  Nothing here has ever asked.

Every ring scan in this work took a vertex as the centre and grouped the rest
by distance from it.  That is a real restriction and an arbitrary one: a bite
turns the graph about a point, and the point is free.  de Grey's two centres
happen to be vertices -- the origin is in Sa, (-2, 0) is in Y -- but nothing
about the mechanism requires it, and choosing only among vertices searches a
set of size n inside a continuum.

For a fixed radius the right question is the classical one: which point of the
plane has the most graph points at exactly that distance from it?  Every such
centre is equidistant from some pair, so it lies on the perpendicular bisector
of that pair at a computable offset -- and a centre with m points on its ring
is produced by all m*(m-1)/2 of their pairs.  So generating the centre of
every pair and counting multiplicities finds the best one, in one pass over
the pairs rather than a search over the plane.

The offset has to exist in the field or the centre is not constructible here:
for a pair at distance squared dd and a target radius squared R, the centre
sits at the midpoint plus sqrt(R/dd - 1/4) times the perpendicular, so
R/dd - 1/4 must be a square in K.

Run at the radii the field can bite, which is where a big ring would actually
be worth something.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G, build_Sa, build_Y
from hn.geometry import DEGREY_FIELD as K, Point

SEED = sys.argv[1] if len(sys.argv) > 1 else "G"
RADII = [Fr(x) for x in (sys.argv[2] if len(sys.argv) > 2
                         else "4,9,16,64").split(",")]
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
builders = {"Sa": build_Sa, "Y": build_Y,
            "G": lambda f: build_G(f, as_graph=False)}
P = builders[SEED](K)
n = len(P)
half = K.rational(Fr(1, 2))
print(f"{SEED}: {n} points, {n*(n-1)//2} pairs  [{time.time()-t0:.0f}s]",
      flush=True)


def rsqrt_rat(q):
    """sqrt of a positive rational, as a field element, when K holds it."""
    if q <= 0:
        return None
    for rr in CLASSES:
        t = q / rr
        a, c = t.numerator, t.denominator
        ra, rc = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc * rc == c:
            s = K.rational(Fr(ra, rc))
            return s if rr == 1 else K.sqrt(rr) * s
    return None


def biteable(D):
    ct = 1 - Fr(1, 2) / D
    s2 = 1 - ct * ct
    for rr in CLASSES:
        q = s2 / rr
        a, c = q.numerator, q.denominator
        ra, rc = int(round(a ** .5)), int(round(c ** .5))
        if ra * ra == a and rc * rc == c:
            return f"cos={ct}, sin={Fr(ra,rc)}" + (f"*sqrt{rr}" if rr != 1
                                                   else "")
    return None


D2 = [(p.x * p.x + p.y * p.y) for p in P]
for R in RADII:
    bt = biteable(R)
    cnt = defaultdict(int)
    cache = {}
    pairs = 0
    for i in range(n):
        A = P[i]
        for j in range(i + 1, n):
            B = P[j]
            dx, dy = B.x - A.x, B.y - A.y
            dd = dx * dx + dy * dy
            if any(dd.c[1:]):
                continue
            q = dd.c[0]
            if not (0 < q <= 4 * R):
                continue
            key = q
            if key not in cache:
                cache[key] = rsqrt_rat(R / q - Fr(1, 4))
            h = cache[key]
            if h is None:
                continue
            pairs += 1
            mx, my = (A.x + B.x) * half, (A.y + B.y) * half
            px, py = -dy * h, dx * h
            for C in (Point(mx + px, my + py), Point(mx - px, my - py)):
                cnt[C] += 1
    if not cnt:
        print(f"R={R}: no constructible centre  [{time.time()-t0:.0f}s]",
              flush=True)
        continue
    best = max(cnt.values())
    # a centre seen t times carries m points with m*(m-1)/2 = t
    m = int((1 + (1 + 8 * best) ** .5) / 2)
    tops = [C for C, t in cnt.items() if t == best]
    top = tops[0]
    inG = top in set(P)
    print(f"R={R} (rho={float(R)**.5:.3f}): {pairs} usable pairs, "
          f"{len(cnt)} candidate centres, best ring {m} points at "
          f"({float(top.x):.4f},{float(top.y):.4f}), that centre is a vertex: "
          f"{inG}; {bt or 'not biteable in K'}  [{time.time()-t0:.0f}s]",
          flush=True)
    with open(f"/tmp/claude-0/-home-user-darwin-50/"
              f"aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
              f"centres_{SEED}_{R}.pkl".replace("/", "_"), "wb") as f:
        pickle.dump(sorted(((t, i) for i, (C, t) in enumerate(cnt.items())),
                           reverse=True)[:50], f)
    hist = defaultdict(int)
    for t in cnt.values():
        hist[int((1 + (1 + 8 * t) ** .5) / 2)] += 1
    print(f"   ring-size histogram: "
          f"{dict(sorted(hist.items(), reverse=True)[:8])}", flush=True)
print("DONE", flush=True)
