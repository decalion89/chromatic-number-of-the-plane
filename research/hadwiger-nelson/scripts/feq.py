"""Forced-equal pairs at k colours, and whether their distance is spindleable.

The ladder verified on Y says exactly what a chi >= 6 graph needs:

    a 5-chromatic unit-distance graph H, and two vertices u, v of it such that
    EVERY proper 5-colouring of H gives u and v the same colour, with
    4|u-v|^2 - 1 a square in the field.

Then rotating H about u by 2*arcsin(1/(2|u-v|)) fixes u, so the two copies
agree there; forced equality carries that agreement to v and to its image;
and those two are one apart.  Contradiction.  That is de Grey's own step, and
it is the whole of it.

This searches for such pairs.  A pair that takes different colours in any
sampled colouring is not forced-equal, so a few dozen genuinely different
colourings kill almost every candidate for free; only the survivors are put to
the solver.  The distance test is applied at the end, not the start, because a
forced-equal pair at an illegal distance is still a result worth having.
"""
import sys, time, random, pickle
from fractions import Fraction as Fr
from math import isqrt
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import Point
from hn.field import Field
from hn.graph import build_graph
from pysat.solvers import Solver

F = Field((3, 5, 7, 11))
which = sys.argv[1] if len(sys.argv) > 1 else "G"
K = int(sys.argv[2]) if len(sys.argv) > 2 else 5
NCOL = int(sys.argv[3]) if len(sys.argv) > 3 else 40
pts = {"Sa": build_Sa, "Y": build_Y, "G": lambda f: build_G(f, as_graph=False)}[which](F)
g = build_graph(pts)
n = g.n
print(f"{which}: n={n}, m={sum(len(a) for a in g.adj)//2}, k={K}", flush=True)

X = lambda v, c: 1 + v * K + c
cnf = [[X(v, c) for c in range(K)] for v in range(n)]
for u, v in g.edges():
    for c in range(K):
        cnf.append([-X(u, c), -X(v, c)])

s = Solver(name="m22", bootstrap_with=cnf)
rng = random.Random(20260922)
cols, tries = [], 0
t0 = time.time()
# Diversity by randomised decision polarity, not by assumptions: a tight graph
# refuses almost every random partial assignment, which is why assumptions
# returned nothing at all here.  Random phases steer the solver to a different
# corner of the solution space without forbidding anything.
while len(cols) < NCOL and tries < NCOL * 6:
    tries += 1
    s.set_phases([(1 if rng.random() < 1.0 / K else -1) * X(v, c)
                  for v in range(n) for c in range(K)])
    if not s.solve():
        continue
    pos = set(l for l in s.get_model() if l > 0)
    cols.append([next(c for c in range(K) if X(v, c) in pos) for v in range(n)])
print(f"  {len(cols)} diverse {K}-colourings in {time.time()-t0:.0f}s", flush=True)

# a pair that ever differs is not forced-equal: intersect the partitions
from collections import defaultdict
groups = {v: () for v in range(n)}
for col in cols:
    for v in range(n):
        groups[v] = groups[v] + (col[v],)
buckets = defaultdict(list)
for v in range(n):
    buckets[groups[v]].append(v)
cand = [(a, b) for vs in buckets.values() if len(vs) > 1
        for i, a in enumerate(vs) for b in vs[i + 1:]]
print(f"  survivors of the colouring filter: {len(cand)} pairs "
      f"of {n*(n-1)//2}  ({100.0*len(cand)/(n*(n-1)/2):.4f}%)", flush=True)

def sqfree(x: Fr):
    num, den = x.numerator, x.denominator
    r = num * den
    if r <= 0: return None
    d = 2
    while d * d <= r:
        while r % (d * d) == 0: r //= d * d
        d += 1
    return r

def legal(d2):
    """Is the spindle at squared distance d2 available in Q(sqrt3,sqrt5,sqrt7,sqrt11)?"""
    if not d2.is_rational(): return None
    q = Fr(d2.c[0])
    if q <= Fr(1, 4): return None
    c = Fr(1) - Fr(1, 2) / q
    s2 = 1 - c * c
    if s2 == 0: return 1
    r = sqfree(s2)
    return r if r in (1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155) else -r

found, t0 = [], time.time()
for i, (a, b) in enumerate(cand):
    s2 = Solver(name="m22", bootstrap_with=cnf)
    diff = s2.solve(assumptions=[X(a, 0), X(b, 1)])
    s2.delete()
    if not diff:
        d2 = (pts[a] - pts[b]).norm2()
        lg = legal(d2)
        found.append((a, b, str(d2), lg))
        print(f"  FORCED EQUAL: v{a} v{b}   d^2 = {d2}   spindle radical {lg}", flush=True)
    if (i + 1) % 200 == 0:
        print(f"    ... {i+1}/{len(cand)} pairs, {len(found)} forced-equal "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(found)} forced-equal pairs at {K} colours in {which}", flush=True)
ok = [f for f in found if f[3] is not None and f[3] > 0]
print(f"  of which spindleable in this field: {len(ok)}", flush=True)
for f in ok[:20]: print("   ", f, flush=True)
pickle.dump(found, open(f"/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/feq_{which}_{K}.pkl", "wb"))
