"""Sa at four, against the per-vertex law -- the contrast that matters.

The 5-chromatic graphs beat the per-vertex null model by a factor of two,
steadily: 0.51, 0.50, 0.51.  The tuned chain beats it by 1.15.  The question is
what Sa does at four, because Sa is the only object here that reaches free@k
EXACTLY zero, and zero is not a factor -- it is a different kind of statement.

It also asks whether zero is a degree effect.  free@k is dominated by the
low-degree tail: a vertex of degree 4 at five colours is free with probability
about 4*(3/4)^4 = 1.27, i.e. essentially always.  So free@k = 0 requires a
graph with no thin boundary at all, and the minimum degree is the number to
look at, not the mean.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, math
from fractions import Fraction as Fr
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.field import Field
from hn.degrey import build_Sa, build_Y
from hn.geometry import Point, rotation_joining, _rot60
from hn.graph import build_graph

F = Field((3, 5, 7, 11))
rot60 = _rot60(F); Sa = build_Sa(F)
def orb(p):
    out, q = [], p
    for _ in range(6):
        if q not in out: out.append(q)
        q = rot60(q)
    for r in [Point(z.x, -z.y) for z in list(out)]:
        if r not in out: out.append(r)
    return out
def grow(b):
    seen, U = set(Sa), list(Sa)
    for w in orb(Sa[b]):
        rot = rotation_joining(Fr(1), F).about(w)
        for p in Sa:
            q = rot(p)
            if q not in seen: seen.add(q); U.append(q)
    return U

def f(d, K):
    return (K - 1) * ((K - 2) / (K - 1)) ** d

CASES = [("Sa", list(Sa), 4), ("Y", list(build_Y(F)), 4),
         ("Sa[25]", grow(25), 4), ("Sa[265]", grow(265), 4)]
print(f"{'graph':10s} {'n':>5s} {'deg':>6s} {'min':>4s} {'sd':>5s} "
      f"{'mean-law':>9s} {'per-vertex':>11s}  low-degree tail")
for name, U, K in CASES:
    g = build_graph(U); n = g.n
    degs = [len(a) for a in g.adj]
    mu = sum(degs) / n
    sd = math.sqrt(sum((x - mu) ** 2 for x in degs) / n)
    pv = sum(f(x, K) for x in degs) / n
    c = Counter(degs)
    tail = " ".join(f"{d}:{c[d]}" for d in sorted(c)[:4])
    print(f"{name:10s} {n:5d} {mu:6.2f} {min(degs):4d} {sd:5.2f} "
          f"{100*f(mu,K):8.3f}% {100*pv:10.3f}%  {tail}", flush=True)
