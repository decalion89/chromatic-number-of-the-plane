"""How many colourings are there, really?  The quantity behind the ladder.

The correlation ratio is an indirect reading of one thing: how large the space
of proper colourings is.  The triangular lattice scores eleven thousand at
three colours because its 3-colouring is unique up to permuting the colours --
there is exactly one -- and G scores one at five because there are so many
that two sampled at random share nothing.

So count them directly.  Sample colourings, rename each one's colours by
order of first appearance so that permutations collapse, and see how many
distinct ones come back.  One means the colouring is forced; twenty-four out
of twenty-four means the space is far larger than the sample and nothing can
be read off a pair of points.

Also reported: how many vertices take the same colour in all samples once
permutations are removed, which is the part of the graph that is pinned
outright.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_Y, build_G
from hn.geometry import DEGREY_FIELD as F, Point
from hn.graph import build_graph
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def lattice(side):
    h = F.sqrt(3) / F.rational(2)
    return [Point(F.rational(i) + F.rational(j) / F.rational(2),
                  h * F.rational(j))
            for i in range(side) for j in range(side)]


def canon(col, k):
    """Rename colours by order of first appearance, killing permutations."""
    seen, out = {}, []
    for c in col:
        if c not in seen:
            seen[c] = len(seen)
        out.append(seen[c])
    return tuple(out)


def count(name, pts, k):
    g = build_graph(pts)
    n = g.n
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: not colourable", flush=True)
        sv.delete()
        return
    rng, seen = random.Random(3011), []
    for s in range(SAMPLES):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        seen.append(canon([next(c for c in range(k) if m[w * k + c] > 0)
                           for w in range(n)], k))
    sv.delete()
    distinct = len(set(seen))
    pinned = sum(1 for v in range(n)
                 if len({t[v] for t in seen}) == 1)
    print(f"  {name} ({n} pts) at {k} colours: {distinct} distinct of "
          f"{SAMPLES} sampled, {pinned} of {n} vertices pinned to one colour "
          f"throughout  [{time.time()-t0:.0f}s]", flush=True)


L = lattice(20)
for k in (3, 4, 5):
    count("triangular lattice", L, k)
Sa = build_Sa(F)
for k in (4, 5):
    count("Sa", Sa, k)
Y = build_Y(F)
for k in (4, 5):
    count("Y", Y, k)
G = build_G(F, as_graph=False)
count("G", G, 5)
