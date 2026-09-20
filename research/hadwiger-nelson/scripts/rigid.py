"""How rigid is ONE copy?  That is the quantity, and it was never cross edges.

Z4 = Sa u rho(Sa) has 793 points and 3954 edges.  Two separate copies of Sa
would have 3948.  So de Grey's forcing at four colours is produced by SIX
cross edges -- and the translate unions built here carried 1552 of them, the
stacks 4306, and forced nothing.  Counting cross edges was measuring the wrong
thing from the start, with the anatomy of Y already written down.

So the quantity to look at is how constrained a single copy is before
anything is glued to it: how many pairs of its points take the same colour in
every proper colouring.  The guess written here first was that six edges can
only close a structure that is nearly closed already.

The measurement says otherwise, and flatly.  Sa at four colours has ZERO pairs
that agree in every sampled colouring -- it is completely free -- and Sa u
rho(Sa), six edges later, has seventy-five, thirty-five of them at a closable
distance, including all six antipodal pairs of the ring.  Nothing was nearly
closed.  Forcing does not accumulate; it arrives.
"""
import sys, time, random
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G, build_Sa, build_Y, build_graph
from hn.geometry import DEGREY_FIELD as F, Point, _rot60
from hn.homcol import agreeing_pairs, closable_distance
from pysat.solvers import Solver

t0 = time.time()


def rigidity(name, pts, k, samples=14):
    g = build_graph(pts)
    E = sorted((min(a, b), max(a, b)) for a, b in g.edges())
    n = len(pts)
    cls = [[1 + w * k + c for c in range(k)] for w in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name}: NOT {k}-colourable  [{time.time()-t0:.0f}s]",
              flush=True)
        sv.delete()
        return
    rng, cols = random.Random(2357), []
    for s in range(samples):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + v)
                       for v in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    sv.delete()
    pairs = agreeing_pairs(cols, colours=k)
    close = 0
    for i, j in pairs:
        v = float((pts[i].x - pts[j].x) ** 2 + (pts[i].y - pts[j].y) ** 2)
        D = Fr(round(v * 55440), 55440)
        if abs(float(D) - v) < 1e-9 and D >= Fr(1, 4) and closable_distance(D):
            close += 1
    print(f"  {name}: {n} points, {len(E)} edges, {k} colours -- "
          f"{len(pairs)} pairs agree in all {samples} samples, {close} of "
          f"them at a closable distance; that is {len(pairs)/n:.3f} per "
          f"point  [{time.time()-t0:.0f}s]", flush=True)


rigidity("Sa at four", build_Sa(F), 4)
rigidity("Y at four", build_Y(F), 4)
rigidity("Sa at five", build_Sa(F), 5)
rigidity("G at five", build_G(F, as_graph=False), 5)

PIV = Point(F.rational(-2), F.zero())
rot = _rot60(F).about(PIV)
G = build_G(F, as_graph=False)
Gs, seen = [], set()
for refl in (False, True):
    for j in range(6):
        for p in G:
            q = p
            for _ in range(j):
                q = rot(q)
            if refl:
                q = Point(q.x, -q.y)
            if q not in seen:
                seen.add(q)
                Gs.append(q)
rigidity("G* at five", Gs, 5)
