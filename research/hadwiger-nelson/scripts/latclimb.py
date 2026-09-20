"""Start the climb from the most correlated object there is, not from de Grey's.

The ladder says correlation at k = chi collapses as chi rises: 11014 for the
triangular lattice at three colours, 24 for Sa at four, 1.1 for G at five.
Every construction in this session started from Sa, at 24.  The lattice sits
a thousand times higher, and it has never been used as a base.

It can be.  Its 3-colouring is unique up to permutation, so two points agree
in every colouring exactly when they lie in the same class -- (0,0) and (3,0)
do -- and that pair is at distance 3, squared 9, with sqrt(4*9-1) = sqrt(35) =
sqrt(5)*sqrt(7), which is in de Grey's field.  So the spindle exists:
rotation_joining(9) about one end brings the other to distance 1 from itself,
and lattice u rho(lattice) cannot be 3-coloured.

The question is not whether that works -- it must, by the lemma -- but what
the result inherits.  If spindling a graph of correlation 11014 gives a
4-chromatic graph of correlation far above Sa's 24, then the base matters and
the climb should be restarted from here.  If it gives 24 again, the ladder's
collapse is a property of the colour count alone and no base helps.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.geometry import DEGREY_FIELD as F, Point, rotation_joining
from hn.graph import build_graph
from hn.homcol import closable_distance
from pysat.solvers import Solver

SAMPLES = 24
t0 = time.time()


def lattice(side):
    h = F.sqrt(3) / F.rational(2)
    return [Point(F.rational(i) + F.rational(j) / F.rational(2),
                  h * F.rational(j))
            for i in range(-side, side + 1) for j in range(-side, side + 1)]


def measure(name, pts, k):
    g = build_graph(pts)
    n = g.n
    E = set((min(a, b), max(a, b)) for a, b in g.edges())
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in sorted(E):
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="g4", bootstrap_with=cls)
    if not sv.solve():
        print(f"  {name} at {k}: NOT {k}-COLOURABLE  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return False
    rng, cols = random.Random(7919), []
    for s in range(SAMPLES):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(n * k)])
        sv.solve()
        m = sv.get_model()
        cols.append([next(c for c in range(k) if m[w * k + c] > 0)
                     for w in range(n)])
    sv.delete()
    C = np.array(cols, dtype=np.int8)
    nd = 0
    for i in range(n - 1):
        nd += int((C[:, i + 1:] != C[:, i:i + 1]).all(axis=0).sum())
    nd -= len(E)
    ch = n * (n - 1) // 2 * ((k - 1) / k) ** SAMPLES
    from itertools import combinations, permutations
    perms = list(permutations(range(k)))
    lut = np.array(perms, dtype=np.int8)
    tot, cnt = 0.0, 0
    for a, b in list(combinations(range(SAMPLES), 2))[:40]:
        tot += min(int((C[a] != lut[i][C[b]]).sum())
                   for i in range(len(perms))) / n
        cnt += 1
    print(f"  {name} at {k}: {n} points, {len(E)} edges, ratio "
          f"{nd/max(ch,1e-9):.1f}, samples disagree on {tot/cnt:.3f} "
          f"(independent {1-1/k:.3f})  [{time.time()-t0:.0f}s]", flush=True)
    return True


assert closable_distance(Fr(9)), "the spindle at distance 3 must exist"
L = lattice(9)
measure("lattice", L, 3)

P = Point(F.rational(3), F.zero())
Q = Point(F.zero(), F.zero())
assert P in set(L) and Q in set(L)
rho = rotation_joining(Fr(9), F).about(Q)
d2 = (P.x - rho(P).x) ** 2 + (P.y - rho(P).y) ** 2
print(f"spindle: rho = {rotation_joining(Fr(9), F).cos} + i "
      f"{rotation_joining(Fr(9), F).sin}; P and rho(P) at squared distance "
      f"{d2}  [{time.time()-t0:.0f}s]", flush=True)

W, ws = list(L), set(L)
for p in L:
    q = rho(p)
    if q not in ws:
        ws.add(q)
        W.append(q)
print(f"W = lattice u rho_9(lattice): {len(W)} points  "
      f"[{time.time()-t0:.0f}s]", flush=True)
measure("W", W, 3)
measure("W", W, 4)
