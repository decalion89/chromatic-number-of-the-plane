"""How rare is a seed whose closure has a capped ring?

Sa is the twelve-element dihedral closure of a 39-point seed, and what makes
the whole construction work is that the closure's D = 4 ring is CAPPED at two
colours.  De Grey found that seed by search.  Nothing here knows how large a
search that was, and the number matters: it is the only guide to how big the
search would be one level up.

So calibrate.  Draw random 39-point seeds from the natural universe -- the
points a short unit walk reaches from the origin -- close each under the same
group, and ask whether ANY ring about the origin is capped below four colours.
The screen is one SAT call per ring: a ring that can show all k colours at once
is not capped, and that is the satisfiable and therefore fast answer, so a
seed that fails costs almost nothing.

Unit rings are excluded: their cap at k-1 is the colourability of their centre
and says nothing.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, deque
from hn.geometry import DEGREY_FIELD as K, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = 4
TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SEED_SIZE = 39
t0 = time.time()
rng = random.Random(17)
half = K.rational(Fr(1, 2))
rot60 = Rotation(half, K.sqrt(3) * half)
ZERO = Point(K.zero(), K.zero())

# The universe has to contain a seed that works, or the rate measures nothing.
# The first attempt drew from the points three unit steps from the origin --
# 253 of them, at radii 0.5 and up -- and only FOUR of de Grey's thirty-nine
# lie there: his seed sits at radii 0.168, 0.264, 0.292, 0.333 and so on,
# well inside the unit disc, where unit walks do not reach.  Zero of four
# hundred was therefore a fact about the wrong space.
#
# Sa's own 397 points are the right universe: they contain S by construction,
# so a capping seed is present and the question is how rare it is among
# subsets of the same size.
from hn.degrey import build_Sa, build_S
pool = build_Sa(K)
_S = set(build_S(K))
assert _S <= set(pool), "the universe must contain de Grey's own seed"
print(f"universe: Sa's {len(pool)} points, containing all {len(_S)} of "
      f"de Grey's seed  [{time.time()-t0:.0f}s]", flush=True)


def closure(seed):
    out, s = [], set()
    for p in seed:
        for base in (p, Point(p.x, -p.y)):
            q = base
            for _ in range(6):
                if q not in s:
                    s.add(q)
                    out.append(q)
                q = rot60(q)
    return out


def capped_rings(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        return None
    n = len(P)
    if ZERO not in P:
        return []
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    dm, D2 = b.dim, b.D * b.D
    zi = P.index(ZERO)
    d = r - r[zi]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), D2)
        if v and v != 1 and closable_distance(v):
            grp[v].append(int(off))
    rings = [(D, m) for D, m in grp.items() if len(m) >= 4]
    if not rings:
        return []
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    nv = n * k
    ind = [nv + 1 + c for c in range(k)]
    enc = CardEnc.atleast(lits=list(ind), bound=k, top_id=nv + k,
                          encoding=EncType.seqcounter)
    sel0 = max([nv + k] + [abs(x) for cl in enc.clauses for x in cl]) + 1
    extra_cls = list(enc.clauses)
    for i, (D, mem) in enumerate(rings):
        for c in range(k):
            extra_cls.append([-(sel0 + i), -ind[c]] +
                             [1 + v * k + c for v in mem])
    sv = Solver(name="cd15", bootstrap_with=cls + extra_cls)
    if not sv.solve():
        sv.delete()
        return "UNCOLOURABLE"
    out = [str(rings[i][0]) for i in range(len(rings))
           if not sv.solve(assumptions=[sel0 + i])]
    sv.delete()
    return out


hits = 0
sizes = []
for t in range(1, TRIALS + 1):
    seed = [ZERO] + rng.sample([p for p in pool if p != ZERO], SEED_SIZE - 1)
    P = closure(seed)
    sizes.append(len(P))
    res = capped_rings(P)
    if res is None:
        continue
    if res == "UNCOLOURABLE":
        print(f"*** trial {t}: {len(P)} points NOT {k}-COLOURABLE ***",
              flush=True)
        hits += 1
        continue
    if res:
        hits += 1
        print(f"*** trial {t}: {len(P)} points, capped rings {res} ***"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if t % 50 == 0:
        print(f"   {t}/{TRIALS} seeds, {hits} with a capped ring, closures "
              f"average {sum(sizes)//len(sizes)} points"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{hits} of {TRIALS} random seeds give a capped ring"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
