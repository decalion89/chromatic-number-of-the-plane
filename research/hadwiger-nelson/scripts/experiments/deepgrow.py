"""Grow towards uncolourability with a pool that deepens instead of running out.

The first version of this greedy scored candidates by how many sampled
colourings their neighbourhood saturates -- a point cannot be coloured when
its neighbours already show all five -- and the criterion works: the first
point it chose killed 193 of 250 colourings.

But its candidate pool was fixed, the first round of unit-circle
intersections of G, and adding every one of those gives the 35132-point
universe, which is 5-colourable.  A greedy that cannot do better than its
pool cannot beat a pool that is already known to colour.  The score was
right and the reservoir was wrong.

So let the pool deepen.  Each point added is intersected against the points
already chosen, and those second-, third- and later-round intersections join
the candidates.  The reachable set is then the full intersection closure
rather than one round of it, and the greedy is choosing where to dig rather
than which of a fixed list to take.

Exactness is kept throughout, and the int64 basis is rebuilt as the point set
grows, with the overflow guard checked every time.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NSAMP = int(sys.argv[2]) if len(sys.argv) > 2 else 200
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
POOLCAP = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
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


def meets(A, B):
    """The two points a unit from both A and B, when the field holds them."""
    D = A.dist2(B)
    if not .05 < float(D) < 3.99:
        return []
    sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
    if sD is None or s4 is None:
        return []
    inv = K.rational(1) / sD
    mx, my = (A.x + B.x) * half, (A.y + B.y) * half
    nx = -(B.y - A.y) * inv * s4 * half
    ny = (B.x - A.x) * inv * s4 * half
    return [Point(mx + nx, my + ny), Point(mx - nx, my - ny)]


chosen = list(build_G(K, as_graph=False))
inset = set(chosen)
print(f"seed {len(chosen)} points  [{time.time()-t0:.0f}s]", flush=True)
rng = random.Random(2971215073)

# A modest first pool from the seed, then deepened as points are added.
pool = []
inpool = set()
src = chosen[::3]
for i, A in enumerate(src):
    for B in src[i + 1:]:
        for q in meets(A, B):
            if q not in inset and q not in inpool:
                inpool.add(q)
                pool.append(q)
    if len(pool) > POOLCAP:
        break
print(f"initial pool {len(pool)}  [{time.time()-t0:.0f}s]", flush=True)


def build(points):
    b = IntBasis.covering(points)
    rows = b.rows(points)
    hr = b.overflow_headroom(rows)
    if hr >= 1.0:
        return None, None, hr
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, rows)))
    return b, E, hr


def sample(points, E, want):
    m = len(points)
    cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        return None
    out = []
    for _ in range(want):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(m * k)])
        sv.solve()
        mo = sv.get_model()
        out.append([next(c for c in range(k) if mo[v * k + c] > 0)
                    for v in range(m)])
    sv.delete()
    return out


for step in range(STEPS):
    allpts = chosen + pool
    b, E, hr = build(allpts)
    if b is None:
        print(f"overflow headroom {hr:.2f} -- stop", flush=True)
        break
    nch = len(chosen)
    nb = defaultdict(set)
    for a, c in E:
        if a < nch and c >= nch:
            nb[c].add(a)
        elif c < nch and a >= nch:
            nb[a].add(c)
        elif a < nch and c < nch:
            pass
    Ecore = [(a, c) for a, c in E if a < nch and c < nch]
    cols = sample(chosen, Ecore, NSAMP)
    if cols is None:
        print(f"\n*** NOT {k}-COLOURABLE at {nch} points "
              f"[{time.time()-t0:.0f}s]", flush=True)
        pickle.dump(chosen, open(SC + "WITNESS_deep.pkl", "wb"))
        break
    best, bestv = None, -1
    for v in range(nch, len(allpts)):
        nbs = nb.get(v)
        if not nbs or len(nbs) < k:
            continue
        cnt = 0
        for c in cols:
            s = set()
            for u in nbs:
                s.add(c[u])
                if len(s) == k:
                    cnt += 1
                    break
        if cnt > bestv:
            best, bestv = v, cnt
    if best is None:
        print("no candidate with k neighbours -- stop", flush=True)
        break
    newp = allpts[best]
    chosen.append(newp)
    inset.add(newp)
    pool.remove(newp)
    # deepen: intersect the new point against a slice of the chosen set
    add = 0
    for A in chosen[::7]:
        if A == newp:
            continue
        for q in meets(A, newp):
            if q not in inset and q not in inpool:
                inpool.add(q)
                pool.append(q)
                add += 1
        if add > 60:
            break
    if step % 5 == 0:
        print(f"   step {step}: {len(chosen)} points, last killed {bestv}/"
              f"{NSAMP}, pool {len(pool)} (+{add})  [{time.time()-t0:.0f}s]",
              flush=True)
    if step % 50 == 0 and step:
        pickle.dump(chosen, open(SC + "deepgrow.pkl", "wb"))
print("DONE", flush=True)
