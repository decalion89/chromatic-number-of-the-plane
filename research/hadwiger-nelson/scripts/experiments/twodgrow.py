"""Grow the TWO-DISTANCE graph towards uncolourability -- the cheap target.

Making the unit-distance graph need six colours is the expensive goal: the
smallest known graph needing five already has about five hundred vertices.
Making the {1, sqrt3} graph need six is far cheaper, and it is just as useful.
If that graph has no proper 5-colouring, then in every 5-colouring of the
points some edge of it is monochromatic; edges at distance 1 cannot be, so
some pair at distance sqrt3 is -- a forced monochromatic pair at a named
distance, with no cap needed anywhere.

The scale difference is already measured: {1, sqrt3} needs five colours on
NINE points, where a unit-distance graph needs about five hundred for the
same.  Two orders of magnitude.

The growth is the same idea as before -- a point cannot be coloured when its
neighbourhood already shows all five colours, so add the one that destroys the
most sampled colourings -- but applied to the two-distance graph, with its
sqrt3 edges included in the neighbourhood and in the colouring constraints.
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
FOCUS = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
CKPTNAME = sys.argv[6] if len(sys.argv) > 6 else "twod.pkl"
DSQ = Fr(sys.argv[7]) if len(sys.argv) > 7 else Fr(3)
BATCH = int(sys.argv[8]) if len(sys.argv) > 8 else 1
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/hn/")
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


# The container restarts, so the run has to survive it: pick up whatever the
# last checkpoint holds instead of starting from the seed again.
import os
CKPT = SC + CKPTNAME
if os.path.exists(CKPT):
    chosen = pickle.load(open(CKPT, "rb"))
    print(f"resumed from checkpoint: {len(chosen)} points", flush=True)
else:
    chosen = list(build_G(K, as_graph=False))
inset = set(chosen)
print(f"seed {len(chosen)} points  [{time.time()-t0:.0f}s]", flush=True)
rng = random.Random(2971215073)
# A 6-chromatic graph has to be dense somewhere, and growth spread over the
# whole of G dilutes every point it adds.  Restricting candidates to a disc
# concentrates the same effort where it can compound.
FX = FY = 0.0
if FOCUS:
    FX = sum(float(q.x) for q in chosen) / len(chosen)
    FY = sum(float(q.y) for q in chosen) / len(chosen)
    print(f"focus disc of radius {FOCUS} at ({FX:.3f},{FY:.3f})", flush=True)


def infocus(q):
    if not FOCUS:
        return True
    return (float(q.x) - FX) ** 2 + (float(q.y) - FY) ** 2 <= FOCUS * FOCUS

# A modest first pool from the seed, then deepened as points are added.
pool = []
inpool = set()
src = chosen[::3]
for i, A in enumerate(src):
    for B in src[i + 1:]:
        for q in meets(A, B):
            if q not in inset and q not in inpool and infocus(q):
                inpool.add(q)
                pool.append(q)
    if len(pool) > POOLCAP:
        break
print(f"initial pool {len(pool)}  [{time.time()-t0:.0f}s]", flush=True)


sys.path.insert(0, SC)
from pairsat import pairs_at


def build(points):
    b = IntBasis.covering(points)
    rows = b.rows(points)
    hr = b.overflow_headroom(rows)
    if hr >= 1.0:
        return None, None, hr
    E1 = set((min(a, c), max(a, c))
             for a, c in fast_edges_complete(b, rows))
    E2 = set(pairs_at(b, rows, int(DSQ * b.D * b.D)))
    return b, sorted(E1 | E2), hr


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
    # Re-scoring after every single addition is the pure greedy and it is slow.
    # Taking the best few at once costs a little accuracy -- the scores shift as
    # soon as the first of them is in -- and multiplies the rate of growth.
    scored = []
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
        scored.append((cnt, v))
    if not scored:
        print("no candidate with k neighbours -- stop", flush=True)
        break
    scored.sort(reverse=True)
    take = scored[:BATCH]
    bestv = take[0][0]
    fresh = []
    for _cnt, v in take:
        q = allpts[v]
        if q in inset:
            continue
        chosen.append(q)
        inset.add(q)
        pool.remove(q)
        fresh.append(q)
    newp = fresh[-1] if fresh else None
    if newp is None:
        break
    # deepen: intersect the new point against a slice of the chosen set
    add = 0
    for A in chosen[::2]:
        if A in fresh:
            continue
        for q in meets(A, newp):
            if q not in inset and q not in inpool and infocus(q):
                inpool.add(q)
                pool.append(q)
                add += 1
        if add > 200:
            break
    if step % 5 == 0:
        print(f"   step {step}: {len(chosen)} points, last killed {bestv}/"
              f"{NSAMP}, pool {len(pool)} (+{add})  [{time.time()-t0:.0f}s]",
              flush=True)
    if step % 10 == 0 and step:
        pickle.dump(chosen, open(SC + "deepgrow.pkl", "wb"))
print("DONE", flush=True)
