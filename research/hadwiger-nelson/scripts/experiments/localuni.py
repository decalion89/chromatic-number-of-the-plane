"""Cap the dense core by thickening the carrier exactly where it lives.

The thirteen-point core needs five colours for its own two distances, so
capping it at four would immediately force a monochromatic pair at sqrt3.  Its
palette is 5 in G and in all three order-twelve closures, but palette only
falls as the carrier grows, and those closures add their points far away --
twelve copies spread round a circle, most of them nowhere near the core.

Constraint is local.  The points that can remove a colouring of the core are
the ones adjacent to it, and the way to make many of those is the unit-circle
intersection: every pair of points less than two apart defines two more, each
a unit from both.  Doing that around the core's own neighbourhood thickens the
graph precisely where it has to bite.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import ring_palette_bound
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
RAD = float(sys.argv[2]) if len(sys.argv) > 2 else 2.2
ROUNDS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
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


W, Wpts = pickle.load(open(SC + "core13.pkl", "rb"))
cx = sum(float(p.x) for p in Wpts) / len(Wpts)
cy = sum(float(p.y) for p in Wpts) / len(Wpts)
print(f"core centroid ({cx:.3f}, {cy:.3f}), thickening within {RAD}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
P = list(build_G(K, as_graph=False))
for rnd in range(ROUNDS):
    have = set(P)
    near = [p for p in P
            if (float(p.x) - cx) ** 2 + (float(p.y) - cy) ** 2 < (RAD + 1) ** 2]
    fresh = []
    for i, A in enumerate(near):
        for B in near[i + 1:]:
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
                if q in have:
                    continue
                if (float(q.x) - cx) ** 2 + (float(q.y) - cy) ** 2 > RAD ** 2:
                    continue
                have.add(q)
                fresh.append(q)
    P = P + fresh
    print(f"round {rnd}: {len(fresh)} new points, {len(P)} total"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    b = IntBasis.covering(P)
    r = b.rows(P)
    hr = b.overflow_headroom(r)
    if hr >= 1.0:
        print(f"   overflow headroom {hr:.2f} -- stopping", flush=True)
        break
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    idx = {p: i for i, p in enumerate(P)}
    T = [idx[p] for p in Wpts]
    pal = ring_palette_bound(cls, n * k, T, k) if ok else None
    print(f"   {n} points, {len(E)} edges, {k}-colourable {ok}, "
          f"palette of the core = {pal}  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("   *** NOT 5-COLOURABLE ***", flush=True)
        pickle.dump(P, open(SC + "WITNESS_localuni.pkl", "wb"))
        break
    if pal is not None and pal <= 4:
        print("   *** CORE CAPPED: a forced pair at sqrt3 follows",
              flush=True)
        pickle.dump((P, T), open(SC + "capped_core.pkl", "wb"))
        break
print("DONE", flush=True)
