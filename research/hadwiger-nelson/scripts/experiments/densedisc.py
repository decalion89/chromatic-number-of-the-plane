"""Build the dense closure of a small disc and measure its two-distance
chromatic number.

The greedy adds one point at a time and its kill rate has settled around ten
per cent: it is still hardening the graph, but slowly, and point-by-point
growth is the wrong tool if what is needed is simply a much denser region.

The target is worth restating because its scale is favourable.  A
unit-distance graph needs about five hundred vertices to force five colours.
The {1, sqrt3} graph forces five on NINE.  If the step from five to six is
anything like proportionate, six might be reachable in hundreds of points
rather than the hundreds of thousands the unit-distance version would need.
And six is enough: a {1, sqrt3} graph with no proper 5-colouring puts a
monochromatic pair at distance sqrt3 in every 5-colouring of the plane's
points, since its distance-1 edges cannot carry one.

So take the points of G inside a small disc, close them under unit-circle
intersection twice over, keep everything that stays inside the disc, and ask
the whole thing for five colours at both distances.  Density where it counts,
rather than volume everywhere.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pairsat import pairs_at
from pysat.solvers import Solver

RAD = float(sys.argv[1]) if len(sys.argv) > 1 else 1.3
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
DSQ = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(3)
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 40000
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


G = build_G(K, as_graph=False)
cx = sum(float(p.x) for p in G) / len(G)
cy = sum(float(p.y) for p in G) / len(G)
P = [p for p in G
     if (float(p.x) - cx) ** 2 + (float(p.y) - cy) ** 2 <= RAD * RAD]
# A disc of radius r has diameter 2r, so it cannot contain a pair at any
# distance above that.  Asking for sqrt3 edges inside a radius-0.5 disc is
# asking for something geometrically impossible, and the answer was zero
# edges -- which the chromatic loop then reported as 5 only because it
# started counting there.
assert 2 * RAD > float(DSQ) ** .5, (
    f"radius {RAD} gives diameter {2*RAD}, below the distance "
    f"{float(DSQ) ** .5:.3f} being asked for")
print(f"disc radius {RAD} at ({cx:.3f},{cy:.3f}): {len(P)} points of G"
      f"  [{time.time()-t0:.0f}s]", flush=True)
have = set(P)
for rnd in range(ROUNDS):
    fresh = []
    cur = list(P)
    for i, A in enumerate(cur):
        for B in cur[i + 1:]:
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
                if (float(q.x) - cx) ** 2 + (float(q.y) - cy) ** 2 > RAD * RAD:
                    continue
                have.add(q)
                fresh.append(q)
        if len(P) + len(fresh) > CAP:
            break
    P = P + fresh
    print(f"round {rnd}: +{len(fresh)} -> {len(P)} points"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if len(P) > CAP:
        break
pickle.dump(P, open(SC + f"disc_{RAD}_{ROUNDS}.pkl", "wb"))
b = IntBasis.covering(P)
r = b.rows(P)
hr = b.overflow_headroom(r)
print(f"overflow headroom {hr:.3f}  [{time.time()-t0:.0f}s]", flush=True)
assert hr < 1.0
E1 = set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r))
E2 = set(pairs_at(b, r, int(DSQ * b.D * b.D)))
n = len(P)
print(f"{n} points, {len(E1)} unit + {len(E2)} at distance^2 {DSQ}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
for name, ed in (("unit only", sorted(E1)),
                 (f"with distance^2 {DSQ}", sorted(E1 | E2))):
    for kk in range(2, 9):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(n)]
        for a, c in ed:
            for col in range(kk):
                cls.append([-(1 + a * kk + col), -(1 + c * kk + col)])
        for col in range(1, kk):
            cls.append([-(1 + col)])
        s = Solver(name="cd15", bootstrap_with=cls)
        ok = s.solve()
        s.delete()
        if ok:
            mark = ("   <<< SIX: a forced monochromatic pair at that distance"
                    if kk >= 6 and "distance" in name else "")
            print(f"   {name}: chi = {kk}{mark}  [{time.time()-t0:.0f}s]",
                  flush=True)
            break
        print(f"   {name}: not {kk}-colourable  [{time.time()-t0:.0f}s]",
              flush=True)
print("DONE", flush=True)
