"""A universe from a DIFFERENT seed, which is the one stone left unturned.

Everything measured in this project descends from de Grey's points.  When
the notes say "75500 configurations across two carriers", the carriers are
both closures of the same seed in the same field, so it is one point of the
space looked at from several angles, not a sweep of the space.

The missing ingredient is forcing at five colours: at four it exists -- 29
non-adjacent pairs in Sa, all at distance 1/3, verified individually -- and
at five whole distance classes come back empty.  Whether that is a fact
about five colours or a fact about THIS seed has never been tested, because
no other seed has ever been built.

Three seeds are tried here, each generating a different quadratic character
in the field, and each closed under circle intersection the same way:

    rhombus + 60 degrees   the classical seed, sqrt 3        (control)
    square  + 45 degrees   brings sqrt 2, absent from the usual field
    pentagon               brings sqrt 5 through a fivefold rotation
                           rather than through a distance

Then each is asked the only question that matters: does it carry a
non-adjacent pair forced to differ at five colours?
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
from collections import defaultdict
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from hn.field import Field
from hn.geometry import Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver
from sqrtK import Roots

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
ROUNDS = int(sys.argv[1]) if len(sys.argv) > 1 else 3
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 2200
RAD = float(sys.argv[3]) if len(sys.argv) > 3 else 2.2
t0 = time.time()


def build(field, seed_pts, name):
    K = field
    R = Roots(K)
    half = K.rational(Fr(1, 2))
    four = K.rational(4)
    cache = {}

    def ratio_root(D):
        k = tuple(D.c)
        if k not in cache:
            cache[k] = R.sqrt((four - D) / D)
        return cache[k]

    def key(p):
        return (tuple(p.x.c), tuple(p.y.c))

    pts = {}
    for p in seed_pts:
        pts[key(p)] = p
    for rnd in range(ROUNDS):
        cur = list(pts.values())
        added = 0
        for i, A in enumerate(cur):
            if len(pts) >= CAP:
                break
            for B in cur[i + 1:]:
                D = A.dist2(B)
                fd = float(D)
                if not (1e-9 < fd < 4.0):
                    continue
                s = ratio_root(D)
                if s is None:
                    continue
                mx, my = (A.x + B.x) * half, (A.y + B.y) * half
                nx = -(B.y - A.y) * s * half
                ny = (B.x - A.x) * s * half
                for P2 in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
                    if float(P2.x) ** 2 + float(P2.y) ** 2 > RAD ** 2:
                        continue
                    k = key(P2)
                    if k not in pts:
                        pts[k] = P2
                        added += 1
                if len(pts) >= CAP:
                    break
        if added == 0 or len(pts) >= CAP:
            break
    P = list(pts.values())
    b = IntBasis.covering(P)
    rows = b.rows(P)
    hr = b.overflow_headroom(rows)
    if hr >= 1.0:
        print(f"{name}: {len(P)} points but the basis overflows "
              f"({hr:.1f})  [{time.time()-t0:.0f}s]", flush=True)
        return
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, rows)))
    n = len(P)
    print(f"{name}: {n} points, {len(E)} edges, mean degree "
          f"{2*len(E)/max(1,n):.2f}  [{time.time()-t0:.0f}s]", flush=True)
    if len(E) < 10:
        return
    adj = defaultdict(set)
    for a, c in E:
        adj[a].add(c)
        adj[c].add(a)
    for KC in (4, 5):
        cls = [[1 + v * KC + c for c in range(KC)] for v in range(n)]
        for a, c in E:
            for col in range(KC):
                cls.append([-(1 + a * KC + col), -(1 + c * KC + col)])
        s = Solver(name="m22", bootstrap_with=cls)
        if not s.solve():
            print(f"   not {KC}-colourable!", flush=True)
            s.delete()
            continue
        forced = 0
        tested = 0
        for u in range(min(n, 250)):
            for v in range(u + 1, n):
                if v in adj[u]:
                    continue
                tested += 1
                if not s.solve(assumptions=[1 + u * KC, 1 + v * KC]):
                    forced += 1
                    if forced <= 2:
                        d = float(P[u].dist2(P[v])) ** 0.5
                        print(f"   k={KC}: FORCED pair {u},{v} at distance "
                              f"{d:.5f}  [{time.time()-t0:.0f}s]", flush=True)
                if tested > 30000:
                    break
            if tested > 30000:
                break
        s.delete()
        print(f"   k={KC}: {forced} forced of {tested} non-adjacent pairs"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    pickle.dump(P, open(SC + f"seed_{name}.pkl", "wb"))


# Small fields strangle the closure: sqrt3 gave 18 points, sqrt2 gave 7,
# and the two together 65.  The intersection of two unit circles needs
# (4-D)/D to be a SQUARE in the field, which almost never happens with one
# or two generators.  So the usual family is not a fashion -- a closure
# needs a field rich in square roots or it dies in the first round.
# The real question is which generators, at equal richness.
def seed4(gens, tag):
    K = Field(gens)
    h = K.rational(Fr(1, 2))
    pts = [Point(K.zero(), K.zero()), Point(K.one(), K.zero())]
    if 3 in gens:
        pts.append(Point(h, K.sqrt(3) * h))
        pts.append(Point(h, -(K.sqrt(3) * h)))
    if 2 in gens:
        r2h = K.sqrt(2) * h
        pts.append(Point(r2h, r2h))
        pts.append(Point(K.one() + r2h, r2h))
    if 5 in gens:
        pts.append(Point(K.rational(Fr(1, 4)) * (K.sqrt(5) - K.one()),
                         K.zero()))
    build(K, pts, tag)


seed4((3, 5, 7, 11), "degrey_field")
seed4((2, 3, 5, 7), "swap11for2")
seed4((3, 5, 7, 13), "swap11for13")
seed4((2, 3, 5, 7, 11), "five_generators")
print("DONE", flush=True)
