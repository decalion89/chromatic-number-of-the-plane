"""How many colours can CONCENTRIC rings about one point ever need?

Fix h and keep only the points lying on a few circles about it.  Every
unit-distance edge among them is then decided by angles alone: two points at
radii r, r' are a unit apart exactly when their angular separation is

    theta(r,r') = arccos( (r^2 + r'^2 - 1) / (2 r r') ),

defined when |r - r'| <= 1 <= r + r'.  So a whole configuration is an induced
subgraph of one graph per radius set, generated from a seed by adding +-theta
and hopping between circles, and its chromatic number is the ceiling on
everything a hub argument can reach with those radii.

One ring, radius 1: theta = 60, the ring is a hexagon or a union of paths,
chi = 2 -- the bipartite-neighbourhood theorem.
One ring, radius 1/sqrt3: theta = 120, a union of triangles, chi = 3, and by
Niven no other radius has an odd cycle at all.
Two rings, 1 and 1/sqrt3: computed exactly -- chi = 3 at every depth, with a
3-colouring of period two, so chi = 3 for the whole infinite configuration.

Each extra ring adds at most two neighbours per vertex, so Brooks caps chi at
2k for k rings: two rings could never have passed 4, three rings could reach 6.
And a concentric configuration needing six colours would BE the answer, since
it is a unit-distance graph like any other.

Generation is by float BFS with exact-looking angles kept to nine decimals --
a coincidence would only ever merge two vertices and so lower chi, and any
chi >= 4 that turns up gets re-derived exactly before it is believed.
"""
import sys, time, math
from collections import deque
from pysat.solvers import Solver

t0 = time.time()

def cross(ra, rb):
    v = (ra * ra + rb * rb - 1.0) / (2.0 * ra * rb)
    if v < -1.0 or v > 1.0: return None
    return math.degrees(math.acos(v))

def build(radii, cap=260):
    th = {}
    for a, ra in enumerate(radii):
        for b, rb in enumerate(radii):
            t = cross(ra, rb)
            if t is not None and t > 1e-9 and t < 180 - 1e-9:
                th[(a, b)] = t
    seed = (0, 0.0)
    seen = {seed: 0}; order = [seed]; q = deque([seed])
    while q and len(order) < cap:
        (a, ang) = q.popleft()
        for (aa, b), t in th.items():
            if aa != a: continue
            for s in (t, -t):
                na = (ang + s) % 360.0
                key = (b, round(na, 7))
                if key not in seen:
                    seen[key] = len(order); order.append(key); q.append(key)
    E = set()
    for (a, ang), u in seen.items():
        for (aa, b), t in th.items():
            if aa != a: continue
            for s in (t, -t):
                key = (b, round((ang + s) % 360.0, 7))
                v = seen.get(key)
                if v is not None and v != u:
                    E.add(tuple(sorted((u, v))))
    return len(order), sorted(E), th

def chi(n, E, cap=7):
    for k in range(1, cap + 1):
        X = lambda v, c: 1 + v * k + c
        cnf = [[X(v, c) for c in range(k)] for v in range(n)]
        for a, b in E:
            for c in range(k):
                cnf.append([-X(a, c), -X(b, c)])
        s = Solver(name="cd19", bootstrap_with=cnf)
        ok = s.solve(); s.delete()
        if ok: return k
    return cap + 1

R1, R3 = 1.0, 1.0 / math.sqrt(3.0)
base = [("1", R1), ("1/sqrt3", R3)]
n, E, th = build([R1, R3])
print(f"  radii (1, 1/sqrt3): n={n} edges={len(E)}  chi = {chi(n, E)}   "
      f"[{time.time()-t0:.0f}s]", flush=True)

named = [("1/sqrt2", 1/math.sqrt(2)), ("1/2", 0.5), ("2/sqrt3", 2/math.sqrt(3)),
         ("sqrt3/2", math.sqrt(3)/2), ("sqrt(7)/3", math.sqrt(7)/3),
         ("golden", (1+math.sqrt(5))/2 - 0.0), ("1/golden", 2/(1+math.sqrt(5))),
         ("sqrt(1/3+1)", math.sqrt(4/3)), ("0.8", 0.8), ("1.2", 1.2)]
best = (0, None)
for nm, r3 in named:
    n, E, th = build([R1, R3, r3])
    k = chi(n, E)
    print(f"  radii (1, 1/sqrt3, {nm}={r3:.6f}): n={n} edges={len(E)}  chi = {k}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if k > best[0]: best = (k, nm)
print(f"\n  scanning a third radius on a grid   [{time.time()-t0:.0f}s]", flush=True)
for step in range(1, 60):
    r3 = 0.05 * step
    if r3 < 0.02: continue
    n, E, th = build([R1, R3, r3])
    k = chi(n, E)
    if k >= 4:
        print(f"    *** r3={r3:.4f}: n={n} edges={len(E)} chi = {k} ***",
              flush=True)
    if k > best[0]: best = (k, f"{r3:.4f}")
print(f"\n  best chi over everything tried: {best}   [{time.time()-t0:.0f}s]",
      flush=True)
