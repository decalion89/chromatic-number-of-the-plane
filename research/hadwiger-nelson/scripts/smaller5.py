"""Is there a 5-chromatic graph here smaller than G?

The census forces the question.  Sixty-degree rotation carries Sa u rho(Sa) to
itself and permutes its three surviving patterns cyclically, so no union that
keeps that symmetry can have one pattern, or two: it has three, or it has NONE
-- and none means the union does not 4-colour at all.  Every rotation about
the origin commutes with the sixty-degree one, so every union of rotated
copies of Sa keeps the symmetry.  The whole family therefore sits on a knife
edge: three patterns, or 5-chromatic.

G is 1581 points, so any union that stops colouring below that is a smaller
5-chromatic unit-distance graph than de Grey's, built from his own pieces.

Lean on purpose.  Only colourability is asked, nothing else, and the
satisfiable answers -- which is what "still colours" means -- are the fast
ones, so the sweep costs almost nothing until it finds something.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K, Rotation, rotation_joining
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = 4
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 1581
t0 = time.time()
Sa = build_Sa(K)
ORDER = [Fr(5, 9), Fr(5, 11), Fr(15, 16), Fr(3, 5), Fr(3, 7), Fr(15, 11),
         Fr(4), Fr(25), Fr(4, 15), Fr(3), Fr(5, 3), Fr(7), Fr(7, 3)]
rots = []
for D in ORDER:
    if not closable_distance(D):
        continue
    r = rotation_joining(D, K)
    rots.append((D, "+", r))
    rots.append((D, "-", Rotation(r.cos, -r.sin)))
print(f"Sa {len(Sa)} points, {len(rots)} rotations, cap {LIMIT} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def colours(P):
    b = IntBasis.covering(P)
    r = b.rows(P)
    if b.overflow_headroom(r) >= 1.0:
        return None
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(P)
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    ok = s.solve()
    s.delete()
    return ok, n, len(E)


# every ordered prefix of the rotation list, cumulative
best = None
for start in range(len(rots)):
    U, seen = list(Sa), set(Sa)
    chain = []
    for step in range(start, len(rots)):
        D, name, rot = rots[step]
        fresh = [q for q in (rot(p) for p in Sa) if q not in seen]
        if not fresh:
            continue
        seen.update(fresh)
        U.extend(fresh)
        chain.append(f"{D}{name}")
        if len(U) > LIMIT:
            break
        res = colours(U)
        if res is None:
            break
        ok, n, ne = res
        if not ok:
            print(f"*** {n} points, {ne} edges, NOT 4-COLOURABLE -- "
                  f"5-chromatic and smaller than G.  chain {chain} ***"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
            best = (n, list(chain))
            break
    if best:
        break
    print(f"  start {start} ({rots[start][0]}{rots[start][1]}): reached "
          f"{len(U)} points, all colour  [{time.time()-t0:.0f}s]", flush=True)
if not best:
    print(f"\nno union under {LIMIT} points stops colouring"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
