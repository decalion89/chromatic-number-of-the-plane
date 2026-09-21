"""Sumsets of unit vectors: dense by construction, and not de Grey's seed.

The bitten lattices came out at 2.5 to 3.1 edges per vertex against Sa's
4.98, which is why they coloured instantly -- the rotated copies land in fresh
territory and add almost nothing.  Density is the missing ingredient, and
there is a standard way to get it that nothing here has used.

Take unit vectors v_1 .. v_m and form all integer combinations sum a_i v_i
with |a_i| <= c.  Every point has an edge to every point differing by one
generator, so the degree is 2m before any coincidence, and every ALGEBRAIC
RELATION between the generators folds the set onto itself and adds more.  With
generic generators the result is a hypercube and bipartite; with related ones
it is not, and de Grey's seed is itself an instance of this shape.

The generators are taken from his field so the arithmetic stays exact: the
sixth root of unity, and the biting rotations rho_D for D = 3, 4, 7, which are
cos 5/6 + i sqrt11/6, cos 7/8 + i sqrt15/8 and cos 13/14 + i 3sqrt3/14.  All
have infinite order, so the point set has many scales -- the property the
lattices lacked.

Reported per configuration: points, edges, edges per vertex, and whether it
colours with five.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from fractions import Fraction as Fr
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60
from hn.graph import build_graph
from pysat.solvers import Solver

k = 5
t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())


def unit(rot):
    return rot(ONE)


GENS = {
    "w":  unit(_rot60(K)),
    "r3": unit(rotation_joining(Fr(3), K)),
    "r4": unit(rotation_joining(Fr(4), K)),
    "r7": unit(rotation_joining(Fr(7), K)),
}
BASE = Point(K.one(), K.zero())


def sumset(names, c, rmax):
    vs = [BASE] + [GENS[n] for n in names]
    out, seen = [], set()
    for a in itertools.product(range(-c, c + 1), repeat=len(vs)):
        p = Point(K.zero(), K.zero())
        for coef, v in zip(a, vs):
            if coef:
                p = Point(p.x + v.x * K.rational(coef),
                          p.y + v.y * K.rational(coef))
        if p in seen:
            continue
        if float(p.norm2()) > rmax * rmax:
            continue
        seen.add(p)
        out.append(p)
    return out


def solve(P):
    g = build_graph(P)
    E = sorted(set((min(a, b), max(a, b)) for a, b in g.edges()))
    cls = [[1 + v * k + c for c in range(k)] for v in range(g.n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return g.n, len(E), ok


CONFIGS = [(("w",), 3, 6), (("w",), 4, 8),
           (("w", "r4"), 2, 5), (("w", "r4"), 3, 6),
           (("w", "r3"), 3, 6), (("w", "r7"), 3, 6),
           (("r3", "r4"), 3, 6), (("w", "r3", "r4"), 2, 5),
           (("w", "r4", "r7"), 2, 5), (("w", "r3", "r4", "r7"), 2, 5),
           (("w", "r3", "r4"), 3, 6)]
for names, c, rmax in CONFIGS:
    P = sumset(names, c, rmax)
    if len(P) < 20 or len(P) > 12000:
        print(f"  gens {('1',)+names} c={c} r<={rmax}: {len(P)} points "
              f"(skipped)  [{time.time()-t0:.0f}s]", flush=True)
        continue
    n, e, ok = solve(P)
    print(f"  gens {('1',)+names} c={c} r<={rmax}: {n} pts, {e} edges, "
          f"{e/n:.2f} per vertex -> "
          f"{'5-colourable' if ok else '*** NOT 5-COLOURABLE ***'}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok:
        print("*** chi >= 6 ***", flush=True)
        break
print("DONE", flush=True)
