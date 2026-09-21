"""The same graph, a richer catalogue: bite Sa with rotations it did not have.

De Grey's field decides which rings are closable and therefore which bites
exist at all -- 52 of the first 279 rational rings.  Adjoining sqrt2 takes that
to 80 and sqrt2 with sqrt13 to 97.  And the seed's points lie in the base
field, so the CLOSURE IS THE SAME GRAPH: identical points, identical unit
edges, identical census of ten.  Only the operations change.

That makes the comparison exact.  Sa's complete bite family in its own field is
104 operations and the best reads three.  In the richer field it is 160, and
the 56 new ones are rings that were simply unavailable before.  If any of them
beats three, the floor was a property of the field rather than of the
construction.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np, math
from fractions import Fraction as Fr
from hn.field import Field, embed
from hn.degrey import build_Sa
from hn.geometry import DEGREY_FIELD as K0, Point, Rotation
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = 4
t0 = time.time()
GENS = tuple(int(x) for x in (sys.argv[1] if len(sys.argv) > 1
                              else "2,3,5,7,11").split(","))
K = Field(GENS)
print(f"field {GENS}, dimension {K.dim}  [{time.time()-t0:.0f}s]", flush=True)
P0 = [Point(embed(p.x, K), embed(p.y, K)) for p in build_Sa(K0)]
TEN = [
    [[0, 1, 2, 3, 4, 5, 6]],
    [[1, 2], [0, 3, 4, 5, 6]], [[2, 3], [0, 1, 4, 5, 6]],
    [[0, 1, 2, 3, 4], [5, 6]], [[3, 4], [0, 1, 2, 5, 6]],
    [[0, 1, 4], [2, 3, 5, 6]], [[0, 2, 3, 4, 5], [1, 6]],
    [[1, 2, 4, 5], [0, 3, 6]], [[4, 5], [0, 1, 2, 3, 6]],
    [[0, 2, 5], [1, 3, 4, 6]],
]
b0 = IntBasis.covering(P0)
r0 = b0.rows(P0)
dm, D2 = b0.dim, b0.D * b0.D
E0 = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b0, r0)))
ZERO = Point(K.zero(), K.zero())
zi = P0.index(ZERO)
d = r0 - r0[zi]
sq = b0._field_square(d[:, :dm]) + b0._field_square(d[:, dm:])
ok = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    ok &= sq[:, j] == 0
ring0 = [int(o) for o in np.nonzero(ok)[0] if Fr(int(sq[o, 0]), D2) == 4]
ring0.sort(key=lambda i: math.atan2(float(P0[i].y), float(P0[i].x)))
SEVEN = [P0[zi]] + [P0[i] for i in ring0]
print(f"Sa embedded: {len(P0)} points, {len(E0)} edges, ring of {len(ring0)}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def field_sqrt(v):
    """The square root of a rational in this field, or None."""
    if v <= 0:
        return None
    for r in K._prod:
        t = Fr(v) / r
        n2, d2 = t.numerator, t.denominator
        rn, rd = int(round(n2 ** .5)), int(round(d2 ** .5))
        if rn * rn == n2 and rd * rd == d2:
            s = K.rational(Fr(rn, rd))
            return s if r == 1 else K.sqrt(r) * s
    return None


def joining(D):
    """cos = 1 - 1/(2D), sin = sqrt(4D-1)/(2D), when the root exists."""
    D = Fr(D)
    s = field_sqrt(4 * D - 1)
    if s is None:
        return None
    return Rotation(K.rational(1 - Fr(1, 2) / D), s * K.rational(Fr(1, 2) / D))


Ds = sorted({Fr(x, y) for y in range(1, 10) for x in range(1, 10 * y + 1)}
            - {Fr(1)})
usable = [D for D in Ds if joining(D) is not None]
old = [D for D in usable
       if all((4 * D - 1) == 0 or True for _ in (0,))]
print(f"{len(usable)} rings closable here  [{time.time()-t0:.0f}s]", flush=True)


def census(U, label):
    b = IntBasis.covering(U)
    r = b.rows(U)
    if b.overflow_headroom(r) >= 1.0:
        return None
    E = sorted(set((min(a, c), max(a, c))
                   for a, c in fast_edges_complete(b, r)))
    n = len(U)
    W = [U.index(q) for q in SEVEN]
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, c in E:
        for col in range(k):
            cls.append([-(1 + a * k + col), -(1 + c * k + col)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    if not sv.solve():
        sv.delete()
        print(f"*** {label}: {n} points NOT 4-COLOURABLE ***", flush=True)
        return 0
    keep = sum(1 for p in TEN
               if sv.solve(assumptions=[1 + W[e] * k + bi
                                        for bi, blk in enumerate(p)
                                        for e in blk]))
    sv.delete()
    return keep


best, tried = 10, 0
for D in usable:
    rot = joining(D)
    for sgn in (+1, -1):
        f = (rot if sgn > 0 else Rotation(rot.cos, -rot.sin)).about(ZERO)
        s2 = set(P0)
        U = list(P0) + [q for q in (f(p) for p in P0) if q not in s2]
        c = census(U, f"D={D} dir {sgn}")
        tried += 1
        if c is None:
            continue
        if c < best:
            best = c
            print(f"*** D={D}, dir {sgn}: {len(U)} points, census {c} of 10 "
                  f"-- new best  [{time.time()-t0:.0f}s]", flush=True)
        if c == 0:
            sys.exit()
        if tried % 20 == 0:
            print(f"   {tried}/{2*len(usable)} bites, best {best}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tried} bites in the richer field, best census {best} of 10"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
