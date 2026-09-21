"""Derive the bite by hand from the pattern census, then check it with a solver.

The joint census of centre and ring says exactly what Sa claims at four
colours: ten patterns of 715, the seven points take at most two colours, and
the centre is never alone.  Sorting the ten by shape:

    1 pattern    everything one colour
    6 patterns   the minority is an ADJACENT ring pair; the centre is with the
                 other four
    3 patterns   the centre is with exactly one ANTIPODAL pair, and the other
                 four ring points carry the second colour

The bite adds six edges, h_j to rho(h_j), and fixes the centre.  So a
4-colouring of Sa u rho(Sa) picks one pattern for each copy, they agree on the
centre's colour a, and corresponding ring points differ.  Write A for the ring
points of Sa coloured a and B for those of rho(Sa) coloured a: h_j = a forces
rho(h_j) != a, so A and B are disjoint, and each is one of the shapes above --
all six, the complement of an adjacent pair, or an antipodal pair.

    A = all six              ->  B empty                      impossible
    A = complement of {j,j+1}->  B inside {j, j+1}, adjacent   impossible
    A = antipodal {j, j+3}   ->  B = {j+1,j+4} or {j+2,j+5}    the only case

And in that case the other FOUR ring points all carry the second colour, so
the other two antipodal pairs are monochromatic as well.  All three are.  That
is de Grey's forced pair, and the naming of one of them is a convenience: the
property is symmetric, as it must be, since rotating by sixty degrees carries
Sa u rho(Sa) to itself and permutes the three pairs cyclically.

Checked twice.  The elimination is finite, so it is enumerated outright; and
the conclusion is put to the solver on the union itself.
"""
import sys, time
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np, math
from itertools import product
from fractions import Fraction as Fr
from hn.degrey import build_Sa, build_Sb
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import forced_same
from pysat.solvers import Solver

t0 = time.time()
# ---- the ten patterns, as the census produced them (centre is position 0) ---
SHAPES = [
    ((0, 1, 2, 3, 4, 5, 6),),
    ((1, 2), (0, 3, 4, 5, 6)), ((2, 3), (0, 1, 4, 5, 6)),
    ((0, 1, 2, 3, 4), (5, 6)), ((3, 4), (0, 1, 2, 5, 6)),
    ((0, 1, 4), (2, 3, 5, 6)), ((0, 2, 3, 4, 5), (1, 6)),
    ((1, 2, 4, 5), (0, 3, 6)), ((4, 5), (0, 1, 2, 3, 6)),
    ((0, 2, 5), (1, 3, 4, 6)),
]
ANTI = [(1, 4), (2, 5), (3, 6)]


def a_set(shape):
    """Ring positions sharing the centre's colour."""
    blk = next(b for b in shape if 0 in b)
    return frozenset(x for x in blk if x)


print("the ten shapes, by the ring points that share the centre's colour:",
      flush=True)
for sh in SHAPES:
    A = a_set(sh)
    kind = ("all six" if len(A) == 6 else
            "antipodal pair" if tuple(sorted(A)) in ANTI else
            f"complement of {tuple(sorted(set(range(1,7)) - A))}")
    print(f"   {sh}  ->  A = {sorted(A)}  [{kind}]", flush=True)

# ---- the elimination, enumerated ------------------------------------------
ok_pairs = []
for sa in SHAPES:
    for sb in SHAPES:
        A, B = a_set(sa), a_set(sb)
        if A & B:                      # h_j = a forces rho(h_j) != a
            continue
        ok_pairs.append((sa, sb, A, B))
print(f"\n{len(ok_pairs)} of {len(SHAPES)**2} shape pairs survive the six "
      f"crossing edges  [{time.time()-t0:.0f}s]", flush=True)
allthree = True
for sa, sb, A, B in ok_pairs:
    # the other class of sa holds the remaining four ring points
    other = set(range(1, 7)) - A
    mono = all((x in A) == (y in A) for x, y in ANTI)
    print(f"   A={sorted(A)}  B={sorted(B)}   all three antipodal pairs "
          f"monochromatic in Sa: {mono}", flush=True)
    allthree &= mono
print(f"\nevery surviving case makes ALL THREE antipodal pairs "
      f"monochromatic: {allthree}", flush=True)

# ---- and the same conclusion, put to the solver on the union ---------------
k = 4
U, seen = [], set()
for p in build_Sa(K) + build_Sb(K):
    if p not in seen:
        seen.add(p)
        U.append(p)
b = IntBasis.covering(U)
r = b.rows(U)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(U)
dm, d2 = b.dim, b.D * b.D
zi = U.index(Point(K.zero(), K.zero()))
d = r - r[zi]
sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
okr = np.ones(len(sq), dtype=bool)
for j in range(1, dm):
    okr &= sq[:, j] == 0
ring = [int(o) for o in np.nonzero(okr)[0] if Fr(int(sq[o, 0]), d2) == 4]
print(f"\nSa u Sb: {n} pts, {len(E)} edges, D=4 ring has {len(ring)} points"
      f"  [{time.time()-t0:.0f}s]", flush=True)
cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
for a, c in E:
    for col in range(k):
        cls.append([-(1 + a * k + col), -(1 + c * k + col)])
sv = Solver(name="cd15", bootstrap_with=cls)
assert sv.solve(), "the union must still colour"
key = {tuple(r[i]): i for i in range(n)}
pairs = []
for i in ring:
    anti = key.get(tuple(2 * r[zi] - r[i]))
    if anti is not None and anti > i:
        pairs.append((i, anti))
print(f"{len(pairs)} antipodal pairs on the ring  [{time.time()-t0:.0f}s]",
      flush=True)
for i, j in pairs:
    f = forced_same(sv, i, j, k)
    print(f"   pair ({i},{j}) at {tuple(float(x) for x in (U[i].x, U[i].y))}"
          f": forced same = {f}  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
