"""The synthesis: rotate a dense sumset patch instead of translating it.

Two directions, each with a clean reason for failing.

  Translation-built sets -- sumsets, Cayley balls -- get dense, up to 6.59
  edges per vertex against the de Grey family's 4.98, but they are chunks of a
  group, and a group is homogeneous.  The lattice version of that is provably
  Cayley-colourable; the general version keeps colouring too.
  De Grey's family is inhomogeneous and carries the weak property at four, but
  it is stuck at 4.98 per vertex and carries nothing at five.

The two failures point the same way: take the DENSITY from the sumset and the
INHOMOGENEITY from the rotational closure.  Build a dense sumset patch, then
close it under the twelve-element dihedral group about the origin.  The result
is not preserved by any translation, so no periodic colouring applies to it,
and it starts from a denser seed than S ever was.

Sa is exactly this construction with a 39-point seed.  The question is what it
gives with a seed built for density instead.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from hn.graph import build_graph
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())
ZERO = Point(K.zero(), K.zero())
W = _rot60(K)
ROT = {D: rotation_joining(Fr(D), K) for D in (2, 3, 4, 7)}


def sumset(rots, c, rmax):
    vs = [ONE, W(ONE)] + [ROT[D](ONE) for D in rots]
    out, seen = [], set()
    for a in itertools.product(range(-c, c + 1), repeat=len(vs)):
        p = ZERO
        for co, v in zip(a, vs):
            if co:
                p = Point(p.x + v.x * K.rational(co),
                          p.y + v.y * K.rational(co))
        if p in seen or float(p.norm2()) > rmax * rmax:
            continue
        seen.add(p)
        out.append(p)
    return out


def dihedral(P):
    rot = W
    out, seen = [], set()
    for refl in (False, True):
        for j in range(6):
            for p in P:
                q = p
                for _ in range(j):
                    q = rot(q)
                if refl:
                    q = Point(q.x, -q.y)
                if q not in seen:
                    seen.add(q)
                    out.append(q)
    return out


def graph(P):
    basis = IntBasis.covering(P)
    rows = basis.rows(P)
    E = sorted(set((min(a, b), max(a, b))
                   for a, b in fast_edges_complete(basis, rows)))
    return basis, rows, E


def kcolour(n, E, k):
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return ok


def weak(basis, rows, E, k, cap=12):
    """Closable distance classes, smallest first -- the gateway test."""
    Eset = set(E)
    dim, D2 = basis.dim, basis.D * basis.D
    byd = defaultdict(list)
    for i in range(len(rows) - 1):
        d = rows[i + 1:] - rows[i]
        sq = basis._field_square(d[:, :dim]) + basis._field_square(d[:, dim:])
        rat = np.ones(len(sq), dtype=bool)
        for m in range(1, dim):
            rat &= sq[:, m] == 0
        for off in np.nonzero(rat)[0]:
            j = i + 1 + int(off)
            if (i, j) not in Eset:
                byd[Fr(int(sq[off, 0]), D2)].append((i, j))
    clo = sorted((d for d in byd if closable_distance(d)),
                 key=lambda d: len(byd[d]))
    hits = []
    for d in clo[:cap]:
        if not kcolour(len(rows), E + byd[d], k):
            hits.append((d, len(byd[d])))
            print(f"      *** weak property on D={d} "
                  f"({len(byd[d])} pairs) ***", flush=True)
    return len(clo), hits


for rots, c, rmax in (((4,), 2, 3), ((3, 4), 2, 3), ((4,), 3, 4),
                      ((3, 4, 7), 2, 3)):
    S = sumset(rots, c, rmax)
    P = dihedral(S)
    if len(P) > 9000:
        print(f"\nseed {rots} c={c}: closure {len(P)} points, too big",
              flush=True)
        continue
    basis, rows, E = graph(P)
    n = len(rows)
    print(f"\nseed {rots} c={c} r<={rmax}: {len(S)} points -> closure {n}, "
          f"{len(E)} edges, {len(E)/n:.2f} per vertex"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    ok5 = kcolour(n, E, 5)
    ok4 = kcolour(n, E, 4) if ok5 else None
    print(f"   5-colourable: {ok5}; 4-colourable: {ok4}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    if not ok5:
        print("*** chi >= 6 ***", flush=True)
        sys.exit(0)
    nclo, hits = weak(basis, rows, E, 5)
    print(f"   weak property at five: {len(hits)} of {min(nclo,12)} smallest "
          f"closable classes  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
