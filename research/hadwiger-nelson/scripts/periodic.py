"""Why every translation-built set here colours, proved rather than sampled.

The k-core said it first: the 8-core of the best ball is denser than anything
else built here, 6.63 edges per vertex, and the 12-core is EMPTY.  No dense
core means the graph is locally tree-like, and that is a structural statement
about what was built, not about the plane.

Here is the structure.  Every point in these constructions -- sumsets, Cayley
balls, confined walks -- is an integer combination of fixed vectors, so the
point set is a chunk of a finitely generated additive group, a copy of Z^r
embedded densely in the plane.  The unit-distance graph on it is a subgraph of
the CAYLEY GRAPH of Z^r whose generators are the unit vectors of the group.

A Cayley graph of Z^r is 5-colourable as soon as there is a homomorphism
phi : Z^r -> Z_5 with phi(u) != 0 for every unit vector u: colour each point
by phi of its coordinate vector, and adjacent points differ by a unit vector,
so their colours differ.  That colours the WHOLE infinite group at once, not a
finite chunk of it, which is stronger than any SAT run can be.

Such a phi exists unless the unit vectors hit every functional.  Up to scaling
there are (5^r - 1)/4 functionals, so beating this needs at least that many
unit vectors -- 3906 of them at r = 6, against the 162 the richest generator
set here produced.  That is the quantitative reason the whole direction is
shaped wrong, and it is worth checking rather than asserting.
"""
import sys, time, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation

t0 = time.time()
K = Field((3, 5, 7, 11))
ONE = Point(K.one(), K.zero())
W = _rot60(K)
ROT = {D: rotation_joining(Fr(D), K) for D in (2, 3, 4, 7)}


def units(rots, emax):
    out, seen = [], set()
    for a in range(6):
        for es in itertools.product(*[range(-emax, emax + 1) for _ in rots]):
            p = ONE
            for _ in range(a):
                p = W(p)
            for D, e in zip(rots, es):
                r = ROT[D]
                use = r if e >= 0 else Rotation(r.cos, -r.sin)
                for _ in range(abs(e)):
                    p = use(p)
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def coords(U):
    """Each unit vector in integer coordinates over a basis of the group.

    The vectors live in a 32-dimensional Q-space (16 field coordinates each
    for x and y), so a basis is read off by Gaussian elimination over Q and
    every u is then an integer vector against it.
    """
    M = np.array([[float(c) for c in list(u.x.c) + list(u.y.c)] for u in U])
    # exact rank and basis via Fractions
    rows = [[Fr(c) for c in list(u.x.c) + list(u.y.c)] for u in U]
    basis, pivots = [], []
    for r in rows:
        v = r[:]
        for b, p in zip(basis, pivots):
            if v[p]:
                f = v[p] / b[p]
                v = [x - f * y for x, y in zip(v, b)]
        nz = [i for i, x in enumerate(v) if x]
        if nz:
            basis.append(v)
            pivots.append(nz[0])
    r = len(basis)
    out = []
    for row in rows:
        v, co = row[:], []
        for b, p in zip(basis, pivots):
            f = v[p] / b[p] if v[p] else Fr(0)
            co.append(f)
            if f:
                v = [x - f * y for x, y in zip(v, b)]
        out.append(co)
    return r, out


for rots, emax in (((4,), 1), ((4,), 2), ((3, 4), 1), ((3, 4, 7), 1)):
    U = units(rots, emax)
    r, co = coords(U)
    dens = [c for c in co]
    lcm = 1
    for c in dens:
        for f in c:
            lcm = lcm * f.denominator // np.gcd(lcm, f.denominator)
    ints = np.array([[int(f * lcm) for f in c] for c in dens], dtype=np.int64)
    print(f"\nrotations {rots} e<={emax}: {len(U)} unit vectors, group rank "
          f"{r}, functionals to check {(5**r - 1)//4}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    found = None
    for phi in itertools.product(range(5), repeat=r):
        if not any(phi):
            continue
        vals = (ints @ np.array(phi, dtype=np.int64)) % 5
        if np.all(vals != 0):
            found = phi
            break
    if found:
        print(f"   phi = {found} kills none of them -> the WHOLE group is "
              f"5-colourable, periodically  [{time.time()-t0:.0f}s]",
              flush=True)
    else:
        print(f"   *** no functional works: the unit vectors hit every one "
              f"***  [{time.time()-t0:.0f}s]", flush=True)
print("\nDONE", flush=True)
