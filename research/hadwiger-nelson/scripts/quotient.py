"""Periodic colourings beyond homomorphisms to Z_5.

The functional test ruled out one family and said so.  A periodic colouring of
a Cayley graph of Z^r needs only a finite-index subgroup H and a 5-colouring
of the quotient Cayley graph Cay(Z^r/H, U mod H); colouring by phi to Z_5 is
the special case where the quotient is Z_5 and the colouring is the identity.

So take the quotients coordinate-wise: reduce mod m and ask whether
Cay((Z_m)^r, U mod m) is 5-colourable.  If it is, LIFT it -- every point of
the infinite group takes the colour of its residue, adjacent points differ by
a unit vector, and a unit vector is never zero mod m as long as none of them
reduces to zero, which is checked.  That is a proof for the whole group, not a
finite sample.

If no m works, the (rho_3, rho_4) group has survived a much larger family than
the first test, and is worth the expensive experiments.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys, time, itertools
sys.path.insert(0, HN_DIR)
import numpy as np
from fractions import Fraction as Fr
from hn.field import Field
from hn.geometry import Point, rotation_joining, _rot60, Rotation
from pysat.solvers import Solver

k = 5
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


def hnf(M):
    """Row echelon over Z by Euclidean elimination: a basis of the row lattice.

    Column by column, the Euclidean algorithm between the pivot row and each
    row below drives the column to a single gcd entry.  What comes out spans
    exactly the same lattice over Z, which a Q-echelon does not.
    """
    M = [r[:] for r in M]
    cols = len(M[0])
    r = 0
    piv = []
    for c in range(cols):
        pick = next((i for i in range(r, len(M)) if M[i][c]), None)
        if pick is None:
            continue
        M[r], M[pick] = M[pick], M[r]
        for i in range(r + 1, len(M)):
            while M[i][c]:
                q = M[r][c] // M[i][c]
                M[r] = [a - q * b for a, b in zip(M[r], M[i])]
                M[r], M[i] = M[i], M[r]
        piv.append(c)
        r += 1
        if r == len(M):
            break
    return [row for row in M[:r] if any(row)], piv


def coords(U):
    """Each unit vector in integer coordinates over a Z-BASIS of the lattice.

    The first version reduced against a Q-basis and scaled by a common
    denominator.  That is not a Z-basis, and it put the vector 1 at
    coordinates (den, 0, ...), which vanishes mod every divisor of den -- why
    every modulus reported "a generator vanishes" and nothing got tested.
    """
    rows = [[Fr(c) for c in list(u.x.c) + list(u.y.c)] for u in U]
    den = 1
    for r in rows:
        for f in r:
            den = den * f.denominator // np.gcd(den, f.denominator)
    M = [[int(f * den) for f in r] for r in rows]
    basis, piv = hnf(M)

    def express(row):
        v, co = row[:], []
        for b, p in zip(basis, piv):
            if v[p] % b[p]:
                raise ValueError("not an integer combination")
            q = v[p] // b[p]
            co.append(q)
            if q:
                v = [x - q * y for x, y in zip(v, b)]
        if any(v):
            raise ValueError("residue left over")
        return co

    out = [express(r) for r in M]
    assert all(any(c) for c in out), "a generator came out as zero"
    # HNF coefficients blow up -- they overflowed int64 -- so they stay Python
    # integers here and are reduced mod m by the caller, where they fit.
    return len(basis), out


def quotient_colours(ints, r, m):
    """Is Cay((Z_m)^r, U mod m) 5-colourable?"""
    U = np.unique(np.array([[c % m for c in co] for co in ints],
                           dtype=np.int64), axis=0)
    if np.any(np.all(U == 0, axis=1)):
        return None                     # a generator dies: self-loop
    n = m ** r
    if n > 200000:
        return "too big"
    pw = np.array([m ** i for i in range(r)], dtype=np.int64)
    allv = np.array(list(itertools.product(range(m), repeat=r)),
                    dtype=np.int64)
    ids = allv @ pw
    E = set()
    for u in U:
        tgt = ((allv + u) % m) @ pw
        for a, b in zip(ids, tgt):
            if a != b:
                E.add((min(int(a), int(b)), max(int(a), int(b))))
    cls = [[1 + v * k + c for c in range(k)] for v in range(n)]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    sv = Solver(name="cd15", bootstrap_with=cls)
    ok = sv.solve()
    sv.delete()
    return (ok, n, len(E))


for rots, emax in (((3, 4), 1), ((3, 4, 7), 1), ((4,), 1)):
    U = units(rots, emax)
    r, ints = coords(U)
    print(f"\nrotations {rots}: {len(U)} unit vectors, group rank {r}"
          f"  [{time.time()-t0:.0f}s]", flush=True)
    for m in (2, 3, 4, 5, 6, 7):
        res = quotient_colours(ints, r, m)
        if res is None:
            print(f"   m={m}: a generator vanishes mod {m}, no good",
                  flush=True)
            continue
        if res == "too big":
            print(f"   m={m}: {m**r} cosets, too big  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            continue
        ok, n, e = res
        print(f"   m={m}: quotient has {n} cosets, {e} edges -> "
              f"{'5-COLOURABLE, so the whole group is' if ok else 'not 5-colourable'}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if ok:
            break
print("\nDONE", flush=True)
