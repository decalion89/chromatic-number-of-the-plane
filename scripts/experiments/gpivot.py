"""Rotations about a pivot, not just about the origin.

Every solve so far has looked for u with |p - u.q| = 1, which is a rotation
about the ORIGIN.  de Grey's final step is not: Ya and Yb are Y turned about
(-2,0).  So a whole family has gone untested, and the fix is a translation --
a cross edge for the rotation about p0 is

    |p - (p0 + u.(q - p0))| = 1,

which is |p' - u.q'| = 1 for p' = p - p0, q' = q - p0.  The same equation on
the translated set.  So: translate G to each of several pivots, solve, and see
whether any pivot offers a rotation that bites harder than the one or two
cross edges the origin gives.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle
from fractions import Fraction as Fr
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from msqrt import madd, msub, mscal, mmul, minv, msqrt
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F
from hn.graph import build_graph

GENS = (3, 5, 7, 11)
DIM = 16
t0 = time.time()
pts = build_G(F, as_graph=False)
P0 = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c))
      for p in pts]
g = build_graph(pts)
deg = {}
for a, b in g.edges():
    deg[a] = deg.get(a, 0) + 1
    deg[b] = deg.get(b, 0) + 1
ONE = (Fr(1),) + (Fr(0),) * (DIM - 1)
print(f"G: {len(P0)} points; pivots by degree  [{time.time()-t0:.0f}s]",
      flush=True)


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


def cmul(z, w):
    return (msub(mmul(z[0], w[0], GENS), mmul(z[1], w[1], GENS)),
            madd(mmul(z[0], w[1], GENS), mmul(z[1], w[0], GENS)))


def conj(z):
    return (z[0], tuple(-c for c in z[1]))


pivots = [v for v, _ in sorted(deg.items(), key=lambda t: -t[1])[:10]]
best_overall = []
for pv in pivots:
    o = P0[pv]
    P = [(msub(p[0], o[0]), msub(p[1], o[1])) for p in P0]
    found = {}
    step = max(1, len(P) // 40)
    for qi in range(0, len(P), step):
        q = P[qi]
        A = norm2(q)
        if all(c == 0 for c in A):
            continue
        Ainv = minv(A, GENS)
        qbar = conj(q)
        for p in P:
            Pn = norm2(p)
            if all(c == 0 for c in Pn):
                continue
            R = mscal(Fr(1, 2), msub(madd(A, Pn), ONE))
            disc = msub(mmul(A, Pn, GENS), mmul(R, R, GENS))
            s = msqrt(disc, GENS)
            if s is None:
                continue
            Pinv = minv(Pn, GENS)
            for sg in (s, tuple(-c for c in s)):
                wx = mmul(msub(mmul(R, p[0], GENS), mmul(sg, p[1], GENS)),
                          Pinv, GENS)
                wy = mmul(madd(mmul(R, p[1], GENS), mmul(sg, p[0], GENS)),
                          Pinv, GENS)
                w = (wx, wy)
                if norm2(w) != A:
                    continue
                d = (msub(p[0], w[0]), msub(p[1], w[1]))
                if norm2(d) != ONE:
                    continue
                u = cmul(w, qbar)
                u = (mmul(u[0], Ainv, GENS), mmul(u[1], Ainv, GENS))
                if norm2(u) != ONE:
                    continue
                found[u] = found.get(u, 0) + 1
    # A rotation that maps the translated set to ITSELF is a symmetry: the
    # union adds no point and no edge, and the multiplicity it scores is just
    # the edges it preserves.  Those have to go before the numbers mean
    # anything.
    # A rotation that maps the translated set to ITSELF is a symmetry: the
    # union adds no point and no edge, and the multiplicity it scores is just
    # the edges it preserves.  Checking that by mapping all 1581 points, for
    # each of sixteen hundred rotations, is millions of degree-16 products --
    # so probe five points first and only confirm the ones that pass.
    Pset = set(P)
    probe = P[::max(1, len(P) // 5)][:5]
    real = {}
    for u, mult in sorted(found.items(), key=lambda t: -t[1])[:40]:
        if all(cmul(u, p) in Pset for p in probe):
            if {cmul(u, p) for p in P} == Pset:
                continue
        real[u] = mult
    sym = 40 - len(real) if len(found) >= 40 else len(found) - len(real)
    top = sorted(real.values(), reverse=True)[:6]
    best_overall.append((max(top) if top else 0, pv, deg.get(pv, 0),
                         len(real), sym))
    print(f"  pivot {pv} (degree {deg.get(pv,0)}): {len(found)} rotations, "
          f"{sym} of them symmetries, {len(real)} genuine; multiplicities "
          f"{top}  [{time.time()-t0:.0f}s]", flush=True)
    if top and top[0] >= 8:
        with open("/tmp/hn/gpivot_%d.pkl" % pv, "wb") as fh:
            pickle.dump((pv, [u for u, m in real.items()
                              if m >= top[0] // 2]), fh)
best_overall.sort(reverse=True)
print(f"\nbest multiplicity over all pivots: {best_overall[:5]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
print("  (multiplicity counts the sampled pairs a rotation serves, which is "
      "a lower bound on the cross edges it makes)", flush=True)
