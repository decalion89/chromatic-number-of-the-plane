"""A 5-chromatic graph that blocks at every modulus up to five.

G is 5-chromatic and its 133 edge directions block at 2, 3 and 4 -- and not at
5, exactly as the barrier theorem says they cannot: his field is multiquadratic,
so every prime above 5 has residue degree at most 2 and a coset colouring mod 5
always exists.

Blocking is monotone, so G u rho_p(G) blocks wherever D(G) u rho.D(G) does, and
rho only has to be a rotation of the plane living in a field where 5 has
residue degree 3.  Adjoin m, a root of x^3 - 10x^2 + 26x - 11, which is
irreducible mod 5.  Every rotation then comes free from the rational
parametrisation of the circle:

    cos = (1 - t^2)/(1 + t^2),   sin = 2t/(1 + t^2),    t in Q(m)(sqrt3,..),

and taking t with a nonzero m-part makes the chord 4t^2/(1+t^2) irrational,
which is what the barrier theorem demands.

No division is needed.  Scaling every direction by one fixed field element is a
Z-module isomorphism of the span, so it leaves blocking alone; scaling the
whole set by the product of the denominators 1 + t_i^2 clears them all at once.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.graph import build_graph
from hn.geometry import DEGREY_FIELD as F
from hn.homcol import edge_vectors, has_homomorphism

t0 = time.time()
CUB = (Fr(11), Fr(-26), Fr(10))        # m^3 = 10 m^2 - 26 m + 11


def cz():
    return (F.zero(), F.zero(), F.zero())


def crat(q):
    return (F.rational(Fr(q)), F.zero(), F.zero())


def cm():
    return (F.zero(), F.rational(1), F.zero())


def cadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def csub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cmul(a, b):
    r = [F.zero()] * 5
    for i in range(3):
        for j in range(3):
            r[i + j] = r[i + j] + a[i] * b[j]
    for d in (4, 3):                      # reduce m^4 then m^3
        c = r[d]
        if c == F.zero():
            continue
        r[d] = F.zero()
        for k in range(3):
            r[d - 3 + k] = r[d - 3 + k] + c * F.rational(CUB[k])
    return tuple(r[:3])


def coords(x):
    """96 rationals: three cubic components, each a 32-rational plane vector."""
    out = []
    for comp in x:
        out.extend(Fr(c) for c in comp.c)
    return out


G = build_graph(build_G(as_graph=False))
D = edge_vectors(G)
print(f"G: {len(list(G.edges()))} edges, {len(D)} directions, dim {len(D[0])}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# rebuild the directions as pairs of field elements, not flattened rationals
pts = build_G(F, as_graph=False)
idx = {p: i for i, p in enumerate(pts)}
vecs = set()
for a, b in G.edges():
    d = pts[b] - pts[a]
    vecs.add((d.x, d.y))
    vecs.add((-d.x, -d.y))
vecs = sorted(vecs, key=lambda v: (float(v[0]), float(v[1])))
print(f"  {len(vecs)} signed direction pairs  [{time.time()-t0:.0f}s]",
      flush=True)

TS = [cm(),                                   # t = m
      cadd(crat(1), cm()),                    # t = 1 + m
      csub(cm(), crat(1)),                    # t = m - 1
      cmul(cm(), cm()),                       # t = m^2
      cadd(crat(2), cm()),                    # t = 2 + m
      cadd(crat(1), cmul(cm(), cm())),        # t = 1 + m^2
      csub(cmul(crat(2), cm()), crat(1)),     # t = 2m - 1
      cadd(cmul(crat(3), cm()), crat(1))]     # t = 3m + 1
PAR = []
for t in TS:
    t2 = cmul(t, t)
    PAR.append((cadd(crat(1), t2), csub(crat(1), t2), cmul(crat(2), t)))
print(f"  {len(PAR)} rotations of irrational chord ready  "
      f"[{time.time()-t0:.0f}s]", flush=True)


def reduce_rank(rows):
    """Express in a basis of the span: a Z-module iso, so blocking is kept."""
    piv, basis = [], []
    for r in rows:
        w = list(r)
        for c, b in zip(piv, basis):
            if w[c]:
                f = w[c]
                w = [x - f * y for x, y in zip(w, b)]
        nz = next((i for i, x in enumerate(w) if x), None)
        if nz is not None:
            f = w[nz]
            basis.append([x / f for x in w])
            piv.append(nz)
    out = []
    for r in rows:
        w, co = list(r), []
        for c, b in zip(piv, basis):
            f = w[c]
            co.append(f)
            if f:
                w = [x - f * y for x, y in zip(w, b)]
        out.append(co)
    return out


def integral(rows):
    dn = 1
    for r in rows:
        for q in r:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for r in rows:
        w = tuple(int(q * dn) for q in r)
        g = 0
        for x in w:
            g = gcd(g, abs(x))
        out.add(tuple(x // g for x in w) if g > 1 else w)
    return sorted(out)


used = []
for k in range(1, len(PAR) + 1):
    used = PAR[:k]
    scale = crat(1)
    for N, _, _ in used:
        scale = cmul(scale, N)
    rows = []
    for vx, vy in vecs:                       # the unrotated copy, scaled
        rows.append(coords(cmul(scale, (vx, F.zero(), F.zero())))
                    + coords(cmul(scale, (vy, F.zero(), F.zero()))))
    for i, (N, P, Q) in enumerate(used):
        other = crat(1)
        for j, (Nj, _, _) in enumerate(used):
            if j != i:
                other = cmul(other, Nj)
        for vx, vy in vecs:
            X = (vx, F.zero(), F.zero())
            Y = (vy, F.zero(), F.zero())
            rx = csub(cmul(P, X), cmul(Q, Y))
            ry = cadd(cmul(Q, X), cmul(P, Y))
            rows.append(coords(cmul(other, rx)) + coords(cmul(other, ry)))
    red = reduce_rank(rows)
    r = max(len(x) for x in red)
    iv = integral([x + [Fr(0)] * (r - len(x)) for x in red])
    blocks = []
    for n in (2, 3, 4, 5):
        phi, _ = has_homomorphism(iv, n)
        if phi is None:
            blocks.append(n)
    print(f"  {k} rotations: {len(iv)} directions, rank {r}, "
          f"blocks at {blocks}  [{time.time()-t0:.0f}s]", flush=True)
    if blocks == [2, 3, 4, 5]:
        print(f"  *** A 5-CHROMATIC GRAPH BLOCKING AT 2, 3, 4 AND 5 ***",
              flush=True)
        break
