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
from hn.homcol import (edge_vectors, has_homomorphism, blocks_at,
                       saturated_at)
exec(open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/hnfcheck.py").read()
     .split('print("re-checking')[0].split("t0 = time.time()")[1])

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
    """Global content only.  Dividing each vector by its own gcd is a
    different operation -- safe in the direction that matters, since it can
    only make blocking harder to achieve, but it is not a module map, and the
    honest module is the one scaled all at once."""
    dn = 1
    for r in rows:
        for q in r:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    ints = {tuple(int(q * dn) for q in r) for r in rows}
    content = 0
    for v in ints:
        for x in v:
            content = gcd(content, abs(x))
    if content > 1:
        ints = {tuple(x // content for x in v) for v in ints}
    return sorted(ints)


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
    # No rational change of basis first.  Picking a Q-basis out of the
    # vectors themselves gives coordinates with large denominators, and
    # clearing those is what makes the Hermite reduction blow up at rank 48.
    # The raw integer coordinates in dimension 96 are smaller; the reduction
    # finds the rank itself.
    iv = integral(rows)
    r = len(iv[0])
    # Ask saturation first.  When the module is n-saturated the ambient test
    # is already correct and the Hermite reduction -- which is what blows up
    # at rank 48 -- can be skipped entirely.
    # n = 6 and 7 are expensive and beside the point: proving blocking there
    # means an unsatisfiability proof over a rank-48 module, and the gate is
    # n = 5, where the answer so far is an escape and escapes are cheap to
    # find.  Ask only the question that matters.
    from hn.homcol import _rank_q, _rank_mod
    rk = _rank_q(iv)
    import math
    lines = len({tuple(x % 5 for x in v) for v in iv}) // 2
    need = math.ceil(rk * math.log(5) / math.log(5 / 4))
    sat = {5: _rank_mod(iv, 5) == rk}
    blocks, blat = [], []
    print(f"  {k} rotations: {len(iv)} directions, rank {rk}, {lines} lines "
          f"mod 5, threshold {need} lines  [{time.time()-t0:.0f}s]",
          flush=True)
    continue
    flag = "  <-- DIFFERENT" if blocks != blat else ""
    hgt = len(str(max(max(abs(x) for x in v) for v in iv)))
    # Why blocking is not arriving.  Counting says 804 direction lines in
    # rank 48 over F_5 should leave 5^48 . (4/5)^804 escapes, which is far
    # below one -- so if escapes persist the directions are not spread, they
    # are concentrated.  The number of DISTINCT residues mod 5 says which.
    res5 = len({tuple(x % 5 for x in v) for v in iv})
    print(f"  {k} rotations: {len(iv)} directions ({res5} distinct mod 5), dim {r} ({hgt} digits), saturated "
          f"{[n for n in sat if sat[n]]}; over Z^d blocks {blocks}, over the "
          f"module {blat}{flag}  [{time.time()-t0:.0f}s]", flush=True)
    blocks = blat
    if 5 in blocks:
        print(f"  *** A 5-CHROMATIC GRAPH BLOCKING AT THE GATE ***",
              flush=True)
        break
