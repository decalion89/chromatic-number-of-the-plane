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
import sys, time
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.degrey import build_G
from hn.graph import build_graph
from hn.geometry import DEGREY_FIELD as F
from hn.homcol import (edge_vectors, has_homomorphism, _rank_mod)
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

# The rank saturates at 96 after five rotations, so the counting threshold
# settles: 5^96 . (4/5)^L escapes, which drops below one at L = 693 lines.
# Each rotation contributes 55 lines mod 5, so eight rotations give 495 and
# blocking cannot arrive; twelve give 715 and it can.  Generate twenty.
TS = []
for a in range(-2, 4):
    for bq in range(0, 4):
        for c in range(0, 3):
            if bq == 0 and c == 0:
                continue
            t = crat(a)
            if bq:
                t = cadd(t, cmul(crat(bq), cm()))
            if c:
                t = cadd(t, cmul(crat(c), cmul(cm(), cm())))
            TS.append(t)
seen_t = set()
uniq_t = []
for t in TS:
    key = tuple(tuple(q for q in comp.c) for comp in t)
    if key in seen_t:
        continue
    seen_t.add(key)
    uniq_t.append(t)
TS = uniq_t[:20]
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
# Order matters.  Just below the threshold the instance is
# satisfiable but the escapes are astronomically rare -- about 5^96.(4/5)^660
# = a thousand of them among 10^67 functionals -- and finding one is a needle
# hunt.  Well above it the instance is over-constrained and refutation is
# comparatively easy.  So start with margin and walk back down.
for k in [16, 14, 13, 12, 11]:
    used = PAR[:k]
    scale = crat(1)
    for N, _, _ in used:
        scale = cmul(scale, N)
    rows = []
    for vx, vy in vecs:
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
    iv = integral(rows)
    d = len(iv[0])
    # Saturation at 5, proved cheaply.  rank_5 <= rank_Q <= d always, so
    # rank_5 = d forces rank_Q = d and the module IS 5-saturated -- no
    # Hermite reduction, no rational elimination, one pass mod 5.
    r5 = _rank_mod(iv, 5)
    sat = (r5 == d)
    # phi(v) depends only on v mod 5, so fold before the solver: 1608 vectors
    # become about 715 lines, and the chain encoding is linear in both.
    small = sorted({tuple(x % 5 for x in v) for v in iv})
    zero = any(not any(v) for v in small)
    lines = (len(small) - (1 if zero else 0)) // 2
    tag = " (5-SATURATED)" if sat else " -- NOT saturated"
    import math
    # 5^d . (4/5)^L < 1  <=>  L > d . ln5 / ln(5/4) = 7.213 d
    need = int(d * math.log(5) / math.log(1.25)) + 1
    print(f"  {k} rotations: {len(iv)} directions, dim {d}, "
          f"rank mod 5 = {r5}{tag}, {lines} lines against the {need} "
          f"threshold  [{time.time()-t0:.0f}s]", flush=True)
    if not sat:
        print("  skipping: the module test would need a Hermite reduction",
              flush=True)
        continue
    phi, why = has_homomorphism(small, 5)
    if phi is None:
        print(f"  *** BLOCKS AT THE GATE with {k} rotations: {why} ***  "
              f"[{time.time()-t0:.0f}s]", flush=True)
        import pickle
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-"
                  "a506-39c59179b415/scratchpad/gate.pkl", "wb") as fh:
            pickle.dump((k, iv), fh)
        break
    print(f"  an escape exists  [{time.time()-t0:.0f}s]", flush=True)
