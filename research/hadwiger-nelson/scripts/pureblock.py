"""Every edge a blocking step.

The rotation construction failed the honest test: its 4-critical core is one
spindle, so the blocking was still decoration.  The core collapses because
each rotated copy is 4-chromatic by itself.

The way to stop that is to leave the cheap steps out entirely.  Build over the
denominator-29 unit steps ALONE.  Then every edge of every subgraph is a
blocking step, so whatever the critical core turns out to be, its directions
are blocking directions -- load bearing by construction, provided the core's
own direction set still blocks.

Triangles survive this: |u| = |v| = |u - v| = 1 forces u.vbar + ubar.v = 1, so
u/v is a primitive sixth root of unity.  A step set closed under zeta_6 -- and
the denominator-29 set is, since zeta_6 is a unit at 29 -- carries triangles.
"""
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver
import cmath

K = CycloField(21)
D = K.degree


def cyc_inv(a):
    rows = [list(K.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        inv = Fraction(1) / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def denom(t):
    d = 1
    for x in t:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


Z6 = K.neg(K.mul(K.zeta(7), K.zeta(7)))
raw = set()
for coeffs in itertools.product(range(-3, 4), repeat=4):
    if not any(coeffs):
        continue
    a, p = K.zero(), K.rational(1)
    for c in coeffs:
        a = K.add(a, tuple(Fraction(c) * x for x in p))
        p = K.mul(p, K.zeta(3))
    try:
        u = K.mul(a, cyc_inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) == K.rational(1) and denom(u) == 29:
        v = u
        for _ in range(6):
            raw.add(v)
            v = K.mul(v, Z6)
steps = sorted(raw)
print(f"{len(steps)} unit steps of denominator exactly 29 "
      f"(closed under zeta_6)", flush=True)


def as_int_vecs(elts):
    den = 1
    for v in elts:
        den = den * denom(v) // gcd(den, denom(v))
    out = set()
    for v in elts:
        w = tuple(int(x * den) for x in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(out)


phi, _ = has_homomorphism(as_int_vecs(steps), 5)
print("  these directions alone: "
      + ("BLOCK" if phi is None else "admit a coset 5-colouring"), flush=True)

tri = sum(1 for u in steps for v in steps
          if u < v and K.norm2(K.sub(u, v)) == K.rational(1))
print(f"  triangles through the origin: {tri} pairs at mutual distance 1",
      flush=True)

t0 = time.time()
pts, seen = [K.zero()], {K.zero()}
frontier = [K.zero()]
for depth in (1, 2):
    nxt = []
    for p in frontier:
        for s in steps:
            q = K.add(p, s)
            if q not in seen:
                seen.add(q); pts.append(q); nxt.append(q)
    frontier = nxt
    zs = [sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
              for i, c in enumerate(q)) for q in pts]
    E = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if abs(abs(zs[i] - zs[j]) - 1) < 1e-7 \
               and K.norm2(K.sub(pts[j], pts[i])) == K.rational(1):
                E.append((i, j))
    k = None
    for kk in (3, 4, 5):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(pts))]
        for a, b in E:
            for c in range(kk):
                cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                k = kk
                break
    print(f"  ball of radius {depth}: {len(pts)} points, {len(E)} edges, "
          f"chi = {k}  [{time.time()-t0:.0f}s]", flush=True)
    if k and k >= 4:
        print("  *** 4-chromatic with every edge a blocking step ***",
              flush=True)
        break
