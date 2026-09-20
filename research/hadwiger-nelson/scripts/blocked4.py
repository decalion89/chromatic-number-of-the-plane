"""A 4-chromatic unit-distance graph with NO coset 5-colouring.

Q(zeta_21, sqrt-11) carries the Moser spindle -- verified, 7 points, 11 edges,
chi 4 -- because it has sqrt(-11). It also contains Q(zeta_7), whose
modulus-one elements of denominator 29 admit no homomorphism to Z/5 at all.

So a graph built from BOTH should be 4-chromatic and blocked: the spindle
supplies the chromatic number, the zeta_7 steps supply the blocking, and
blocking passes upward through a field extension for free -- a phi on the
degree-24 module would restrict to one on the Q(zeta_7) submodule, and none
exists.

That would be the first graph in this package with a real chromatic number
AND no coset colouring. Everything blocked so far has been 3-chromatic, and
everything with a chromatic number has admitted a coset colouring.
"""
import sys, time, itertools
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from field24 import (K, D, Ext, e_zero, e_of, e_add, e_sub, e_neg, e_mul,
                     e_conj, e_norm2, ONE, S11, Z6, Z7, RHO)
from hn.homcol import has_homomorphism
from pysat.solvers import Solver


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


# blocking steps: modulus-one elements of Q(zeta_7) with denominator 29
blk = []
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
        blk.append(e_of(u))
print(f"{len(blk)} blocking steps of denominator 29 in the degree-24 field",
      flush=True)

rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
spindle = list(rh) + [e_mul(RHO, p) for p in rh[1:]]
seen, P = set(), []
for p in spindle:
    if p not in seen:
        seen.add(p)
        P.append(p)
print(f"spindle: {len(P)} points", flush=True)

# Blocking is inherited only if EVERY direction of the Q(zeta_7) submodule is
# present -- forty of the hundred and forty-four was not enough, and the graph
# duly admitted a coset colouring. Attach all of them, and to ONE vertex only:
# that already makes each a genuine edge vector, and keeps the point count
# down where the degree-24 arithmetic can afford it.
anchor = P[0]
for u in blk:
    for w in (e_add(anchor, u), e_sub(anchor, u)):
        if w not in seen:
            seen.add(w)
            P.append(w)
print(f"  with blocking steps attached: {len(P)} points", flush=True)

t0 = time.time()
# A float pre-filter: exact norms in a degree-24 algebra are expensive, and
# almost every pair is nowhere near distance one.
import cmath
zs = []
for q in P:
    za = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.a))
    zb = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.b))
    zs.append(za + zb * cmath.sqrt(-11))
E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
     if abs(abs(zs[i] - zs[j]) - 1) < 1e-6
     and e_norm2(e_sub(P[j], P[i])) == 1]
print(f"  {len(E)} edges  [{time.time()-t0:.0f}s]", flush=True)

# the edge-vector module, as 24 rational coordinates
vecs = set()
for i, j in E:
    d = e_sub(P[j], P[i])
    vecs.add(tuple(d.a) + tuple(d.b))
den = 1
for v in vecs:
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
iv = set()
for v in vecs:
    w = tuple(int(x * den) for x in v)
    g = 0
    for t in w:
        g = gcd(g, abs(t))
    iv.add(tuple(t // g for t in w) if g > 1 else w)
phi, why = has_homomorphism(sorted(iv), 5)
print(f"  {len(iv)} edge directions; coset 5-colouring: "
      + ("EXISTS" if phi else "NONE -- BLOCKED"), flush=True)

for k in (3, 4, 5):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, b in E:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  {k}-colourable: {r}", flush=True)
    if r:
        print(f"\n  => chi = {k}" + (" AND BLOCKED" if not phi else ""),
              flush=True)
        break
