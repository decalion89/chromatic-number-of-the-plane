"""A ball over every blocking step at once: blocked by construction, and folded.

The earlier ball failed because fourteen Z-independent steps build a tree.
The fix is not fewer steps but ALL of them: blocking and folding are both
monotone in the step set, so taking every modulus-one element of Q(zeta_7)
with denominator at most 29 gives the blocking directions AND the relations
u1 + u2 = u3 + u4 that close the short cycles.

Edges are looked up rather than searched: a pair is joined when their
difference is one of the steps, which is O(n |S|) instead of O(n^2).  That
finds a SUBGRAPH of the true unit-distance graph -- differences of denominator
29^2 are missed -- so every chromatic lower bound it yields is honest.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, json
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

F = CycloField(7)
D = F.degree
one = F.rational(1)


def inverse(a):
    rows = [list(F.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
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


def denom(u):
    d = 1
    for x in u:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


steps = set()
for coeffs in itertools.product(range(-3, 4), repeat=4):
    a = tuple(Fraction(c) for c in coeffs) + (Fraction(0), Fraction(0))
    if not any(a):
        continue
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) != one:
        continue
    d = denom(u)
    if d % 5 or d > 29:
        steps.add(u)
        steps.add(F.neg(u) if hasattr(F, "neg") else tuple(-x for x in u))
steps = sorted(s for s in steps if denom(s) <= 29 and denom(s) % 5)
print(f"{len(steps)} steps, denominators <= 29", flush=True)

pts = {F.zero()}
for rnd in range(2):
    pts = {F.add(p, s) for p in pts for s in steps} | pts
    print(f"  round {rnd+1}: {len(pts)} points", flush=True)

P = sorted(pts)
idx = {p: i for i, p in enumerate(P)}
edges = set()
for i, p in enumerate(P):
    for s in steps:
        j = idx.get(F.add(p, s))
        if j is not None and j > i:
            edges.add((i, j))
edges = sorted(edges)
print(f"graph: {len(P)} vertices, {len(edges)} edges", flush=True)

ev = {tuple(t for t in F.sub(P[j], P[i])) for i, j in edges}
den = 1
for v in ev:
    for x in v:
        den = den * x.denominator // gcd(den, x.denominator)
iv = set()
for v in ev:
    w = tuple(int(x * den) for x in v)
    g = 0
    for t in w:
        g = gcd(g, abs(t))
    iv.add(tuple(t // g for t in w) if g > 1 else w)
phi, why = has_homomorphism(sorted(iv), 5)
print(f"{len(iv)} edge vectors; coset 5-colouring: "
      + ("EXISTS" if phi else "NONE -- blocked"), flush=True)

for k in (3, 4, 5):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
    for a, b in edges:
        for c in range(k):
            cls.append([-(1 + a * k + c), -(1 + b * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        r = s.solve()
    print(f"  {k}-colourable: {r}" + ("" if r else f"   *** chi > {k} ***"),
          flush=True)
    if r:
        break
