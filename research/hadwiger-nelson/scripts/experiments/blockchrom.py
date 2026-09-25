"""A 4-chromatic graph over a field that BLOCKS: new, if it works.

de Grey's G is 5-chromatic and admits coset 5-colourings, so its colour
classes slide and every rigidity measurement reads flat.  A graph over a
blocking field has none -- which is the rigidity this whole search has been
short of.

Q(zeta_21) has both halves: zeta_6, so the Eisenstein lattice and its
triangles live there, and zeta_7, whose modulus-one elements of denominator 29
block every homomorphism to Z/5.  Blocking is inherited upward, so anything
built there blocks.

The mechanism to raise the chromatic number is already verified in this
package: the Eisenstein lattice is uniquely 3-colourable, and unioning it with
a rotated copy destroys that colouring and makes it 4-chromatic.  The rotation
does not have to be Moser's -- any modulus-one element that is not a root of
unity moves the lattice off itself.

So: an Eisenstein patch in Q(zeta_21), plus its image under a BLOCKING step,
and the chromatic number of the union.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from fractions import Fraction
from math import gcd
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

F = CycloField(21)
D = F.degree
one = F.rational(1)
z3 = F.zeta(7)
z6 = F.neg(F.mul(z3, z3))
z7 = F.zeta(3)


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


# blocking steps: modulus-one elements of Q(zeta_7), denominator 29
blk = []
for coeffs in itertools.product(range(-3, 4), repeat=4):
    if not any(coeffs):
        continue
    a, p = F.zero(), one
    for c in coeffs:
        a = F.add(a, tuple(Fraction(c) * x for x in p))
        p = F.mul(p, z7)
    try:
        u = F.mul(a, inverse(F.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if F.norm2(u) == one and denom(u) == 29:
        blk.append(u)
print(f"{len(blk)} blocking steps of denominator 29", flush=True)


def patch(r):
    pts = set()
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            v = F.zero()
            for _ in range(abs(a)):
                v = F.add(v, one if a > 0 else F.neg(one))
            for _ in range(abs(b)):
                v = F.add(v, z6 if b > 0 else F.neg(z6))
            pts.add(v)
    return sorted(pts)


def chi_of(P, hi=5):
    idx = {p: i for i, p in enumerate(P)}
    E = []
    for i, p in enumerate(P):
        for j in range(i + 1, len(P)):
            if F.norm2(F.sub(P[j], p)) == one:
                E.append((i, j))
    for k in range(2, hi + 1):
        cls = [[1 + v * k + c for c in range(k)] for v in range(len(P))]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k, len(E)
    return hi + 1, len(E)


for r in (3, 4):
    base = patch(r)
    k0, e0 = chi_of(base)
    print(f"\nEisenstein patch radius {r}: {len(base)} points, {e0} edges, "
          f"chi = {k0}", flush=True)
    for u in blk[:3]:
        both = sorted(set(base) | {F.mul(u, p) for p in base})
        t = time.time()
        k1, e1 = chi_of(both)
        print(f"  + rotated by a blocking step: {len(both)} points, {e1} "
              f"edges, chi = {k1}"
              + ("   *** 4-CHROMATIC ***" if k1 >= 4 else "")
              + f"  [{time.time()-t:.0f}s]", flush=True)
        if k1 >= 4:
            break
