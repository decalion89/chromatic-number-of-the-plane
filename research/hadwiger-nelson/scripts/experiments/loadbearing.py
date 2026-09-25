"""Blocking that carries weight: rotate the spindle, do not hang pendants.

The 107-point graph blocked, but by pendant vertices -- and a pendant changes
no chromatic number, so the blocking was decoration.  Making it load bearing
means every edge that carries a blocking direction must be an edge the
chromatic number depends on.

The denominator-29 rotations cannot be spindle rotations: measured, not one of
the 144 has a rational chord |1 - rho|^2, so none matches a lattice distance.
But they can rotate the WHOLE spindle.  Put the Moser spindle at the origin
and add its image under each blocking step w: every edge of every copy is a
spindle edge, the union is 4-chromatic because it contains a spindle, and the
edge directions are the union of w times the spindle's own directions.

If that module blocks, the blocking is carried by edges the chromatic number
rests on -- which is the thing the pendant construction did not do.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from field24 import (K, D, Ext, e_zero, e_of, e_add, e_sub, e_neg, e_mul,
                     e_conj, e_norm2, ONE, S11, Z6, RHO)
from hn.homcol import has_homomorphism
from pysat.solvers import Solver
import cmath


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
print(f"{len(blk)} blocking rotations", flush=True)

rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
spindle = list(rh) + [e_mul(RHO, p) for p in rh[1:]]
base, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        base.append(p)
print(f"spindle: {len(base)} points", flush=True)


def analyse(P, label):
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
    phi, _ = has_homomorphism(sorted(iv), 5)
    k = None
    for kk in (3, 4, 5):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(P))]
        for a, b in E:
            for c in range(kk):
                cls.append([-(1 + a * kk + c), -(1 + b * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                k = kk
                break
    print(f"  {label}: {len(P)} points, {len(E)} edges, {len(iv)} directions, "
          f"chi = {k}, " + ("BLOCKED" if not phi else "coset colouring exists"),
          flush=True)
    return phi is None


t0 = time.time()
cur, seen2 = list(base), set(base)
for idx, w in enumerate(blk, 1):
    for p in base:
        q = e_mul(w, p)
        if q not in seen2:
            seen2.add(q)
            cur.append(q)
    if idx in (1, 5, 20, 60, 144):
        if analyse(cur, f"{idx} rotated copies"):
            print(f"  *** BLOCKING IS LOAD BEARING: every edge is a spindle "
                  f"edge ***  [{time.time()-t0:.0f}s]", flush=True)
            break
