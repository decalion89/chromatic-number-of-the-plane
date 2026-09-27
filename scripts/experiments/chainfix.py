"""The closing condition, with the target corrected -- and a valuation filter.

|1 + a + b|^2 expands to 3 + 2[Re(a) + Re(b) + Re(a.bbar)], so closing at 1/3
needs Re(a) + Re(b) + Re(a.bbar) = -4/3, i.e. the quantity the code actually
forms,

    (a + abar) + (b + bbar) + (a.bbar + abar.b)  =  -8/3.

The earlier sweep compared it to -4/3 and so searched for |1 + a + b|^2 = 5/3,
a distance of sqrt5 between the chain's ends rather than 1.  The k = 2 control
was right -- 2 + (a + abar) = 1/3 really is a + abar = -5/3 -- which is
exactly why the fault survived: the control it was there to catch passed.
Recorded rather than quietly fixed.

The valuation filter comes from the arithmetic at 3.  Writing m, n for the
valuations of a, b at a prime q above 3, unit steps give v_qbar = -v_q, so

    v_q(x) + v_qbar(x) = min(0,m,n) + min(0,-m,-n) = -2,

and the solutions include (m,n) = (0, +-2) and (+-2, 0).  Over Q(zeta_33) the
Moser rotation has valuation 2 and the zeta_6 and Q(zeta_11) units have 0, so
pairing one of each is where a chain can be.  That cuts the sweep by a factor
of three and aims it at the only place it can succeed.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, time, cmath
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.cyclotomic import CycloField
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

K = CycloField(33)
D = K.degree
ONE = K.rational(1)


def inv(a):
    rows = [list(K.mul(a, tuple(Fraction(1 if j == i else 0) for j in range(D))))
            for i in range(D)]
    M = [[rows[j][i] for j in range(D)] + [Fraction(1 if i == 0 else 0)]
         for i in range(D)]
    for c in range(D):
        p = next(r for r in range(c, D) if M[r][c])
        M[c], M[p] = M[p], M[c]
        s = Fraction(1) / M[c][c]
        M[c] = [v * s for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return tuple(M[i][D] for i in range(D))


def den(t):
    d = 1
    for x in t:
        d = d * x.denominator // gcd(d, x.denominator)
    return d


Z11, Z3 = K.zeta(3), K.zeta(11)
Z6 = K.neg(K.mul(Z3, Z3))
g, p = K.zero(), ONE
for k in range(1, 11):
    p = K.mul(p, Z11) if k > 1 else Z11
    g = K.add(g, p if k in {1, 3, 4, 5, 9} else K.neg(p))
assert K.mul(g, g) == K.rational(-11)
RHO = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))
assert K.norm2(RHO) == ONE

t0 = time.time()
extra = set()
for coeffs in itertools.product(range(-1, 2), repeat=6):
    if not any(coeffs):
        continue
    a, q = K.zero(), ONE
    for c in coeffs:
        if c:
            a = K.add(a, tuple(Fraction(c) * x for x in q))
        q = K.mul(q, Z11)
    try:
        u = K.mul(a, inv(K.conj(a)))
    except (StopIteration, ZeroDivisionError):
        continue
    if K.norm2(u) == ONE:
        extra.add(u)
print(f"{len(extra)} unit steps from Q(zeta_11)  [{time.time()-t0:.0f}s]",
      flush=True)

zp, z = [], ONE
for _ in range(6):
    zp.append(z)
    z = K.mul(z, Z6)
V0 = set()
for c in zp:
    for e in {ONE} | extra:
        V0.add(K.mul(c, e))
V2 = set()
for c in V0:
    V2.add(K.mul(c, RHO))
    V2.add(K.mul(c, K.conj(RHO)))
V0 = sorted(V0)
V2 = sorted(V2)
print(f"valuation 0: {len(V0)} steps; valuation +-2: {len(V2)} steps  "
      f"[{time.time()-t0:.0f}s]", flush=True)

TWO = K.rational(Fraction(-5, 3))
EIGHT3 = K.rational(Fraction(-8, 3))
ctrl = [u for u in V2 if K.add(u, K.conj(u)) == TWO]
print(f"k = 2 control: {len(ctrl)} Moser rotations among the valuation-2 steps",
      flush=True)

BASIS = [ONE, Z3, g, K.mul(Z3, g)]


def rank_of(vs):
    rows = [[Fraction(c) for c in v] for v in vs]
    r = 0
    for c in range(D):
        pr = next((t for t in range(r, len(rows)) if rows[t][c]), None)
        if pr is None:
            continue
        rows[r], rows[pr] = rows[pr], rows[r]
        f = rows[r][c]
        rows[r] = [x / f for x in rows[r]]
        for t in range(len(rows)):
            if t != r and rows[t][c]:
                k2 = rows[t][c]
                rows[t] = [x - k2 * y for x, y in zip(rows[t], rows[r])]
        r += 1
    return r


BR = rank_of(BASIS)
found = []
for i, a in enumerate(V0):
    w = K.add(ONE, K.conj(a))
    ra = K.add(a, K.conj(a))
    target = K.sub(EIGHT3, ra)
    for b in V2:
        t = K.mul(b, w)
        if K.add(t, K.conj(t)) == target:
            big = rank_of(BASIS + [a]) > BR or rank_of(BASIS + [b]) > BR
            found.append((a, b, big))
            print(f"  *** CHAIN (target -8/3): step {i}, "
                  + ("OUTSIDE the spindle subfield" if big else
                     "inside the spindle subfield")
                  + f"  [{time.time()-t0:.0f}s]", flush=True)
    if i % 50 == 0:
        print(f"  ... {i}/{len(V0)}, {len(found)} chains "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(found)} chains, {sum(1 for f in found if f[2])} leaving the "
      f"subfield  [{time.time()-t0:.0f}s]", flush=True)
import pickle
with open("/tmp/hn/chains.pkl",
          "wb") as fh:
    pickle.dump(found, fh)
