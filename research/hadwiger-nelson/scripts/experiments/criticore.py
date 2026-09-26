"""Does the blocking survive the critical core?

"No pendants" is not yet "load bearing".  The 367-point graph is blocked and
every edge is a spindle edge, but its 4-chromaticity might still rest on one
copy while the other 143 only decorate the module.  The honest test is to
strip the graph to a 4-critical subgraph -- every remaining vertex needed for
chi = 4 -- and ask whether THAT is still blocked.  If it is, the blocking is
carried by edges the chromatic number depends on, which is the definition.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, itertools
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/hn")
from field24 import (K, D, Ext, e_zero, e_of, e_add, e_sub, e_mul,
                     e_conj, e_norm2, ONE, S11, Z6, RHO)
from hn.homcol import has_homomorphism
from pysat.solvers import Solver
import cmath

exec(open("/tmp/hn/loadbearing.py")
     .read().split("t0 = time.time()")[0].split('print(f"{len(blk)}')[0]
     + "\npass\n")

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

rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
sp = list(rh) + [e_mul(RHO, p) for p in rh[1:]]
base, seen = [], set()
for p in sp:
    if p not in seen:
        seen.add(p); base.append(p)
P, seen2 = list(base), set(base)
for w in blk:
    for p in base:
        q = e_mul(w, p)
        if q not in seen2:
            seen2.add(q); P.append(q)
zs = []
for q in P:
    za = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.a))
    zb = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.b))
    zs.append(za + zb * cmath.sqrt(-11))
E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
     if abs(abs(zs[i] - zs[j]) - 1) < 1e-6 and e_norm2(e_sub(P[j], P[i])) == 1]
print(f"{len(P)} points, {len(E)} edges", flush=True)


def three_col(keep):
    idx = {v: i for i, v in enumerate(sorted(keep))}
    cls = [[1 + i * 3 + c for c in range(3)] for i in range(len(idx))]
    for a, b in E:
        if a in keep and b in keep:
            for c in range(3):
                cls.append([-(1 + idx[a] * 3 + c), -(1 + idx[b] * 3 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


def blocked(keep):
    vecs = set()
    for a, b in E:
        if a in keep and b in keep:
            d = e_sub(P[b], P[a]); vecs.add(tuple(d.a) + tuple(d.b))
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
    return phi is None, len(iv)


keep = set(range(len(P)))
assert not three_col(keep), "not 4-chromatic?"
t0 = time.time()
deg = sorted(range(len(P)),
             key=lambda v: sum(1 for a, b in E if a == v or b == v))
for v in deg:
    keep.discard(v)
    if three_col(keep):
        keep.add(v)
n_edges = sum(1 for a, b in E if a in keep and b in keep)
bl, nd = blocked(keep)
print(f"4-critical core: {len(keep)} points, {n_edges} edges, {nd} directions",
      flush=True)
print("  core is " + ("BLOCKED -- the blocking carries chromatic weight"
                      if bl else "NOT blocked: the core has a coset colouring"),
      f" [{time.time()-t0:.0f}s]", flush=True)
