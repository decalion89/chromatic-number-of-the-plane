"""Force the critical core off a single copy.

The 367-point union of rotated spindles is blocked with every edge a spindle
edge, but its 4-critical core is ONE spindle: each copy is 4-chromatic alone,
so the core keeps one and drops the rest.  And a spindle can never block --
its field has degree 4, below the residue-degree bound proved here.

So remove the option.  Delete a hitting set for the spindles: at least one
vertex from every rotated copy.  No copy is then 4-chromatic on its own, and
if what is left is STILL not 3-colourable, its core has to span copies -- a
4-critical graph that is not a spindle, over a field big enough to block.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, itertools, time
from fractions import Fraction
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from field24 import (K, D, e_zero, e_of, e_add, e_sub, e_mul, e_norm2,
                     ONE, Z6, RHO)
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
        s = Fraction(1) / M[c][c]
        M[c] = [v * s for v in M[c]]
        for r in range(D):
            if r != c and M[r][c]:
                f = M[r][c]
                M[r] = [xx - f * y for xx, y in zip(M[r], M[c])]
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

rh = [e_zero(), ONE, Z6, e_add(ONE, Z6)]
base, seen = [], set()
for p in list(rh) + [e_mul(RHO, q) for q in rh[1:]]:
    if p not in seen:
        seen.add(p); base.append(p)
P, idx = list(base), {p: i for i, p in enumerate(base)}
copies = [[idx[p] for p in base]]
for w in blk:
    cp = []
    for p in base:
        q = e_mul(w, p)
        if q not in idx:
            idx[q] = len(P); P.append(q)
        cp.append(idx[q])
    copies.append(cp)
zs = []
for q in P:
    za = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.a))
    zb = sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 21)
             for i, c in enumerate(q.b))
    zs.append(za + zb * cmath.sqrt(-11))
E = [(i, j) for i in range(len(P)) for j in range(i + 1, len(P))
     if abs(abs(zs[i] - zs[j]) - 1) < 1e-7 and e_norm2(e_sub(P[j], P[i])) == 1]
print(f"{len(P)} points, {len(E)} edges, {len(copies)} spindle copies",
      flush=True)


def three_col(keep):
    ix = {v: i for i, v in enumerate(sorted(keep))}
    cls = [[1 + i * 3 + c for c in range(3)] for i in range(len(ix))]
    for a, b in E:
        if a in ix and b in ix:
            for c in range(3):
                cls.append([-(1 + ix[a] * 3 + c), -(1 + ix[b] * 3 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as s:
        return s.solve()


t0 = time.time()
alive = [set(c) for c in copies]
hit = set()
while alive:
    cnt = {}
    for c in alive:
        for v in c:
            cnt[v] = cnt.get(v, 0) + 1
    v = max(cnt, key=cnt.get)
    hit.add(v)
    alive = [c for c in alive if v not in c]
keep = set(range(len(P))) - hit
print(f"hitting set: {len(hit)} vertices removed, {len(keep)} left  "
      f"[{time.time()-t0:.0f}s]", flush=True)
col = three_col(keep)
print(f"  still 3-colourable after every spindle is broken: {col}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
if not col:
    print("  *** 4-CHROMATIC WITH NO SPINDLE IN IT ***", flush=True)
    order = sorted(keep, key=lambda v: sum(1 for a, b in E
                                           if (a == v or b == v)))
    core = set(keep)
    for v in order:
        core.discard(v)
        if three_col(core):
            core.add(v)
    vecs = []
    for a, b in E:
        if a in core and b in core:
            d = e_sub(P[b], P[a]); vecs.append(tuple(d.a) + tuple(d.b))
    den = 1
    for v in vecs:
        for x in v:
            den = den * x.denominator // gcd(den, x.denominator)
    iv = set()
    for v in vecs:
        w = tuple(int(x * den) for x in v)
        gg = 0
        for t in w:
            gg = gcd(gg, abs(t))
        iv.add(tuple(t // gg for t in w) if gg > 1 else w)
    phi, _ = has_homomorphism(sorted(iv), 5)
    print(f"  4-critical core: {len(core)} points, {len(iv)} directions, "
          + ("BLOCKED" if phi is None else "coset colouring exists")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
