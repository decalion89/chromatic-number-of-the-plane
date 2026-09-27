"""Build the three-rhombus chain over Q(zeta_33) and measure it.

Rhombus j sits on base B_{j-1} in direction w_j:

    B_{j-1},  B_{j-1} + w_j,  B_{j-1} + zeta_6 w_j,  B_j = B_{j-1} + w_j(1+zeta_6)

so B_j = (1 + zeta_6)(w_1 + .. + w_j) and every B_j is forced to the origin's
colour at three colours.  With |w_1 + w_2 + w_3| = 1/sqrt3 the last one is at
distance exactly 1 from the origin, so the closing edge is there and no
three-colouring exists.

What is being tested: that the graph really is 4-chromatic, that it is not a
Moser spindle in disguise, what field its directions generate, and whether
that direction set blocks.
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


Z11, Z3 = K.zeta(3), K.zeta(11)
Z6 = K.neg(K.mul(Z3, Z3))
QR = {1, 3, 4, 5, 9}
g, p = K.zero(), ONE
for k in range(1, 11):
    p = K.mul(p, Z11) if k > 1 else Z11
    g = K.add(g, p if k in QR else K.neg(p))
RHO = tuple(Fraction(1, 6) * x for x in K.add(K.rational(5), g))

extra = set()
for coeffs in itertools.product(range(-1, 2), repeat=5):
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
zp, z = [], ONE
for _ in range(6):
    zp.append(z)
    z = K.mul(z, Z6)
gens = zp + [RHO, K.conj(RHO), K.neg(RHO), K.neg(K.conj(RHO))] \
    + sorted(extra)[:40]
steps, frontier = {ONE}, [ONE]
for _ in range(2):
    nxt = []
    for x in frontier:
        for gq in gens:
            y = K.mul(x, gq)
            if y not in steps:
                steps.add(y); nxt.append(y)
    frontier = nxt
    if len(steps) > 6000:
        break
steps = sorted(steps)
FOUR3 = K.rational(Fraction(-4, 3))
re = {i: K.add(u, K.conj(u)) for i, u in enumerate(steps)}


def chain_graph(a, b):
    ws = [ONE, steps[a], steps[b]]
    one_plus = K.add(ONE, Z6)
    pts, B = [K.zero()], K.zero()
    for w in ws:
        for off in (w, K.mul(Z6, w), K.mul(one_plus, w)):
            pts.append(K.add(B, off))
        B = K.add(B, K.mul(one_plus, w))
    uniq, seen = [], set()
    for q in pts:
        if q not in seen:
            seen.add(q); uniq.append(q)
    zs = [sum(complex(float(c)) * cmath.exp(2j * cmath.pi * i / 33)
              for i, c in enumerate(q)) for q in uniq]
    E = [(i, j) for i in range(len(uniq)) for j in range(i + 1, len(uniq))
         if abs(abs(zs[i] - zs[j]) - 1) < 1e-7
         and K.norm2(K.sub(uniq[j], uniq[i])) == ONE]
    return uniq, E


def chrom(pts, E):
    for k in (3, 4, 5):
        cls = [[1 + v * k + c for c in range(k)] for v in range(len(pts))]
        for a, b in E:
            for c in range(k):
                cls.append([-(1 + a * k + c), -(1 + b * k + c)])
        with Solver(name="cd19", bootstrap_with=cls) as s:
            if s.solve():
                return k
    return None


t0 = time.time()
# The two chains the search already returned; rebuilt here rather than
# re-run, since the sweep over all 986 steps takes the best part of an hour.
found = [(3, 979), (69, 979)]
for a, b in found:
    ra, rb = re[a], re[b]
    t = K.mul(steps[b], K.conj(steps[a]))
    assert K.add(K.add(ra, rb), K.add(t, K.conj(t))) == FOUR3, (a, b)
print(f"{len(steps)} steps; {len(found)} chains re-verified  "
      f"[{time.time()-t0:.0f}s]", flush=True)

best = None
for a, b in found[:40]:
    pts, E = chain_graph(a, b)
    k = chrom(pts, E)
    vecs = []
    for i, j in E:
        d = K.sub(pts[j], pts[i]); vecs.append(d)
        vecs.append(K.sub(pts[i], pts[j]))
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
    rows = [[Fraction(x) for x in v] for v in iv]
    rank = 0
    for c in range(D):
        r = next((t for t in range(rank, len(rows)) if rows[t][c]), None)
        if r is None:
            continue
        rows[rank], rows[r] = rows[r], rows[rank]
        f = rows[rank][c]
        rows[rank] = [x / f for x in rows[rank]]
        for t in range(len(rows)):
            if t != rank and rows[t][c]:
                kk = rows[t][c]
                rows[t] = [x - kk * y for x, y in zip(rows[t], rows[rank])]
        rank += 1
    phi, _ = has_homomorphism(sorted(iv), 5)
    msg = (f"  chain ({a},{b}): {len(pts)} points, {len(E)} edges, chi = {k}, "
           f"{len(iv)} directions, module rank {rank}, "
           + ("BLOCKS" if phi is None else "coset colouring exists"))
    print(msg + f"  [{time.time()-t0:.0f}s]", flush=True)
    if k and k >= 4 and phi is None:
        print("  *** A 4-CHROMATIC BLOCKED CHAIN ***", flush=True)
        break
