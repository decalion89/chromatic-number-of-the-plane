"""Chains whose steps LEAVE the spindle subfield.

The first two chains found over Q(zeta_33) are 10 points, 19 edges, chi = 4 --
genuinely new 4-chromatic unit-distance graphs, three chained rhombi rather
than two paired ones.  But their 14 directions have module rank 4: the steps
fell inside Q(zeta_3, sqrt-11), the spindle'"'"'s own field, which the
residue-degree theorem says can never block.

So the search is narrowed to what matters.  A step lies in the spindle
subfield exactly when it is a rational combination of 1, zeta_3, sqrt-11 and
zeta_3.sqrt-11; requiring at least one of the two free steps to lie OUTSIDE
that span is the condition for the chain to generate a field with a chance.

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
import sys, itertools, time, cmath
from fractions import Fraction
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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


def in_spindle_field(x):
    """Is x in the Q-span of 1, zeta_3, sqrt-11, zeta_3 sqrt-11?"""
    basis = [ONE, Z3, g, K.mul(Z3, g)]
    rows = [[Fraction(c) for c in bb] + [Fraction(0)] for bb in basis]
    rows.append([Fraction(c) for c in x] + [Fraction(1)])
    M = [[rows[r][c] for r in range(5)] for c in range(D + 1)]
    rank, piv = 0, []
    for c in range(5):
        r = next((t for t in range(rank, len(M)) if M[t][c]), None)
        if r is None:
            continue
        M[rank], M[r] = M[r], M[rank]
        f = M[rank][c]
        M[rank] = [v / f for v in M[rank]]
        for t in range(len(M)):
            if t != rank and M[t][c]:
                kk = M[t][c]
                M[t] = [p - kk * q for p, q in zip(M[t], M[rank])]
        piv.append(c)
        rank += 1
    return 4 in piv or rank <= 4


outside = [i for i, u in enumerate(steps) if not in_spindle_field(u)]
print(f"{len(steps)} steps, {len(outside)} outside the spindle subfield",
      flush=True)
found = []
for ii, i in enumerate(outside):
    ci, ri = K.conj(steps[i]), re[i]
    for j in range(len(steps)):
        if j == i:
            continue
        t = K.mul(steps[j], ci)
        if K.add(K.add(ri, re[j]), K.add(t, K.conj(t))) == FOUR3:
            found.append((min(i, j), max(i, j)))
            print(f"  *** chain leaving the subfield: {i}, {j}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
    if ii % 25 == 0:
        print(f"  ... {ii}/{len(outside)}, {len(found)} so far "
              f"[{time.time()-t0:.0f}s]", flush=True)
found = sorted(set(found))
print(f"{len(found)} chains with a step outside  [{time.time()-t0:.0f}s]",
      flush=True)

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
