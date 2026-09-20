"""Build the chain over a field chosen to block, and measure it.

The tower is K = L(W)(Y) with

    L = Q[T]/(T^3 - 9T^2 + 14T + 8),   m = T,
    W^2 = V(m) = -39m^4 + 540m^3 - 666m^2 - 324m + 297,
    Y^2 = -3,

so F = L(W) is totally real of degree 6 -- V is positive at all three
embeddings -- and K = F(Y) is its CM quadratic extension, degree 12, with
complex conjugation Y -> -Y.  The cubic was chosen so that every prime of F
above 5 has residue degree at least 3, which is exactly the condition the
residue-degree theorem requires of a graph that blocks.

Inside it the chain is explicit:

    r  = (m^2 - 6m - 3)/(m^2 + 3)        |r| <= 2 automatically
    y1 = 3 + m(r - 1)                    y1^2 = 3(4 - r^2)
    y2 = W/(m^2 + 3)                     y2^2 = 3(8 - 12r - 9r^2)
    a  = (r + y1 Y/3)/2                  |a| = 1
    b  = (-8/3 - r + y2 Y/9) / (2(1 + abar))
    zeta_6 = (1 + Y)/2

and |1 + a + b|^2 = 1/3, so the three rhombi on 1, a, b close.
"""
import sys, time, cmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.homcol import has_homomorphism
from pysat.solvers import Solver

CUB = (8, 14, -9)          # T^3 - 9T^2 + 14T + 8  (ascending constant first)


def l_zero():
    return (Fr(0),) * 3


def l_rat(q):
    return (Fr(q), Fr(0), Fr(0))


L_T = (Fr(0), Fr(1), Fr(0))


def l_add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def l_sub(x, y):
    return tuple(a - b for a, b in zip(x, y))


def l_scal(c, x):
    return tuple(Fr(c) * a for a in x)


def l_mul(x, y):
    r = [Fr(0)] * 5
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                r[i + j] += a * b
    for k in (4, 3):
        c = r[k]
        if c:
            r[k] = Fr(0)
            r[k - 3] -= c * CUB[0]
            r[k - 2] -= c * CUB[1]
            r[k - 1] -= c * CUB[2]
    return tuple(r[:3])


def l_inv(x):
    M = [[l_mul(x, (Fr(1) if k == j else Fr(0),
                    Fr(1) if k == j and False else Fr(0),
                    Fr(0)))[0] for k in range(3)] for j in range(3)]
    cols = []
    for j in range(3):
        e = tuple(Fr(1) if k == j else Fr(0) for k in range(3))
        cols.append(l_mul(x, e))
    A = [[cols[j][i] for j in range(3)] + [Fr(1) if i == 0 else Fr(0)]
         for i in range(3)]
    for c in range(3):
        p = next(r for r in range(c, 3) if A[r][c])
        A[c], A[p] = A[p], A[c]
        s = Fr(1) / A[c][c]
        A[c] = [v * s for v in A[c]]
        for r in range(3):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [u - f * v for u, v in zip(A[r], A[c])]
    return tuple(A[i][3] for i in range(3))


MM = l_mul(L_T, L_T)
V = l_add(l_add(l_scal(-39, l_mul(MM, MM)), l_scal(540, l_mul(MM, L_T))),
          l_add(l_scal(-666, MM), l_add(l_scal(-324, L_T), l_rat(297))))


def m_of(x):
    return (x, l_zero())


def m_add(x, y):
    return (l_add(x[0], y[0]), l_add(x[1], y[1]))


def m_sub(x, y):
    return (l_sub(x[0], y[0]), l_sub(x[1], y[1]))


def m_mul(x, y):
    return (l_add(l_mul(x[0], y[0]), l_mul(l_mul(x[1], y[1]), V)),
            l_add(l_mul(x[0], y[1]), l_mul(x[1], y[0])))


def m_inv(x):
    d = l_sub(l_mul(x[0], x[0]), l_mul(l_mul(x[1], x[1]), V))
    di = l_inv(d)
    return (l_mul(x[0], di), l_mul(l_scal(-1, x[1]), di))


M_ONE = m_of(l_rat(1))
M_W = (l_zero(), l_rat(1))


def k_of(x):
    return (x, (l_zero(), l_zero()))


def k_add(x, y):
    return (m_add(x[0], y[0]), m_add(x[1], y[1]))


def k_sub(x, y):
    return (m_sub(x[0], y[0]), m_sub(x[1], y[1]))


def k_mul(x, y):
    ac = m_mul(x[0], y[0])
    bd = m_mul(x[1], y[1])
    return (m_sub(ac, (l_scal(3, bd[0]), l_scal(3, bd[1]))),
            m_add(m_mul(x[0], y[1]), m_mul(x[1], y[0])))


def k_conj(x):
    return (x[0], (l_scal(-1, x[1][0]), l_scal(-1, x[1][1])))


def k_inv(x):
    d = m_add(m_mul(x[0], x[0]), (l_scal(3, m_mul(x[1], x[1])[0]),
                                  l_scal(3, m_mul(x[1], x[1])[1])))
    di = m_inv(d)
    return (m_mul(x[0], di), m_mul((l_scal(-1, x[1][0]),
                                    l_scal(-1, x[1][1])), di))


def k_norm2(x):
    return k_mul(x, k_conj(x))[0]


def k_rat(q):
    return k_of(m_of(l_rat(q)))


K_ONE = k_rat(1)
K_ZERO = k_rat(0)
K_Y = (m_of(l_zero()), M_ONE)

mm = k_of(m_of(L_T))
mm2 = k_mul(mm, mm)
den = k_add(mm2, k_rat(3))
r = k_mul(k_sub(k_sub(mm2, k_mul(k_rat(6), mm)), k_rat(3)), k_inv(den))
y1 = k_add(k_rat(3), k_mul(mm, k_sub(r, K_ONE)))
assert k_mul(y1, y1) == k_mul(k_rat(3), k_sub(k_rat(4), k_mul(r, r)))
y2 = k_mul(k_of(M_W), k_inv(den))
lhs = k_mul(y2, y2)
rhs = k_mul(k_rat(3), k_sub(k_sub(k_rat(8), k_mul(k_rat(12), r)),
                            k_mul(k_rat(9), k_mul(r, r))))
assert lhs == rhs
HALF = k_rat(Fr(1, 2))
a = k_mul(HALF, k_add(r, k_mul(y1, k_mul(k_rat(Fr(1, 3)), K_Y))))
assert k_norm2(a) == M_ONE, "a must be a unit step"
Tv = k_sub(k_rat(Fr(-8, 3)), r)
b = k_mul(k_add(Tv, k_mul(y2, k_mul(k_rat(Fr(1, 9)), K_Y))),
          k_inv(k_mul(k_rat(2), k_add(K_ONE, k_conj(a)))))
print(f"|a|^2 = 1 ? {k_norm2(a) == M_ONE};  |b|^2 = 1 ? "
      f"{k_norm2(b) == M_ONE}", flush=True)
s3 = k_add(k_add(K_ONE, a), b)
print(f"|1 + a + b|^2 = 1/3 ? "
      f"{k_norm2(s3) == m_of(l_rat(Fr(1, 3)))}", flush=True)
Z6 = k_mul(HALF, k_add(K_ONE, K_Y))
assert k_norm2(Z6) == M_ONE
ONE_PLUS = k_add(K_ONE, Z6)
assert k_norm2(ONE_PLUS) == m_of(l_rat(3))

# Every step splits.  w' + w'' = w with both of modulus one needs the roots of
# X^2 - wX + w/wbar, whose discriminant is -3 w/wbar = (Y w)^2 -- always a
# square, because wbar = 1/w.  The roots are w.zeta_6 and w.zeta_6bar, which
# is just zeta_6 + zeta_6bar = 1 seen from the other side.  So a chain of k
# steps becomes one of k+1 by replacing any step with its two halves: the sum
# is unchanged, the closing edge survives, and the graph gains three points
# and new directions every time.
import random

Z6B = k_conj(Z6)


def build(ws):
    pts, B = [K_ZERO], K_ZERO
    for w in ws:
        for off in (w, k_mul(Z6, w), k_mul(ONE_PLUS, w)):
            pts.append(k_add(B, off))
        B = k_add(B, k_mul(ONE_PLUS, w))
    uniq, seen = [], set()
    for q in pts:
        if q not in seen:
            seen.add(q)
            uniq.append(q)
    E = [(i, j) for i in range(len(uniq)) for j in range(i + 1, len(uniq))
         if k_norm2(k_sub(uniq[j], uniq[i])) == M_ONE]
    return uniq, E


def flat(x):
    out = []
    for mpart in x:
        for lpart in mpart:
            out.extend(lpart)
    return tuple(out)


def directions(uniq, E, keep=None):
    vecs = []
    for x, y in E:
        if keep is not None and (x not in keep or y not in keep):
            continue
        vecs.append(flat(k_sub(uniq[y], uniq[x])))
        vecs.append(flat(k_sub(uniq[x], uniq[y])))
    dn = 1
    for v in vecs:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    iv = set()
    for v in vecs:
        w = tuple(int(q * dn) for q in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        iv.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(iv)


def rank_of(iv):
    rows = [[Fr(x) for x in v] for v in iv]
    rk = 0
    for c in range(12):
        pr = next((t for t in range(rk, len(rows)) if rows[t][c]), None)
        if pr is None:
            continue
        rows[rk], rows[pr] = rows[pr], rows[rk]
        f = rows[rk][c]
        rows[rk] = [x / f for x in rows[rk]]
        for t in range(len(rows)):
            if t != rk and rows[t][c]:
                k2 = rows[t][c]
                rows[t] = [x - k2 * y for x, y in zip(rows[t], rows[rk])]
        rk += 1
    return rk


def chrom(uniq, E, cap=4):
    for kk in (3, 4, 5):
        cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(uniq))]
        for x, y in E:
            for c in range(kk):
                cls.append([-(1 + x * kk + c), -(1 + y * kk + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            if sv.solve():
                return kk
    return None


# Splitting a step adds no direction: the rhombus on w has edges w, zeta_6 w
# and the diagonal w(1 - zeta_6) = w.zeta_6bar, so a step contributes exactly
# its zeta_6-ORBIT -- three projective directions -- and the halves
# w.zeta_6, w.zeta_6bar lie in the same orbit.  Measured: 20 directions at
# k = 3 and still 20 at k = 32.
#
# New directions need new orbits, and the chain gives that freely.  Only the
# closing sum is constrained, so w_1, .., w_{k-2} may be ANY unit steps; the
# last two have to satisfy w + w' = z with z = s - sum of the rest, which is
# solvable in K exactly when z is a sum of two of them.  So: choose a handful
# of orbits at random, then look up whether the remainder splits.
import random

Z6B = k_conj(Z6)


def hilbert90(x):
    return k_mul(x, k_inv(k_conj(x)))


rng = random.Random(2024)
seeds = [K_ONE, a, b, Z6, k_add(K_ONE, mm), k_add(K_ONE, k_mul(mm, mm)),
         k_add(mm, K_Y), k_add(K_ONE, k_mul(k_rat(2), mm)),
         k_sub(mm, k_rat(2)), k_add(k_of(M_W), K_ONE),
         k_add(k_of(M_W), mm), k_sub(k_of(M_W), k_rat(3))]
raw = set()
for x in seeds:
    for y in seeds:
        for z0 in (k_add(x, y), k_sub(x, y), k_mul(x, y)):
            if z0 == K_ZERO:
                continue
            try:
                u = hilbert90(z0)
            except Exception:
                continue
            if k_norm2(u) == M_ONE:
                v = u
                for _ in range(6):
                    raw.add(v)
                    v = k_mul(v, Z6)
steps_all = sorted(raw)
sset = set(steps_all)
orbit = {}
for u in steps_all:
    v, best = u, None
    for _ in range(6):
        if best is None or v < best:
            best = v
        v = k_mul(v, Z6)
    orbit[u] = best
n_orb = len(set(orbit.values()))
print(f"{len(steps_all)} unit steps of the tower in {n_orb} zeta_6-orbits",
      flush=True)

# Blocking is monotone in the direction set, and the free steps of a chain
# contribute their zeta_6-orbits to it whatever the last two turn out to be.
# So the order of work is the other way round from a random completion search:
# first ask whether the ORBITS AVAILABLE IN THIS FIELD can block at all, and
# only then look for an arrangement of them that closes.  If they cannot, no
# chain over this field ever blocks and the field has to change; if they can,
# a blocking subset is the target the closing search has to hit.
# The 450 directions of this field's 75 zeta_6-orbits block, at rank 12.  So
# the field has the material; what is needed is a CHAIN that uses enough of
# it, since blocking is monotone and every step of a chain is load bearing --
# remove a rhombus and the forcing breaks.
#
# The move that adds orbits without touching the closing sum: replace a step w
# by three unit steps summing to w.  Any u leaves z = w - u to be split as
# v1 + v2, which is a hash lookup once all pairwise sums of the enumerated
# steps are tabulated.  Each success adds up to three orbits to the direction
# set, and the sum -- hence the closing edge, hence chi = 4 -- is untouched.
import time

t0 = time.time()
pairsum = {}
for i1, u in enumerate(steps_all):
    for v in steps_all:
        pairsum.setdefault(k_add(u, v), (u, v))
print(f"{len(pairsum)} distinct pairwise sums tabulated "
      f"[{time.time()-t0:.0f}s]", flush=True)


def orbits_of(ws):
    return {orbit[w] for w in ws}


def dirs_of(ws):
    raw = []
    for w in ws:
        v = w
        for _ in range(6):
            raw.append(flat(v))
            v = k_mul(v, Z6)
    dn2 = 1
    for v in raw:
        for q in v:
            dn2 = dn2 * q.denominator // gcd(dn2, q.denominator)
    out = set()
    for v in raw:
        w2 = tuple(int(q * dn2) for q in v)
        gg2 = 0
        for t in w2:
            gg2 = gcd(gg2, abs(t))
        out.add(tuple(t // gg2 for t in w2) if gg2 > 1 else w2)
    return sorted(out)


ws = [K_ONE, a, b]
assert k_norm2(k_add(k_add(ws[0], ws[1]), ws[2])) == m_of(l_rat(Fr(1, 3)))
rounds = 0
while rounds < 40:
    grew = False
    for idx in range(len(ws)):
        w = ws[idx]
        for u in steps_all:
            z = k_sub(w, u)
            hit = pairsum.get(z)
            if hit is None:
                continue
            new = {orbit[u], orbit[hit[0]], orbit[hit[1]]}
            if new <= orbits_of(ws):
                continue
            ws = ws[:idx] + [u, hit[0], hit[1]] + ws[idx + 1:]
            grew = True
            break
        if grew:
            break
    if not grew:
        print("  no replacement adds a new orbit", flush=True)
        break
    rounds += 1
    tot = k_add(k_add(ws[0], ws[1]), ws[2]) if len(ws) == 3 else None
    acc = ws[0]
    for w in ws[1:]:
        acc = k_add(acc, w)
    assert k_norm2(acc) == m_of(l_rat(Fr(1, 3))), "the chain stopped closing"
    dd = dirs_of(ws)
    phi, _ = has_homomorphism(dd, 5)
    print(f"  round {rounds}: {len(ws)} steps, {len(orbits_of(ws))} orbits, "
          f"{len(dd)} directions, rank {rank_of(dd)}, "
          + ("*** BLOCKS ***" if phi is None else "coset colouring")
          + f"  [{time.time()-t0:.0f}s]", flush=True)
    if phi is None:
        uniq, E = build(ws)
        chi = chrom(uniq, E)
        full = directions(uniq, E)
        ph2, _ = has_homomorphism(full, 5)
        print(f"  GRAPH: {len(uniq)} points, {len(E)} edges, chi = {chi}, "
              f"{len(full)} directions, "
              + ("BLOCKED" if ph2 is None else "coset colouring"), flush=True)
        break
