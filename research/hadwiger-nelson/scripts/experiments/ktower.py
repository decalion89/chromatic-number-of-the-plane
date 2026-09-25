"""Solve for the rotations of K that bite Sa, instead of sampling for them.

Sampling found none: all 3030 units of K leave Sa and its image sharing only
the origin, with zero cross edges in exact arithmetic.  That is not surprising
-- a cross edge is an exact algebraic condition, and de Grey did not stumble
on 2 arcsin(1/4) either.

Set it up properly.  A cross edge is p, q in Sa with |p - u.q| = 1 and
|u| = 1.  Put w = u.q, so |w|^2 = |q|^2 =: A and |w - p|^2 = 1 with
|p|^2 =: P.  Expanding,

    w.pbar + wbar.p = A + P - 1,     |w.pbar|^2 = A.P,

so w.pbar = R + sqrt(R^2 - A.P) with R = (A + P - 1)/2, and the whole question
is whether that square root lies in K.  Both R and A.P are real, the radicand
is negative, and K = F(sqrt-3) with F = Q(m, sqrt33) -- so

    the rotation exists over K  iff  (R^2 - A.P) / (-3)  is a square in F.

F is Q(m)(sqrt33), so x = a + b.sqrt33 is a square iff a^2 - 33b^2 is a square
in Q(m) and then (a +- that root)/2 is too.  Square roots in the cubic Q(m)
are read off its three real embeddings and confirmed exactly.

Sa is dihedrally symmetric and u.Sa = (zeta_6.u).Sa, so q may be taken one per
zeta_6-orbit.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import isqrt
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.degrey import S_POINTS


def flat(x):
    out = []
    for part in x:
        out.extend(part)
    return tuple(out)


def kinvF(x):
    """Inverse of a REAL element of K -- one F-inversion, not twelve."""
    a, b = k_to_f(x)
    d = lsub(lmul(a, a), lmul(L33, lmul(b, b)))
    di = linv(d)
    return f_to_k((lmul(a, di), lneg(lmul(b, di))))

t0 = time.time()
SQ3 = ksub(kmul(krat(2), Z6), K1)
assert kmul(SQ3, SQ3) == krat(-3)
SQ33 = kmul(krat(-1), kmul(SQ3, S))
assert kmul(SQ33, SQ33) == krat(33), "sqrt33 = -sqrt-3.sqrt-11"
assert kconj(SQ33) == SQ33, "and it is real"

# the three real embeddings of m, to fifty digits -- float has sixteen and
# the coefficients to be reconstructed are far longer than that
from decimal import Decimal as Dc, getcontext
getcontext().prec = 60
ROOTS = []
for start in ("0.5", "4.0", "9.0"):
    z = Dc(start)
    for _ in range(300):
        f = z ** 3 - 10 * z ** 2 + 26 * z - 11
        d = 3 * z ** 2 - 20 * z + 26
        z -= f / d
    ROOTS.append(z)
print(f"cubic roots: {[str(r)[:12] for r in ROOTS]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
assert all(abs(r ** 3 - 10 * r ** 2 + 26 * r - 11) < Dc(10) ** -40
           for r in ROOTS), "not all three roots are real"
MS = ROOTS


def lneg(x):
    return tuple(-c for c in x)


def lval(y, i):
    m = MS[i]
    return Dc(y[0].numerator) / Dc(y[0].denominator) + \
        Dc(y[1].numerator) / Dc(y[1].denominator) * m + \
        Dc(y[2].numerator) / Dc(y[2].denominator) * m * m


def lnorm(y):
    """N(y) as a rational -- a square in Q is necessary for y to be a square."""
    cols = [lmul(y, e) for e in (L1, LM, lmul(LM, LM))]
    a = [[cols[j][i] for j in range(3)] for i in range(3)]
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def _ratsq(q):
    n, d = q.numerator, q.denominator
    if n < 0:
        return False
    rn, rd = isqrt(n), isqrt(d)
    return rn * rn == n and rd * rd == d


def lsqrt(y):
    """A square root of y in Q(m), or None."""
    if y == LZ:
        return LZ
    if not _ratsq(lnorm(y)):
        return None
    vs = [lval(y, i) for i in range(3)]
    if any(v < 0 for v in vs):
        return None
    rs = [v.sqrt() for v in vs]
    V = [[Dc(1), MS[i], MS[i] * MS[i]] for i in range(3)]

    def det3(W):
        return (W[0][0] * (W[1][1] * W[2][2] - W[1][2] * W[2][1])
                - W[0][1] * (W[1][0] * W[2][2] - W[1][2] * W[2][0])
                + W[0][2] * (W[1][0] * W[2][1] - W[1][1] * W[2][0]))

    dv = det3(V)
    for sg in range(8):
        bb = [rs[i] * (1 if (sg >> i) & 1 == 0 else -1) for i in range(3)]
        co = []
        for c in range(3):
            W = [row[:] for row in V]
            for i in range(3):
                W[i][c] = bb[i]
            co.append(det3(W) / dv)
        cand = tuple(Fr(str(c)).limit_denominator(10 ** 20) for c in co)
        if lmul(cand, cand) == y:
            return cand
    return None


L33 = (Fr(33), Fr(0), Fr(0))
LHALF = (Fr(1, 2), Fr(0), Fr(0))
LTWO = (Fr(2), Fr(0), Fr(0))


def fadd(x, y):
    return (ladd(x[0], y[0]), ladd(x[1], y[1]))


def fsub(x, y):
    return (lsub(x[0], y[0]), lsub(x[1], y[1]))


def fmul(x, y):
    a, b = x
    c, d = y
    return (ladd(lmul(a, c), lmul(L33, lmul(b, d))),
            ladd(lmul(a, d), lmul(b, c)))


def fscal(q, x):
    r = (Fr(q), Fr(0), Fr(0))
    return (lmul(r, x[0]), lmul(r, x[1]))


def fsqrt(x):
    """A square root of a + b.sqrt33 in F = Q(m, sqrt33), or None."""
    a, b = x
    if b == LZ:
        c = lsqrt(a)
        return (c, LZ) if c is not None else None
    delta = lsub(lmul(a, a), lmul(L33, lmul(b, b)))
    sq = lsqrt(delta)
    if sq is None:
        return None
    for t in (lmul(LHALF, ladd(a, sq)), lmul(LHALF, lsub(a, sq))):
        c = lsqrt(t)
        if c is None or c == LZ:
            continue
        d = lmul(b, linv(lmul(LTWO, c)))
        if fmul((c, d), (c, d)) == x:
            return (c, d)
    return None


def k_to_f(x):
    """A real element of K is (a0, 0, c0, 2c0); it equals a0 - c0.sqrt33."""
    assert x[1] == LZ and x[3] == lmul(LTWO, x[2]), "not real"
    return (x[0], lneg(x[2]))


def f_to_k(x):
    a, b = x
    return (a, LZ, lneg(b), lneg(lmul(LTWO, b)))


for _t in (krat(7), knorm2(kadd(K1, Z6)), krat(Fr(1, 3))):
    assert f_to_k(k_to_f(_t)) == _t, "F embedding is wrong"
assert f_to_k(fmul(k_to_f(SQ33), k_to_f(SQ33))) == krat(33)
print(f"F = Q(m, sqrt33) arithmetic verified  [{time.time()-t0:.0f}s]",
      flush=True)


# --- the next floor: solve for what bites U1, not what bites Sa ----------
#
# Sa u u.Sa forces nothing, at any K-closable distance or at de Grey's own.
# But his tower is not one union either -- S, then Sa, then Y, then G -- so the
# honest continuation is to take the richest union K gives, U1 with its 156
# cross edges, and ask the same question of it.  The equation is the same; only
# the point set changes, and U1 is not zeta_6-invariant so every q counts.
import pickle
from hn.degrey import S_POINTS


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    return kadd(kadd(krat(p), kmul(krat(r), SQ3)),
                kmul(ksub(krat(t), kmul(krat(q), SQ3)), S))


Sa, seenp = [], set()
for z in [to_k(x, y) for x, y in S_POINTS]:
    for base in (z, kconj(z)):
        w = base
        for _ in range(6):
            if w not in seenp:
                seenp.add(w)
                Sa.append(w)
            w = kmul(Z6, w)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/ksolve.pkl", "rb") as fh:
    raw = pickle.load(fh)
ROT = [tuple(tuple(v[3 * i:3 * i + 3]) for i in range(4)) for v in raw]


def edges_of(pts):
    zs = [zof(q) for q in pts]
    cell = {}
    for i, z in enumerate(zs):
        cell.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)
    out = []
    for i, z in enumerate(zs):
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i or abs(abs(z - zs[j]) - 1) > 1e-6:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


best = None
for u in ROT:
    pts, seen2 = list(Sa), set(Sa)
    for q in Sa:
        z = kmul(u, q)
        if z not in seen2:
            seen2.add(z)
            pts.append(z)
    E = edges_of(pts)
    if best is None or len(E) > best[0]:
        best = (len(E), u, pts)
print(f"U1: {len(best[2])} points, {best[0]} edges  "
      f"[{time.time()-t0:.0f}s]", flush=True)
U1 = best[2]

FONE = k_to_f(K1)
found, tried = {}, 0
step = max(1, len(U1) // 120)
for qi in range(0, len(U1), step):
    q = U1[qi]
    if q == kz():
        continue
    A = k_to_f(knorm2(q))
    qbar = kconj(q)
    for p in U1:
        if p == kz():
            continue
        P = k_to_f(knorm2(p))
        R = fscal(Fr(1, 2), fsub(fadd(A, P), FONE))
        disc = fsub(fmul(R, R), fmul(A, P))
        tried += 1
        r = fsqrt(fscal(Fr(-1, 3), disc))
        if r is None:
            continue
        SQd = kmul(f_to_k(r), SQ3)
        Pk = f_to_k(P)
        for sgn in (K1, kmul(krat(-1), K1)):
            wp = kadd(f_to_k(R), kmul(sgn, SQd))
            wq = kmul(wp, kmul(p, kinvF(Pk)))
            u = kmul(wq, kmul(qbar, kinvF(f_to_k(A))))
            if knorm2(u) != K1:
                continue
            if knorm2(ksub(p, kmul(u, q))) != K1:
                continue
            found.setdefault(flat(u), u)
    if (qi // step) % 20 == 0:
        print(f"  ... q {qi}/{len(U1)}, {tried} pairs, {len(found)} rotations"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(found)} rotations bite U1, from {tried} pairs  "
      f"[{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/ktower.pkl", "wb") as fh:
    pickle.dump(([flat(q) for q in U1], [flat(u) for u in found.values()]), fh)
