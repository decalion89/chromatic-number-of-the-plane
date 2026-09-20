"""Which rotations of de Grey's own field bite G?

The step that worked at four colours was a union whose rotation BITES: Sa and
rho_4(Sa) share one point and are joined by six edges, and that is the whole
difference between no forced pair and one.  The unions of G tried earlier used
CLOSING rotations instead -- chosen because G realises the distance, not
because they produce a cross edge -- and they produced almost none (+5 on a
ring of eleven).  So the analogous experiment at five colours has not actually
been run.

Run it.  Same equation: a cross edge is p, q in G with |p - u.q| = 1, |u| = 1,
and with w = u.q, A = |q|^2, P = |p|^2 it forces

    w.pbar = R +- i.sqrt(A.P - R^2),    R = (A + P - 1)/2,

so the rotation exists over F(i) exactly when sqrt(A.P - R^2) lies in
F = Q(sqrt3, sqrt5, sqrt7, sqrt11).  That is a multiquadratic square test,
exact and recursive -- no floats.
"""
import sys, time, pickle
from fractions import Fraction as Fr
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
from msqrt import madd, msub, mscal, mmul, minv, msqrt
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as F

GENS = (3, 5, 7, 11)
D = 16
t0 = time.time()
pts = build_G(F, as_graph=False)
P = [(tuple(Fr(c) for c in p.x.c), tuple(Fr(c) for c in p.y.c)) for p in pts]
print(f"G: {len(P)} points  [{time.time()-t0:.0f}s]", flush=True)
ZERO = (Fr(0),) * D
ONE = (Fr(1),) + (Fr(0),) * (D - 1)


def norm2(z):
    return madd(mmul(z[0], z[0], GENS), mmul(z[1], z[1], GENS))


def cmul(z, w):
    return (msub(mmul(z[0], w[0], GENS), mmul(z[1], w[1], GENS)),
            madd(mmul(z[0], w[1], GENS), mmul(z[1], w[0], GENS)))


def conj(z):
    return (z[0], tuple(-c for c in z[1]))


step = max(1, len(P) // 90)
found, tried = {}, 0
for qi in range(0, len(P), step):
    q = P[qi]
    A = norm2(q)
    if all(c == 0 for c in A):
        continue
    Ainv = minv(A, GENS)
    qbar = conj(q)
    for p in P:
        Pn = norm2(p)
        if all(c == 0 for c in Pn):
            continue
        R = mscal(Fr(1, 2), msub(madd(A, Pn), ONE))
        disc = msub(mmul(A, Pn, GENS), mmul(R, R, GENS))
        tried += 1
        s = msqrt(disc, GENS)
        if s is None:
            continue
        Pinv = minv(Pn, GENS)
        for sg in (s, tuple(-c for c in s)):
            wx = mmul(msub(mmul(R, p[0], GENS), mmul(sg, p[1], GENS)),
                      Pinv, GENS)
            wy = mmul(madd(mmul(R, p[1], GENS), mmul(sg, p[0], GENS)),
                      Pinv, GENS)
            w = (wx, wy)
            if norm2(w) != A:
                continue
            d = (msub(p[0], w[0]), msub(p[1], w[1]))
            if norm2(d) != ONE:
                continue
            u = cmul(w, qbar)
            u = (mmul(u[0], Ainv, GENS), mmul(u[1], Ainv, GENS))
            if norm2(u) != ONE:
                continue
            found.setdefault(u, u)
    if (qi // step) % 15 == 0:
        print(f"  ... q {qi}/{len(P)}, {tried} pairs, {len(found)} rotations"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{len(found)} rotations of de Grey's field bite G, from {tried} "
      f"pairs  [{time.time()-t0:.0f}s]", flush=True)
with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-"
          "39c59179b415/scratchpad/gsolve.pkl", "wb") as fh:
    pickle.dump(list(found.values()), fh)
