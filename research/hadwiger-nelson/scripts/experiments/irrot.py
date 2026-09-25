"""Rotate by something of irrational chord -- which is what the theorem asks.

The closure under the six rotations of K reached 594 directions and still did
not block at 5.  The reason is the barrier theorem, applied to my own
construction: every one of those chords is RATIONAL (1, 1/3, 3, 11/3, 11/9,
25/9), so every rotation lies in Q(sqrt-3, sqrt-11) -- degree 4, the Moser
spindle's own field, which provably never blocks.  The construction was using
a quarter of the field it was built in.

What is needed is a rotation whose chord is irrational, i.e. one that involves
the cubic part m.  Hilbert 90 supplies them: u = a/abar for any a in K.  Used
as a ROTATION, X u u.X brings in a whole new orbit of directions at once,
where a translation brought one.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import gcd
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from hn.homcol import has_homomorphism

t0 = time.time()
M = (LM, LZ, LZ, LZ)


def kinv(x):
    basis = []
    for i in range(4):
        for j in range(3):
            e = [list(LZ) for _ in range(4)]
            e[i][j] = Fr(1)
            basis.append(tuple(tuple(t) for t in e))
    cols = [kmul(x, b) for b in basis]

    def fl(v):
        out = []
        for part in v:
            out.extend(part)
        return out

    A = [[fl(cols[j])[i] for j in range(12)] + [Fr(1 if i == 0 else 0)]
         for i in range(12)]
    for c in range(12):
        p = next(t for t in range(c, 12) if A[t][c])
        A[c], A[p] = A[p], A[c]
        sc = Fr(1) / A[c][c]
        A[c] = [v * sc for v in A[c]]
        for t in range(12):
            if t != c and A[t][c]:
                f = A[t][c]
                A[t] = [u - f * v for u, v in zip(A[t], A[c])]
    co = [A[i][12] for i in range(12)]
    return tuple(tuple(co[3 * i:3 * i + 3]) for i in range(4))


def flat(x):
    out = []
    for part in x:
        out.extend(part)
    return tuple(out)


def irrational_chord(u):
    """chord = 2 - (u + ubar); irrational means it moves with m."""
    t = kadd(u, kconj(u))
    # the chord is 2 - t, and t is rational exactly when its only nonzero
    # component is the constant one of the cubic part
    return bool(t[0][1]) or bool(t[0][2])


seeds = [K1, W, S, M, kadd(K1, M), ksub(M, K1), kadd(M, W), kadd(M, S),
         kadd(kmul(krat(2), M), K1), kmul(M, M), kadd(kmul(M, M), K1),
         kadd(kmul(M, W), K1), kadd(kmul(M, S), K1)]
cands = []
for a in seeds:
    for b in seeds:
        for z in (kadd(a, b), ksub(a, b), kadd(a, kmul(krat(2), b)),
                  kadd(kmul(krat(3), a), b)):
            if z == kz():
                continue
            try:
                u = kmul(z, kinv(kconj(z)))
            except (StopIteration, ZeroDivisionError):
                continue
            if knorm2(u) != K1:
                continue
            if not irrational_chord(u):
                continue
            if u not in [c for c in cands]:
                cands.append(u)
print(f"{len(cands)} rotations of IRRATIONAL chord  [{time.time()-t0:.0f}s]",
      flush=True)

rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
X = close(seed, [r for _, r in rots])
print(f"closure over the rational chords: {len(X)} points", flush=True)

base = set()
for p in X:
    for q in X:
        if p != q and knorm2(ksub(q, p)) == K1:
            base.add(ksub(q, p))
print(f"  {len(base)} edge vectors", flush=True)


def dirs_of(vecs):
    raw = [flat(v) for v in vecs]
    dn = 1
    for v in raw:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for v in raw:
        w = tuple(int(q * dn) for q in v)
        g = 0
        for t in w:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w) if g > 1 else w)
    return sorted(out)


cur = list(base)
for k, u in enumerate(cands[:20], 1):
    cur = cur + [kmul(u, v) for v in base]
    iv = dirs_of(cur)
    blocks = [n for n in (2, 3, 4, 5) if has_homomorphism(iv, n)[0] is None]
    print(f"  after {k} irrational rotations: {len(iv)} directions, "
          f"blocks at {blocks}  [{time.time()-t0:.0f}s]", flush=True)
    if blocks == [2, 3, 4, 5]:
        print("  *** BLOCKS AT 2, 3, 4 AND 5 ***", flush=True)
        break
