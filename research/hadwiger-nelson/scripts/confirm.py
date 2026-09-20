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
import sys, time
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
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

# Blocking depends only on the DIRECTIONS, so there is no need to carry the
# 597-point closure, whose points accumulate hundred-digit denominators under
# repeated multiplication.  The spindle's own directions -- the two rhombi and
# the closing edge -- have small ones, and rotating those is the same
# question at a thousandth of the cost.  Any graph realising the resulting
# set, for instance the spindle together with its rotated copies, blocks
# whenever the set does.
rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
base = []
for p in seed:
    for q in seed:
        if p != q and knorm2(ksub(q, p)) == K1:
            base.append(ksub(q, p))
print(f"spindle: {len(seed)} points, {len(base)} edge vectors", flush=True)

MOD = 60


def den_of(v):
    d = 1
    for q in v:
        d = d * q.denominator // gcd(d, q.denominator)
    return d


def scaled(v, d):
    return tuple((q.numerator * (d // q.denominator)) % MOD for q in v)


BD = den_of([q for v in base for q in flat(v)])
D0 = sorted({t for t in (scaled(flat(v), BD) for v in base) if any(t)})
print(f"  base scaled by {BD}: {len(D0)} directions mod {MOD}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

UD, kept = {}, 0
for i2, u in enumerate(cands):
    du = den_of(flat(u))
    if gcd(du, MOD) != 1:
        continue
    kept += 1
    orb = set()
    ok = True
    for v in base:
        w = flat(kmul(u, v))
        dw = den_of(w)
        if dw % BD:
            dw = dw * BD // gcd(dw, BD)
        if gcd(dw // gcd(dw, BD), MOD) != 1:
            ok = False
            break
        orb.add(scaled(w, dw))
    if ok and orb:
        UD[i2] = sorted(t for t in orb if any(t))
print(f"  {kept} rotations with denominator prime to {MOD}, {len(UD)} usable "
      f"orbits  [{time.time()-t0:.0f}s]", flush=True)

# Confirm without the modular shortcut: build the actual graph -- the spindle
# together with its images under the chosen rotations -- and test its exact
# direction set at every modulus.  The mod-60 reduction is sound for these
# four, but after two scaling faults it is worth spending the time.
chosen, cur = [], set(D0)
for rnd in range(1, 200):
    iv = sorted(cur)
    phi, _ = has_homomorphism(iv, 5)
    mod = 5
    if phi is None:
        blocks = [n for n in (2, 3, 4, 5)
                  if has_homomorphism(iv, n)[0] is None]
        if blocks == [2, 3, 4, 5]:
            print(f"  mod-60 screen: {len(chosen)} rotations, {len(iv)} "
                  f"directions, blocks at 2,3,4,5  [{time.time()-t0:.0f}s]",
                  flush=True)
            break
        mod = next(n for n in (2, 3, 4) if n not in blocks)
        phi, _ = has_homomorphism(iv, mod)
    best, bk = None, -1
    for i2, orb in UD.items():
        if i2 in chosen:
            continue
        k = sum(1 for d in orb
                if sum(a * b for a, b in zip(phi, d)) % mod == 0)
        if k > bk:
            best, bk = i2, k
    if best is None or bk <= 0:
        print("  stuck", flush=True)
        break
    chosen.append(best)
    cur |= set(UD[best])

pts, seenp = list(seed), set(seed)
for i2 in chosen:
    u = cands[i2]
    for p in list(seed):
        q = kmul(u, p)
        if q not in seenp:
            seenp.add(q)
            pts.append(q)
print(f"\ngraph: spindle plus {len(chosen)} rotated copies -> {len(pts)} "
      f"points  [{time.time()-t0:.0f}s]", flush=True)
E = [(i2, j2) for i2 in range(len(pts)) for j2 in range(i2 + 1, len(pts))
     if knorm2(ksub(pts[j2], pts[i2])) == K1]
chi = None
for k in (3, 4, 5, 6):
    cls = [[1 + v * k + c for c in range(k)] for v in range(len(pts))]
    for a2, b2 in E:
        for c in range(k):
            cls.append([-(1 + a2 * k + c), -(1 + b2 * k + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if sv.solve():
            chi = k
            break
vecs = []
for a2, b2 in E:
    vecs.append(flat(ksub(pts[b2], pts[a2])))
    vecs.append(flat(ksub(pts[a2], pts[b2])))
DN = 1
for v in vecs:
    for q in v:
        DN = DN * q.denominator // gcd(DN, q.denominator)
ivs = set()
for v in vecs:
    w = tuple(int(q * DN) for q in v)
    g = 0
    for t in w:
        g = gcd(g, abs(t))
    ivs.add(tuple(t // g for t in w) if g > 1 else w)
ivs = sorted(ivs)
exact = [n for n in (2, 3, 4, 5) if has_homomorphism(ivs, n)[0] is None]
print(f"  {len(E)} edges, chi = {chi}, {len(ivs)} exact directions, "
      f"blocks at {exact}  [{time.time()-t0:.0f}s]", flush=True)
