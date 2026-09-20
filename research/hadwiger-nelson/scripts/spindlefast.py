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
# The necklace over K.  The chain of length two -- w_1 = 1, w_2 = -rho -- has
# |w_1 + w_2|^2 = |1 - rho|^2 = 1/3, so it closes: it IS the Moser spindle,
# 7 points and 11 edges, which is 3k+1 and 5k+1 at k = 2.  The shortest
# necklace there is, and 4-critical by the necklace theorem.
#
# Grow it with the same move as before: replace a step w by three unit steps
# summing to w, which leaves the closing sum untouched, and reject any
# replacement that lets a shorter necklace close inside.  Stop when the
# directions block at 2, 3, 4 and 5.
#
# Unlike the tower, this field carries the Moser rotation, so de Grey's step
# from four colours to five is available over it.
import pickle

MRHO = kmul(krat(-1), RHO)
assert knorm2(kadd(K1, MRHO)) == krat(Fr(1, 3)), "the 2-chain must close"

steps = set()
for c, r in rots:
    v = r
    for _ in range(6):
        steps.add(v)
        v = kmul(v, Z6)
for u in cands:
    v = u
    for _ in range(6):
        steps.add(v)
        v = kmul(v, Z6)
for u in list(steps):
    steps.add(kconj(u))
    steps.add(kmul(krat(-1), u))
# 3030 steps give nine million pairwise sums, which is past what memory will
# hold.  The replacement only needs enough choice to keep growing, so the
# table is built on a slice -- the rest stay available as the single step u.
steps = sorted(steps)
sset = set(steps)
TABLE = steps[::max(1, len(steps) // 420)]
print(f"{len(steps)} unit steps in K, {len(TABLE)} in the sum table  "
      f"[{time.time()-t0:.0f}s]", flush=True)

pairsum = {}
for u in TABLE:
    for v in TABLE:
        pairsum.setdefault(kadd(u, v), (u, v))
print(f"{len(pairsum)} pairwise sums  [{time.time()-t0:.0f}s]", flush=True)

ONE_PLUS = kadd(K1, Z6)
assert knorm2(ONE_PLUS) == krat(3)


def parts_of(seq):
    out, acc = [kz()], kz()
    for w in seq:
        acc = kadd(acc, kmul(ONE_PLUS, w))
        out.append(acc)
    return out


def dirs_of(seq):
    raw = []
    for w in seq:
        v = w
        for _ in range(6):
            raw.append(flat(v))
            v = kmul(v, Z6)
    dn = 1
    for v in raw:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for v in raw:
        w2 = tuple(int(q * dn) for q in v)
        g = 0
        for t in w2:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w2) if g > 1 else w2)
    return sorted(out)


# A step contributes exactly its zeta_6-orbit, so a replacement that stays
# inside the orbits already present adds no direction at all -- measured:
# the necklace grew to k = 76 with the direction count stuck at 18.  Require
# a new orbit.
def orbit_of(w):
    best, v = None, w
    for _ in range(6):
        if best is None or v < best:
            best = v
        v = kmul(v, Z6)
    return best


def orbits_of(seq):
    return {orbit_of(w) for w in seq}


ws = [K1, MRHO]
PART = parts_of(ws)
rounds = 0
while rounds < 80:
    dd = dirs_of(ws)
    blocks = [n for n in (2, 3, 4, 5)
              if has_homomorphism(dd, n)[0] is None]
    print(f"  k = {len(ws):3d}: {len(orbits_of(ws)):3d} orbits, {len(dd):4d} "
          f"directions, blocks at {blocks}  [{time.time()-t0:.0f}s]",
          flush=True)
    if blocks == [2, 3, 4, 5]:
        print("  *** NECKLACE BLOCKING AT 2, 3, 4 AND 5 OVER K ***",
              flush=True)
        with open("/tmp/claude-0/-home-user-darwin-50/aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/kneck.pkl", "wb") as fh:
            pickle.dump(ws, fh)
        break
    grew = False
    for idx in range(len(ws)):
        w = ws[idx]
        for u in steps:
            hit = pairsum.get(ksub(w, u))
            if hit is None:
                continue
            fresh = {orbit_of(u), orbit_of(hit[0]), orbit_of(hit[1])}
            if fresh <= orbits_of(ws):
                continue
            cand = ws[:idx] + [u, hit[0], hit[1]] + ws[idx + 1:]
            n1 = kadd(PART[idx], kmul(ONE_PLUS, u))
            n2 = kadd(n1, kmul(ONE_PLUS, hit[0]))
            bad = False
            for q in PART:
                for nn in (n1, n2):
                    if knorm2(ksub(nn, q)) == K1:
                        bad = True
                        break
                if bad:
                    break
            if not bad and knorm2(ksub(n2, n1)) == K1:
                bad = True
            if bad:
                continue
            ws = cand
            grew = True
            break
        if grew:
            break
    if not grew:
        print("  no replacement left", flush=True)
        break
    rounds += 1
    PART = parts_of(ws)
    acc = ws[0]
    for w in ws[1:]:
        acc = kadd(acc, w)
    assert knorm2(acc) == krat(Fr(1, 3)), "the necklace stopped closing"


# Build it and check every claim: the counts 3k+1 and 5k+1 exactly, which is
# the necklace theorem's hypothesis; the chromatic number; the 4-critical
# core; and the blocking of the graph's own directions at all four moduli.
pts, B = [kz()], kz()
for w in ws:
    for off in (w, kmul(Z6, w), kmul(ONE_PLUS, w)):
        pts.append(kadd(B, off))
    B = kadd(B, kmul(ONE_PLUS, w))
uniq, seenq = [], set()
for q in pts:
    if q not in seenq:
        seenq.add(q)
        uniq.append(q)
E = [(i4, j4) for i4 in range(len(uniq)) for j4 in range(i4 + 1, len(uniq))
     if knorm2(ksub(uniq[j4], uniq[i4])) == K1]
k = len(ws)
print(f"\nnecklace of {k} rhombi: {len(uniq)} points (3k+1 = {3*k+1}), "
      f"{len(E)} edges (5k+1 = {5*k+1})  [{time.time()-t0:.0f}s]", flush=True)
chi = None
for kk in (3, 4, 5):
    cls = [[1 + v * kk + c for c in range(kk)] for v in range(len(uniq))]
    for a4, b4 in E:
        for c in range(kk):
            cls.append([-(1 + a4 * kk + c), -(1 + b4 * kk + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        if sv.solve():
            chi = kk
            break


def exact_of(keep):
    vv = []
    for a4, b4 in E:
        if keep is None or (a4 in keep and b4 in keep):
            vv.append(flat(ksub(uniq[b4], uniq[a4])))
            vv.append(flat(ksub(uniq[a4], uniq[b4])))
    dn = 1
    for v in vv:
        for q in v:
            dn = dn * q.denominator // gcd(dn, q.denominator)
    out = set()
    for v in vv:
        w2 = tuple(int(q * dn) for q in v)
        g = 0
        for t in w2:
            g = gcd(g, abs(t))
        out.add(tuple(t // g for t in w2) if g > 1 else w2)
    return sorted(out)


full = exact_of(None)
fb = [n for n in (2, 3, 4, 5) if has_homomorphism(full, n)[0] is None]
print(f"  chi = {chi}, {len(full)} directions, blocks at {fb}  "
      f"[{time.time()-t0:.0f}s]", flush=True)


# de Grey's move: Sa does not force at four colours, Sa u Sb does, where Sb
# is Sa rotated about a point IN the graph.  The necklace has the structure
# for it -- consecutive tips are at sqrt3, since B_j - B_{j-1} = w(1 + zeta_6)
# -- and K carries the Moser rotation, which sends a point at sqrt3 to
# distance exactly 1.
#
# So: for every vertex p and every rotation the field provides, build
# N u rho_p(N) and measure the chromatic number.  A 5 anywhere is the jump,
# and it would be a 5-chromatic graph blocking at 2, 3, 4 and 5 -- which
# de Grey's G, being multiquadratic, can never be.
def about(rho, piv, q):
    return kadd(piv, kmul(rho, ksub(q, piv)))


allrots = [r for _, r in rots] + [RHO, kconj(RHO)]
print(f"\ntrying {len(uniq)} centres x {len(allrots)} rotations  "
      f"[{time.time()-t0:.0f}s]", flush=True)
best = None
for ci, piv in enumerate(uniq):
    for ri, rho in enumerate(allrots):
        pts2 = list(uniq)
        seen2 = set(uniq)
        for q in uniq:
            z = about(rho, piv, q)
            if z not in seen2:
                seen2.add(z)
                pts2.append(z)
        if len(pts2) == len(uniq):
            continue
        # 133 centres times 8 rotations, each an O(n^2) exact pass over
        # 266 points, is 37 million degree-12 norms.  Float first, exact only
        # on the survivors -- the same prefilter the graph builder uses.
        zs2 = [zof(q) for q in pts2]
        cell = {}
        for i6, z in enumerate(zs2):
            cell.setdefault((int(z.real // 1), int(z.imag // 1)),
                            []).append(i6)
        E2 = []
        for i6, z in enumerate(zs2):
            cx, cy = int(z.real // 1), int(z.imag // 1)
            for da in (-1, 0, 1):
                for db in (-1, 0, 1):
                    for j6 in cell.get((cx + da, cy + db), ()):
                        if j6 <= i6:
                            continue
                        if abs(abs(z - zs2[j6]) - 1) > 1e-7:
                            continue
                        if knorm2(ksub(pts2[j6], pts2[i6])) == K1:
                            E2.append((i6, j6))
        cls = [[1 + v * 4 + c for c in range(4)] for v in range(len(pts2))]
        for a6, b6 in E2:
            for c in range(4):
                cls.append([-(1 + a6 * 4 + c), -(1 + b6 * 4 + c)])
        with Solver(name="cd19", bootstrap_with=cls) as sv:
            four = sv.solve()
        if not four:
            print(f"  *** chi >= 5 at centre {ci}, rotation {ri}: "
                  f"{len(pts2)} points, {len(E2)} edges  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            best = (ci, ri, pts2, E2)
            break
        if best is None or len(E2) > best[0]:
            best = (len(E2), ci, ri)
    if best is not None and len(best) == 4:
        break
    if ci % 10 == 0:
        print(f"  ... centre {ci}/{len(uniq)}  [{time.time()-t0:.0f}s]",
              flush=True)
if best is None or len(best) != 4:
    print(f"  no union reached five colours  [{time.time()-t0:.0f}s]",
          flush=True)
