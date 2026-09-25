"""The Sa -> Y step, over K.

X has no forced pair, exactly as Sa has none; de Grey's forcing appears only
after the union Y = Sa u rho(Sa), rho the rotation about the origin closing
distance 2.  K cannot close 2 (that needs sqrt-15), but it closes 1, 1/3, 3
and 7, and X realises all four.  The new one is

    rho_7 = (13 + 3 sqrt-3) / 14,      |1 - rho_7|^2 = 1/7,

so a point at distance sqrt7 from the pivot lands one away from its image.
rho_3 = RHO is the Moser rotation, already in hand.

For each (distance, pivot) the union X u rho_p(X) is built and scanned for a
pair forced monochromatic at a distance K can close.  One such pair ends it:
spindling it gives a 5-chromatic graph over K, and K's directions block at
every modulus up to five, which no multiquadratic field can do.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
src = open("/tmp/claude-0/-home-user-darwin-50/"
           "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/"
           "degreyfield.py").read()
exec(src[:src.index("t0 = time.time()\nrh = ")])
from pysat.solvers import Solver

t0 = time.time()
SQ3 = ksub(kmul(krat(2), Z6), K1)                  # zeta_6 = (1 + sqrt-3)/2
assert kmul(SQ3, SQ3) == krat(-3)


def squarefree(n):
    if n == 0:
        return 0
    s, d = (-1 if n < 0 else 1), abs(n)
    f = 2
    while f * f <= d:
        e = 0
        while d % f == 0:
            d //= f
            e += 1
        if e % 2:
            s *= f
        f += 1
    return s * d


def closable(D):
    r = Fr(1) - 4 * Fr(D)
    if r == 0:
        return False
    return squarefree(r.numerator * r.denominator) in (1, -3, -11, 33)


# the rotation that closes distance sqrt(D):  (2 - 1/D + sqrt(1-4D)/D) / 2
ROT = {Fr(3): RHO,
       Fr(7): kmul(krat(Fr(1, 14)), kadd(krat(13), kmul(krat(3), SQ3)))}
for D, r in ROT.items():
    assert knorm2(r) == K1, D
    assert knorm2(ksub(K1, r)) == krat(Fr(1, 1) / D), D
print(f"closing rotations verified for D = {sorted(map(str, ROT))}",
      flush=True)

rh = [kz(), K1, Z6, kadd(K1, Z6)]
spindle = list(rh) + [kmul(RHO, q) for q in rh[1:]]
seed, seen = [], set()
for p in spindle:
    if p not in seen:
        seen.add(p)
        seed.append(p)
X = close(seed, [r for _, r in rots])
zsX = [zof(p) for p in X]
print(f"X: {len(X)} points  [{time.time()-t0:.0f}s]", flush=True)


def edges_of(pts, zs):
    cell = {}
    for i, z in enumerate(zs):
        cell.setdefault((int(z.real // 1), int(z.imag // 1)), []).append(i)
    out = []
    for i, z in enumerate(zs):
        cx, cy = int(z.real // 1), int(z.imag // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i:
                        continue
                    if abs(abs(z - zs[j]) - 1) > 1e-7:
                        continue
                    if knorm2(ksub(pts[j], pts[i])) == K1:
                        out.append((i, j))
    return out


def scan(pts, zs, name):
    E = edges_of(pts, zs)
    n = len(pts)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name}: {n} points, {len(E)} edges, NOT 4-COLOURABLE "
              f"*** [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return "chi5", n, len(E)
    hits = []
    for i in range(n):
        for j in range(i + 1, n):
            v = abs(zs[i] - zs[j]) ** 2
            if v > 60:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable(D):
                continue
            if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)]):
                hits.append((i, j, D))
    sv.delete()
    if hits:
        print(f"  *** {name}: {n} points, {len(E)} edges, "
              f"{len(hits)} FORCED PAIRS, first at D = {hits[0][2]} "
              f"*** [{time.time()-t0:.0f}s]", flush=True)
        return "forced", hits, (pts, E)
    return None, n, len(E)


def about(rho, piv, q):
    return kadd(piv, kmul(rho, ksub(q, piv)))


for D in (Fr(7), Fr(3)):
    ring = Counter()
    for i in range(len(X)):
        for j in range(len(X)):
            if i != j and abs(abs(zsX[i] - zsX[j]) ** 2 - float(D)) < 1e-7:
                ring[i] += 1
    order = [i for i, _ in ring.most_common()]
    print(f"\nD = {D}: {len(order)} pivots carry the distance, "
          f"largest ring {ring[order[0]] if order else 0}  "
          f"[{time.time()-t0:.0f}s]", flush=True)
    for rho in (ROT[D], kconj(ROT[D])):
        for k, pi_ in enumerate(order[:60]):
            piv = X[pi_]
            pts, seen2 = list(X), set(X)
            for q in X:
                z = about(rho, piv, q)
                if z not in seen2:
                    seen2.add(z)
                    pts.append(z)
            zs = [zof(q) for q in pts]
            tag, a, b = scan(pts, zs, f"D={D} pivot {pi_}")
            if tag:
                print("  STOP", flush=True)
                sys.exit(0)
            if k % 10 == 0:
                print(f"  ... pivot {k}/60 ring {ring[pi_]}: {a} points, "
                      f"{b} edges  [{time.time()-t0:.0f}s]", flush=True)
