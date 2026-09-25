"""de Grey's Sa lives inside K -- so start from his seed, not from mine.

Write his points as complex numbers.  Each has x = p + q sqrt33 and
y = r sqrt3 + t sqrt11, so

    z = x + iy = (p + r.sqrt-3) + (t - q.sqrt-3).sqrt-11,

using i sqrt3 = sqrt-3, i sqrt11 = sqrt-11 and sqrt33 = -sqrt-3.sqrt-11.
Every one of those lies in Q(sqrt-3, sqrt-11), which is a subfield of
K = Q(m, sqrt-3, sqrt-11).  So Sa ports verbatim; it is Sb and the final pair
of rotations that do not, needing sqrt-15 and sqrt-7.

Sa is far denser than the closure I had been using -- ten neighbours a vertex
against five -- and density is what the forcing lives on.  So: take Sa, turn
it about each of its own vertices by a rotation K does have, and look for a
pair forced monochromatic at a distance K can close.
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
from hn.degrey import S_POINTS, build_Sa
from pysat.solvers import Solver

t0 = time.time()
SQ3 = ksub(kmul(krat(2), Z6), K1)
assert kmul(SQ3, SQ3) == krat(-3)
S11 = S
assert kmul(S11, S11) == krat(-11)


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


def to_k(xs, ys):
    p, q = Fr(xs.get(1, 0)), Fr(xs.get(33, 0))
    r, t = Fr(ys.get(3, 0)), Fr(ys.get(11, 0))
    re = kadd(krat(p), kmul(krat(r), SQ3))
    im = ksub(krat(t), kmul(krat(q), SQ3))
    return kadd(re, kmul(im, S11))


Sk = [to_k(x, y) for x, y in S_POINTS]
# Sa: every image under the 60-degree rotation and reflection in the x-axis
Sa, seen = [], set()
for z in Sk:
    for base in (z, kconj(z)):
        q = base
        for _ in range(6):
            if q not in seen:
                seen.add(q)
                Sa.append(q)
            q = kmul(Z6, q)
print(f"Sa over K: {len(Sa)} points (de Grey's builder gives "
      f"{len(build_Sa())})  [{time.time()-t0:.0f}s]", flush=True)
zsA = [zof(p) for p in Sa]


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


EA = edges_of(Sa, zsA)
print(f"  {len(EA)} edges, mean degree {2*len(EA)/len(Sa):.1f}  "
      f"[{time.time()-t0:.0f}s]", flush=True)

ROT = {Fr(3): RHO,
       Fr(7): kmul(krat(Fr(1, 14)), kadd(krat(13), kmul(krat(3), SQ3)))}
for D, r in ROT.items():
    assert knorm2(r) == K1 and knorm2(ksub(K1, r)) == krat(Fr(1) / D)
print(f"  closing rotations for D = 3 and 7 verified", flush=True)


def scan(pts, zs, name, quiet=True):
    E = edges_of(pts, zs)
    n = len(pts)
    cls = [[1 + v * 4 + c for c in range(4)] for v in range(n)]
    for a, b in E:
        for c in range(4):
            cls.append([-(1 + a * 4 + c), -(1 + b * 4 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name}: {n} points, {len(E)} edges, "
              f"NOT 4-COLOURABLE *** [{time.time()-t0:.0f}s]", flush=True)
        sv.delete()
        return "chi5", (pts, E)
    cand = []
    for i in range(n):
        for j in range(i + 1, n):
            v = abs(zs[i] - zs[j]) ** 2
            if v > 30:
                continue
            D = Fr(round(v * 1584), 1584)
            if abs(float(D) - v) > 1e-7 or D == 1 or not closable(D):
                continue
            cand.append((i, j, D))
    hits = [(i, j, D) for i, j, D in cand
            if not sv.solve(assumptions=[1 + i * 4, -(1 + j * 4)])]
    sv.delete()
    if hits:
        print(f"  *** {name}: {n} points, {len(E)} edges, {len(hits)} "
              f"FORCED PAIRS (D = {sorted({str(h[2]) for h in hits})}) *** "
              f"[{time.time()-t0:.0f}s]", flush=True)
        return "forced", (pts, E, hits)
    if not quiet:
        print(f"  {name}: {n} points, {len(E)} edges, {len(cand)} candidate "
              f"pairs, none forced  [{time.time()-t0:.0f}s]", flush=True)
    return None, (n, len(E), len(cand))


print("\nSa itself:", flush=True)
scan(Sa, zsA, "Sa", quiet=False)


def about(rho, piv, q):
    return kadd(piv, kmul(rho, ksub(q, piv)))


for D in (Fr(7), Fr(3)):
    ring = Counter()
    for i in range(len(Sa)):
        for j in range(len(Sa)):
            if i != j and abs(abs(zsA[i] - zsA[j]) ** 2 - float(D)) < 1e-7:
                ring[i] += 1
    order = [i for i, _ in ring.most_common()]
    print(f"\nD = {D}: largest ring {ring[order[0]] if order else 0}, "
          f"{len(order)} pivots  [{time.time()-t0:.0f}s]", flush=True)
    for rho in (ROT[D], kconj(ROT[D])):
        for k, pi_ in enumerate(order[:40]):
            piv = Sa[pi_]
            pts, seen2 = list(Sa), set(Sa)
            for q in Sa:
                z = about(rho, piv, q)
                if z not in seen2:
                    seen2.add(z)
                    pts.append(z)
            zs = [zof(q) for q in pts]
            tag, info = scan(pts, zs, f"D={D} pivot {pi_} ring {ring[pi_]}")
            if tag:
                sys.exit(0)
            if k % 5 == 0:
                print(f"  ... pivot {k}/40 (ring {ring[pi_]}): {info[0]} pts, "
                      f"{info[1]} edges, {info[2]} pairs  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
