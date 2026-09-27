"""Where a union puts its forced pair -- and the pruning that follows.

Y = Sa u rho(Sa) about the ORIGIN, rho closing D = 4.  Its one forced pair is
(2,0) and (-2,0): both on the ring of radius sqrt(D) = 2 about the pivot, and
antipodal on it, so at distance 2 sqrt(D) and squared distance 4D.  That is
not a coincidence of that example, it is where the rigidity is: the union
pins the ring, because the ring is exactly the set the rotation moves by one.

So the chain doubles, D -> 4D, and the step from D needs sqrt(16D - 1):

    D    =    1     4    16     64      256
    16D-1=   15    63   255   1023     4095
    rad  =   15     7   255   1023      455

de Grey has 3, 5, 7, 11, so he closes 1, 4 and 16 and stops -- 255 = 3.5.17
asks for sqrt17, which he has not got.  His construction is exactly as long as
his field allows.

The pruning: scan the ring, not the graph.  Sixty pairs an union instead of
eighty thousand.  Validated by rerunning his own step and recovering (2,0),
(-2,0), then applied at five colours to unions of G.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import isqrt
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD, Rotation, Point
from hn.homcol import closable_distance, closing_radicand
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD


def rot_closing(D):
    D = Fr(D)
    if not closable_distance(D):
        return None
    sq = closing_radicand(D)
    c = 1 - Fr(1, 2) / D
    k2 = (1 - c * c) / sq
    num, den = k2.numerator, k2.denominator
    t = isqrt(num * den)
    if t * t != num * den:
        return None
    r = Rotation(F.rational(c), F.sqrt(sq) * F.rational(Fr(t, den)))
    assert r.cos * r.cos + r.sin * r.sin == F.rational(1)
    return r


def edges_of(pts, zs):
    cell = {}
    for i, (a, b) in enumerate(zs):
        cell.setdefault((int(a // 1), int(b // 1)), []).append(i)
    out = []
    for i, (a, b) in enumerate(zs):
        cx, cy = int(a // 1), int(b // 1)
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for j in cell.get((cx + da, cy + db), ()):
                    if j <= i:
                        continue
                    if abs((a - zs[j][0]) ** 2 + (b - zs[j][1]) ** 2 - 1) > 1e-7:
                        continue
                    d = pts[i] - pts[j]
                    if d.x * d.x + d.y * d.y == F.rational(1):
                        out.append((i, j))
    return out


def ring_scan(base, pivot_pt, D, K, name):
    """Union base with its rotated copy about pivot_pt, scan the ring."""
    rot = rot_closing(D)
    turn = rot.about(pivot_pt)
    pts, seen = list(base), set(base)
    for p in base:
        q = turn(p)
        if q not in seen:
            seen.add(q)
            pts.append(q)
    zs = [(float(p.x), float(p.y)) for p in pts]
    E = edges_of(pts, zs)
    n = len(pts)
    px, py = float(pivot_pt.x), float(pivot_pt.y)
    ring = [i for i, (a, b) in enumerate(zs)
            if abs((a - px) ** 2 + (b - py) ** 2 - float(D)) < 1e-7
            and (pts[i] - pivot_pt).x ** 2 + (pts[i] - pivot_pt).y ** 2
            == F.rational(D)]
    cls = [[1 + v * K + c for c in range(K)] for v in range(n)]
    for a, b in E:
        for c in range(K):
            cls.append([-(1 + a * K + c), -(1 + b * K + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    if not sv.solve():
        print(f"  *** {name}: {n} pts, {len(E)} edges, NOT {K}-COLOURABLE ***",
              flush=True)
        sv.delete()
        return "chi", None
    hits = []
    for a in range(len(ring)):
        for b in range(a + 1, len(ring)):
            i, j = ring[a], ring[b]
            if not sv.solve(assumptions=[1 + i * K, -(1 + j * K)]):
                d = pts[i] - pts[j]
                dd = float(d.x) ** 2 + float(d.y) ** 2
                hits.append((i, j, dd))
                print(f"  *** {name}: FORCED at {K} colours, ring pair "
                      f"{i},{j}, d^2 ~ {dd:.4f} (4D = {4*float(D)})  "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    sv.delete()
    return ("forced" if hits else None), (n, len(E), len(ring), len(hits))


# --- validation: rerun de Grey's own step and recover his pair -------------
Sa = build_Sa(F)
print(f"validation: Sa {len(Sa)} points, pivot the origin, D = 4  "
      f"[{time.time()-t0:.0f}s]", flush=True)
tag, info = ring_scan(Sa, Point(F.zero(), F.zero()), Fr(4), 4, "Sa u rho_4(Sa)")
print(f"  -> {tag}, {info}  [{time.time()-t0:.0f}s]", flush=True)

# --- the real run: unions of G at five colours -----------------------------
G = build_G(F, as_graph=False)
zs0 = [(float(p.x), float(p.y)) for p in G]
print(f"\nG: {len(G)} vertices  [{time.time()-t0:.0f}s]", flush=True)
DS = [D for D in (1, 2, 3, 4, 7, 9, 14, 16, 19, 25, 34, 37, 44)
      if rot_closing(D) is not None]
for D in DS:
    fd = float(D)
    ring = Counter()
    for i, (ai, bi) in enumerate(zs0):
        for j, (aj, bj) in enumerate(zs0):
            if i != j and abs((ai - aj) ** 2 + (bi - bj) ** 2 - fd) < 1e-7:
                ring[i] += 1
    if not ring:
        continue
    order = ring.most_common(80)
    print(f"\nD = {D} (next step needs sqrt{16*D-1} -> "
          f"{closing_radicand(4*D)}): {len(ring)} pivots, largest ring "
          f"{order[0][1]}  [{time.time()-t0:.0f}s]", flush=True)
    for k, (pi_, sz) in enumerate(order):
        tag, info = ring_scan(G, G[pi_], Fr(D), 5, f"D={D} p={pi_}")
        if tag:
            print(f"  STOP: {tag}", flush=True)
            sys.exit(0)
        if k % 20 == 0:
            print(f"  ... pivot {k}/{len(order)} (ring {sz}): {info}  "
                  f"[{time.time()-t0:.0f}s]", flush=True)
