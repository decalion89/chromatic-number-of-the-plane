"""More hexagons: every orientation then activates several classes at once.

An auxiliary sits between two hexagons and is confined exactly when their
orientation bits disagree with its neighbours' parities, so with h hexagons
each orientation activates one parity class per PAIR -- C(h,2) classes at
once, not one. Two hexagons gave the densest class 18 points and six edges,
far too thin for anything 4-chromatic; three and four give the confined set
room.

The angles are still solved rather than scanned, from
cos(alpha - beta - phi) = (1 - r^2 - s^2)/(2 r s) with r, s the moduli of
differences within a hexagon, which are 1, sqrt(3) or 2. Ten values survive
per pair, and combinations of those are the search space.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, math, itertools
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]

TWO = 2 * math.pi
# de Grey's own offsets, read off the witness: 0, th/2, th, 60-th, 60-th/2
# with th = 2 arcsin(1/(2 sqrt 3)) = 33.5573 deg, the chord-1 angle on the
# circle of radius sqrt(3). The half-angles come from edges BETWEEN different
# hexagon pairs, which do not exist when there are only two hexagons -- which
# is exactly why the two-hexagon solve missed them.
TH = 2 * math.asin(1 / (2 * math.sqrt(3)))
ANGLES = [0.0, TH / 2, TH, math.radians(60) - TH, math.radians(60) - TH / 2,
          math.radians(28.95502), math.radians(31.04498)]
H = int(sys.argv[1]) if len(sys.argv) > 1 else 3


def pts_of(phis):
    hexs = [[(math.cos(p + TWO * k / 6), math.sin(p + TWO * k / 6))
             for k in range(6)] for p in phis]
    aux = []          # (point, hexagon pair, parity of index sum)
    for a, b in itertools.combinations(range(len(phis)), 2):
        for i, u in enumerate(hexs[a]):
            for j, v in enumerate(hexs[b]):
                q = (u[0] + v[0], u[1] + v[1])
                if q[0] ** 2 + q[1] ** 2 > 1e-12:
                    aux.append((q, (a, b), (i + j) % 2))
    return aux


def chrom(P):
    n = len(P)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            d = (P[a][0] - P[b][0]) ** 2 + (P[a][1] - P[b][1]) ** 2
            if abs(d - 1.0) < 1e-7:
                adj[a].add(b)
                adj[b].add(a)
    m = sum(len(x) for x in adj) // 2
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    pos = {v: i for i, v in enumerate(order)}
    nb = [[pos[w] for w in adj[order[i]]] for i in range(n)]
    for chi in range(1, 6):
        col = [-1] * n

        def go(v):
            if v == n:
                return True
            used = {col[w] for w in nb[v] if w < v and col[w] >= 0}
            for c in range(chi):
                if c not in used:
                    col[v] = c
                    if go(v + 1):
                        return True
                    col[v] = -1
            return False
        if go(0):
            return chi, n, m
    return 6, n, m


t0, best = time.time(), None
for combo in itertools.product(ANGLES, repeat=H - 1):
    phis = (0.0,) + combo
    aux = pts_of(phis)
    worst = None
    for bits in range(1 << H):
        active = [q for q, (a, b), par in aux
                  if par != (((bits >> a) & 1) ^ ((bits >> b) & 1))]
        c = chrom(active)
        if worst is None or c[0] < worst[0]:
            worst = c
    if best is None or worst[0] > best[0][0]:
        best = (worst, phis)
        print(f"  angles {[round(math.degrees(p),4) for p in phis]}: worst "
              f"orientation gives chi {worst[0]} over {worst[1]} points, "
              f"{worst[2]} edges  [{time.time()-t0:.0f}s]", flush=True)
print(f"H={H} best worst-case chi: {best[0][0]} at angles "
      f"{[round(math.degrees(p),4) for p in best[1]]}  "
      f"(need 4)  [{time.time()-t0:.0f}s]", flush=True)
