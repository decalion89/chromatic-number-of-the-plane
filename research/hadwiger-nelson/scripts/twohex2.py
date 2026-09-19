"""Solve for the angle instead of scanning it.

The auxiliary one away from circle points u and v is the pivot reflected
across the chord, which is just u + v. So with two hexagons on the pivot's
unit circle the whole auxiliary set is a MINKOWSKI SUM H1 + H2, and two
auxiliaries are one apart exactly when

    |(u - u') + (v - v')| = 1.

Differences within a hexagon have modulus 0, 1, sqrt(3) or 2, so with r and s
those moduli and alpha, beta their directions the condition is

    cos(alpha - beta - phi) = (1 - r^2 - s^2) / (2 r s),

one equation with a discrete solution set. Scanning phi finds none of it --
the good angles are measure zero. Solving gives them all.

Within a parity class i+j is fixed mod 2, so the grid edges (u'=u or v'=v)
never appear there; the class's edges come entirely from these solved angles.
"""
import sys, time, math, collections, itertools
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
from hn.coloring import is_k_colorable
from hn.graph import build_graph

TWO = 2 * math.pi


def hexes(phi):
    H1 = [(math.cos(TWO * k / 6), math.sin(TWO * k / 6)) for k in range(6)]
    H2 = [(math.cos(phi + TWO * k / 6), math.sin(phi + TWO * k / 6))
          for k in range(6)]
    return H1, H2


# every phi at which some pair of auxiliaries becomes one apart
cands = set()
for i, ip in itertools.product(range(6), repeat=2):
    du = (math.cos(TWO * i / 6) - math.cos(TWO * ip / 6),
          math.sin(TWO * i / 6) - math.sin(TWO * ip / 6))
    r = math.hypot(*du)
    if r < 1e-9:
        continue
    al = math.atan2(du[1], du[0])
    for j, jp in itertools.product(range(6), repeat=2):
        # v_j - v_jp at phi = 0, its direction rotates with phi
        dv0 = (math.cos(TWO * j / 6) - math.cos(TWO * jp / 6),
               math.sin(TWO * j / 6) - math.sin(TWO * jp / 6))
        s = math.hypot(*dv0)
        if s < 1e-9:
            continue
        be = math.atan2(dv0[1], dv0[0])
        c = (1 - r * r - s * s) / (2 * r * s)
        if abs(c) > 1:
            continue
        for sign in (1, -1):
            phi = (al - be - sign * math.acos(c)) % (TWO / 6)
            cands.add(round(phi, 9))
cands = sorted(cands)
print(f"{len(cands)} candidate angles, solved not scanned", flush=True)

t0, hits = time.time(), []
for phi in cands:
    H1, H2 = hexes(phi)
    cls = {0: [], 1: []}
    for i, u in enumerate(H1):
        for j, v in enumerate(H2):
            q = (u[0] + v[0], u[1] + v[1])
            if q[0] ** 2 + q[1] ** 2 > 1e-12:
                cls[(i + j) % 2].append(q)
    info = []
    for c in (0, 1):
        P = cls[c]
        n = len(P)
        adj = [set() for _ in range(n)]
        for a in range(n):
            for b in range(a + 1, n):
                d = (P[a][0] - P[b][0]) ** 2 + (P[a][1] - P[b][1]) ** 2
                if abs(d - 1.0) < 1e-7:
                    adj[a].add(b)
                    adj[b].add(a)
        m = sum(len(x) for x in adj) // 2
        # chromatic number, exactly, on at most 18 vertices
        chi = 1
        while chi < 6:
            ok = [None] * n

            def go(v):
                if v == n:
                    return True
                for col in range(chi):
                    if all(ok[w] != col for w in adj[v] if w < v):
                        ok[v] = col
                        if go(v + 1):
                            return True
                        ok[v] = None
                return False
            if go(0):
                break
            chi += 1
        info.append((chi, n, m))
    lo = min(x[0] for x in info)
    if True:
        hits.append((phi, info))
        print(f"  phi = {math.degrees(phi):9.5f} deg: classes {info}  "
              f"[{time.time()-t0:.0f}s]", flush=True)
print(f"{len(hits)} angles give BOTH parity classes chromatic number >= 3  "
      f"[{time.time()-t0:.0f}s]", flush=True)
mx = max((min(x[0] for x in info), phi, info) for phi, info in hits) if hits else None
print(f"best: {mx}", flush=True)
