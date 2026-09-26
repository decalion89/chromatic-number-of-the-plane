"""G u rho_p(G): the move that produced forcing at four colours, tried at five.

Sa carries no forced pair; Y = Sa u rho(Sa) carries exactly one.  So the union
is where forcing is born, and the honest next rung is to apply it to G itself.

Two passes.  First the cheap one: is any union outright 6-chromatic?  That is
one SAT call per union.  Then, on the unions that survive, the pair scan --
a pair forced at five colours, at a distance the field closes, ends it.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
from fractions import Fraction as Fr
from math import isqrt
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD, Rotation
from hn.homcol import closable_distance, closing_radicand
from pysat.solvers import Solver

t0 = time.time()
F = DEGREY_FIELD


def rot_closing(D):
    D = Fr(D)
    sq = closing_radicand(D)
    if not closable_distance(D):
        return None
    c = 1 - Fr(1, 2) / D
    k2 = (1 - c * c) / sq
    num, den = k2.numerator, k2.denominator
    t = isqrt(num * den)
    if t * t != num * den:
        return None
    r = Rotation(F.rational(c), F.sqrt(sq) * F.rational(Fr(t, den)))
    assert r.cos * r.cos + r.sin * r.sin == F.rational(1)
    return r


G = build_G(F, as_graph=False)
n0 = len(G)
zs0 = [(float(p.x), float(p.y)) for p in G]
print(f"G: {n0} vertices  [{time.time()-t0:.0f}s]", flush=True)


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


DS = [Fr(3), Fr(4), Fr(7), Fr(9), Fr(16), Fr(2), Fr(14), Fr(19), Fr(25)]
DS = [D for D in DS if rot_closing(D) is not None]
print(f"  closing rotations available for D = {[str(D) for D in DS]}",
      flush=True)

rings = {}
for D in DS:
    r = Counter()
    fd = float(D)
    for i in range(n0):
        ai, bi = zs0[i]
        for j in range(n0):
            if i != j and abs((ai - zs0[j][0]) ** 2
                              + (bi - zs0[j][1]) ** 2 - fd) < 1e-7:
                r[i] += 1
    rings[D] = r
    print(f"  D = {D}: {len(r)} pivots, largest ring "
          f"{r.most_common(1)[0][1] if r else 0}  [{time.time()-t0:.0f}s]",
          flush=True)

jobs = []
for D in DS:
    for pi_, sz in rings[D].most_common(10):
        jobs.append((D, pi_, sz))
print(f"\n{len(jobs)} unions to try  [{time.time()-t0:.0f}s]", flush=True)

survivors = []
for D, pi_, sz in jobs:
    turn = rot_closing(D).about(G[pi_])
    pts, seen = list(G), set(G)
    for p in G:
        q = turn(p)
        if q not in seen:
            seen.add(q)
            pts.append(q)
    zs = [(float(p.x), float(p.y)) for p in pts]
    E = edges_of(pts, zs)
    cls = [[1 + v * 5 + c for c in range(5)] for v in range(len(pts))]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    with Solver(name="cd19", bootstrap_with=cls) as sv:
        ok = sv.solve()
    cross = len(E) - 2 * 7877
    print(f"  D={D} pivot {pi_} ring {sz}: {len(pts)} pts, {len(E)} edges "
          f"({cross:+d} cross), 5-colourable {ok}  [{time.time()-t0:.0f}s]",
          flush=True)
    if not ok:
        print(f"  *** SIX COLOURS: D={D} pivot {pi_} ***", flush=True)
        import pickle
        with open("/tmp/hn/six.pkl",
                  "wb") as fh:
            pickle.dump((len(pts), E), fh)
        sys.exit(0)
    survivors.append((D, pi_, sz, len(pts), len(E), cross))

survivors.sort(key=lambda t: -t[5])
print(f"\nrichest unions by cross edges: {survivors[:6]}  "
      f"[{time.time()-t0:.0f}s]", flush=True)
for D, pi_, sz, np_, ne, cross in survivors[:6]:
    turn = rot_closing(D).about(G[pi_])
    pts, seen = list(G), set(G)
    for p in G:
        q = turn(p)
        if q not in seen:
            seen.add(q)
            pts.append(q)
    zs = [(float(p.x), float(p.y)) for p in pts]
    E = edges_of(pts, zs)
    cls = [[1 + v * 5 + c for c in range(5)] for v in range(len(pts))]
    for a, b in E:
        for c in range(5):
            cls.append([-(1 + a * 5 + c), -(1 + b * 5 + c)])
    sv = Solver(name="cd19", bootstrap_with=cls)
    sv.solve()
    hits, hard, seenp = [], [], 0
    for i in range(len(pts)):
        ai, bi = zs[i]
        for j in range(i + 1, len(pts)):
            v = (ai - zs[j][0]) ** 2 + (bi - zs[j][1]) ** 2
            if v > 36.0:
                continue
            DD = Fr(round(v * 1584), 1584)
            if abs(float(DD) - v) > 1e-7 or DD == 1 \
                    or not closable_distance(DD):
                continue
            seenp += 1
            sv.conf_budget(40000)
            r = sv.solve_limited(assumptions=[1 + i * 5, -(1 + j * 5)])
            if r is False:
                hits.append((i, j, DD))
                print(f"  *** FORCED AT FIVE: D={D} pivot {pi_}, "
                      f"pair {i},{j} at {DD} ***  [{time.time()-t0:.0f}s]",
                      flush=True)
            elif r is None:
                hard.append((i, j, DD))
    sv.delete()
    print(f"  D={D} pivot {pi_}: {seenp} pairs, {len(hits)} forced, "
          f"{len(hard)} hard  [{time.time()-t0:.0f}s]", flush=True)
