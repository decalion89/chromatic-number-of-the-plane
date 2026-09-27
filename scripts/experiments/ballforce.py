"""Forced pairs on balls: the same trick, aimed at what actually matters.

A capped ring gives a disjunction and needs a bite to sharpen.  A FORCED PAIR
needs nothing: the spindle converts it straight into one more colour.  Every
forced-pair census in this work ran on whole graphs, which is why it only ever
ran on G.

But forcing travels upward exactly as capping does.  "u and v share a colour in
every proper k-colouring" is "the graph plus the edge uv does not colour", and
adding vertices only removes colourings -- so if a BALL around u and v plus
that edge already fails to colour, the whole graph fails, and so does every
supergraph of it.

So the census can be local.  Take a pair at a closable distance, grow a ball
around it, add the edge, ask for a colouring.  Unsatisfiable on the ball proves
the pair forced in everything containing it; satisfiable says only that the
ball is too small.  Seconds per pair instead of minutes, which is the
difference between censusing one graph and censusing a universe.

Calibrated first at four colours, where de Grey's pair is known to be forced.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_Sa, build_G, build_Y
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SEED = sys.argv[2] if len(sys.argv) > 2 else "Y"
RADIUS = float(sys.argv[3]) if len(sys.argv) > 3 else 4.0
MAXBALL = int(sys.argv[4]) if len(sys.argv) > 4 else 2500
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))


def rsqrt(e):
    if any(x for x in e.c[1:]):
        return None
    v = e.c[0]
    if v <= 0:
        return None
    for rr in CLASSES:
        q = v / rr
        nn, dd = q.numerator, q.denominator
        rn, rd = int(round(nn ** .5)), int(round(dd ** .5))
        if rn * rn == nn and rd * rd == dd:
            s = K.rational(Fr(rn, rd))
            return s if rr == 1 else K.sqrt(rr) * s
    return None


builders = {"Sa": build_Sa, "Y": build_Y,
            "G": lambda f: build_G(f, as_graph=False)}
P = list(builders[SEED](K))
have, fresh = set(P), []
n0 = len(P)
for i in range(n0):
    A = P[i]
    for j in range(i + 1, n0):
        B = P[j]
        D = A.dist2(B)
        if not .05 < float(D) < 3.99:
            continue
        sD, s4 = rsqrt(D), rsqrt(K.rational(4) - D)
        if sD is None or s4 is None:
            continue
        inv = K.rational(1) / sD
        mx, my = (A.x + B.x) * half, (A.y + B.y) * half
        nx = -(B.y - A.y) * inv * s4 * half
        ny = (B.x - A.x) * inv * s4 * half
        for q in (Point(mx + nx, my + ny), Point(mx - nx, my - ny)):
            if q in have or float(q.x) ** 2 + float(q.y) ** 2 > RADIUS ** 2:
                continue
            have.add(q)
            fresh.append(q)
P = P + fresh
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
Eset = set(E)
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
deg = {v: len(adj[v]) for v in range(n)}
print(f"{SEED} universe: {n} points, {len(E)} edges, k={k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

dm, D2 = b.dim, b.D * b.D
hot = sorted(range(n), key=lambda v: -deg[v])[:25]
pairs = []
for u in hot:
    d = r - r[u]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    for off in np.nonzero(ok)[0]:
        v = int(off)
        if v == u or (min(u, v), max(u, v)) in Eset:
            continue
        D = Fr(int(sq[v, 0]), D2)
        if D and closable_distance(D) and deg[v] >= 8:
            pairs.append((D, min(u, v), max(u, v)))
pairs = sorted(set(pairs), key=lambda t: (-min(deg[t[1]], deg[t[2]]), t[0]))
print(f"{len(pairs)} candidate pairs from {len(hot)} busiest vertices"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def ball(seeds, hops):
    cur = set(seeds)
    for _ in range(hops):
        nxt = set(cur)
        for v in cur:
            nxt |= adj[v]
        cur = nxt
    return sorted(cur)


def forced_on(sub, u, v):
    idx = {x: i for i, x in enumerate(sub)}
    m = len(sub)
    ss = set(sub)
    cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        if a in ss and c in ss:
            for col in range(k):
                cls.append([-(1 + idx[a] * k + col), -(1 + idx[c] * k + col)])
    for col in range(k):          # the extra edge: u and v must differ
        cls.append([-(1 + idx[u] * k + col), -(1 + idx[v] * k + col)])
    s = Solver(name="cd15", bootstrap_with=cls)
    out = s.solve()
    s.delete()
    return not out


tested = hits = 0
for D, u, v in pairs[:400]:
    for hops in (1, 2, 3):
        sub = ball([u, v], hops)
        if len(sub) > MAXBALL:
            break
        tested += 1
        if forced_on(sub, u, v):
            hits += 1
            print(f"*** FORCED on a ball: pair ({u},{v}) at D={D}, ball "
                  f"{len(sub)} points, {hops} hops -- forced in every "
                  f"supergraph  [{time.time()-t0:.0f}s]", flush=True)
            break
    if tested % 50 == 0 and tested:
        print(f"   {tested} ball tests, {hits} forced"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tested} ball tests, {hits} forced pairs  [{time.time()-t0:.0f}s]",
      flush=True)
print("DONE", flush=True)
