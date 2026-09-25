"""Cap-test on balls, not on the whole graph: the sound direction is the cheap
one.

Screening a ring on the full universe costs minutes because the whole graph
must be coloured.  It does not have to be.  A cap is "this ring cannot show all
k colours", and that property travels UPWARD: any colouring of a supergraph
restricts to one of a subgraph, so if a BALL around the ring already cannot
show them, neither can anything containing it.

So the search becomes: take the ring, grow a ball around it inside the
universe, and test.  A positive transfers to the full graph and to every
supergraph of it; a negative on a ball says nothing and simply means growing
the ball.  Each test is tiny, so hundreds of rings can be screened in the time
one full-graph screen takes.

This is the same asymmetry that makes racing a solver worthwhile, used on the
other axis.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict, deque
from hn.degrey import build_Sa, build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from hn.homcol import closable_distance
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
SEED = sys.argv[2] if len(sys.argv) > 2 else "Sa"
RADIUS = float(sys.argv[3]) if len(sys.argv) > 3 else 4.5
HOPS = [1, 2, 3, 4]
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
half = K.rational(Fr(1, 2))
ZERO = Point(K.zero(), K.zero())


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


base = build_Sa(K) if SEED == "Sa" else build_G(K, as_graph=False)
P, have, fresh = list(base), set(base), []
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
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
print(f"{SEED} universe: {n} points, {len(E)} edges, k={k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

dm, D2 = b.dim, b.D * b.D
deg = {v: len(adj[v]) for v in range(n)}
centres = sorted(range(n), key=lambda v: -deg[v])[:40]
print(f"testing {len(centres)} centres by degree, top {deg[centres[0]]}"
      f"  [{time.time()-t0:.0f}s]", flush=True)


def ball(seeds, hops):
    cur = set(seeds)
    for _ in range(hops):
        nxt = set(cur)
        for v in cur:
            nxt |= adj[v]
        cur = nxt
    return sorted(cur)


def capped_on(sub, ring):
    idx = {v: i for i, v in enumerate(sub)}
    m = len(sub)
    ss = set(sub)
    cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        if a in ss and c in ss:
            for col in range(k):
                cls.append([-(1 + idx[a] * k + col), -(1 + idx[c] * k + col)])
    nv = m * k
    ind = [nv + 1 + c for c in range(k)]
    enc = CardEnc.atleast(lits=list(ind), bound=k, top_id=nv + k,
                          encoding=EncType.seqcounter)
    aux = [[-ind[c]] + [1 + idx[v] * k + c for v in ring] for c in range(k)]
    s = Solver(name="cd15", bootstrap_with=cls + aux + list(enc.clauses))
    out = s.solve()
    s.delete()
    return not out


tested = hits = 0
for ci in centres:
    d = r - r[ci]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    grp = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        v = Fr(int(sq[off, 0]), D2)
        if v and v != 1 and closable_distance(v):
            grp[v].append(int(off))
    for D, mem in grp.items():
        if len(mem) < k:
            continue
        for hops in HOPS:
            sub = ball(mem + [ci], hops)
            if len(sub) > 4000:
                break
            tested += 1
            if capped_on(sub, mem):
                hits += 1
                print(f"*** CAPPED on a ball: centre {ci}, ring D={D} "
                      f"({len(mem)} pts), ball {len(sub)} points, {hops} hops "
                      f"-- transfers to every supergraph ***"
                      f"  [{time.time()-t0:.0f}s]", flush=True)
                break
    if ci == centres[len(centres) // 4]:
        print(f"   a quarter through, {tested} ball tests, {hits} caps"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n{tested} ball tests over {len(centres)} centres, {hits} capped"
      f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
