"""Ring alone, or ring with its centre?  The two are not the same cap.

The binary cap test asks whether a RING can show all k colours.  de Grey's
lemma is about the ring TOGETHER WITH its centre, and bounds that set by two
out of four.  Sixty-eight rings in Sa's universe came back capped at four
colours; a first pass at the centre+ring sets found none capped at all, which
would mean adding the single centre point destroys every one of them.

That is worth getting exactly right, because it decides what the instrument
has been measuring.  So: the same balls, the same rings, both palettes side by
side, and the depth recorded at whatever hop count the ball limit allows
rather than only at three.

A palette bound proved on a ball holds in every supergraph, both versions
alike: adding vertices only removes colourings, so the largest palette a set
can display only falls.
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
from hn.homcol import ring_palette_bound
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 4
SEED = sys.argv[2] if len(sys.argv) > 2 else "Sa"
RADIUS = float(sys.argv[3]) if len(sys.argv) > 3 else 4.5
NCENTRES = int(sys.argv[4]) if len(sys.argv) > 4 else 40
MAXBALL = int(sys.argv[5]) if len(sys.argv) > 5 else 3000
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
n = len(P)
adj = defaultdict(set)
for a, c in E:
    adj[a].add(c)
    adj[c].add(a)
deg = {v: len(adj[v]) for v in range(n)}
print(f"{SEED} universe: {n} points, {len(E)} edges, k={k}"
      f"  [{time.time()-t0:.0f}s]", flush=True)

dm, D2 = b.dim, b.D * b.D
centres = sorted(range(n), key=lambda v: -deg[v])[:NCENTRES]


def rings_of(cidx):
    d = r - r[cidx]
    sq = b._field_square(d[:, :dm]) + b._field_square(d[:, dm:])
    ok = np.ones(len(sq), dtype=bool)
    for j in range(1, dm):
        ok &= sq[:, j] == 0
    byD = defaultdict(list)
    for off in np.nonzero(ok)[0]:
        if off == cidx:
            continue
        D = Fr(int(sq[off, 0]), D2)
        if 0 < D <= 4:
            byD[D].append(int(off))
    return {D: v for D, v in byD.items() if 4 <= len(v) <= 40}


def ball(seeds, hops):
    cur = set(seeds)
    for _ in range(hops):
        nxt = set(cur)
        for v in cur:
            nxt |= adj[v]
        cur = nxt
    return sorted(cur)


def make_base(sub):
    idx = {x: i for i, x in enumerate(sub)}
    m, ss = len(sub), set(sub)
    base = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        if a in ss and c in ss:
            for col in range(k):
                base.append([-(1 + idx[a] * k + col), -(1 + idx[c] * k + col)])
    return base, idx, m


def depth(base, idx, m, watch):
    """Largest palette the set can show, via the library's own encoding."""
    return ring_palette_bound(base, m * k, [idx[w] for w in watch], k)


joint = defaultdict(int)
tested = 0
for ci in centres:
    for D, ring in sorted(rings_of(ci).items()):
        sub, used_hops = None, 0
        for hops in (3, 2, 1):            # deepest ball the limit allows
            cand = ball([ci] + ring, hops)
            if len(cand) <= MAXBALL:
                sub, used_hops = cand, hops
                break
        if sub is None:
            continue
        base, idx, m = make_base(sub)
        dr = depth(base, idx, m, ring)
        dc = depth(base, idx, m, [ci] + ring)
        joint[(dr, dc)] += 1
        tested += 1
        if dc <= k - 2:
            print(f"*** centre+ring at depth {dc} (de Grey's is 2 of 4): "
                  f"centre {ci}, D={D}, {len(ring)} pts, ball {len(sub)}, "
                  f"{used_hops} hops  [{time.time()-t0:.0f}s]", flush=True)
    if ci == centres[len(centres) // 2]:
        print(f"   halfway, {tested} rings, joint(ring,centre+ring)="
              f"{dict(sorted(joint.items()))}  [{time.time()-t0:.0f}s]",
              flush=True)
print(f"\n{tested} rings at k={k}.  joint (ring depth, centre+ring depth):",
      flush=True)
for key, cnt in sorted(joint.items()):
    print(f"   ring {key[0]}, centre+ring {key[1]}: {cnt}", flush=True)
capped_r = sum(c for key, c in joint.items() if key[0] < k)
capped_c = sum(c for key, c in joint.items() if key[1] < k)
print(f"rings capped: {capped_r}   centre+rings capped: {capped_c}   "
      f"at de Grey depth: "
      f"{sum(c for key, c in joint.items() if key[1] <= k-2)}", flush=True)
print("DONE", flush=True)
