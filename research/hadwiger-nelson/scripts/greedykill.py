"""Grow a graph towards uncolourability by killing colourings, not by volume.

Every universe built here grew by brute force -- take all the unit-circle
intersections inside a radius and hope.  Thirty-five thousand points later
the answer was still a colouring, because most of those points constrain
nothing: they sit where the colouring has room.

There is a sharper criterion, and it is exact.  A new point p cannot be
coloured in a colouring c precisely when its neighbourhood already shows all
five colours under c.  So adding p destroys exactly those sampled colourings,
and the point worth adding is the one that destroys the most.

That makes growth a greedy optimisation over a cheap score instead of a
sweep: keep a sample of proper colourings, score every candidate intersection
point by how many of them its neighbourhood saturates, add the best, drop the
colourings it killed, refill the sample from the enlarged graph, repeat.  The
sample shrinks and is replenished, and the graph grows only where colourings
actually live.

If the sample ever cannot be refilled, the graph has no proper 5-colouring.
"""
import sys, time, pickle, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
sys.path.insert(0, "/tmp/claude-0/-home-user-darwin-50/"
                   "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad")
import numpy as np
from fractions import Fraction as Fr
from collections import defaultdict
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete
from pysat.solvers import Solver

k = int(sys.argv[1]) if len(sys.argv) > 1 else 5
NSAMP = int(sys.argv[2]) if len(sys.argv) > 2 else 300
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 4000
RAD = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
t0 = time.time()
CLASSES = [1, 3, 5, 7, 11, 15, 21, 33, 35, 55, 77, 105, 165, 231, 385, 1155]
SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
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


P = list(build_G(K, as_graph=False))
have = set(P)
cx = sum(float(p.x) for p in P) / len(P)
cy = sum(float(p.y) for p in P) / len(P)
print(f"seed {len(P)} points, centroid ({cx:.3f},{cy:.3f})"
      f"  [{time.time()-t0:.0f}s]", flush=True)

# Candidate pool: unit-circle intersections of nearby pairs, generated once.
cand = []
for i, A in enumerate(P):
    for B in P[i + 1:]:
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
            if q in have:
                continue
            if (float(q.x) - cx) ** 2 + (float(q.y) - cy) ** 2 > RAD ** 2:
                continue
            have.add(q)
            cand.append(q)
print(f"{len(cand)} candidate points  [{time.time()-t0:.0f}s]", flush=True)
ALL = P + cand
b = IntBasis.covering(ALL)
rows = b.rows(ALL)
assert b.overflow_headroom(rows) < 1.0
E = sorted(set((min(a, c), max(a, c))
               for a, c in fast_edges_complete(b, rows)))
nb = defaultdict(set)
for a, c in E:
    nb[a].add(c)
    nb[c].add(a)
n0 = len(P)
N = len(ALL)
print(f"{N} points, {len(E)} edges in the closure  [{time.time()-t0:.0f}s]",
      flush=True)

chosen = set(range(n0))
rng = random.Random(433494437)


def make_solver(S):
    idx = {v: i for i, v in enumerate(sorted(S))}
    m = len(idx)
    cls = [[1 + i * k + c for c in range(k)] for i in range(m)]
    for a, c in E:
        if a in S and c in S:
            for col in range(k):
                cls.append([-(1 + idx[a] * k + col), -(1 + idx[c] * k + col)])
    return Solver(name="cd15", bootstrap_with=cls), idx, m


def sample(S, want):
    sv, idx, m = make_solver(S)
    out = []
    if not sv.solve():
        sv.delete()
        return None
    for _ in range(want):
        sv.set_phases([(1 if rng.random() < .5 else -1) * (1 + w)
                       for w in range(m * k)])
        sv.solve()
        mo = sv.get_model()
        col = {}
        for v, i in idx.items():
            col[v] = next(c for c in range(k) if mo[i * k + c] > 0)
        out.append(col)
    sv.delete()
    return out


samples = sample(chosen, NSAMP)
print(f"{len(samples)} colourings sampled  [{time.time()-t0:.0f}s]",
      flush=True)
for step in range(STEPS):
    best, bestv = None, -1
    for v in range(n0, N):
        if v in chosen:
            continue
        nbs = [u for u in nb[v] if u in chosen]
        if len(nbs) < k:
            continue
        cnt = 0
        for col in samples:
            s = set()
            for u in nbs:
                s.add(col[u])
                if len(s) == k:
                    cnt += 1
                    break
        if cnt > bestv:
            best, bestv = v, cnt
    if best is None:
        print("no candidate has k neighbours inside -- stop", flush=True)
        break
    chosen.add(best)
    samples = [c for c in samples
               if len({c[u] for u in nb[best] if u in chosen and u != best})
               < k]
    if len(samples) < NSAMP // 3:
        fresh = sample(chosen, NSAMP)
        if fresh is None:
            print(f"\n*** NOT {k}-COLOURABLE at {len(chosen)} points "
                  f"[{time.time()-t0:.0f}s]", flush=True)
            pickle.dump([ALL[v] for v in sorted(chosen)],
                        open(SC + "WITNESS_greedy.pkl", "wb"))
            break
        samples = fresh
    if step % 25 == 0:
        print(f"   step {step}: {len(chosen)} points, added one killing "
              f"{bestv} of the sample, {len(samples)} colourings left"
              f"  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
