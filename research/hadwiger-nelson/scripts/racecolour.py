"""Race the solver with a local search, because waiting settles nothing.

A long CDCL run does not distinguish "unsatisfiable" from "satisfiable and
hard", and this work has the receipts: a ball instance here climbed from 21 to
6893 seconds and was SATISFIABLE, and a local-search plateau was twice read as
evidence of unsatisfiability and twice falsified.  Cost is not evidence.

So resolve it from the other side.  TabuCol lands on a colouring or it does
not; landing PROVES satisfiability outright and ends the wait, while failing
proves nothing at all and is reported as proving nothing.  The asymmetry is the
whole point and it is why this is worth running in parallel rather than
instead.

Same universe as the solver is chewing on: G closed once under the
unit-circle-intersection generator.
"""
import sys, time, random
sys.path.insert(0, "/home/user/darwin-50/research/hadwiger-nelson")
import numpy as np
from fractions import Fraction as Fr
from hn.degrey import build_G
from hn.geometry import DEGREY_FIELD as K, Point
from hn.fast import IntBasis, fast_edges_complete

k = 5
RADIUS = 4.5
CAP = 26000
t0 = time.time()
rng = random.Random(9)
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


P = build_G(K, as_graph=False)
have = set(P)
fresh = []
n0 = len(P)
for i in range(n0):
    A = P[i]
    for j in range(i + 1, n0):
        B = P[j]
        D = A.dist2(B)
        f = float(D)
        if not .05 < f < 3.99:
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
            if float(q.x) ** 2 + float(q.y) ** 2 > RADIUS * RADIUS:
                continue
            have.add(q)
            fresh.append(q)
        if n0 + len(fresh) > CAP:
            break
    if n0 + len(fresh) > CAP:
        break
P = P + fresh
b = IntBasis.covering(P)
r = b.rows(P)
assert b.overflow_headroom(r) < 1.0
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, r)))
n = len(P)
print(f"{n} points, {len(E)} edges ({len(E)/n:.2f}/v)"
      f"  [{time.time()-t0:.0f}s]", flush=True)

adj = [[] for _ in range(n)]
for a, c in E:
    adj[a].append(c)
    adj[c].append(a)
adj = [np.array(a, dtype=np.int32) for a in adj]

best_overall = None
for restart in range(1, 200):
    col = np.array([rng.randrange(k) for _ in range(n)], dtype=np.int8)
    cnt = np.zeros((n, k), dtype=np.int32)
    for v in range(n):
        if len(adj[v]):
            np.add.at(cnt[v], col[adj[v]], 1)
    confl = int(sum(cnt[v][col[v]] for v in range(n)) // 2)
    tabu = np.zeros((n, k), dtype=np.int64)
    it = 0
    while confl > 0 and it < 60000:
        it += 1
        bad = [v for v in range(n) if cnt[v][col[v]]]
        if not bad:
            break
        bestd, moves = 10 ** 9, []
        for v in bad[:400]:
            cv = col[v]
            base = cnt[v][cv]
            for c in range(k):
                if c == cv:
                    continue
                d = cnt[v][c] - base
                if tabu[v][c] > it and not (confl + d == 0):
                    continue
                if d < bestd:
                    bestd, moves = d, [(v, c)]
                elif d == bestd:
                    moves.append((v, c))
        if not moves:
            break
        v, c = moves[rng.randrange(len(moves))]
        cv = col[v]
        confl += cnt[v][c] - cnt[v][cv]
        if len(adj[v]):
            np.add.at(cnt[:, cv], adj[v], -1)
            np.add.at(cnt[:, c], adj[v], 1)
        col[v] = c
        tabu[v][cv] = it + int(0.6 * confl) + rng.randrange(10)
        if it % 5000 == 0:
            print(f"   restart {restart}, iter {it}, conflicts {confl}"
                  f"  [{time.time()-t0:.0f}s]", flush=True)
    if best_overall is None or confl < best_overall:
        best_overall = confl
    if confl == 0:
        bad = [(a, c) for a, c in E if col[a] == col[c]]
        assert not bad, "claimed colouring has a monochromatic edge"
        print(f"*** SATISFIABLE: proper {k}-colouring found on restart "
              f"{restart}, verified over {len(E)} edges -- the solver's long "
              f"run means nothing ***  [{time.time()-t0:.0f}s]", flush=True)
        sys.exit()
    print(f"  restart {restart}: stuck at {confl} conflicts (best {best_overall})"
          f"  [{time.time()-t0:.0f}s]", flush=True)
print(f"\nno colouring found -- which proves NOTHING, by this work's own "
      f"record  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
