"""The independence ratio of the Erdos lattice unit-distance graphs.

They were dropped as colouring candidates -- every primitive one is
5-colourable by cosets -- but their ratio was never measured, and that is a
different question with a different use.  For any finite unit-distance graph
H, m_1 <= alpha(H)/n(H), so a low ratio bounds the density of measurable
distance-avoiding sets, and a ratio under 0.25 proves measurable chi >= 5.

These graphs reach degree 16, 32, 64 where this project's best cores manage
12, and degree is the lever on independence.  Croft's construction puts an
absolute floor of 0.229 under every unit-distance graph, the published
record is 0.2544 and the threshold is 0.25, so the whole game sits in a band
nine per cent wide.

Greedy understates alpha badly on dense graphs -- by 21 per cent on the last
candidates -- so every ratio here is greedy followed by (1,2)-swap local
search, and only the local-search value is reported.
"""
import sys, time, random
from itertools import product
from math import gcd

t0 = time.time()


def conn_set(N):
    out = set()
    a = 0
    while a * a <= N:
        b2 = N - a * a
        b = int(b2 ** 0.5)
        for bb in (b - 1, b, b + 1):
            if bb >= 0 and a * a + bb * bb == N:
                for sa, sb in product((1, -1), repeat=2):
                    out.add((sa * a, sb * bb))
                    out.add((sb * bb, sa * a))
        a += 1
    out.discard((0, 0))
    g = 0
    for x, y in out:
        g = gcd(g, gcd(abs(x), abs(y)))
    if g > 1:
        out = {(x // g, y // g) for x, y in out}
    return sorted(out)


def ratio(N, side):
    S = conn_set(N)
    if not S:
        return None
    R = max(max(abs(x), abs(y)) for x, y in S)
    if side < 2 * R + 4:
        side = 2 * R + 4
    pts = [(x, y) for x in range(side) for y in range(side)]
    idx = {p: i for i, p in enumerate(pts)}
    n = len(pts)
    if n > 60000:
        return ("too big", len(S), n)
    adj = [set() for _ in range(n)]
    for (x, y) in pts:
        i = idx[(x, y)]
        for (dx, dy) in S:
            q = (x + dx, y + dy)
            j = idx.get(q)
            if j is not None:
                adj[i].add(j)
                adj[j].add(i)
    m = sum(len(a) for a in adj) // 2
    random.seed(5)
    best = 0
    for restart in range(6):
        order = sorted(range(n), key=lambda v: (len(adj[v]), random.random()))
        cur, banned = set(), set()
        for v in order:
            if v not in banned:
                cur.add(v)
                banned |= adj[v] | {v}
        cnt = [0] * n
        for v in cur:
            for u in adj[v]:
                cnt[u] += 1
        for it in range(40000):
            free = [v for v in range(n) if v not in cur and cnt[v] == 0]
            if free:
                v = random.choice(free)
                cur.add(v)
                for u in adj[v]:
                    cnt[u] += 1
                continue
            ones = [v for v in range(n) if v not in cur and cnt[v] == 1]
            if not ones:
                break
            v = random.choice(ones)
            w = next(u for u in adj[v] if u in cur)
            cur.discard(w)
            for u in adj[w]:
                cnt[u] -= 1
            cur.add(v)
            for u in adj[v]:
                cnt[u] += 1
        best = max(best, len(cur))
    return (best / n, len(S), n, 2 * m / n)


print(f"{'N':>8s} {'degree':>7s} {'points':>7s} {'mean deg':>9s} "
      f"{'ratio':>8s}   (floor 0.229, record 0.2544, threshold 0.25)",
      flush=True)
best = (1.0, None)
for N in (25, 65, 325, 1105, 5525, 27625, 145, 725, 3625, 169, 2197,
          85, 425, 2125, 221, 1445, 4225):
    r = ratio(N, 40)
    if r is None or r[0] == "too big":
        continue
    rat, deg, n, md = r
    mark = "  <<<" if rat < 0.25 else ""
    print(f"{N:8d} {deg:7d} {n:7d} {md:9.2f} {rat:8.4f}{mark}"
          f"   [{time.time()-t0:.0f}s]", flush=True)
    if rat < best[0]:
        best = (rat, N)
print(f"\nlowest: {best[0]:.4f} at N = {best[1]}", flush=True)
print("DONE", flush=True)
