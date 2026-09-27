"""Push alpha UPWARD by local search, to kill candidates cheaply.

m_1 <= alpha(H)/n(H) for any finite unit-distance graph H, by averaging a
measurable distance-avoiding set over random isometries.  Four measurable
classes cover the plane only if 4*m_1 >= 1, so a graph with
alpha(H) < n(H)/4 proves the measurable chromatic number is at least 5 --
a third route to a theorem known two other ways, and a better bound than the
published m_1 <= 0.2544.

Greedy bounds alpha from below, so a greedy ratio under 0.25 is necessary
and proves nothing.  The expensive half is the UNSAT that bounds alpha from
above.  The cheap half is this: push the lower bound up with local search,
and if it crosses n/4 the candidate is dead in seconds.

The move is the standard (1,2)-swap for independent sets: drop one chosen
vertex, add two of its neighbours that are free once it goes.  Plus random
restarts and plateau walks.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/hn/")
NAME = sys.argv[1]
KCORE = int(sys.argv[2]) if len(sys.argv) > 2 else 8
ITERS = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
t0 = time.time()
P = pickle.load(open(SC + NAME, "rb"))
b = IntBasis.covering(P)
rows = b.rows(P)
E = sorted(set((min(a, c), max(a, c)) for a, c in fast_edges_complete(b, rows)))
n0 = len(P)
adj0 = [set() for _ in range(n0)]
for a, c in E:
    adj0[a].add(c)
    adj0[c].add(a)
keep = set(range(n0))
while True:
    drop = {v for v in keep if len(adj0[v] & keep) < KCORE}
    if not drop:
        break
    keep -= drop
ks = sorted(keep)
idx = {v: i for i, v in enumerate(ks)}
n = len(ks)
adj = [set() for _ in range(n)]
for v in ks:
    for u in adj0[v] & keep:
        adj[idx[v]].add(idx[u])
m = sum(len(a) for a in adj) // 2
thresh = n / 4.0
print(f"{NAME} {KCORE}-core: {n} points, {m} edges, mean degree "
      f"{2*m/n:.2f}; threshold n/4 = {thresh:.1f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
random.seed(11)
best = 0
for restart in range(40):
    order = sorted(range(n), key=lambda v: (len(adj[v]), random.random()))
    cur, banned = set(), set()
    for v in order:
        if v not in banned:
            cur.add(v)
            banned |= adj[v] | {v}
    # count how many chosen neighbours each free vertex has
    cnt = [0] * n
    for v in cur:
        for u in adj[v]:
            cnt[u] += 1
    for it in range(ITERS // 40):
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
    if len(cur) > best:
        best = len(cur)
        print(f"   restart {restart}: independent set {best}, ratio "
              f"{best/n:.4f}{'   >= n/4, candidate dead' if best >= thresh else ''}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
        if best >= thresh:
            break
print(f"\nbest independent set {best} of {n}, ratio {best/n:.4f}"
      f"  [{time.time()-t0:.0f}s]", flush=True)
if best < thresh:
    print(f"   still below n/4 = {thresh:.1f}: alpha is somewhere in "
          f"[{best}, ?] and an UNSAT at {int(thresh)+1} would settle it",
          flush=True)
else:
    print(f"   alpha >= {best} >= n/4, so this graph gives m_1 >= "
          f"{best/n:.4f} and cannot beat 0.25", flush=True)
print("DONE", flush=True)
