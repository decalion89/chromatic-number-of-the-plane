"""Sculpt a unit-distance graph towards a low independence ratio.

m_1 <= alpha(H)/n(H) for any finite unit-distance graph, so a graph with
ratio below 1/4 proves the measurable chromatic number is at least 5 by
density alone.  The record is 0.2544, the threshold 0.25, and Croft's
construction puts an absolute floor of 0.229 under every candidate.  This
corpus sits at 0.2637.

Removing a vertex that lies in a maximum independent set lowers the ratio
whenever alpha < n, since (alpha-1)/(n-1) < alpha/n.  Alpha only falls if
EVERY maximum independent set loses a member, so the move is to strike at
the vertices that appear in many large independent sets: sample independent
sets by local search, count participation, and remove the most popular.

Greedy is not trusted anywhere here -- it understated alpha by 21 per cent
on these graphs -- so every ratio quoted is from (1,2)-swap local search.
"""
import os
HN_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import sys, time, pickle, random
from collections import Counter
sys.path[:0] = [HN_DIR, os.path.join(HN_DIR, "scripts")]
from hn.fast import IntBasis, fast_edges_complete

SC = ("/tmp/claude-0/-home-user-darwin-50/"
      "aceaa9ec-f432-5848-a506-39c59179b415/scratchpad/")
NAME = sys.argv[1] if len(sys.argv) > 1 else "fiveD.pkl"
KCORE = int(sys.argv[2]) if len(sys.argv) > 2 else 8
STEPS = int(sys.argv[3]) if len(sys.argv) > 3 else 400
BATCH = int(sys.argv[4]) if len(sys.argv) > 4 else 10
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
adj = {v: (adj0[v] & keep) for v in keep}
random.seed(17)


def big_independent(alive, tries=8, iters=15000):
    """Largest independent set found, plus a participation count."""
    part = Counter()
    best, bestset = 0, None
    al = sorted(alive)
    for t in range(tries):
        order = sorted(al, key=lambda v: (len(adj[v] & alive), random.random()))
        cur, banned = set(), set()
        for v in order:
            if v not in banned:
                cur.add(v)
                banned |= (adj[v] & alive) | {v}
        cnt = Counter()
        for v in cur:
            for u in adj[v] & alive:
                cnt[u] += 1
        for it in range(iters):
            free = [v for v in al if v not in cur and cnt[v] == 0]
            if free:
                v = random.choice(free)
                cur.add(v)
                for u in adj[v] & alive:
                    cnt[u] += 1
                continue
            ones = [v for v in al if v not in cur and cnt[v] == 1]
            if not ones:
                break
            v = random.choice(ones)
            w = next(u for u in adj[v] & alive if u in cur)
            cur.discard(w)
            for u in adj[w] & alive:
                cnt[u] -= 1
            cur.add(v)
            for u in adj[v] & alive:
                cnt[u] += 1
        part.update(cur)
        if len(cur) > best:
            best, bestset = len(cur), set(cur)
    return best, part


alive = set(keep)
a0, part = big_independent(alive)
print(f"{NAME} {KCORE}-core: {len(alive)} points; alpha >= {a0}, ratio "
      f"{a0/len(alive):.4f}  [{time.time()-t0:.0f}s]", flush=True)
best = a0 / len(alive)
for step in range(STEPS):
    if len(alive) < 300:
        break
    victims = [v for v, _ in part.most_common(BATCH)]
    alive -= set(victims)
    a, part = big_independent(alive)
    r = a / len(alive)
    if r < best:
        best = r
        pickle.dump([P[v] for v in sorted(alive)],
                    open(SC + f"sculpt_{NAME}", "wb"))
    if step % 5 == 0 or r < 0.25:
        print(f"   step {step}: {len(alive)} points, alpha >= {a}, ratio "
              f"{r:.4f}{'   <<< BELOW 0.25' if r < 0.25 else ''}"
              f"  [{time.time()-t0:.0f}s]", flush=True)
    if r < 0.25:
        print("   *** below the threshold -- needs an exact alpha to count "
              "***", flush=True)
        break
print(f"\nbest ratio reached {best:.4f} (threshold 0.25, record 0.2544, "
      f"floor 0.229)  [{time.time()-t0:.0f}s]", flush=True)
print("DONE", flush=True)
